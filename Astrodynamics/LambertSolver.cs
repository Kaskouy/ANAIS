using SFS.World;
using SFS.WorldBase;
using System;
using System.Diagnostics;

public class LambertSolver
{
    // Input variables
    // ---------------
    private double mu;     // gravitational parameter

    private double R1;
    private double arg1;
    private double R2;
    private double arg2;

    private double T1;        // date at P1
    private double T2;        // date at P2
    private int direction;    // direction of the orbit (1: anti-clockwise; -1: clockwise)
    private uint nbTurns = 0; // the number of full turns we have to take into account (elliptic trajectories only, obviously)

    // Internal variables
    // ------------------
    private double delta;  // semi-angle (= (arg2 - arg1) / 2.0)
    private double C;      // chord; distance from P1 to P2
    private double S;      // semi-perimeter
    private double lambda; // lambda parameter (lambda = sqrt(R1*R2)*cos(delta))
    private double Tau;    // A characteristic duration used in the equation (Tau = sqrt(S^3/2mu))
    private double deltaT; // An adimensional parameter, proportional (but not equal!) to (T2 - T1)

    // A set of variables used for initial guess for low times of flight
    private bool lowToF = false;
    private double TM0 = 0.0;
    private double x_min = 0.0;
    private double TMin = 0.0;

    // Constants and internal stuff
    // ---------
    private const double C_MIN_ANGLE_DIFFERENCE = 0.0003; // 1 minute of arc (the value is in radians) - the calculation is not done if the arguments are close by that value
    private const double C_BATTIN_THRESHOLD = 0.001; // A threshold to decide wether or not we switch to Battin's method (used for nearly parabolic trajectories)
    private const uint C_NB_MAX_ITERATIONS = 12; // The maximum number of Newton/Halley/Householder iterations we allow (3 are enough in most cases - a few cases can require more)
    private const double C_PRECISION_THRESHOLD = 0.000000001; // The precision required to end Newton/Halley iterations

    /* Gives the value of gamma(n) as defined in equation 7.46 from Battin's paper - only values from T_Gamma[1] to T_Gamma[5] are useful in practice
    // The gamma(n) values are obtained through the following formula:
    // - gamma =    n * (n-3)/ (2*n+1) / (2*n+3) for n even
    // - gamma = (n+2)* (n+5)/ (2*n+1) / (2*n+3) for n odd*/
    private static readonly double[] T_Gamma = {0.0, 6.0/5.0, -2.0/35.0, 40.0/63.0, 4.0/99.0, 70.0/143.0, 6.0/65.0, 36.0/85.0, 40.0/323.0};

    private enum E_TRANSFER_TYPE
    {
        HIGH_ARC,
        LOW_ARC
    };


    // The main function: CalculateTrajectory
    // --------------------------------------
    //
    // DESCRIPTION:
    // That function allows to calculate a trajectory solving what's usually called the Lambert problem: you start from point P1 at date T1
    // and you want to be at position T2 at date T2. The Lambert solver calculates the trajectory that fulfills those constraints.
    // It's also possible to specify the number of full turns to perform (0 for a direct trajectory).
    // You can also specify the direction: clockwise (-1) or anti-clockwise (+1). There's a solution in both cases. A negative direction
    // does NOT mean going from P2 to P1, this is really about the direction of the orbit.
    // --------------------------------------
    //
    // RESTRICTIONS:
    // The Lambert problem usually always have a solution. There are a few restrictions though:
    // - T2 must be greater than T1: you can't come back to the past!
    // - The arguments of P1 and P2 mustn't be too close. The solver would give a correct result, but there would be a risk that you end up
    //   with a rectilinear orbit, which you probably don't want to deal with.
    // - When a strictly positive number of turns is required, the problem usually has two solutions, but the time difference has to be
    //   high enough. No solution exists if the time difference is too low. With no additional turn, a solution is guaranteed
    //
    // The function will return null in those 2 situations.
    // --------------------------------------
    //
    // RETURN VALUE:
    // Two orbits are returned, with all orbital parameters set. The start and end times are set to the input parameters T1 and T2, so that
    // when the orbits are displayed, they appear to be bounded to P1 and P2 (the starting and ending points).
    // Note: if nbTurns is 0, then both orbits are equal (since there's only one solution). If nbTurns > 0, two solutions usually exist,
    //       but this might not be the case if the time difference is too low. In this case, both orbits are null.
    // --------------------------------------
    public static (Orbit orbit1, Orbit orbit2) CalculateTrajectory(Planet planet, Double2 _P1, Double2 _P2, double _T1, double _T2, uint _nbTurns, int _direction)
    {
        double R1 = _P1.magnitude;
        double R2 = _P2.magnitude;
        double arg1 = _P1.AngleRadians;
        double arg2 = _P2.AngleRadians;
        
        // CONTROL PART
        // ------------
        if (!(_T2 > _T1))
        {
            // The arrival date must ABSOLUTELY be after the departure date!
            LOG(LOG_LEVEL.ERROR, "Impossible to initialize Lambert solver : T2 <= T1");
            return (null,null);
        }

        // The angular difference between the 2 positions
        double deltaTheta = arg2 - arg1;

        // deltaTheta must be between 0.0 and 2.0*PI (assuming positive direction)
        if(deltaTheta < 0.0) deltaTheta += 2.0 * Math.PI;

        if ((deltaTheta < C_MIN_ANGLE_DIFFERENCE) || (2.0 * Math.PI - deltaTheta < C_MIN_ANGLE_DIFFERENCE))
        {
            // Rejected if the 2 arguments are too close: this is to avoid rectilinear orbits (the solver would give a correct answer, but you probably don't want to deal with that...)
            LOG(LOG_LEVEL.ERROR, "Impossible to initialize Lambert solver : the chosen points have their argument too close");
            return (null, null);
        }

        // if direction is negative, deltaTheta must be between -2.0*PI and 0.0
        if (_direction < 0) { deltaTheta -= 2.0 * Math.PI; }


        // SOLVE THE LAMBERT PROBLEM
        // -------------------------
        Orbit orbit1 = null;
        Orbit orbit2 = null;

        LambertSolver theSolver = new LambertSolver(planet.mass, R1, R2, arg1, arg2, deltaTheta, _T1, _T2, _direction, _nbTurns);

        if(_nbTurns == 0)
        {
            // Resolution for single turn
            // --------------------------
            orbit1 = orbit2 = theSolver.SolveLambertProblem(planet); // we return twice the same orbit since there's only one solution
        }
        else
        {
            // Resolution for multi-turn
            // -------------------------
            if (theSolver.checkIfSolutionExists())
            {
                orbit1 = theSolver.SolveLambertProblem(planet, E_TRANSFER_TYPE.LOW_ARC);
                orbit2 = theSolver.SolveLambertProblem(planet, E_TRANSFER_TYPE.HIGH_ARC);
            }
            else
            {
                LOG(LOG_LEVEL.WARNING, "No solution found to exist for the Lambert problem\n" + theSolver.ToString());
            }
        }

        return (orbit1, orbit2);
    }


    // Constructor
    // -----------
    private LambertSolver(double _mu, double _R1, double _R2, double _arg1, double _arg2, double _deltaTheta, double _T1, double _T2, int _direction, uint _nbTurns)
    {
        // SET ALL PARAMETERS
        // ------------------
        mu = _mu;
        R1 = _R1;
        R2 = _R2;
        arg1 = _arg1;
        arg2 = _arg2;
        T1 = _T1;
        T2 = _T2;
        nbTurns = _nbTurns;

        if (_direction < 0) { direction = -1; } // Note that direction = 0 can't work - Direction is assumed positive in this case
        else { direction = 1; }

        // Calculate all the necessary parameters to solve the Lambert problem
        // -------------------------------------------------------------------
        delta = _deltaTheta / 2.0;
        C = Math.Sqrt(Math.Max(0.0, R1 * R1 + R2 * R2 - 2.0 * R1 * R2 * Math.Cos(_deltaTheta))); // Chord
        S = (R1 + R2 + C) / 2.0;
        lambda = Math.Sqrt(R1 * R2) * Math.Cos(delta) / S;
        Tau = Math.Sqrt(S / (2.0 * mu)) * S;
        deltaT = (T2 - T1) / Tau;

        // Correct lambda if necessary (because of numerical precision problems...)
        if (lambda < -1.0) lambda = -1.0;
        else if (lambda > 1.0) lambda = 1.0;

        LOG(LOG_LEVEL.INFO, "Lambert problem resolution for: R1 = " + R1 + "; R2 = " + R2 + "; delta = " + delta + "; T2-T1 = " + (T2 - T1) + "; nbTurns = " + nbTurns);
    }


    // -------------------------------------
    //          INTERNAL METHODS
    // -------------------------------------
    // Modify what follows at your own risk!

    // The function that runs the main procedure to solve the Lambert problem
    private Orbit SolveLambertProblem(Planet planet, E_TRANSFER_TYPE transfer_type = E_TRANSFER_TYPE.LOW_ARC)
    {
        // Get a first estimation
        // ----------------------
        double x = GetInitialGuess(transfer_type);
        double x_next = 0.0;

        uint iteration_number = 0;
        double relative_precision = 0.0;


        // Perform Newton/Householder iterations
        // -------------------------------------
        do
        {
            iteration_number++;

            // calculate next x value from current approximation
            if ((nbTurns == 0) && (Math.Abs(x - 1.0) < C_BATTIN_THRESHOLD))
            {
                x_next = Apply_BattinMethod_NewtonIteration(x);
            }
            else
            {
                x_next = Apply_HouseholderIteration(x);
            }

            // Evaluate the relative difference between old and new value
            if (Math.Abs(x) > 1.0)
            {
                relative_precision = Math.Abs((x_next - x) / x);
            }
            else
            {
                relative_precision = Math.Abs(x_next - x);
            }

            // memorize new value (for next iteration or for the end of the procedure)
            x = x_next;

        } while ((relative_precision > C_PRECISION_THRESHOLD) && (iteration_number < C_NB_MAX_ITERATIONS));

        LOG(LOG_LEVEL.DEBUG, "Found solution x = " + x + "; relative_precision = " + relative_precision + "; nb iterations = " + iteration_number);

        // Check success of procedure 
        if(relative_precision > C_PRECISION_THRESHOLD)
        {
            LOG(LOG_LEVEL.ERROR, "Failed to find a solution by excessive number of Householder iterations to the Lambert problem\n" + ToString());
            LOG(LOG_LEVEL.ERROR, "Found solution x = " + x + "; relative_precision = " + relative_precision + "; nb iterations = " + iteration_number);
            return null;
        }


        // Extract the orbitals parameters
        // -------------------------------
        ComputeOrbitParameters(x, out double slr, out double ecc, out double argOfPeriapsis, out double periapsisPassageTime);


        // Create the orbit with all parameters set and return it
        // ------------------------------------------------------
        Orbit orbit = Orbit_Utils.CreateOrbit(slr, ecc, argOfPeriapsis, direction, planet, PathType.Encounter, null);

        if (orbit != null)
        {
            orbit.periapsisPassageTime = periapsisPassageTime;
            orbit.orbitStartTime = T1;
            orbit.orbitEndTime = T2;
        }

        return orbit;
    }

    // Evaluates the existence of solutions, especially in multi-turn configurations (a solution is guaranteed in single turn)
    private bool checkIfSolutionExists()
    {
        bool exists = false;

        if (nbTurns == 0) 
        { 
            // A unique solution is guaranteed in the single turn case
            exists = true;
        }
        else
        {
            // Multi-turn case : either no solution exists, or 2 solutions exist
            // ---------------
            // Above this value, a couple of solution is guaranteed
            TM0 = nbTurns * Math.PI + Math.Acos(lambda) + lambda * Math.Sqrt(1 - lambda*lambda);

            if(deltaT >= TM0)
            {
                exists = true;
            }
            else
            {
                // Threshold determined from both theoretical study and experimentation - below this value, there's no solution
                double T_inf = TM0 - 43.0 / (15.0 * (nbTurns + 1) * Math.PI); 

                if(deltaT < T_inf)
                {
                    exists = false;
                }
                else
                {
                    // deltaT is between T_inf and TM0 - it's necessary to calculate the minimum time of flight accurately
                    (x_min, TMin) = calculateMinimumTimeOfFlight();
                    lowToF = true; // allows the use of custom approximations adapted for low time of flights for initial guesses

                    if(deltaT < TMin)
                    {
                        exists = false;
                    }
                    else
                    {
                        exists = true;
                    }
                }
            }
        }

        return exists;
    }

    // Calculates the minimum time of flight and the x for which it happens
    // WARNING: Made for multi-turn configurations; result is unpredictable (or predictably funky!) for single turn configurations.
    private (double x_min, double TMin) calculateMinimumTimeOfFlight()
    {
        double x = 0.0;
        double x_next = 0.0;
        double minToF = 0.0;

        uint iteration_number = 0;
        double relative_precision = 0.0;


        // Perform Halley iterations
        // -------------------------
        do
        {
            iteration_number++;

            CalculateValueAndDerivatives(x, out double F0, out double F1, out double F2, out double F3);

            // calculate next x value from current approximation - that's Halley formula applied to the derivative, to find the x that cancels the derivative
            x_next = x - 2.0 * F1 * F2 / (2.0 * F2 * F2 - F1 * F3);

            // Evaluate the relative difference between old and new value
            relative_precision = Math.Abs(x_next - x);

            // memorize new value (for next iteration or for the end of the procedure)
            x = x_next;
            minToF = F0;

        } while ((relative_precision > C_PRECISION_THRESHOLD) && (iteration_number < C_NB_MAX_ITERATIONS));

        // Check success of procedure 
        if (relative_precision > C_PRECISION_THRESHOLD)
        {
            LOG(LOG_LEVEL.ERROR, "Failed to evaluate minimum time of flight by excessive number of Halley iterations for the Lambert problem\n" + ToString());
            minToF = double.PositiveInfinity; // to force the resolution algorithm to fail safely
        }

        LOG(LOG_LEVEL.DEBUG, "Minimum time of flight : found solution x = " + x + "; relative_precision = " + relative_precision + "; nb iterations = " + iteration_number);

        return (x, minToF);
    }


    // Initial guess based on equation (30) from Dario Izzo's paper
    // WARNING: the equation given for range [T1, T0] is wrong in the paper, it's been corrected
    private double GetInitialGuess(E_TRANSFER_TYPE transferType = E_TRANSFER_TYPE.LOW_ARC)
    {
        double x = 0;

        if(nbTurns == 0)
        {
            // Single turn case
            // ----------------
            double T_0 = Math.Acos(lambda) + lambda * (1.0 - lambda * lambda);
            double T_1 = 2.0 / 3.0 * (1.0 - lambda * lambda * lambda);

            if (deltaT >= T_0)
            {
                x = Math.Pow(T_0 / deltaT, 2.0 / 3.0) - 1.0;
            }
            else if (deltaT < T_1)
            {
                x = 2.5 * T_1 * (T_1 - deltaT) / (deltaT * (1 - Math.Pow(lambda, 5.0))) + 1;
            }
            else
            {
                double exponent = Math.Log(2.0) / Math.Log(T_1 / T_0);
                x = Math.Pow(deltaT / T_0, exponent) - 1.0;
            }
        }
        else
        {
            // Multi-turn case
            // ---------------
            if (lowToF)
            {
                // Use the experimental low time of flight approximation
                x = Math.Sqrt((deltaT - TMin) / (TM0 - TMin));

                if (transferType == E_TRANSFER_TYPE.HIGH_ARC)
                {
                    x = 1.0 - x;
                }
                else
                {
                    x = 1.0 + x;
                }

                x *= x_min;
            }
            else
            {
                // Use the usual approximation as indicated in Dario Izzo's paper
                if (transferType == E_TRANSFER_TYPE.HIGH_ARC)
                {
                    x = 1 - 2 / (1 + Math.Pow((nbTurns + 1) * Math.PI / (8.0 * deltaT), 2.0 / 3.0));
                }
                else
                {
                    x = 1 - 2 / (1 + Math.Pow((8.0 * deltaT) / (nbTurns * Math.PI), 2.0 / 3.0));
                }
            }
        }

        LOG(LOG_LEVEL.DEBUG, "Initial guess : x = " + x);

        return x;
    }

    // Calculates the value and the first three derivatives of Lancaster equation for a given x
    private void CalculateValueAndDerivatives(double x, out double F0, out double F1, out double F2, out double F3)
    {
        double _1_minus_x2 = 1 - x * x;
        double lambda2 = lambda * lambda;
        double y = Math.Sqrt(Math.Max(0.0, 1.0 - lambda2 * _1_minus_x2));
        double _1_minus_lambda2 = 1.0 - lambda2;
        double lambda3 = lambda2 * lambda;
        double cos_psi = x * y + lambda * _1_minus_x2;
        double psi;

        if (cos_psi <= 1.0)
        {
            psi = Math.Acos(cos_psi);
        }
        else
        {
            psi = Math.Log(cos_psi + Math.Sqrt(cos_psi * cos_psi - 1.0)); // Inverse hyperbolic cosine
        }

        // Calculation of the function value and its derivatives according to equations (18) and (22) from Dario Izzo's paper
        F0 = ((psi + nbTurns * Math.PI) / Math.Sqrt(Math.Abs(_1_minus_x2)) - x + lambda * y) / _1_minus_x2;
        F1 = (3.0 * x * F0 - 2.0 + 2.0 * lambda3 * x / y) / _1_minus_x2;
        F2 = (3.0 * F0 + 5.0 * x * F1 + 2.0 * _1_minus_lambda2 * Math.Pow(lambda / y, 3.0)) / _1_minus_x2;
        F3 = (7.0 * x * F2 + 8.0 * F1 - 6.0 * _1_minus_lambda2 * x * Math.Pow(lambda / y, 5.0)) / _1_minus_x2;
    }

    // Applies an iteration of the Householder method according to the generic method defined in Dario Izzo's paper
    private double Apply_HouseholderIteration(double x)
    {
        CalculateValueAndDerivatives(x, out double F0, out double F1, out double F2, out double F3);

        F0 = F0 - deltaT;

        double denominator = F1 * (F1*F1 - F0*F2) + F3*F0*F0 / 6.0;
        double new_x = x - F0 * (F1*F1 - F0*F2/2.0) / denominator;

        return new_x;
    }


    // Calculates the function and its derivatives using the hypergeometric function as described in Battin's paper
    private double Apply_BattinMethod_NewtonIteration(double x)
    {
        double new_x = 0.0;

        double _1_minus_x2 = 1.0 - x * x;
        double y = Math.Sqrt(1.0 - lambda * lambda * _1_minus_x2);
        double eta = y - lambda * x;
        double eta3 = eta * eta * eta;
        double S1 = (1.0 - lambda - eta * x) / 2.0;

        double Q = 0.0;
        double Q_derivative = 0.0;
        GetQ_and_Qderivative(S1, y, eta, ref Q, ref Q_derivative);

        double F0 = (Q * eta3 + 4.0 * lambda * eta) / 2.0;
        double F1 = (Q_derivative * eta3 - lambda * eta * (3.0 * Q * eta * eta + 4.0 * lambda) / y) / 2.0;

        F0 = F0 - deltaT;
        new_x = x - F0 / F1;

        return new_x;
    }


    // Calculate the result of the hypergeometric function as a continued fraction
    // with the algorithm given in 7.47 in Battin's paper 
    private double GetHypergeometricSerieValue(double z, uint firstGammaIndice)
    {
        double delta = 1.0;
        double u = 1.0;
        double Sigma = 1.0;
        uint gammaIndice = firstGammaIndice;

        double rel_diff = 1.0;
        const double precision = 0.000000000001;
	    
        // usually iterates 5 times
	    while((rel_diff > precision) && (gammaIndice < 9))
	    {
		    delta = 1.0 / (1.0 - T_Gamma[gammaIndice] * delta* z);
            u = u* (delta - 1.0);
		    Sigma = Sigma + u;
		    rel_diff = Math.Abs(u / Sigma);
            gammaIndice++;
	    }

        return Sigma;
    }


    // Calculates dQ/dx combining equations 7.49 and 7.50 from Battin's paper
    // Calculates Q from gamma1 and G as described in the same paragraph of Battin's paper
    private void GetQ_and_Qderivative(double S1, double y, double eta, ref double Q, ref double Q_der)
    {
	    double gamma1 = T_Gamma[1] ;
        double G = GetHypergeometricSerieValue(S1, 2);

        double dQ_over_dS1 = 2.0 * (2.0 - gamma1 * G) / (1.0 - S1) / (1.0 - gamma1 * G * S1);

        Q_der = (gamma1 * G - 2.0) * eta * eta / y / (1 - S1) / (1 - gamma1 * S1 * G);
        Q = 4.0 / 3.0 / (1.0 - gamma1 * G * S1);
    }


    // Extracts all orbital parameters from the x value that satisfies the equation from the Lambert problem
    private void ComputeOrbitParameters(double x, out double slr, out double ecc, out double argOfPeriapsis, out double periapsisPassageTime)
    {
        double _1_minus_x2 = 1 - x * x;
        double y = Math.Sqrt(Math.Max(0.0, 1 - lambda * lambda * _1_minus_x2));

        double rho; // rho = (R2 - R1)/C ; -1 <= rho <=1

        if (C > 0.0)
        {
            rho = (R2 - R1) / C;
        }
        else
        {
            // rho is not defined because P1 = P2 and C = 0
            // Several solutions exist, so we arbitrarily choose the solution that minimizes eccentricity
            // (Note that it can't happen as long as the angle between the 2 points has been limited to a minimum)
            rho = 0.0;
        }

        // A few constants that will be useful for future operations
        double _1_minus_rho2 = 1 - rho*rho;
        double _y_plus_lx2 = y + lambda * x; _y_plus_lx2 *= _y_plus_lx2;
        double q = _1_minus_rho2 * _y_plus_lx2;

        // Kept for memory, but useless for future calculations (Warning: sma is potentially infinite, in case of parabolic trajectory)
        // specificEnergy = -mu * _1_minus_x2 / S;
        // sma = S / (2.0 * _1_minus_x2); 

        // semi-latus-rectum
        // -----------------
        slr = q * S / 2.0;

        // eccentricity
        // ------------
        ecc = 1.0 - q * _1_minus_x2; // That's actually ecc^2
        if (ecc < 0.0) ecc = 0.0; // Make sure the square root works, in case of numerical imprecisions...

        ecc = Math.Sqrt(ecc);


        if (1.0 + ecc > 1.0) // eccentricity != 0.0 - tested that way (instead of ecc > 0.0) so that ecc is considered as zero if it's negligible (< 10^(-16))
        {
            // Argument of periapsis
            // ---------------------
            argOfPeriapsis = arg1 + delta - Math.Atan2(rho * Math.Sqrt(_1_minus_rho2) * _y_plus_lx2, q * (x * y + lambda * _1_minus_x2) - 2 * lambda) * direction;

            while (argOfPeriapsis < -Math.PI) argOfPeriapsis += 2.0 * Math.PI;
            while (argOfPeriapsis >  Math.PI) argOfPeriapsis -= 2.0 * Math.PI;

            // Make sure we work with the closest point to periapsis for better precision
            double S_term = rho * (x + lambda * y);
            double R_ref;
            double T_ref;
            double T_diff;
            if (R1 < R2) { S_term -= x - lambda * y; T_ref = T1; R_ref = R1; }
            else         { S_term += x - lambda * y; T_ref = T2; R_ref = R2; }

            // Periapsis passage time
            // ----------------------
            periapsisPassageTime = 0.0;

            if (Math.Abs(1.0 - x) > 0.001)
            {
                // Calculate the eccentric term of the Kepler equation: e*sin(E) / e*sinh(H)
                double eccentricTerm = S_term * Math.Sqrt(Math.Abs(_1_minus_x2));

                // Calculate eccentric anomaly
                double E;
                if (x < 1.0)
                {
                    double theSinus = eccentricTerm / ecc;
                    if (theSinus < -1.0) theSinus = -1.0;
                    else if(theSinus > 1.0) theSinus = 1.0;
                    E = Math.Asin(theSinus);

                    // The quantity tested gives the sign of cos(E1), which is necessary to extract the eccentric anomaly
                    if (1.0 - 2.0 * R_ref * _1_minus_x2 / S < 0.0)
                    {
                        E = Math.PI - E;
                        if (E > Math.PI) E -= 2.0 * Math.PI; // The angle must be between -PI and +PI
                    }
                }
                else
                {
                    E = eccentricTerm / ecc;
                    E = Math.Log(E + Math.Sqrt(E * E + 1.0)); // Inverse hyperbolic sine
                }

                // Calculate mean anomaly
                double M = E - eccentricTerm;
                if (x > 1.0) M = -M;

                // Finally, calculate the periapsis passage time from the Kepler equation
                T_diff = 0.5 * Tau * M / Math.Pow(Math.Abs(_1_minus_x2), 1.5);
            }
            else
            {
                // Quasi-parabolic orbits need a special formula!
                double u = q * _1_minus_x2;
                double v = Math.Pow(S_term / ecc,2.0) * _1_minus_x2;

                T_diff = q * S_term * ((((63*u + 70)*u + 80)*u + 96)*u + 128)/256.0 + Math.Pow(S_term/ecc, 3.0) * ((((19845*v + 26950)*v + 39600)*v + 66528)*v + 147840)/887040.0;
                T_diff *= 0.5 * Tau;
            }

            periapsisPassageTime = T_ref - T_diff;
        }
        else // circular orbit
        {
            // arg of periapsis
            // ----------------
            argOfPeriapsis = 0.0; // default value

            // Periapsis passage time
            // ----------------------
            double period = Math.PI * Tau / Math.Pow(_1_minus_x2, 1.5);
            double trueAnomaly1 = arg1;

            // We work with true anomaly between 0 and 2Pi
            while (trueAnomaly1 > 2.0 * Math.PI) trueAnomaly1 -= 2.0 * Math.PI;
            while (trueAnomaly1 < 0.0) trueAnomaly1 += 2.0 * Math.PI;

            periapsisPassageTime = T1 - direction * trueAnomaly1 * period;
        }
    }

    public override string ToString()
    {
        string str_LambertParameters;

        str_LambertParameters = " Parameters for Lambert solver : \n" +
                                "  R1 = " + R1 + "; arg1 = " + arg1 + "; T1 = " + T1 + "\n" +
                                "  R2 = " + R2 + "; arg2 = " + arg2 + "; T2 = " + T2 + "\n" +
                                "  dir = " + direction + "; nbTurns = " + nbTurns + "\n" +
                                " Internal variables :\n" +
                                "  delta = " + delta + "; C = " + C + "; S = " + S + "\n" +
                                "  lambda = " + lambda + "; Tau = " + Tau + "; deltaT = " + deltaT + "\n" +
                                "  TM0 = " + TM0;

        if (lowToF)
        {
            str_LambertParameters += "\n Low time of flight parameters :\n" +
                                     "  x_min = " + x_min + "; TMin = " + TMin +
                                     "\n -------------------------------------";
        }

        return str_LambertParameters;
    }

    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.LAMBERT_SOLVER, level, message);
    }

}
