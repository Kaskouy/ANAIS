using SFS.World;
using SFS.WorldBase;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;



public static class LambertSolver
{
    // Input variables
    // ---------------
    private static double mu;     // gravitational parameter

    private static double R1;
    private static double arg1;
    private static double R2;
    private static double arg2;

    private static double T1;     // date at P1
    private static double T2;     // date at P2
    private static int direction; // direction of the orbit (1: anti-clockwise; -1: clockwise)

    // Internal variables
    // ------------------
    private static double delta;  // semi-angle (= (arg2 - arg1) / 2.0)
    private static double C;      // chord; distance from P1 to P2
    private static double S;      // semi-perimeter
    private static double lambda; // lambda parameter (lambda = sqrt(R1*R2)*cos(delta))
    private static double deltaT; // An adimensional parameter, proportional (but not equal!) to (T2 - T1)

    // Constants and internal shit
    // ---------
    private const double C_MIN_ANGLE_DIFFERENCE = 0.0003; // 1 minute of arc (the value is in radians) - the calculation is not done if the arguments are close by that value
    private const double C_BATTIN_THRESHOLD = 0.001; // A threshold to decide wether or not we switch to Battin's method
    private const uint C_NB_MAX_ITERATIONS = 10; // The maximum number of Newton/Halley iterations we allow (3 are enough in most cases - a few very specific cases can require more)
    private const double C_PRECISION_THRESHOLD = 0.0000000001; // The precision required to end Newton/Halley iterations

    /* Gives the value of gamma(n) as defined in equation 7.46 from Battin's paper - only values from T_Gamma[1] to T_Gamma[5] are useful in practice
    // The gamma(n) values are obtained through the following formula:
    // - gamma =    n * (n-3)/ (2*n+1) / (2*n+3) for n even
    // - gamma = (n+2)* (n+5)/ (2*n+1) / (2*n+3) for n odd*/
    private static readonly double[] T_Gamma = {0.0, 6.0/5.0, -2.0/35.0, 40.0/63.0, 4.0/99.0, 70.0/143.0, 6.0/65.0, 36.0/85.0, 40.0/323.0};


    // The main function: CalculateTrajectory
    // --------------------------------------
    //
    // DESCRIPTION:
    // That function allows to calculate a trajectory solving what's usually called the Lambert problem: you start from point P1 at date T1
    // and you want to be at position T2 at date T2. The Lambert solver calculates the trajectory that fulfills those constraints.
    // You can also specify the direction: clockwise (-1) or anti-clockwise (+1). There's a solution in both cases. A negative direction
    // does NOT mean going from P2 to P1, this is really about the direction of the orbit.
    // --------------------------------------
    //
    // RESTRICTIONS:
    // The Lambert problem usually always have a solution. There are a few restrictions though:
    // - T2 must be greater than T1: you can't come back to the past!
    // - The arguments of P1 and P2 mustn't be too close. The solver would give a correct result, but there would be a risk that you end up
    //   with a rectilinear orbit, which you probably don't want to deal with.
    //
    // The function will return null in those 2 situations.
    // --------------------------------------
    //
    // RETURN VALUE:
    // The orbit returned (if not null) has all its orbital parameters set, including its periapsis passage time.
    // --------------------------------------
    public static Orbit CalculateTrajectory(Planet planet, Double2 _P1, Double2 _P2, double _T1, double _T2, int _direction)
    {
        // CONTROL PART
        // ------------
        if (!(_T2 > _T1))
        {
            // The arrival date must ABSOLUTELY be after the departure date!
            return null;
        }

        arg1 = _P1.AngleRadians;
        arg2 = _P2.AngleRadians;

        // The angular difference between the 2 positions
        double deltaTheta = arg2 - arg1;

        // deltaTheta must be between 0.0 and 2.0*PI (assuming positive direction)
        while (deltaTheta > 2.0 * Math.PI) { deltaTheta -= 2.0 * Math.PI; }
        while (deltaTheta < 0.0) { deltaTheta += 2.0 * Math.PI; }

        if ((deltaTheta < C_MIN_ANGLE_DIFFERENCE) || (2.0 * Math.PI - deltaTheta < C_MIN_ANGLE_DIFFERENCE))
        {
            // Rejected if the 2 arguments are too close: this is to avoid rectilinear orbits (the solver would give a correct answer, but you probably don't want to deal with that...)
            return null;
        }

        // if direction is negative, deltaTheta must be between -2.0*PI and 0.0
        if (_direction < 0) { deltaTheta -= 2.0 * Math.PI; }


        // SET ALL PARAMETERS
        // ------------------
        mu = planet.mass;
        R1 = _P1.magnitude;
        R2 = _P2.magnitude;
        T1 = _T1;
        T2 = _T2;

        if (_direction < 0) { direction = -1; } // Note that direction = 0 can't work - Direction is assumed positive in this case
        else { direction = 1; }

        delta = deltaTheta / 2.0;
        C = Math.Sqrt(Math.Max(0.0, R1 * R1 + R2 * R2 - 2.0 * R1 * R2 * Math.Cos(deltaTheta))); // Chord
        S = (R1 + R2 + C) / 2.0;
        lambda = Math.Sqrt(R1 * R2) * Math.Cos(delta) / S;
        deltaT = Math.Sqrt(2.0 * mu / (S * S * S)) * (T2 - T1);

        // Correct lambda if necessary (because of numerical precision problems...)
        if (lambda < -1.0) lambda = -1.0;
        else if(lambda > 1.0) lambda = 1.0;

        // SOLVE THE LAMBERT PROBLEM
        // -------------------------
        double specificEnergy, slr, ecc, argOfPeriapsis, periapsisPassageTime;
        SolveLambertProblem(out specificEnergy, out slr, out ecc, out argOfPeriapsis, out periapsisPassageTime);

        Orbit orbit = Orbit_Utils.CreateOrbit(slr, ecc, argOfPeriapsis, direction, planet, PathType.Encounter, null);

        if (orbit != null)
        {
            orbit.periapsisPassageTime = periapsisPassageTime;
            orbit.orbitStartTime = T1; // A bit hacky...
            orbit.orbitEndTime = T2;
        }

        return orbit;
    }


    // -------------------------------------
    //          INTERNAL METHODS
    // -------------------------------------
    // Modify what follows at your own risk!

    // The function that runs the main procedure to solve the Lambert problem
    private static void SolveLambertProblem(out double specificEnergy, out double slr, out double ecc, out double argOfPeriapsis, out double periapsisPassageTime)
    {
        // Get a first estimation
        // ----------------------
        double x = GetInitialGuess();
        double x_next = 0.0;

        uint iteration_number = 0;
        double relative_precision = 0.0;


        // Perform Newton/Halley iterations
        // --------------------------------
        do
        {
            iteration_number++;

            // calculate next x value from current approximation
            if (Math.Abs(x - 1.0) < C_BATTIN_THRESHOLD)
            {
                x_next = Apply_BattinMethod_NewtonIteration(x);
            }
            else
            {
                x_next = Apply_GenericMethod_HalleyIteration(x);
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

        
		// Extract the orbitals parameters
        // -------------------------------
        bool periapsisPassageTimeValid = ComputeOrbitParameters(x, out specificEnergy, out slr, out ecc, out argOfPeriapsis, out periapsisPassageTime);

        // In case the orbit is nearly parabolic we have to calculate the periapsis passage time by ourselves
        if(!periapsisPassageTimeValid)
        {
            double timeFromPeri = 0.0;
            new KeplerSolver(mu, slr / (1.0 + ecc), specificEnergy).GetTimeAtPosition(R1, direction*(arg1 - argOfPeriapsis), ref timeFromPeri);
            periapsisPassageTime = T1 - timeFromPeri;
        }
    }


    // Initial guess based on equation (30) from Dario Izzo's paper
    // WARNING: the equation given for range [T1, T0] is wrong in the paper, it's been corrected
    private static double GetInitialGuess()
    {
        double x = 0;

        double T_0 = Math.Acos(lambda) + lambda * (1.0 - lambda * lambda);
        double T_1 = 2.0 / 3.0 * (1.0 - lambda * lambda * lambda);

        if(deltaT >= T_0)
        {
            x = Math.Pow(T_0 / deltaT, 2.0 / 3.0) - 1.0;
        }
        else if(deltaT < T_1)
        {
            x = 2.5 * T_1 * (T_1 - deltaT) / (deltaT * (1 - Math.Pow(lambda, 5.0))) + 1;
        }
        else
        {
            double exponent = Math.Log(2.0) / Math.Log(T_1 / T_0);
            x = Math.Pow(deltaT / T_0, exponent) - 1.0;
        }

        return x;
    }


    // Calculates the function and its derivatives according to the generic method defined in Dario Izzo's paper
    private static double Apply_GenericMethod_HalleyIteration(double x)
    {
        double new_x = 0.0;

        double _1_minus_x2 = 1 - x*x;
        double y = Math.Sqrt(Math.Max(0.0, 1 - lambda * lambda * _1_minus_x2));
        double cos_psi = x * y + lambda * _1_minus_x2;
        double psi;

        if (cos_psi <= 1.0)
        {
            psi = Math.Acos(cos_psi);
        }
        else
        {
            psi = Math.Log(cos_psi + Math.Sqrt(cos_psi* cos_psi - 1.0)); // Inverse hyperbolic cosine
        }

        // Calculation of the function value and its derivatives according to equations (18) and (22) from Dario Izzo's paper
        double F0 = (psi / Math.Sqrt(Math.Abs(_1_minus_x2)) - x + lambda * y) / _1_minus_x2;
        double F1 = (3.0 * x * F0 - 2.0 + 2.0 * lambda * lambda * lambda * x / y) / _1_minus_x2;
        double F2 = (3.0 * F0 + 5.0 * x * F1 + 2.0 * (1 - lambda * lambda) * Math.Pow(lambda/y, 3.0)) / _1_minus_x2;

        F0 = F0 - deltaT;
        new_x = x - 2.0 * F1 * F0 / (2.0 * F1 * F1 - F0 * F2);

        return new_x;
    }


    // Calculates the function and its derivatives using the hypergeometric function as described in Battin's paper
    private static double Apply_BattinMethod_NewtonIteration(double x)
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
    private static double GetHypergeometricSerieValue(double z, uint firstGammaIndice)
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
    private static void GetQ_and_Qderivative(double S1, double y, double eta, ref double Q, ref double Q_der)
    {
	    double gamma1 = T_Gamma[1] ;
        double G = GetHypergeometricSerieValue(S1, 2);

        double dQ_over_dS1 = 2.0 * (2.0 - gamma1 * G) / (1.0 - S1) / (1.0 - gamma1 * G * S1);

        Q_der = (gamma1 * G - 2.0) * eta * eta / y / (1 - S1) / (1 - gamma1 * S1 * G);
        Q = 4.0 / 3.0 / (1.0 - gamma1 * G * S1);
    }


    // Extracts all orbital parameters from the x value that satisfies the equation from the Lambert problem
    private static bool ComputeOrbitParameters(double x, out double specificEnergy, out double slr, out double ecc, out double argOfPeriapsis, out double periapsisPassageTime)
    {
        bool periapsisPassageTimeValid = false;
        
        double _1_minus_x2 = 1 - x * x;
        double y = Math.Sqrt(Math.Max(0.0, 1 - lambda * lambda * _1_minus_x2));
        double cos_phi = x * y - lambda * _1_minus_x2;
        double cos_psi = x * y + lambda * _1_minus_x2;

        double rho; // rho = (R2 - R1)/C ; -1 <= rho <=1

        if (C > 0.0)
        {
            rho = (R2 - R1) / C * direction;
        }
        else
        {
            // rho is not defined because P1 = P2 and C = 0
            // Several solutions exist, so we arbitrarily choose the solution that minimizes eccentricity
            // (Note that it can't happen as long as the angle between the 2 points has been limited to a minimum)
            rho = 0.0;
        }

        double rho2 = rho * rho;
        double D = S * Math.Pow(y + lambda * x, 2.0) / 2.0; // A characteristic distance used in p and arg of periapsis calculations


        // Specific energy
        // ---------------
        specificEnergy = -mu * _1_minus_x2 / S;

        // semi-latus-rectum
        // -----------------
        slr = D * (1 - rho2);

        // eccentricity
        // ------------
        ecc = cos_phi * cos_phi * (1 - rho2) + rho2; // That's actually ecc^2
        if (ecc < 0.0) ecc = 0.0; // Make sure the square root works, in case of numerical imprecisions...

        ecc = Math.Sqrt(ecc);


        if (1.0 + ecc > 1.0) // eccentricity != 0.0 - tested that way (instead of ecc > 0.0) so that ecc is considered as zero if it's negligible (< 10^(-16))
        {
            // Argument of periapsis
            // ---------------------
            argOfPeriapsis = arg1 + delta - Math.Atan2(D * rho * Math.Sqrt(Math.Max(0.0, 1 - rho2)), slr * cos_psi - lambda * S);

            while (argOfPeriapsis < -Math.PI) argOfPeriapsis += 2.0 * Math.PI;
            while (argOfPeriapsis >  Math.PI) argOfPeriapsis -= 2.0 * Math.PI;

            // Periapsis passage time
            // ----------------------
            periapsisPassageTime = 0.0;

            if (Math.Abs(1.0 - x) > 0.001) // because it doesn't work for FUCKING parabolic orbits - I HATE PARABOLAS!!!
            {
                // Calculate the eccentric term of the Kepler equation: e*sin(E) / e*sinh(H)
                double eccentricTerm = (rho * (lambda * y + x) - direction * (x - lambda * y)) * Math.Sqrt(Math.Abs(_1_minus_x2));

                // Calculate eccentric anomaly
                double E1;
                if (x < 1.0)
                {
                    double theSinus = eccentricTerm / ecc;
                    if (theSinus < -1.0) theSinus = -1.0;
                    else if(theSinus > 1.0) theSinus = 1.0;
                    E1 = Math.Asin(theSinus);

                    // The quantity tested gives the sign of cos(E1), which is necessary to extract the eccentric anomaly
                    if (1.0 - 2.0 * R1 * _1_minus_x2 / S < 0.0)
                    {
                        E1 = Math.PI - E1;
                        if (E1 > Math.PI) E1 -= 2.0 * Math.PI; // The angle must be between -PI and +PI
                    }
                }
                else
                {
                    E1 = eccentricTerm / ecc;
                    E1 = Math.Log(E1 + Math.Sqrt(E1 * E1 + 1.0)); // Inverse hyperbolic sine
                }

                // Calculate mean anomaly
                double M1 = E1 - eccentricTerm;
                if (x > 1.0) M1 = -M1;

                // Finally, calculate the periapsis passage time from the Kepler equation
                periapsisPassageTime = T1 - mu / Math.Pow(2.0 * Math.Abs(specificEnergy), 1.5) * M1 * direction;
                periapsisPassageTimeValid = true;
            }
        }
        else // circular orbit
        {
            // arg of periapsis
            // ----------------
            argOfPeriapsis = 0.0; // default value

            // Periapsis passage time
            // ----------------------
            double period = mu / Math.Pow(-2.0 * specificEnergy, 1.5);
            double trueAnomaly1 = arg1;

            // We work with true anomaly between 0 and 2Pi
            while (trueAnomaly1 > 2.0 * Math.PI) trueAnomaly1 -= 2.0 * Math.PI;
            while (trueAnomaly1 < 0.0) trueAnomaly1 += 2.0 * Math.PI;

            periapsisPassageTime = T1 - direction * trueAnomaly1 * period;
            periapsisPassageTimeValid = true;
        }

        return periapsisPassageTimeValid;
    }

}
