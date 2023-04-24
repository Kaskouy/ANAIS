using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

public static class KeplerSolver
{
    // The precision for a value calculated with the quasi-parabolic equation will be like C_QUASI_PARABOLIC_THRESHOLD to the power C_NB_TERMS_QUASI_PARABOLIC_EQUATION
    // Raising the threshold allows to use it in more situations, but it will force to raise the number of terms to keep an equal precision, which also means more calculations.
    private const double C_QUASI_PARABOLIC_THRESHOLD = 0.0001;
    private const uint C_NB_TERMS_QUASI_PARABOLIC_EQUATION = 4;

    // Newton's method ensures quadratic convergence, meaning that the number of correct digits roughly doubles on each iteration
    // Halley's method ensures cubic convergence, so the precision threshold can be higher (allows to save iterations)
    private const double C_PRECISION_NEWTON = 0.000000004; // approximately 1/2^28 - 52 bits (+1) are significative for a double
    private const double C_PRECISION_HALLEY = 0.000001; // approximately 1/2^20

    // The number max of iterations for the resolution algorithms - the max ever needed in practice is 3
    private const uint NB_MAX_ITERATIONS = 5;

    // INTERNAL VARIABLES
    // ------------------
    private static double mu;    // gravitational parameter
    private static double peri;  // periapsis
    private static double K;     // specific energy (= v^2/2 - mu/r)
    private static double ecc;   // eccentricity
    private static double gamma; // coeff gamma = (g_ecc - 1) / (g_ecc + 1) - used in quasi-parabolic calculations
    private static double omega; // coeff omega = K / (mu + peri * K) - used in quasi-parabolic calculations


	// ------------------------
	// METHOD GetTimeAtPosition
	// ------------------------
	// This method calculates the date at which an object will reach a given position (specified through radius and true anomaly).
	// It works with any kind of orbits, including rectilinear ones. It automatically chooses the best suited algorithm for the situation.
	// ------------------------
	// Parameters:
	// - p_mu           : gravitational parameter of the orbited body
	// - p_peri         : periapsis
	// - p_K            : specific energy (= v^2/2 - mu/r)
	// - radius         : radius at current position
	// - trueAnomaly    : argument measured from periapsis at current position (arg - arg_of_periapsis)
	// - time_from_peri : the date calculated. Add the time at periapsis to get the absolute time.
	// ------------------------
	public static void GetTimeAtPosition(double p_mu, double p_peri, double p_K, double radius, double trueAnomaly, ref double time_from_peri)
	{
		// Calculate internal variables
		mu = p_mu;
		peri = p_peri;
		K = p_K;
		ecc = 1.0 + 2.0 * peri * K / mu;
		omega = K / (mu + peri * K);
		// gamma will be calculated in the quasi-parabolic method if needed

		// Make sure we work with an angle included in [-Pi, +Pi]
		double trueAnomaly_normalized = trueAnomaly;
		while (trueAnomaly_normalized >  Math.PI) { trueAnomaly_normalized -= 2.0 * Math.PI; }
		while (trueAnomaly_normalized < -Math.PI) { trueAnomaly_normalized += 2.0 * Math.PI; }

		if (isQuasiParabolicEquationSuitableForTimeCalculation(radius))
		{
			// CALCULATION METHOD: QUASI-PARABOLIC
			// -----------------------------------
			GetTimeAtPos_QuasiParabolic(radius, trueAnomaly_normalized, ref time_from_peri);
		}
		else if (K < 0.0)
		{
			// CALCULATION METHOD: ELLIPTIC
			// ----------------------------
			GetTimeAtPos_Elliptic(radius, trueAnomaly_normalized, ref time_from_peri);
		}
		else // K > 0.0 - Note that K can't be null because this case automatically falls into the quasi-parabolic case
		{
			// CALCULATION METHOD: HYPERBOLIC
			// ------------------------------
			GetTimeAtPos_Hyperbolic(radius, trueAnomaly_normalized, ref time_from_peri);
		}
	}

	// ------------------------
	// METHOD GetPositionAtTime
	// ------------------------
	// This method calculates the position (specified through radius and true anomaly) where will be a given object at a given time.
	// It works with any kind of orbits, including rectilinear ones. It automatically chooses the best suited algorithm for the situation.
	// If the orbit is elliptic, it also works if the time given is several orbital periods ahead in time.
	// ------------------------
	// Parameters:
	// - p_mu           : gravitational parameter of the orbited body
	// - p_peri         : periapsis
	// - p_K            : specific energy (= v^2/2 - mu/r)
	// - time_from_peri : the date at which position must be found. You may want to pass time - time_at_periapsis.
	// - radius         : radius at current position
	// - trueAnomaly    : argument measured from periapsis at current position (in range [-Pi, Pi]). Add arg_of_periapsis to get the argument.
	// ------------------------
	public static void GetPositionAtTime(double p_mu, double p_peri, double p_K, double time_from_peri, ref double radius, ref double trueAnomaly)
	{
		mu = p_mu;
		peri = p_peri;
		K = p_K;
		ecc = 1.0 + 2.0 * peri * K / mu;
		omega = K / (mu + peri * K);
		gamma = peri * omega;

		/*LogFile.Log("--- GetPositionAtTime  ---");
		LogFile.Log("  mu              = " + mu);
		LogFile.Log("  periapsis       = " + peri);
		LogFile.Log("  specific energy = " + K);
		LogFile.Log("  eccentricity    = " + ecc);
		LogFile.Log("  omega           = " + omega);
		LogFile.Log("  gamma           = " + gamma);
		LogFile.Log("  ----------------------");*/


		double date_currentTurn;

		// Correct the date if the orbit is periodic (Make sure we work in the "current" period)
		// -----------------------------------------
		if (K < 0.0) // If orbit is elliptic
		{
			double period = Math.PI * mu / ((-K) * Math.Sqrt(-2.0 * K));
			double date_in_periods = time_from_peri / period;
			double nbTurns = Math.Floor(date_in_periods + 0.5); // nbTurns is actually the number of times the object will pass (or has passed if negative) the apoapsis before specified date
			date_currentTurn = time_from_peri - nbTurns * period;
		}
        else
        {
			date_currentTurn = time_from_peri;
		}


		/*LogFile.Log("  time from periapsis = " + date_currentTurn);
		LogFile.Log("  ----------------------");*/

		if (isQuasiParabolicEquationSuitableForPositionCalculation(date_currentTurn))
		{
			/*LogFile.Log(" CALCULATION METHOD USED: QUASI-PARABOLIC");*/
			// CALCULATION METHOD: QUASI-PARABOLIC
			// -----------------------------------
			GetPositionAtTime_QuasiParabolic(date_currentTurn, ref radius, ref trueAnomaly);
		}
		else if (K < 0.0)
		{
			/*LogFile.Log(" CALCULATION METHOD USED: ELLIPTIC");*/
			// CALCULATION METHOD: ELLIPTIC
			// ----------------------------
			GetPositionAtTime_Elliptic(date_currentTurn, ref radius, ref trueAnomaly);
		}
		else
		{
			/*LogFile.Log(" CALCULATION METHOD USED: HYPERBOLIC");*/
			// CALCULATION METHOD: HYPERBOLIC
			// ------------------------------
			GetPositionAtTime_Hyperbolic(date_currentTurn, ref radius, ref trueAnomaly);
		}

		/*LogFile.Log("   Found radius       = " + radius);
		LogFile.Log("   Found true anomaly = " + trueAnomaly);
		LogFile.Log("  ----------------------\n");*/
	}


	// ----------------------------------------------------------------------
	// --------------------       INTERNAL METHODS       --------------------
	// ----------------------------------------------------------------------

	// METHODS FOR TIME CALCULATION FROM POSITION DATA (the easy part)
	// -----------------------------------------------
	private static bool isQuasiParabolicEquationSuitableForTimeCalculation(double radius)
	{
		if ((ecc > 0.99) && (ecc < 1.01))
		{
			if (Math.Abs(omega) * (radius - peri) < C_QUASI_PARABOLIC_THRESHOLD * (1.0 + omega * radius))
			{
				return true;
			}
		}

		return false;
	}

	private static void GetTimeAtPos_QuasiParabolic(double radius, double trueAnomaly, ref double date)
	{
		double gamma = peri * omega; // also equal to (ecc-1)/(ecc+1) but this is faster to calculate that way

		// Q is the equivalent of the eccentric anomaly in the quasi-parabolic equation
		double Q2;
		if (Math.Abs(trueAnomaly) > 0.0026) // threshold evaluated experimentally to get the best precision
		{
			Q2 = (radius - peri) / (1.0 + omega * radius);
		}
		else
		{
			// approximation for small angles
			Q2 = ecc * peri * trueAnomaly * trueAnomaly / (1.0 + ecc * Math.Cos(trueAnomaly)) / 2.0 / (1.0 + omega * radius);
		}
		double Q = Math.Sqrt(Q2);
		if (trueAnomaly < 0.0) { Q = -Q; }

		// Applying the Quasi-Parabolic equation : sqrt((mu+peri*K)/2) * T = peri*Q + Q^3 * SUM{(i + (i+1)*gamma)/(2i+1) * (omega*Q2)^(i-1)} (i=1..4)
		// Note: we made sure that omega*Q2 < C_QUASI_PARABOLIC_THRESHOLD by calling isQuasiParabolicEquationSuitableForTimeCalculation before, so factor is very small
		double factor = omega * Q2;
		double coeff = 1.0;

		date = 0.0;

		for (uint i = 1; i <= C_NB_TERMS_QUASI_PARABOLIC_EQUATION; i++)
		{
			date += (i + (i + 1) * gamma) / (2 * i + 1) * coeff;  // (i + (i+1)*gamma)/(2i+1) * (omega*Q2)^(i-1)
			coeff *= factor;
		}

		date *= Q2;
		date += peri;
		date *= Q;

		date *= Math.Sqrt(2.0 / (mu + peri * K));
	}

	private static void GetTimeAtPos_Elliptic(double radius, double trueAnomaly, ref double date)
	{
		double E, M; // eccentric anomaly, mean anomaly
		double cos_E;

		if (ecc < 0.8)
		{
			// General formula : works as long as eccentricity != 1.0 (rectilinear orbits have their eccentricity exactly equal to 1 but still can be elliptic!)
			cos_E = (ecc + Math.Cos(trueAnomaly)) / (1.0 + ecc * Math.Cos(trueAnomaly));
		}
		else
		{
			// Alternative formula that works for rectilinear orbits (ecc = 1) (That one doesn't work for e=0 though)
			cos_E = (1.0 + 2.0 * K / mu * radius) / ecc;
		}

		// Make sure cosinus is valid
		if (cos_E < -1.0) { cos_E = -1.0; }
		if (cos_E >  1.0) { cos_E =  1.0; }

		// calculate eccentric anomaly
		E = Math.Acos(cos_E);
		if (trueAnomaly < 0.0) { E = -E; }

		//LogFile.Log(" GetTimeAtPos_Elliptic: trueAnomaly = " + trueAnomaly + " - radius = " + radius + " - ecc = " + ecc + " - cos_E = " + cos_E + " E = " + E);

		// Fonctionne!
		//E = 2.0 * Math.Atan(Math.Sqrt((1.0 - ecc) / (1.0 + ecc)) * Math.Tan(trueAnomaly / 2.0));

		// Fonctionne aussi
		//E = Math.Atan2(Math.Sqrt(1.0 - ecc * ecc) * Math.Sin(trueAnomaly), ecc + Math.Cos(trueAnomaly));

		// then mean anomaly
		M = E - ecc * Math.Sin(E);

		// Finally get the date (that formula also works for rectilinear orbits)
		date = mu * M / (-2.0 * K * Math.Sqrt(-2.0 * K));
	}

	private static void GetTimeAtPos_Hyperbolic(double radius, double trueAnomaly, ref double date)
	{
		double H, M; // eccentric anomaly, mean anomaly
		double cosh_H; // hyperbolic cosine of H

		cosh_H = (1.0 + 2.0 * K / mu * radius) / ecc;

		// Make sure hyperbolic cosinus is valid
		if (cosh_H < 1.0) { cosh_H = 1.0; }

		// calculate eccentric anomaly
		H = Math.Log(cosh_H + Math.Sqrt(cosh_H * cosh_H - 1.0)); // calculate Acosh
		if (trueAnomaly < 0.0) { H = -H; }

		// then mean anomaly
		M = ecc * Math.Sinh(H) - H;

		// Finally get the date (that formula also works for rectilinear orbits)
		date = mu * M / (2.0 * K * Math.Sqrt(2.0 * K));
	}


	// METHODS FOR POSITION CALCULATION FROM TIME DATA (the hard part)
	// -----------------------------------------------
	private static bool isQuasiParabolicEquationSuitableForPositionCalculation(double time)
	{
		if ((ecc > 0.99) && (ecc < 1.01))
		{
			double sqrt_threshold = Math.Sqrt(C_QUASI_PARABOLIC_THRESHOLD);
			double C_THRESHOLD = Math.Abs(gamma) * sqrt_threshold + sqrt_threshold*sqrt_threshold*sqrt_threshold / 3.0;

			if (Math.Abs(omega) * Math.Sqrt(Math.Abs(K) / 2.0) * Math.Abs(time) <= C_THRESHOLD)
			{
				return true;
			}
		}

		return false;
	}

	private static void GetPositionAtTime_QuasiParabolic(double time, ref double radius, ref double trueAnomaly)
	{
		double YVal = Math.Sqrt((mu + peri * K) / 2.0) * time;

		double relative_precision;
		double Q_n, Q_n_plus_1;

		// We'll solve numerically with the Newton method the Quasi-Parabolic equation : sqrt((mu+peri*K)/2) * T = peri*Q + Q^3 * SUM{(i + (i+1)*gamma)/(2i+1) * (omega*Q2)^(i-1)} (i=1..4)
		// Note: we made sure that omega*Q2 < C_QUASI_PARABOLIC_THRESHOLD by calling isQuasiParabolicEquationSuitableForPositionCalculation before.

		// Calculate a first approximation
		Q_n_plus_1 = GetFirstApproximation_QuasiParabolic(YVal);

		if (Math.Abs(Q_n_plus_1) > 0.0)
		{
			uint iteration_number = 0;
			double F0, F1;

			// Perform Newton iterations (1 or 2 are enough in practice, 3 may be needed in worst cases)
			// -------------------------
			do
			{
				iteration_number++;

				// start from the previous estimation (the first approximation if it's the first iteration)
				Q_n = Q_n_plus_1;

				// Calculating the function's value (F0) and its derivative (F1) in Q_n
				F0 = F1 = 0.0;

				double coeff = 1.0;
				double factor = omega * Q_n * Q_n;

				for (uint i = 1; i <= C_NB_TERMS_QUASI_PARABOLIC_EQUATION; i++)
				{
					double factor_i = (i + (i + 1) * gamma) * coeff; // (i + (i+1)*gamma) * (omega*Q2)^(i-1)

					F0 += factor_i / (2 * i + 1);
					F1 += factor_i;
					coeff *= factor;
				}

				F0 *= Q_n * Q_n * Q_n;
				F0 += peri * Q_n;

				F1 *= Q_n * Q_n;
				F1 += peri;

				// Calculate the following approximation
				Q_n_plus_1 = Q_n - (F0 - YVal) / F1;
	
				// relative precision
				relative_precision = Math.Abs((Q_n_plus_1 - Q_n) / Q_n_plus_1);

			} while ((relative_precision > C_PRECISION_NEWTON) && (iteration_number < NB_MAX_ITERATIONS));
		}

		// The solution has been found, calculate the corresponding radius
		// ---------------------------------------------------------------
		double Q2 = Q_n_plus_1 * Q_n_plus_1;
		radius = (peri + Q2) / (1.0 - omega * Q2);

		// Radius found, calculate true anomaly
		// ------------------------------------
		if (radius > 0.0)
		{
			if (Q2 > 0.0001 * peri)
			{
				double cos_arg = (peri * (1.0 + ecc) / radius - 1.0) / ecc;

				if (cos_arg >  1.0) { cos_arg =  1.0; }
				if (cos_arg < -1.0) { cos_arg = -1.0; }

				trueAnomaly = Math.Acos(cos_arg);
			}
			else
			{
				// Approximation for small angles
				trueAnomaly = 2.0 * Math.Sqrt(Q2 / (peri + Q2));
			}

			if (time < 0) { trueAnomaly = -trueAnomaly; }
		}
		else
		{
			// radius is null, true anomaly is not defined
			trueAnomaly = 0.0;
		}
	}

	private static void GetPositionAtTime_Elliptic(double time, ref double radius, ref double trueAnomaly)
	{
		double meanAnomaly = -2.0 * K * Math.Sqrt(-2.0 * K) / mu * time;

		double relative_precision;
		double E_n, E_n_plus_1;

		// We'll solve the Kepler equation (M = E - e.sin(E)) with the Halley method

		// Calculate a first approximation
		E_n_plus_1 = GetFirstApproximation_Elliptic(meanAnomaly);

		if (Math.Abs(E_n_plus_1) > 0.0)
		{
			uint iteration_number = 0;
			double F0, F1, F2;

			// Perform Halley iterations (1 or 2 are enough in practice, 3 may be needed in worst cases)
			// -------------------------
			do
			{
				iteration_number++;

				// start from the previous estimation (the first approximation if it's the first iteration)
				E_n = E_n_plus_1;

				// Calculating the function's value (F0), its derivative (F1) and second derivative (F2) in E_n
				F2 = ecc * Math.Sin(E_n);
				F1 = 1.0 - ecc * Math.Cos(E_n);
				F0 = E_n - F2 - meanAnomaly;

				// Calculate the next approximation
				E_n_plus_1 = E_n - 2.0 * F1 * F0 / (2.0 * F1 * F1 - F0 * F2);
	
				// relative precision
				relative_precision = Math.Abs((E_n_plus_1 - E_n) / E_n_plus_1);

			} while ((relative_precision > C_PRECISION_HALLEY) && (iteration_number < NB_MAX_ITERATIONS));
		}

		// The solution has been found, calculate the corresponding radius
		// ---------------------------------------------------------------
		double cos_eccentricAnomaly = Math.Cos(E_n_plus_1);
		radius = mu / (2.0 * K) * (ecc * cos_eccentricAnomaly - 1.0);

		if (radius > 0.0)
		{
			double cos_arg = (cos_eccentricAnomaly - ecc) / (1.0 - ecc * cos_eccentricAnomaly);

			if (cos_arg >  1.0) { cos_arg =  1.0; }
			if (cos_arg < -1.0) { cos_arg = -1.0; }

			trueAnomaly = Math.Acos(cos_arg);
			if (time < 0) { trueAnomaly = -trueAnomaly; }
		}
		else
		{
			// radius is null, true anomaly is not defined
			trueAnomaly = 0.0;
		}
	}

	private static void GetPositionAtTime_Hyperbolic(double time, ref double radius, ref double trueAnomaly)
	{
		double meanAnomaly = 2.0 * K * Math.Sqrt(2.0 * K) / mu * time;

		double relative_precision;
		double H_n, H_n_plus_1;

		// We'll solve the Kepler equation (M = e.sinh(H) - H) with the Halley method

		// Calculate a first approximation
		H_n_plus_1 = GetFirstApproximation_Hyperbolic(meanAnomaly);

		if (Math.Abs(H_n_plus_1) > 0.0)
		{
			uint iteration_number = 0;
			double F0, F1, F2;

			// Perform Halley iterations (1 or 2 are enough in practice, 3 may be needed in worst cases)
			// -------------------------
			do
			{
				iteration_number++;

				// start from the previous estimation (the first approximation if it's the first iteration)
				H_n = H_n_plus_1;

				// Calculating the function's value (F0), its derivative (F1) and second derivative (F2) in H_n
				F2 = ecc * Math.Sinh(H_n);
				F1 = ecc * Math.Cosh(H_n) - 1.0;
				F0 = F2 - H_n - meanAnomaly;

				// Calculate the next approximation
				H_n_plus_1 = H_n - 2.0 * F1 * F0 / (2.0 * F1 * F1 - F0 * F2);
	
				// relative precision
				relative_precision = Math.Abs((H_n_plus_1 - H_n) / H_n_plus_1);

			} while ((relative_precision > C_PRECISION_HALLEY) && (iteration_number < NB_MAX_ITERATIONS));
		}

		// The solution has been found, calculate the corresponding radius
		// ---------------------------------------------------------------
		double cosh_eccentricAnomaly = Math.Cosh(H_n_plus_1);
		radius = mu / (2.0 * K) * (ecc * cosh_eccentricAnomaly - 1.0);

		if (radius > 0.0)
		{
			double cos_arg = (cosh_eccentricAnomaly - ecc) / (1.0 - ecc * cosh_eccentricAnomaly);

			if (cos_arg >  1.0) { cos_arg =  1.0; }
			if (cos_arg < -1.0) { cos_arg = -1.0; }

			trueAnomaly = Math.Acos(cos_arg);
			if (time < 0) { trueAnomaly = -trueAnomaly; }
		}
		else
		{
			// radius is null, true anomaly is not defined
			trueAnomaly = 0.0;
		}
	}

	private static double GetFirstApproximation_QuasiParabolic(double YVal)
	{
		// Calculate a first approximation by solving the corresponding grade 3 equation (by neglecting the terms of greater degree)
		// The equation is under the form x^3 + p.x + q = 0 and is solved by the Cardan method
		double grade3_coeff = (1.0 + 2.0 * gamma) / 3.0;
		double p = peri / grade3_coeff;
		double q = -YVal / grade3_coeff;

		double sqrt_delta = Math.Sqrt(q*q + 4.0 * p*p*p / 27.0);

		double Q0 = CubicRoot((-q - sqrt_delta) / 2.0) + CubicRoot((-q + sqrt_delta) / 2.0);

		return Q0;
	}

	private static double GetFirstApproximation_Elliptic(double meanAnomaly)
	{
		// The Kepler equation is approximated with 5 linear equations defined on 5 intervals (separated by 4 bounds)
		// Each bound is computed as "A * eccentricity + B"
		// This allows to have a very precise approximation
		// ----------------------------------------------------------------------------------------------------------
		const double C_COMPUTING_STARTING_VAL_BOUND1_VAL_A = -0.684853256372279547; // -(sqrt(3) - PI/3)
		const double C_COMPUTING_STARTING_VAL_BOUND1_VAL_B = 0.684853256372279547; // sqrt(3) - PI/3
		const double C_COMPUTING_STARTING_VAL_BOUND2_VAL_A = -1.0;
		const double C_COMPUTING_STARTING_VAL_BOUND2_VAL_B = 1.315146743627720452; // 2 - sqrt(3) + PI/3
		const double C_COMPUTING_STARTING_VAL_BOUND3_VAL_A = -1.0;
		const double C_COMPUTING_STARTING_VAL_BOUND3_VAL_B = 1.826445909962072785; // 2*PI/3 - 2 + sqrt(3)
		const double C_COMPUTING_STARTING_VAL_BOUND4_VAL_A = -0.684853256372279547; // -(sqrt(3) - PI/3)
		const double C_COMPUTING_STARTING_VAL_BOUND4_VAL_B = 2.456739397217513691; // 4*PI/3 - sqrt(3)

		// 2 constants used for the linear approximation
		const double C_COMPUTING_STARTING_VAL_CONSTANT_2 = 0.342426628186139773; // sqrt(3)/2 - PI/6
		const double C_COMPUTING_STARTING_VAL_CONSTANT_4 = 1.913222954981036392; // sqrt(3)/2 + PI/3

		// For higher eccentricities, a more accurate method is used to get a first approximation if mean anomaly is low
		const double C_HIGH_ECCENTRICITY_THRESHOLD = 0.8;

		double E0;
		double Y = Math.Abs(meanAnomaly);

		if (Y < C_COMPUTING_STARTING_VAL_BOUND2_VAL_A * ecc + C_COMPUTING_STARTING_VAL_BOUND2_VAL_B)
		{
			if (ecc > C_HIGH_ECCENTRICITY_THRESHOLD)
			{
				// High eccentricities require a more accurate approximation
				// --> Solve Y = (1-_ecc) * x + 1/6 . x^3 with Cardan method
				double p = 6.0 * (1.0 - ecc) / ecc;
				double q = -6.0 * Y / ecc;

				// Get the unique solution of the equation under the form x^3 + p*x + q = 0
				double sqrt_delta = Math.Sqrt(q*q + 4.0 * p*p*p / 27.0);
				E0 = CubicRoot((-q - sqrt_delta) / 2.0) + CubicRoot((-q + sqrt_delta) / 2.0);
			}
			else
			{
				if (Y < C_COMPUTING_STARTING_VAL_BOUND1_VAL_A * ecc + C_COMPUTING_STARTING_VAL_BOUND1_VAL_B)
				{
					// Y is between 0 and bound 1
					E0 = Y / (1 - ecc);

				}
				else
				{
					// Y is between bound 1 and bound 2
					E0 = (Y + ecc * C_COMPUTING_STARTING_VAL_CONSTANT_2) / (1 - ecc / 2.0);
				}
			}
		}
		else
		{
			if (Y < C_COMPUTING_STARTING_VAL_BOUND3_VAL_A * ecc + C_COMPUTING_STARTING_VAL_BOUND3_VAL_B)
			{
				// Y is between bound 2 and bound 3
				E0 = Y + ecc;
			}
			else if (Y < C_COMPUTING_STARTING_VAL_BOUND4_VAL_A * ecc + C_COMPUTING_STARTING_VAL_BOUND4_VAL_B)
			{
				// Y is between bound 3 and bound 4
				E0 = (Y + ecc * C_COMPUTING_STARTING_VAL_CONSTANT_4) / (1 + ecc / 2.0);
			}
			else
			{
				// Y is between bound 4 and bound 5
				E0 = (Y + ecc * Math.PI) / (1 + ecc);
			}
		}

		// Negate found value if Y value is negative
		if (meanAnomaly < 0.0) { E0 = -E0; }

		return E0;
	}

	private static double GetFirstApproximation_Hyperbolic(double meanAnomaly)
	{
		double H0;
		double Y = Math.Abs(meanAnomaly);

		// Compute approximation
		if (Y < 1.0)
		{
			// Solve the cubic approximation with the Cardan method
			double sqrt_delta = Math.Sqrt(9.0 * Y*Y + 8.0 * Math.Pow(ecc - 1.0, 3.0) / ecc);

			H0 = CubicRoot((3.0 * Y - sqrt_delta) / ecc) + CubicRoot((3.0 * Y + sqrt_delta) / ecc);
		}
		else
		{
			// Use M = _ecc/2 * (exp(x) - 1) - (exp(x/2) - 1), which works great with high values
			H0 = 2.0 * Math.Log((1.0 + Math.Sqrt(1 + ecc * (2.0 * Y + ecc - 2.0))) / ecc);
		}

		// negate result is original value was negative
		if (meanAnomaly < 0.0) { H0 = -H0; }

		return H0;
	}

	private static double CubicRoot(double n)
	{
		double cbrt = Math.Pow(Math.Abs(n), 1.0 / 3.0);
		if(n<0.0) { cbrt = -cbrt; }

		return cbrt;
	}
}



