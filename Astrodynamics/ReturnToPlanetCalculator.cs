using SFS.World;
using SFS.WorldBase;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.Tracing;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

class ReturnToPlanetCalculator
{
    private static double R_start = 0.0; // The radius at which we start from
    private static double arg_start = 0.0; // The argument at which we start from
    private static double R_end = 0.0;   // The radius we target
    private static double R_max = 0.0; // The SOI limit, to avoid having trajectories that go past that limit
    

    private static double Vorb = 0.0;

    private static double lambda = 0.0; // lambda = R_start / R_end
    private static double sign = 1.0; // 1.0 if lambda > 1, -1.0 if lambda < 1
    private static double x0 = 0.0; // the minimal or maximal solution depending on cases

    private static double wr0 = 0.0; // initial radial speed (divided by orbital velocity, in absolute value)
    private static double wt0 = 0.0; // initial tangential speed (divided by orbital velocity, in absolute value)

    private static double sign_wr = 0.0;
    private static double sign_wt = 0.0;

    private static Planet planet = null;
    private static double currentTime = 0.0;

    private static bool includeInsertionCost;


    // CONSTANTS TO DEFINE THE BEHAVIOUR OF THE ALGORITHM
    // --------------------------------------------------

    // EPSILON_X: an x value that is closer to x0 than by this value will be equal to x0
    //            This is to avoid singularity problems around x0
    private const double EPSILON_X = 0.00000001;

    // LARGE_EPSILON_X: a slightly larger epsilon threshold to make sure that we test a value
    //                  that's close to x0 while still not too close
    private const double LARGE_EPSILON_X = 1.0625 * EPSILON_X;

    // MIN_DELTA_W_THRESHOLD: If a delta_w (delta_v divided by orbital speed) is less than this value, it's considered to be 0
    //                        This is to avoid divisions by 0 in the search process that would occur if there was a solution using 0 m/s of delta_V.
    //                        NOTE: it can be simply interpreted as the minimal fraction of orbital speed (at targeted radius) that we allow (delta_V/orbital_speed must be higher than this threshold)
    private const double MIN_DELTA_W_THRESHOLD = 0.0001;

    // MAX_RADIUS_FACTOR: When we have to set a limit for the apoapsis, this limit will be MAX_RADIUS_FACTOR times the targeted radius.
    //                    NOTE: The SOI limit will be used instead if it appears to be lower.
    private const double MAX_RADIUS_FACTOR = 10000.0;

    // NEWTON_PRECISION_THRESHOLD: This is the precision threshold that will be applied when using the Newton method
    private const double NEWTON_PRECISION_THRESHOLD = 0.00000001;

    private const double MAX_ECCENTRICITY_FINAL_ORBIT = 0.01;

    public static bool calculateTrajectory(Planet planet, Location location, double targetRadius, bool includeInsertionCost, bool suggestOrbitInsertion, out Orbit orbit)
    {
        orbit = null;

        // Initialize parameters
        // ---------------------
        InitParameters(planet, location, targetRadius, includeInsertionCost);


        // Manage the case when the ship is practically at the target altitude
        // -------------------------------------------------------------------
        if(suggestOrbitInsertion && (Math.Abs(lambda - 1.0) < MAX_ECCENTRICITY_FINAL_ORBIT))
        {
            if(includeInsertionCost)
            {
                orbit = generateFinalOrbit();
                return (orbit != null);
            }
        }
        
        // Solve the case without insertion first (aka fly-by mode)
        // --------------------------------------
        bool success = SolveForX_WithoutInsertion(out double x_no_insert);

        if (!success)
        {
            LOG(LOG_LEVEL.ERROR, "deltaV1 squared value could not be calculated for x0");
            return false;
        }

        // This will be the default value
        double x = x_no_insert;

        // Solve the case with insertion (rendez-vous mode)
        // -----------------------------
        if (includeInsertionCost)
        {
            // Check the minimum deltaW1 value (obtained for the x_no_insert value calculated above) - if it's null (or practically), then this is considered the optimal x in Rendez-vous mode
            success = getDeltaV1Squared_value(x_no_insert, out double dw1sqr, out double derivative, out double derivative2);
            LOG(LOG_LEVEL.DEBUG, "  DeltaV (no insertion cost) = " + (Math.Sqrt(dw1sqr) * Vorb));

            if (!success)
            {
                LOG(LOG_LEVEL.ERROR, "deltaV1 squared value could not be calculated for x_no_insert = " + x_no_insert + "; lambda = " + lambda + "; x0 = " + x0);
                return false;
            }

            if (dw1sqr * lambda > MIN_DELTA_W_THRESHOLD * MIN_DELTA_W_THRESHOLD) // NOTE : the multiplication by lambda is important: dw1 is dv1 expressed as a fraction of orbital speed at start radius
                                                                                 //        but we want to reason with the orbital speed of the end radius. It needs a multiplication by lambda for that.
            {
                // Search 
                success = SolveForX(out x, x_no_insert);

                if (!success)
                {
                    LOG(LOG_LEVEL.ERROR, "Could not determine a minimum");
                    return false;
                }

                getDeltaV1Squared_value(x, out dw1sqr, out derivative, out derivative2);
                LOG(LOG_LEVEL.DEBUG, "  DeltaV (with insertion) = " + (Math.Sqrt(dw1sqr) * Vorb));
            }
            else
            {
                LOG(LOG_LEVEL.DEBUG, "deltaV1 squared value is considered null for x = " + x_no_insert + " (no insertion case); that value will be used for the resolution");
            }
        }
        

        orbit = DetermineOrbit(x);
        return true;
    }


    // ------------------------------------------------------------------------
    // -----                       InitParameters                         -----
    // ------------------------------------------------------------------------
    // Allows to initialize all internal parameters for the calculation
    private static void InitParameters(Planet planet, Location location, double targetRadius, bool includeInsertionCost)
    {
        // Init all parameters
        R_start = location.Radius;
        arg_start = location.position.AngleRadians; // argument, between -Pi and +Pi
        R_end = targetRadius;

        // Determine the maximum radius in case we have to set a limit to the apoapsis
        R_max = Math.Max(MAX_RADIUS_FACTOR * targetRadius, 2.0 * R_start);
        if (R_max > planet.SOI) { R_max = planet.SOI; }

        lambda = R_start / targetRadius;

        if (lambda > 1.0) sign = 1.0;
        else sign = -1.0;

        // orbital velocity at the ship's position (used to normalize the equations...)
        Vorb = Math.Sqrt(planet.mass / R_start);

        // Velocity at the starting point in polar coordinates, divided by orbital velocity at that level to normalize them.
        wr0 = (location.velocity.x * Math.Cos(arg_start) + location.velocity.y * Math.Sin(arg_start)) / Vorb;
        wt0 = (-location.velocity.x * Math.Sin(arg_start) + location.velocity.y * Math.Cos(arg_start)) / Vorb;

        // memorize the sign of the velocity components
        if (wr0 > 0) sign_wr = 1.0;
        else sign_wr = -1.0;

        if (wt0 > 0) sign_wt = 1.0;
        else sign_wt = -1.0;

        // Then take the absolute value
        wr0 = Math.Abs(wr0);
        wt0 = Math.Abs(wt0);

        x0 = Math.Sqrt(2.0 / (1.0 + lambda));

        ReturnToPlanetCalculator.planet = planet;
        currentTime = location.time;

        ReturnToPlanetCalculator.includeInsertionCost = includeInsertionCost;
    }


    // ------------------------------------------------------------------------
    // -----                 SolveForX_WithoutInsertion                   -----
    // ------------------------------------------------------------------------
    // Allows to calculate the x value that minimizes the transfer with no insertion cost (aka fly-by mode)
    // The way it's calculated is easier and less risky since there's no square root in the equation, so the
    // derivative has no definition problem. In particular, there's no definition problem if the minimum value
    // happens to be 0. Even in the case of a transfer with insertion cost included, it's a good idea to
    // calculate this one first, to prevent some problems that could happen with the insertion took into account.
    private static bool SolveForX_WithoutInsertion(out double x)
    {
        x = x0;

        double x0_sqr = x0 * x0;
        double lambda_sqr = lambda * lambda;

        LOG(LOG_LEVEL.DEBUG, " Solving the \"no insertion\" case - lambda = " + lambda + "; x0 = " + x0 + "; wt0 = " + wt0 + "; wr0 = " + wr0);

        if (lambda < 1.0)
        {
            double C = 2.0 * (1.0 - lambda) * wr0 * wr0 / (lambda_sqr * lambda_sqr);
            double sqrt_C = Math.Sqrt(C);

            // Calculate the minimum candidate value x_min1 based on the study of P1 and P2 on the upper bound
            // --------------------------------------------
            double x_min1 = 0.0;
            
            if(x0 < wt0 / lambda_sqr)
            {
                double A = wt0 / lambda_sqr - x0;
                double B = sqrt_C / x0;

                if(B > 0.000001 * A)
                {
                    A *= A;
                    B *=B;

                    x_min1 = x0 * Math.Sqrt(1.0 - B / (A + B));
                }
                else
                {
                    // If B is really lower than A (or A and B are both 0!), there would be no measurable difference so we assume the answer is x0
                    x_min1 = x0;
                }
                
                LOG(LOG_LEVEL.DEBUG, "   Val min based on the study of P1 and P2 on the upper bound (x0 < wt0/lambda_sqr): x_min1 = " + x_min1);
            }
            else if(x0 > wt0 / lambda_sqr)
            {
                x_min1 = wt0 / lambda_sqr;
                double denominator = x0 * Math.Sqrt(x0_sqr - x_min1 * x_min1);

                if(sqrt_C > denominator)
                {
                    x_min1 = 0.0;
                }
                else if(sqrt_C > 0.000000000001 * denominator)
                {
                    x_min1 = (1.0 - sqrt_C / denominator) * x_min1;
                }
                
                LOG(LOG_LEVEL.DEBUG, "   Val min based on the study of P1 and P2 on the upper bound (x0 > wt0/lambda_sqr): x_min1 = " + x_min1);
            }

            // Calculate the minimum candidate value x_min2 based on which value P1 and P2 are valued to sqrt(C)
            // --------------------------------------------
            double x_min2 = 0.0;

            {
                // Calculate x1
                double x1 = sqrt_C - C / x0_sqr;

                if (x1 > 0.0) { x1 = Math.Sqrt(x1); }
                else { x1 = 0.0; }

                x1 = wt0 / lambda_sqr - x1;

                if (x1 < 0.0) { x1 = 0.0; }

                // Calculate x2
                double x2 = x0_sqr - sqrt_C;

                if (x2 > 0.0) { x2 = Math.Sqrt(x2); }
                else { x2 = 0.0; }

                // take the lowest value of the 2 as the minimum value (we know that x has to be between x1 and x2)
                x_min2 = Math.Min(x1, x2);
                LOG(LOG_LEVEL.DEBUG, "   Val min based on where P1 and P2 are valued to sqrt_C: x_min2 = " + x_min2);
            }
            

            // Calculate initial guess (x is intentionally initialized to the lowest bound - if new_x was equal to it it's the solution,
            //                          but the derivative would be null and Newton would fail, so this allows to skip the loop) 
            x = Math.Min(x0, wt0 / lambda_sqr);
            double new_x = Math.Max(x_min1, x_min2); // take the best of the 2 minimum as the initial estimation

            int i_iter = 0;

            while((i_iter < 32) && (Math.Abs(new_x - x) > 0.00000001))
            {
                i_iter++;
                x = new_x;

                // Calculating the derivative in x...
                double P1 = (wt0 / lambda_sqr - x);
                P1 = P1 * P1 + C / x0_sqr;

                double P2 = x0_sqr - x * x;

                double P1_prime = -2.0 * (wt0 / lambda_sqr - x);
                double P2_prime = -2.0 * x;

                double val = P1 * P2 - C;
                double P_prime = P1 * P2_prime + P2 * P1_prime;
                double P_sec_halved = -P2 + P1_prime * P2_prime - P1; // second derivative halved - implicitely takes into account that the second derivative is just -2 for P1 and P2

                LOG(LOG_LEVEL.DEBUG, "  Halley iteration " + i_iter + ": x = " + x + "; val = " + val + "; derivative = " + P_prime + "; second derivative (halved) = " + P_sec_halved);

                // Apply Halley formula
                new_x = x - val * P_prime / (P_prime * P_prime - val * P_sec_halved);
            }

            if(Math.Abs(new_x - x) < 0.00000001)
            {
                x = new_x;
                LOG(LOG_LEVEL.DEBUG, " Success: x = " + x);
                return true;
            }
            else
            {
                LOG(LOG_LEVEL.ERROR, " Failed to calculate Return to planet without insertion");
                return false;
            }
        }
        else // lambda > 1.0
        {
            double C = 2.0 * Math.Abs(lambda - 1.0) * wr0 * wr0 / (lambda_sqr * lambda_sqr);
            double sqrt_C = Math.Sqrt(C);
            double x_min = Math.Max(x0, wt0 / lambda_sqr);

            // Test the maximal bound if the ship is getting further from the planet (the apoapsis shouldn't go passed the SOI in this case)
            // ----------------------
            double x_max = Math.Sqrt(2.0 / (lambda + R_start / R_max));

            if (sign_wr > 0.0)
            {
                LOG(LOG_LEVEL.DEBUG, " Testing x_max = " + x_max + "...");
                if (!(x_max > x0))
                {
                    // This can happen if the player is on a very high position, higher than what we consider the maximal radius
                    LOG(LOG_LEVEL.DEBUG, "  x_max is actually lower than x0 --> use x0!");
                    x = x0;
                    return true;
                }
                else if(!(x_max > wt0 / lambda_sqr))
                {
                    // This can happen if the player passed the target altitude, and there's no optimal solution for a fly-by in the future -> we force the use of x_max
                    LOG(LOG_LEVEL.DEBUG, "  x_max is actually lower than wt0 / lambda_sqr --> use x_max!");
                    x = x_max;
                    return true;
                }
                else
                {
                    // Test if x_max is lower than the solution - if it's the case then use x_max
                    double P1 = (x_max - wt0 / lambda_sqr);
                    P1 = P1 * P1 - C / x0_sqr;

                    double P2 = x_max * x_max - x0_sqr;

                    if(!(P1*P2 > C))
                    {
                        LOG(LOG_LEVEL.DEBUG, "  the solution would be past x_max --> use x_max!");
                        x = x_max;
                        return true;
                    }
                }
            }

            // Calculate x1
            double x1 = Math.Sqrt(sqrt_C + C / x0) + wt0 / lambda_sqr;

            // Calculate x2
            double x2 = Math.Sqrt(x0_sqr + sqrt_C);

            // Calculate initial guess (x is intentionally initialized to the lowest bound - if new_x was equal to it it's the solution,
            //                          but the derivative would be null and Newton would fail, so this allows to skip the loop) 
            x = x_min;
            double new_x = Math.Max(x1, x2);

            int i_iter = 0;

            while ((i_iter < 32) && (Math.Abs(new_x - x) > 0.00000001))
            {
                i_iter++;
                x = new_x;

                // Calculating the derivative in x...
                double P1 = (x - wt0 / lambda_sqr);
                P1 = P1 * P1 - C / x0_sqr;

                double P2 = x * x - x0_sqr;

                double P1_prime = 2.0 * (x - wt0 / lambda_sqr);
                double P2_prime = 2.0 * x;

                double val = P1 * P2 - C;
                double P_prime = P1 * P2_prime + P2 * P1_prime;
                double P_sec_halved = P2 + P1_prime * P2_prime + P1; // second derivative halved - implicitely takes into account that the second derivative is just 2 for P1 and P2

                LOG(LOG_LEVEL.DEBUG, "  Halley iteration " + i_iter + ": x = " + x + "; val = " + val + "; derivative = " + P_prime + "; second derivative (halved) = " + P_sec_halved);

                // Apply Halley formula
                new_x = x - val * P_prime / (P_prime * P_prime - val * P_sec_halved);
            }

            if (Math.Abs(new_x - x) < 0.00000001)
            {
                x = new_x;
                LOG(LOG_LEVEL.DEBUG, " Success: x = " + x);
                return true;
            }
            else
            {
                LOG(LOG_LEVEL.ERROR, " Failed to calculate Return to planet without insertion");
                return false;
            }
        }
    }


    // ------------------------------------------------------------------------
    // -----                   getDeltaV1Squared_value                    -----
    // ------------------------------------------------------------------------
    // Allows to calculate the "deltaV1Squared" value and its derivatives for the given value of x.
    // deltaV1Squared is the deltaV (squared) the player will have to apply if he performs the corresponding transfer.
    // NOTE 1 : It doesn't take into account the insertion cost in orbit, this is only part of the job
    // NOTE 2 : The derivatives are undefined for x = x0, so those values will be infinite (+/-) if x is too close to x0.
    //          "too close" means if x differs from x0 by less than EPSILON_X.
    private static bool getDeltaV1Squared_value(double x, out double val, out double derivative, out double derivative2)
    {
        // Check that value is in range
        if( ((lambda > 1) && (x < x0))                  ||   // If lambda > 1, x must be in [x0, +infinity[
            ((lambda < 1) && ((x < 0.0) || (x > x0))) )      // If lambda < 1, x must be in [0.0, x0] - Note that if lambda is exactly 1, it doesn't matter.
        {
            LOG(LOG_LEVEL.ERROR, "Attempt to get deltaV value for incorrect argument x = " + x + " (lambda = " + lambda + ", x0 = " + x0 + ")");
            val = 0.0;
            derivative = 0.0;
            derivative2 = 0.0;
            return false; // return false to let the user know that an error occured
        }

        double lambda_sqr = lambda * lambda;
        double x_sqr = x * x;
        double x0_sqr = x0 * x0;

        // Calculate the deltaV1 value, in normalized form
        // ---------------------------
        double val1 = Math.Sqrt(Math.Abs((lambda_sqr - 1) * (x_sqr - x0_sqr))) - wr0; // absolute value in the sqrt should not be needed, it's just in case
        double val2 = x - wt0;

        val = val1 * val1 + val2 * val2;

        // calculate the derivatives
        // -------------------------
        if(Math.Abs(x-x0) < EPSILON_X)
        {
            // We are on the singularity (x0), the derivatives are not defined
            // Return +/-Infinity for the derivatives, It's the responsibility of the caller to take care about that case!
            if (lambda < 1.0) derivative = Double.PositiveInfinity; 
            else derivative = Double.NegativeInfinity;

            derivative2 = Double.PositiveInfinity;
        }
        else
        {
            // We are far enough from the singularity, derivatives are defined.
            // "quotient" is to simplify the following expressions and calculate it only once.
            double quotient = Math.Sqrt(Math.Abs((lambda_sqr - 1) / (x_sqr - x0_sqr))) * wr0; // absolute value in the sqrt should not be needed, it's just in case

            // calculate the derivative
            derivative = quotient * x;
            if (lambda > 1.0) derivative = -derivative;

            derivative += lambda_sqr * x - wt0;
            derivative *= 2.0;

            // calculate the second derivative
            derivative2 = 2.0 * (lambda_sqr + quotient * x0_sqr / Math.Abs(x0_sqr - x_sqr)); // This absolute value is intentional however...
        }

        return true;
    }


    // ------------------------------------------------------------------------
    // -----                      calculateEtaValues                      -----
    // ------------------------------------------------------------------------
    // Allows to calculate the "Eta" value and its derivative for the given value of x.
    // Eta is the derivative of the whole transfer cost multiplied by 2 times the DeltaV1 cost ("DeltaV1" is sqrt_f0 in the function).
    // Since we want to minimize the transfer cost, it's equivalent to searching the value that makes the derivative (so eta) null.
    // Eta's derivative is also calculated to allow the use of the Newton method.
    // NOTE 1 : Eta values are undefined for x0, so it's important that the caller makes sure not to find himseld in that case.
    // NOTE 2 : the DeltaV1 is supposed to be not null, it's the user responsibility to make sure we can't fall into those cases.
    private static bool calculateEtaValues(double x, out double val, out double derivative)
    {
        // Init values
        val = 0.0;
        derivative = 0.0;

        // Check that value is in range
        if (((lambda > 1) && (x < x0 - EPSILON_X)) ||               // If lambda > 1, x must be in ]x0, +infinity[
            ((lambda < 1) && ((x < 0.0) || (x > x0 + EPSILON_X))))  // If lambda < 1, x must be in [0.0, x0[ - Note that if lambda is exactly 1, it doesn't matter.
        {
            LOG(LOG_LEVEL.ERROR, "Attempt to get deltaV value for incorrect argument x = " + x + " (lambda = " + lambda + ", x0 = " + x0 + ", wr0 = " + wr0 + ", wt0 = " + wt0 + ")");
            return false; // return false to let the user know that an error occured
        }

        // We need the dv1_squared value and its derivative (noted f0, f1 and f2 here)
        double f0, f1, f2;

        bool success = getDeltaV1Squared_value(x, out f0, out f1, out f2);
        if (!success) return false;

        LOG(LOG_LEVEL.DEBUG, "Calculate Eta values for x = " + x + " (lambda = " + lambda + "; x0 = " + x0 + ", wr0 = " + wr0 + ", wt0 = " + wt0 + ")");
        LOG(LOG_LEVEL.DEBUG, "  f0 = " + f0 + "; f1 = " + f1 + "; f2 = " + f2);

        val = f1;
        derivative = f2;

        if(includeInsertionCost)
        {
            double sqrt_f0 = Math.Sqrt(f0);

            if(lambda < 1.0)
            {
                val -= 2.0 * lambda * sqrt_f0;
                derivative -= lambda * f1 / sqrt_f0;
            }
            else
            {
                val += 2.0 * lambda * sqrt_f0;
                derivative += lambda * f1 / sqrt_f0;
            }
        }

        LOG(LOG_LEVEL.DEBUG, "  eta = " + val + "; derivative = " + derivative);

        return true;

    }


    private static bool SolveForX(out double x, double x_no_insert)
    {
        x = 0.0;

        // We search the values of x that will make eta equal to 0
        double eta, eta_prime;
        double x_start = 0.0;
        bool success;

        if(lambda > 1.0)
        {
            // Test the minimal value
            // ----------------------
            LOG(LOG_LEVEL.DEBUG, "Testing x0 + eps = " + (x0 + LARGE_EPSILON_X) + "...");
            success = calculateEtaValues(x0 + LARGE_EPSILON_X, out eta, out eta_prime);

            if (!success) return false;
            
            if(eta > 0.0)
            {
                // calculate an approximate value of x based on an approximation of the derivative of eta close to x0
                double num = Math.Sqrt(2.0 * (lambda * lambda - 1.0) * x0) * wr0;
                double eps = num / (eta + num / Math.Sqrt(LARGE_EPSILON_X));
                x = x0 + eps * eps;

                LOG(LOG_LEVEL.DEBUG, "x0 + eps appears to give a positive value of eta - take approximate value of x = " + x);
                return true;
            }

            // Test the maximal bound if the ship is getting further from the planet (the apoapsis shouldn't go passed the SOI in this case)
            // ----------------------
            double x_max = Math.Sqrt(2.0 / (lambda + R_start / R_max));

            if (sign_wr > 0.0)
            {
                LOG(LOG_LEVEL.DEBUG, "Testing x_max = " + x_max + "...");
                if(x_max < x0 + EPSILON_X)
                {
                    // This can happen if the player is on a very high position, higher than what we consider the maximal radius (Note: VERY unlikely now that Rmax is at least something times start radius...)
                    LOG(LOG_LEVEL.DEBUG, "  x_max is actually lower than x0 --> use x0!");
                    x = x0;
                    return true;
                }
                else
                {
                    success = calculateEtaValues(x_max, out eta, out eta_prime);

                    if (!success) return false;

                    if (eta < 0.0)
                    {
                        // Note : eta is a function that increases. Since x_max is the highest possible value, if eta is still negative then it's the minimum
                        LOG(LOG_LEVEL.DEBUG, "x_max appears to give a negative value of eta - x_max is considered to be the minimum");
                        x = x_max;
                        return true;
                    }
                }
            }

            // Determine a starting value for the Newton method
            // --------------------------

            // The solution must be between x0 + LARGE_EPSILON_X and x_no_insert --> we start from the right bound
            // and move gradually closer to x0 while calculating eta until we find an x value for which eta is negative.
            x_start = x_no_insert;

            do
            {
                // We move the starting value by this value to the left (getting closer to x0 by a factor 4)
                double delta_x = (x_start - x0) / 4.0;

                if (delta_x < LARGE_EPSILON_X)
                {
                    // Make sure we don't closer to x0 than by LARGE_EPSILON_X
                    x_start = x0 + LARGE_EPSILON_X;
                    break;
                }
                else
                {
                    // calculate the new start value and its associated eta value
                    x_start = x0 + delta_x;

                    success = calculateEtaValues(x_start, out eta, out eta_prime);
                    if (!success)
                    {
                        LOG(LOG_LEVEL.ERROR, " eta could not be calculated for x_start = " + x_start);
                        return false;
                    }
                    else
                    {
                        LOG(LOG_LEVEL.DEBUG, "Testing x_start = " + x_start + " : eta = " + eta);
                    }
                }
            } while (eta > 0.0);
        }
        else // lambda < 1.0
        {
            // Test the maximal value
            // ----------------------
            LOG(LOG_LEVEL.DEBUG, "Testing x0 - eps = " + (x0 - LARGE_EPSILON_X) + "...");
            success = calculateEtaValues(x0 - LARGE_EPSILON_X, out eta, out eta_prime);

            if (!success) return false;

            if (eta < 0.0)
            {
                // calculate an approximate value of x based on an approximation of the derivative of eta close to x0
                double num = Math.Sqrt(2.0 * (1.0 - lambda * lambda) * x0) * wr0;
                double eps = num / (-eta + num / Math.Sqrt(LARGE_EPSILON_X));
                x = x0 - eps * eps;

                LOG(LOG_LEVEL.DEBUG, "x0 - eps appears to give a negative value of eta - take approximate value of x = " + x);
                return true;
            }

            // Determine a starting value for the Newton method
            // --------------------------

            // The solution must be between x_no_insert and x0 - LARGE_EPSILON_X --> we start from the left bound
            // and move gradually closer to x0 while calculating eta until we find an x value for which eta is positive.
            x_start = x_no_insert;

            do
            {
                // We move the starting value by this value to the right (getting closer to x0 by a factor 4)
                double delta_x = (x0 - x_start) / 4.0;

                if(delta_x < LARGE_EPSILON_X)
                {
                    // Make sure we don't closer to x0 than by LARGE_EPSILON_X
                    x_start = x0 - LARGE_EPSILON_X;
                    break;
                }
                else
                {
                    // calculate the new start value and its associated eta value
                    x_start = x0 - delta_x;
                    
                    success = calculateEtaValues(x_start, out eta, out eta_prime);
                    if (!success)
                    {
                        LOG(LOG_LEVEL.ERROR, " eta could not be calculated for x_start = " + x_start);
                        return false;
                    }
                    else
                    {
                        LOG(LOG_LEVEL.DEBUG, "Testing x_start = " + x_start + "; eta = " + eta);
                    }
                }
            } while (eta < 0.0);
        }


        // Apply the Newton method
        // -----------------------
        const int nb_max_iterations = 20;
        int i_Newton = 0;

        double new_x = x_start;
        double relativeDifference = Double.PositiveInfinity;

        do
        {
            i_Newton++;

            // Set x on the value found on the previous iteration
            x = new_x;

            LOG(LOG_LEVEL.DEBUG, "Newton method (iteration " + i_Newton + "): testing x = " + x + "...");
            success = calculateEtaValues(x, out eta, out eta_prime);
            if (!success) return false;

            // calculate new x estimate
            new_x = x - eta / eta_prime;

            // Test that the value isn't out of bounds...
            if(lambda < 1.0)
            {
                if ((new_x < 0.0) || (new_x > x0 - EPSILON_X))
                {
                    LOG(LOG_LEVEL.DEBUG, "Newton method (iteration " + i_Newton + "): new_x = " + new_x + "is out of bounds - give up!");
                    return false;
                }
            }
            else
            {
                if (new_x < x0 + EPSILON_X)
                {
                    LOG(LOG_LEVEL.DEBUG, "Newton method (iteration " + i_Newton + "): new_x = " + new_x + "is out of bounds - give up!");
                    return false;
                }
            }

            // Calculate relative difference
            relativeDifference = Math.Abs(new_x - x);

            if (new_x > NEWTON_PRECISION_THRESHOLD) relativeDifference /= new_x;
            else relativeDifference /= NEWTON_PRECISION_THRESHOLD;

        } while ((i_Newton < nb_max_iterations) && (relativeDifference > NEWTON_PRECISION_THRESHOLD));

        if(relativeDifference > NEWTON_PRECISION_THRESHOLD)
        {
            LOG(LOG_LEVEL.ERROR, "Problem with Newton method: too many iterations - Exit");
            return false;
        }

        // SUCCESS
        x = new_x;
        LOG(LOG_LEVEL.INFO, "Successfully found value x = " + x + " for return to mother planet calculation");
        return true;
    }

    private static Orbit DetermineOrbit(double x) 
    {
        double slr = R_start * x * x;
        double ecc;
        double periapsis;

        if(lambda > 1.0)
        {
            // R_end is a periapsis...
            ecc = slr / R_end - 1.0;
            periapsis = R_end;
        }
        else
        {
            // R_end is an apoapsis...
            ecc = 1.0 - slr / R_end;
            periapsis = slr / (1.0 + ecc);
        }

        // specific energy - needed for the Kepler solver and time determination
        double specificEnergy = planet.mass * (ecc * ecc - 1.0) / slr / 2.0;

        double argOfPeriapsis = 0.0; // default value (used for circular orbits)
        double trueAnomaly = 0.0;

        // calculate argument of periapsis
        if (1.0 + ecc > 1.0) // eccentricity != 0.0 - tested that way (instead of ecc > 0.0) so that ecc is considered as zero if it's negligible (< 10^(-16))
        {
            double theCosinus = (slr / R_start - 1.0) / ecc;

            // to fix for numerical imprecisions...
            if (theCosinus > 1.0) theCosinus = 1.0;
            if (theCosinus < -1.0) theCosinus = -1.0;

            trueAnomaly = Math.Acos(theCosinus);

            if(sign_wt * sign_wr < 0.0)
            {
                trueAnomaly = -trueAnomaly;
            }

            argOfPeriapsis = arg_start - trueAnomaly;
        }

        // Determine path type
        PathType pathType;

        if(includeInsertionCost) { pathType = PathType.Encounter; } // if insertion is required, always encounter
        else if((sign_wr < 0) && (slr > planet.SOI*(1.0 - ecc))) { pathType = PathType.Escape; } // if insertion is not required and apoapsis above SOI, Escape, but ONLY if the ship approaches (escape trajectories are forbidden otherwise)
        else { pathType = PathType.Eternal; } // Otherwise eternal orbit

        Orbit orbit = Orbit_Utils.CreateOrbit(slr, ecc, argOfPeriapsis, (int)sign_wt, planet, pathType, null);

        if (orbit != null)
        {
            // Set periapsis passage time
            double timeFromPeri = 0.0;
            new KeplerSolver(planet.mass, periapsis, specificEnergy).GetTimeAtPosition(R_start, trueAnomaly, ref timeFromPeri);
            orbit.periapsisPassageTime = currentTime - timeFromPeri * (double)orbit.direction - 10.0 * orbit.period;

            double argDestination = argOfPeriapsis;
            if (lambda < 1.0) argDestination += Math.PI; // If we start below the target radius, the target radius is an apoapsis... so argument is arg of periapsis + PI

            // Set start and end time
            orbit.orbitStartTime = currentTime;

            if(pathType == PathType.Encounter)
            {
                orbit.orbitEndTime = orbit.GetNextAnglePassTime(currentTime, argDestination);
            }
            else if(pathType == PathType.Escape)
            {
                double theCosinus = (slr / planet.SOI - 1.0) / ecc;

                // to fix for numerical imprecisions...
                if (theCosinus > 1.0) theCosinus = 1.0;
                if (theCosinus < -1.0) theCosinus = -1.0;

                // get true anomaly of escape point
                double trueAnomaly_escape = Math.Acos(theCosinus) * orbit.direction;

                // get passage time at said true anomaly
                orbit.orbitEndTime = orbit.GetNextAnglePassTime(currentTime, trueAnomaly_escape + argOfPeriapsis);
            }

            LOG(LOG_LEVEL.DEBUG, "Created orbit: slr = " + orbit.slr + "; ecc = " + orbit.ecc + "; arg of periapsis = " + orbit.arg + "; dir = " + orbit.direction);
            LOG(LOG_LEVEL.DEBUG, "               orbitStartTime = " + orbit.orbitStartTime + "; orbitEndTime = " + orbit.orbitEndTime + "; period = " + orbit.period);
        }

        return orbit;
    }


    private static Orbit generateFinalOrbit() 
    {
        double ecc = Math.Abs(1.0 - lambda);
        double slr = (2.0 - lambda) * R_start;
        double periapsis = slr / (1.0 + ecc);

        double argOfPeriapsis = 0.0; // default value (used for circular orbits)
        double trueAnomaly = 0.0;

        // calculate argument of periapsis
        if (1.0 + ecc > 1.0) // eccentricity != 0.0 - tested that way (instead of ecc > 0.0) so that ecc is considered as zero if it's negligible (< 10^(-16))
        {
            double theCosinus = (slr / R_start - 1.0) / ecc;

            // to fix for numerical imprecisions...
            if (theCosinus > 1.0) theCosinus = 1.0;
            if (theCosinus < -1.0) theCosinus = -1.0;

            trueAnomaly = Math.Acos(theCosinus);

            if (sign_wt * sign_wr < 0.0)
            {
                trueAnomaly = -trueAnomaly;
            }

            argOfPeriapsis = arg_start - trueAnomaly;
        }

        Orbit orbit = Orbit_Utils.CreateOrbit(slr, ecc, argOfPeriapsis, (int)sign_wt, planet, PathType.Eternal, null);

        if (orbit != null)
        {
            // specific energy - needed for the Kepler solver and time determination
            double specificEnergy = planet.mass * (ecc * ecc - 1.0) / slr / 2.0;

            // Set periapsis passage time
            double timeFromPeri = 0.0;
            new KeplerSolver(planet.mass, periapsis, specificEnergy).GetTimeAtPosition(R_start, trueAnomaly, ref timeFromPeri);
            orbit.periapsisPassageTime = currentTime - timeFromPeri * (double)orbit.direction - 10.0 * orbit.period;

            // Set start
            orbit.orbitStartTime = currentTime;

            LOG(LOG_LEVEL.DEBUG, "Created orbit: slr = " + orbit.slr + "; ecc = " + orbit.ecc + "; arg of periapsis = " + orbit.arg + "; dir = " + orbit.direction);
            LOG(LOG_LEVEL.DEBUG, "               orbitStartTime = " + orbit.orbitStartTime + "; period = " + orbit.period);
        }

        return orbit;
    }


    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.RETURN_TRAJECTORY, level, message);
    }
}
