using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

using SFS.World;
using SFS.WorldBase;

class EjectionTrajectoryCalculator
{

    private static bool calculateSolution(double A, double B, double C, double D, double initialGuess, out double x)
    {
        x = initialGuess;
        bool ok = false;

        // exit if initial guess is out of bounds (stop at 0.9999 instead of 1 because the derivative couldn't be calculated for -1 and 1)
        if ((initialGuess < -0.9999) || (initialGuess > 0.9999))
        {
            LOG(LOG_LEVEL.DEBUG, "calculateSolution: Error, x_ini out of bounds : " + initialGuess);
            return false;
        }

        uint nbIterations = 0;

        while((nbIterations < 10) && !ok)
        {
            nbIterations++;

            double x2 = x * x;
            double theSqrt = Math.Sqrt(1 - x2);
            double value = A*x2 + B*x + C + D*theSqrt;
            double derivative = 2.0 * A * x + B - D*x/theSqrt;

            // Can't apply Newton method with derivative equal to 0
            if (Math.Abs(derivative) < 0.0001)
            {
                LOG(LOG_LEVEL.DEBUG, "calculateSolution: Error, derivative is null: " + derivative);
                return false;
            }

            // Apply Newton method to get new estimation
            double new_x = x - value / derivative;
            LOG(LOG_LEVEL.DEBUG, "calculateSolution: iteration " + nbIterations + ": new_x = " + new_x);

            // Check that the new value is valid
            if ((new_x < -0.9999) || (new_x > 0.9999))
            {
                LOG(LOG_LEVEL.DEBUG, "calculateSolution: Error, x out of bounds: " + new_x);
                return false;
            }

            // set new values before looping
            if (Math.Abs(new_x - x) < 0.000001)
            {
                LOG(LOG_LEVEL.INFO, "calculateSolution: x = " + x + " is accepted as a solution");
                ok = true;
            }
            x = new_x;
        }

        if(!ok)
        {
            LOG(LOG_LEVEL.INFO, "calculateSolution: Error, the solution could not be calculated, too many iterations");
        }

        return ok;
    }

    private static Orbit CalculateOrbit(Planet planet, double alpha, double heading, double sign, double x, double startRadius, double startArg, double startTime)
    {
        double x2 = x * x;
        double p = planet.SOI * alpha * x2;
        double ecc = Math.Sqrt(1.0 + alpha * (alpha - 2.0) * x2);
        double phi = heading - Math.Atan2((alpha - 1.0) * x, -sign * Math.Sqrt(1 - x2));
        int direction = Math.Sign(x);
        double exitAngle = phi + Math.Sign(sign * x) * Math.Acos((p / planet.SOI - 1.0) / ecc);

        Orbit orbit = Orbit_Utils.CreateOrbit(p, ecc, phi, direction, planet, PathType.Escape, planet.parentBody);

        if (orbit != null)
        {
            // Set periapsis passage time
            double timeFromPeri = 0.0;
            new KeplerSolver(planet.mass, orbit.periapsis, orbit.getSpecificEnergy()).GetTimeAtPosition(startRadius, startArg - orbit.arg, ref timeFromPeri);
            orbit.periapsisPassageTime = startTime - timeFromPeri * (double)direction - 10.0 * orbit.period;

            // Set start and end time
            orbit.orbitStartTime = startTime;
            orbit.orbitEndTime = orbit.GetNextAnglePassTime(startTime, exitAngle);

            LOG(LOG_LEVEL.INFO, "CalculateOrbit: SUCCESS: p = " + orbit.slr + ", ecc = " + orbit.ecc + ", phi = " + orbit.arg + ", peri = " + orbit.periapsis + ", peri pass time = " + orbit.periapsisPassageTime + ", start time = " + orbit.orbitStartTime + ", end time = " + orbit.orbitEndTime);

            return orbit;
        }
        else
        {
            // Orbit is degenerated -> give up
            LOG(LOG_LEVEL.INFO, "CalculateOrbit: FAILURE: orbit is degenerated: p = " + p);
            return null;
        }
    }

    public static Orbit calculateEjectionTrajectories(Planet planet, Location startLocation, double startTime, Double2 exitVelocity, int direction, bool exit = true)
    {
        // List of variables
        double startRadius = startLocation.position.magnitude;
        double startArgument = startLocation.position.AngleRadians;
        double escapeVelocity = exitVelocity.magnitude;
        double heading = exitVelocity.AngleRadians;

        double delta = startArgument - heading;

        double sign = 1.0;
        if (exit == false) sign = -1.0;

        double alpha = planet.SOI * escapeVelocity * escapeVelocity / planet.mass;

        // List of coefficients: equation to solve is A.x^2 + B.x + C + D * sqrt(1 - x^2) = 0.0
        double A = planet.SOI * alpha / startRadius;
        double B = (alpha - 1.0) * Math.Sin(delta);
        double C = -1.0;
        double D = sign * Math.Cos(delta);

        LOG(LOG_LEVEL.INFO, "  EJECTION TRAJECTORY : New calculation");
        LOG(LOG_LEVEL.DEBUG, "A = " + A + ", B = " + B + ", C = " + C + ", D = " + D);

        // Solve A*x^2 + B*x + C2 = 0 with C2 = max possible value of C + D * sqrt(1 - x2)
        double C2 = C + Math.Max(0.0, D);

        double discri = B * B - 4.0 * A * C2;
        LOG(LOG_LEVEL.DEBUG, "C2 = " + C2 + ", Delta = " + discri);

        if(discri > 0.0)
        {
            // make sure direction is 1 or -1
            int dir = 1;
            if(direction < 0) dir = -1;
            
            // Choose the most likely solution depending on direction
            double x_ini = (-B + dir * Math.Sqrt(discri)) / (2.0 * A);

            LOG(LOG_LEVEL.DEBUG, "x_ini = " + x_ini);

            if (Math.Sign(x_ini) != dir)
            {
                LOG(LOG_LEVEL.INFO, "ERROR; x_ini has the wrong sign");
                return null; // Something is fucked
            }

            bool ok = calculateSolution(A, B, C, D, x_ini, out double x);

            if(ok)
            {
                return CalculateOrbit(planet, alpha, heading, sign, x, startRadius, startArgument, startTime);
            }
        }
        else
        {
            LOG(LOG_LEVEL.INFO, "ERROR! Delta < 0");
        }

        return null;
    }

    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.EJECTION_TRAJECTORY, level, message);
    }
}
