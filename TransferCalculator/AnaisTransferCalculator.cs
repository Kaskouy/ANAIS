using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

using SFS.World;
using SFS.WorldBase;
using SFS.World.Maps;

class AnaisTransferCalculator
{
    
    // FUNCTIONS TO CALCULATE THE TIME RANGE FOR THE RESEARCH OF AN OPTIMAL TRANSFER
    // -----------------------------------------------------------------------------
    private static double getClosestPassageTimeAtArg(Orbit orbit, double time_ref, double arg)
    {
        if (orbit == null) return -1.0;

        double time = orbit.GetNextAnglePassTime(time_ref, arg);

        if(orbit.period > 0.0)
        {
            if(time > time_ref + orbit.period / 2.0)
            {
                time -= orbit.period;
            }
        }

        return time;
    }


    public static bool computeTimeRangeAroundTimeOfArrival(Orbit playerOrbit, Orbit targetOrbit, double timeNow, double arrivalDate, double playerArg, out double startTime, out double endTime)
    {
        const double C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG = 0.001; // radians - approximately 3 minutes of arc; we'll always make sure to never calculate a rectilinear transfer
        const double C_MIN_DELAY_BEYOND_START_TIME = 60.0;

        startTime = targetOrbit.GetLastAnglePassTime(arrivalDate, playerArg + C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG * targetOrbit.direction);
        endTime   = targetOrbit.GetNextAnglePassTime(arrivalDate, playerArg - C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG * targetOrbit.direction);

        if(endTime < timeNow)
        {
            // If end time is in the past (which basically means that the most favorable encounter is past)
            // we search on the next period so that at least we can propose something to the player
            startTime += targetOrbit.period;
            endTime += targetOrbit.period;
        }

        LOG(LOG_LEVEL.DEBUG, "   Calculate startTime at arg = " + (playerArg + C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG * targetOrbit.direction) + "; start time = now + " + (startTime - WorldTime.main.worldTime));
        LOG(LOG_LEVEL.DEBUG, "   Calculate endTime at arg = " + (playerArg - C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG * targetOrbit.direction) + "; end time = now + " + (endTime - WorldTime.main.worldTime));

        double minStartTime = timeNow + C_MIN_DELAY_BEYOND_START_TIME;
        if (startTime < minStartTime) startTime = minStartTime;

        if(endTime > targetOrbit.orbitEndTime) endTime = targetOrbit.orbitEndTime;

        return (startTime < endTime);
    }


    private static bool computeTimeRangeForAnaisTransferCalculation(Orbit playerOrbit, Orbit targetOrbit, double timeNow, out double startTime, out double endTime)
    {
        startTime = endTime = 0.0;

        Location playerLocationNow = playerOrbit.GetLocation(timeNow);
        double playerArg = playerLocationNow.position.AngleRadians;

        // Check if orbits intersect
        double angleA, angleB;
        bool intersect = Orbit_Utils.GetIntersectionAngles(playerOrbit, targetOrbit, out angleA, out angleB);

        if (intersect)
        {
            // The orbits intersect - calculate a time range based on the passage times of player and target at the intersections
            // ------------------------------------------------------------------------------------------------------------------
            // take as reference date the passage time at the next crossing points
            double timeOfPlayerPassageA = playerOrbit.GetNextAnglePassTime(timeNow, angleA);
            double timeOfPlayerPassageB = playerOrbit.GetNextAnglePassTime(timeNow, angleB);

            LOG(LOG_LEVEL.DEBUG, "  computeTimeRangeForAnaisTransferCalculation: Two intersections were found with the target orbit");
            LOG(LOG_LEVEL.DEBUG, "  computeTimeRangeForAnaisTransferCalculation:   Intersection A at argument " + angleA + "; time of next player arrival = " + timeOfPlayerPassageA);
            LOG(LOG_LEVEL.DEBUG, "  computeTimeRangeForAnaisTransferCalculation:   Intersection B at argument " + angleB + "; time of next player arrival = " + timeOfPlayerPassageB);

            double bestTimeOfPassageTarget;

            if((timeOfPlayerPassageA > timeNow) && (timeOfPlayerPassageB > timeNow))
            {
                LOG(LOG_LEVEL.DEBUG, " --- CALCULATION OF START TIME AND END TIME ---");
                LOG(LOG_LEVEL.DEBUG, " Player period = " + playerOrbit.period + "; Target period = " + targetOrbit.period);
                LOG(LOG_LEVEL.DEBUG, " Time of player passage at intersection A = now + " + (timeOfPlayerPassageA - timeNow));
                LOG(LOG_LEVEL.DEBUG, " Time of player passage at intersection B = now + " + (timeOfPlayerPassageB - timeNow));

                // Get the time of passage of target for each intersection (the closest to the time of passage of player)
                double timeOfTargetPassageA = getClosestPassageTimeAtArg(targetOrbit, timeOfPlayerPassageA, angleA);
                double timeOfTargetPassageB = getClosestPassageTimeAtArg(targetOrbit, timeOfPlayerPassageB, angleB);

                // Calculate the difference of time of passage at each intersection
                double deltaTimeA = timeOfTargetPassageA - timeOfPlayerPassageA;
                double deltaTimeB = timeOfTargetPassageB - timeOfPlayerPassageB;

                LOG(LOG_LEVEL.DEBUG, " Time of closest target passage in A = now + " + (timeOfTargetPassageA - timeNow) + "; deltaTimeA = " + deltaTimeA);
                LOG(LOG_LEVEL.DEBUG, " Time of closest target passage in B = now + " + (timeOfTargetPassageB - timeNow) + "; deltaTimeB = " + deltaTimeB);

                if (Math.Abs(deltaTimeA) < Math.Abs(deltaTimeB))
                {
                    bestTimeOfPassageTarget = timeOfTargetPassageA;
                }
                else
                {
                    bestTimeOfPassageTarget = timeOfTargetPassageB;
                }
            }
            else if(timeOfPlayerPassageA > timeNow)
            {
                // Passage of player at intersection B is in the past (if the orbit is hyperbolic, so it only happens once) --> search around passage at intersection A
                bestTimeOfPassageTarget = getClosestPassageTimeAtArg(targetOrbit, timeOfPlayerPassageA, angleA);
            }
            else if(timeOfPlayerPassageB > timeNow)
            {
                // Passage of player at intersection A is in the past (if the orbit is hyperbolic, so it only happens once) --> search around passage at intersection B
                bestTimeOfPassageTarget = getClosestPassageTimeAtArg(targetOrbit, timeOfPlayerPassageB, angleB);
            }
            else
            {
                // both passage times are in the past (hyperbolic orbit)
                return false;
            }

            LOG(LOG_LEVEL.DEBUG, "Best passage time of target = " + bestTimeOfPassageTarget);

            // Compute a time range around the best time of passage
            return computeTimeRangeAroundTimeOfArrival(playerOrbit, targetOrbit, timeNow, bestTimeOfPassageTarget, playerArg, out startTime, out endTime);
        }
        else
        {
            LOG(LOG_LEVEL.DEBUG, "  computeTimeRangeForAnaisTransferCalculation: Simulate a Hohmann transfer");

            // Simulate an Hohmann transfer and search an optimal transfer base on that one
            // ----------------------------------------------------------------------------
            Orbit_Utils.CalculateHohmannTransfer(playerOrbit, targetOrbit, playerLocationNow, out Orbit hohmannTransfer, out double arrivalArg);

            if (hohmannTransfer == null)
            {
                LOG(LOG_LEVEL.DEBUG, "ERROR: Hohmann transfer could not be calculated");
                return false;
            }

            // Hohmann transfer valid; use the arrival time as target date
            double playerArrivalDate = hohmannTransfer.orbitEndTime;
            LOG(LOG_LEVEL.DEBUG, "  calculateReferenceDateAndPos: A Hohmann transfer was found! RefArg = " + arrivalArg + "; time of player arrival at refArg = " + playerArrivalDate);

            // Calculate the closest time at which the target will pass at ref point
            double timeOfTargetPassageAtRefPoint = getClosestPassageTimeAtArg(targetOrbit, playerArrivalDate, arrivalArg);

            LOG(LOG_LEVEL.DEBUG, "  Passage time of player at ref point = " + playerArrivalDate);
            LOG(LOG_LEVEL.DEBUG, "  Closest passage time of target at ref point = " + timeOfTargetPassageAtRefPoint + " (orbital period = " + targetOrbit.period + ")");

            // Compute a time range around the best time of passage
            return computeTimeRangeAroundTimeOfArrival(playerOrbit, targetOrbit, timeNow, timeOfTargetPassageAtRefPoint, playerArg, out startTime, out endTime);
        }
    }


    // FUNCTION THAT COMPUTES A LIST OF TRANSFERS EVENLY DISPATCHED IN TERMS OF ARRIVAL TIME
    // -------------------------------------------------------------------------------------
    public static List<AnaisTransfer> calculateTransferList(Orbit playerOrbit, Orbit targetOrbit, Planet targetPlanet, double start_time, double arrival_time_start_window, double arrival_time_end_window)
    {
        const uint C_TIME_INTERVALS_PER_PERIOD = 12; // On a full period of target object, the time window would be splitted in that many periods
        
        List<AnaisTransfer> anaisTransferList = new List<AnaisTransfer>();

        double time_window_length = arrival_time_end_window - arrival_time_start_window;
        double time_window_length_in_fraction_of_period = time_window_length / targetOrbit.period;

        uint nbTimeIntervals = Math.Max(2, (uint)(C_TIME_INTERVALS_PER_PERIOD * time_window_length_in_fraction_of_period));

        double time_interval = time_window_length / nbTimeIntervals;

        LOG(LOG_LEVEL.DEBUG, "    Calculating the transfer list");

        for(int i = 0; i < nbTimeIntervals + 1; i++)
        {
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet);

            // calculate a transfer for that arrival time
            double arrival_time = arrival_time_start_window + i * time_interval;
            anaisTransfer.calculateTransfer(start_time, arrival_time);

            if(anaisTransfer.transferOrbit != null)
            {
                // Add it to the list
                LOG(LOG_LEVEL.DEBUG, "      Transfer calculated; arrival time = " + anaisTransfer.transferOrbit.orbitEndTime + "; total ΔV = " + anaisTransfer.total_dv);
                anaisTransferList.Add(anaisTransfer);
            }
            else
            {
                LOG(LOG_LEVEL.DEBUG, "    transfer is null");
            }
        }

        // Insert a transfer close to the first bound
        // ------------------------------------------
        {
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet);

            double arrival_time = arrival_time_start_window + time_interval / 8.0;
            anaisTransfer.calculateTransfer(start_time, arrival_time);

            if (anaisTransfer.transferOrbit != null)
            {
                // Add it to the list in second position
                LOG(LOG_LEVEL.DEBUG, "      Transfer calculated; arrival time = " + anaisTransfer.transferOrbit.orbitEndTime + "; total ΔV = " + anaisTransfer.total_dv);
                anaisTransferList.Insert(1, anaisTransfer);
            }
        }

        // Insert a transfer close to the last bound
        // ------------------------------------------
        {
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet);

            double arrival_time = arrival_time_end_window - time_interval / 8.0;
            anaisTransfer.calculateTransfer(start_time, arrival_time);

            if (anaisTransfer.transferOrbit != null)
            {
                // Add it to the list before last position
                LOG(LOG_LEVEL.DEBUG, "      Transfer calculated; arrival time = " + anaisTransfer.transferOrbit.orbitEndTime + "; total ΔV = " + anaisTransfer.total_dv);
                anaisTransferList.Insert(anaisTransferList.Count - 1, anaisTransfer);
            }
        }

        return anaisTransferList;
    }


    // CALCULATES THE OPTIMAL TRANSFER WITHIN THE SPECIFIED TIME RANGE
    // ---------------------------------------------------------------
    public static AnaisTransfer calculateAnaisTransferForArrivalTimeRange(Orbit playerOrbit, Orbit targetOrbit, Planet targetPlanet, double departureTime, double start_arrivalTime, double end_arrivalTime)
    {
        // Controls
        // --------
        if (!(start_arrivalTime < end_arrivalTime)) return null; // time range must be defined consistently

        if (!(departureTime < start_arrivalTime)) return null; // arrival must be past departure no matter what happens


        // Calculate several transfers with arrival times evenly dispatched on the time range specified
        // --------------------------------------------------------------------------------------------
        LOG(LOG_LEVEL.DEBUG, "  Building transfer list");

        List<AnaisTransfer> anaisTransferList = calculateTransferList(playerOrbit, targetOrbit, targetPlanet, departureTime, start_arrivalTime, end_arrivalTime);

        if (anaisTransferList.Count == 0) return null; // safety check

        foreach (AnaisTransfer anaisTransfer in anaisTransferList)
        {
            if (anaisTransfer.transferOrbit == null)
            {
                LOG(LOG_LEVEL.DEBUG, "  Transfer is null - give up!");
                return null; // give-up
            }
        }

        LOG(LOG_LEVEL.DEBUG, "  Transfer list calculated\n");

        // Search for a minimum
        // --------------------
        double dvMini = anaisTransferList[0].total_dv;
        int i_dvMini = 0;

        for (int i = 1; i < anaisTransferList.Count; i++)
        {
            LOG(LOG_LEVEL.DEBUG, "AnaisTransfer[" + i + "] : t = " + anaisTransferList[i].arrivalTime + "; dv = " + anaisTransferList[i].total_dv);

            if (anaisTransferList[i].total_dv < dvMini)
            {
                dvMini = anaisTransferList[i].total_dv;
                i_dvMini = i;
            }
        }

        if ((i_dvMini == 0) || (i_dvMini == anaisTransferList.Count - 1))
        {
            LOG(LOG_LEVEL.DEBUG, "  No local minimum found - Give up");
            return null; // No minimum located
        }
        else
        {
            LOG(LOG_LEVEL.DEBUG, "  Minimum found at position " + i_dvMini + "; total ΔV = " + anaisTransferList[i_dvMini].total_dv);
        }

        // Search for a more accurate minimum around the located minimum
        // -------------------------------------------------------------
        NumericalMinimumCalculator minCalculator = NumericalMinimumCalculator.createInstance(anaisTransferList[i_dvMini - 1].arrivalTime, anaisTransferList[i_dvMini - 1].total_dv,
                                                                                             anaisTransferList[i_dvMini].arrivalTime    , anaisTransferList[i_dvMini].total_dv,
                                                                                             anaisTransferList[i_dvMini + 1].arrivalTime, anaisTransferList[i_dvMini + 1].total_dv);

        if(minCalculator == null)
        {
            LOG(LOG_LEVEL.ERROR, "  Minimum calculator could not be instantiated - no transfer returned");
            return null;
        }

        AnaisTransfer bestAnaisTransfer = anaisTransferList[i_dvMini];
        int nbIterations = 0;
        bool getOutNow = false;

        LOG(LOG_LEVEL.DEBUG, "  Now perform iterations to locate precisely the minimum!");

        do
        {
            nbIterations++;

            // Calculate a new arrival time
            double newArrivalTime = minCalculator.calculateNewGuess();

            // Calculate the transfer at new date
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet);
            anaisTransfer.calculateTransfer(departureTime, newArrivalTime);

            if (anaisTransfer.transferOrbit == null)
            {
                LOG(LOG_LEVEL.DEBUG, "  Error - Anais transfer calculation failed for t = " + newArrivalTime + "; Give-up") ;
                return null;
            }
            
            // memorize the best transfer
            if (anaisTransfer.total_dv < bestAnaisTransfer.total_dv) bestAnaisTransfer = anaisTransfer;

            // Insert appropriately the new transfer in the calculator data
            bool isOK = minCalculator.insertNewPoint(anaisTransfer.arrivalTime, anaisTransfer.total_dv);

            if (!isOK)
            {
                LOG(LOG_LEVEL.DEBUG, "  Error - Failed to insert new transfer in the calculator");
                return null;
            }

            // Check if the latest transfer found is good enough
            const double C_MAX_DV_DIFFERENCE = 0.1;
            bool conditionOverSpeed = (minCalculator.getMaxYdifference() < C_MAX_DV_DIFFERENCE);

            double C_DELTA_TIME_TOLERANCE = targetOrbit.period / 3600.0; // 21600 = 360 * 60 : 1/21600 corresponds to one minute of angle on a circular orbit
            bool conditionOverTime = (minCalculator.getXwidth() < C_DELTA_TIME_TOLERANCE);

            getOutNow = (conditionOverSpeed || conditionOverTime); 

        } while ((nbIterations < 50) && !getOutNow);

        if (getOutNow)
        {
            LOG(LOG_LEVEL.DEBUG, "  Calculation SUCCESSFUL (" + nbIterations + " iterations) - return found transfer");
            return bestAnaisTransfer;
        }
        else
        {
            LOG(LOG_LEVEL.DEBUG, "  Calculation FAILED by excessive number of iterations");
            return null;
        }
    }


    // MAIN FUNCTION
    // -------------
    public static AnaisTransfer calculateTransferToTarget(Orbit playerOrbit, Orbit targetOrbit, Planet targetPlanet, double timeNow)
    {
        LOG(LOG_LEVEL.DEBUG, "-------------------------------------------");
        LOG(LOG_LEVEL.DEBUG, "---   NEW ANAIS TRANSFER CALCULATION    ---");
        LOG(LOG_LEVEL.DEBUG, "-------------------------------------------");

        // Deal with the return to mother planet case first
        if ((targetPlanet != null) && ((playerOrbit.Planet == targetPlanet) || (playerOrbit.Planet.parentBody == targetPlanet)))
        {
            // Return to mother planet
            AnaisTransfer anaisReturnTransfer = new AnaisTransfer(playerOrbit, null, targetPlanet); // No destination orbit in this case

            anaisReturnTransfer.calculateReturnToPlanet(timeNow);

            return anaisReturnTransfer;
        }

        // We can't deal with a non periodic orbit for now.
        if(!(targetOrbit.period > 0))
        {
            return null;
        }

        Orbit startOrbit, endOrbit;

        if(playerOrbit.Planet == targetOrbit.Planet)
        {
            // player and target orbit the same body
            startOrbit = playerOrbit;
            endOrbit = targetOrbit;
        }
        else if((playerOrbit.Planet.parentBody != null) && (playerOrbit.Planet.parentBody == targetOrbit.Planet))
        {
            // player is in the SOI of the child planet of body orbited by target body (interplanetary transfer)
            startOrbit = playerOrbit.Planet.orbit;
            endOrbit = targetOrbit;
        }
        else
        {
            // case not treated
            return null;
        }


        double startTime, endTime;
        bool timeRangeValid = computeTimeRangeForAnaisTransferCalculation(startOrbit, endOrbit, timeNow, out startTime, out endTime);

        if (!timeRangeValid)
        {
            LOG(LOG_LEVEL.DEBUG, "ERROR: Time range is invalid - exit");
            return null;
        }

        // Calculate the optimal transfer on that time range
        // -------------------------------------------------
        AnaisTransfer newAnaisTransfer = calculateAnaisTransferForArrivalTimeRange(playerOrbit, targetOrbit, targetPlanet, timeNow, startTime, endTime);
        return newAnaisTransfer;
    }


    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.ANAIS_TRANSFER_CALCULATOR, level, message);
    }
}
