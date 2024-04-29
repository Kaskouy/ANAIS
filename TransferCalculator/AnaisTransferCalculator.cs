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
    // C_MINIMUM_PERIOD_FOR_ACCURATE_CALCULATION: If the player orbit's period is lower than this, ejection trajectory calculation will be skipped for interplanetary transfers
    const double C_MINIMUM_PERIOD_FOR_ACCURATE_CALCULATION = 0.6; // 0.6 seconds - if player orbit's period is lower than this (in real time - factoring time warp factor), only the main transfer will be shown in case of interplanetary transfer


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


    public static bool computeTimeRangeAroundTimeOfArrival(Orbit playerOrbit, Orbit targetOrbit, double timeNow, double arrivalDate, double playerArg, double minimalStartTime, out double startTime, out double endTime)
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

        LOG(LOG_LEVEL.DEBUG, "   Calculate startTime at arg = " + (playerArg + C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG * targetOrbit.direction) + "; start time = now + " + (startTime - timeNow));
        LOG(LOG_LEVEL.DEBUG, "   Calculate endTime at arg = " + (playerArg - C_MIN_DIFF_ANGLE_WITH_PLAYER_ARG * targetOrbit.direction) + "; end time = now + " + (endTime - timeNow));

        // Make sure the start time isn't below a minimal value to avoid excessively energetic transfers
        if(startTime < minimalStartTime) startTime = minimalStartTime;

        // Ultimate security over start time
        double ultimateMinStartTime = timeNow + C_MIN_DELAY_BEYOND_START_TIME;
        if (startTime < ultimateMinStartTime) startTime = ultimateMinStartTime;

        if(endTime > targetOrbit.orbitEndTime) endTime = targetOrbit.orbitEndTime;

        LOG(LOG_LEVEL.DEBUG, "   Final search window: startTime = now + " + (startTime - timeNow) + " - endTime = now + " + (endTime - timeNow));

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
            double minimumTransferDelay = 0.0;

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
                    minimumTransferDelay = (timeOfPlayerPassageA - timeNow) / 4.0;
                }
                else
                {
                    bestTimeOfPassageTarget = timeOfTargetPassageB;
                    minimumTransferDelay = (timeOfPlayerPassageB - timeNow) / 4.0;
                }
            }
            else if(timeOfPlayerPassageA > timeNow)
            {
                // Passage of player at intersection B is in the past (if the orbit is hyperbolic, so it only happens once) --> search around passage at intersection A
                bestTimeOfPassageTarget = getClosestPassageTimeAtArg(targetOrbit, timeOfPlayerPassageA, angleA);
                minimumTransferDelay = (timeOfPlayerPassageA - timeNow) / 4.0;
            }
            else if(timeOfPlayerPassageB > timeNow)
            {
                // Passage of player at intersection A is in the past (if the orbit is hyperbolic, so it only happens once) --> search around passage at intersection B
                bestTimeOfPassageTarget = getClosestPassageTimeAtArg(targetOrbit, timeOfPlayerPassageB, angleB);
                minimumTransferDelay = (timeOfPlayerPassageB - timeNow) / 4.0;
            }
            else
            {
                // both passage times are in the past (hyperbolic orbit)
                return false;
            }

            LOG(LOG_LEVEL.DEBUG, "Best passage time of target = " + bestTimeOfPassageTarget);

            // Calculate the minimal start time (to avoid excessively energetic transfers that are difficult to calculate and aren't viable anyway)
            double minimal_start_time = timeNow + minimumTransferDelay;

            // Compute a time range around the best time of passage
            return computeTimeRangeAroundTimeOfArrival(playerOrbit, targetOrbit, timeNow, bestTimeOfPassageTarget, playerArg, minimal_start_time, out startTime, out endTime);
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

            // Calculate the minimal start time (to avoid excessively energetic transfers that are difficult to calculate and aren't viable anyway)
            double minimal_start_time = timeNow + (hohmannTransfer.orbitEndTime - timeNow) / 4.0;

            // Compute a time range around the best time of passage
            return computeTimeRangeAroundTimeOfArrival(playerOrbit, targetOrbit, timeNow, timeOfTargetPassageAtRefPoint, playerArg, minimal_start_time, out startTime, out endTime);
        }
    }

    public static double GetMaxTransferStartingVelocity(Orbit playerOrbit, Orbit targetOrbit, double timeNow)
    {
        // The maximum velocity will be evaluated to that many times orbital velocity
        const double C_MAX_VELOCITY_COEFFICIENT = 2.0;
        
        double maxVelocity = 0.0;

        if(playerOrbit.Planet == targetOrbit.Planet)
        {
            // local transfer
            Location location = playerOrbit.GetLocation(timeNow);
            maxVelocity = C_MAX_VELOCITY_COEFFICIENT * Math.Sqrt(playerOrbit.Planet.mass / location.position.magnitude);
        }
        else if(playerOrbit.Planet.parentBody == targetOrbit.Planet)
        {
            // interplanetary transfer
            Location location = playerOrbit.Planet.orbit.GetLocation(timeNow);
            maxVelocity = C_MAX_VELOCITY_COEFFICIENT * Math.Sqrt(playerOrbit.Planet.parentBody.mass / location.position.magnitude);
        }
        else
        {
            LOG(LOG_LEVEL.ERROR, "GetMaxTransferStartingVelocity: unknown transfer configuration");
        }

        LOG(LOG_LEVEL.DEBUG, "    Maximum starting transfer velocity evaluated to: " + maxVelocity);

        return maxVelocity;
    }

    // FUNCTION THAT COMPUTES A LIST OF TRANSFERS EVENLY DISPATCHED IN TERMS OF ARRIVAL TIME
    // -------------------------------------------------------------------------------------
    public static List<AnaisTransfer> calculateTransferList(Orbit playerOrbit, Orbit targetOrbit, Planet targetPlanet, double targetAltitude, double start_time, double arrival_time_start_window, double arrival_time_end_window, ANAIS_Settings.E_TRANSFER_TYPE transferType)
    {
        const uint C_TIME_INTERVALS_PER_PERIOD = 12; // On a full period of target object, the time window would be splitted in that many periods
        
        List<AnaisTransfer> anaisTransferList = new List<AnaisTransfer>();

        double time_window_length = arrival_time_end_window - arrival_time_start_window;
        double time_window_length_in_fraction_of_period = time_window_length / targetOrbit.period;

        uint nbTimeIntervals = Math.Max(2, (uint)(C_TIME_INTERVALS_PER_PERIOD * time_window_length_in_fraction_of_period));

        double time_interval = time_window_length / nbTimeIntervals;

        double maxTransferVelocity = GetMaxTransferStartingVelocity(playerOrbit, targetOrbit, start_time);
        bool maxVelocityExcessed = false; // to allow to enter the loop

        LOG(LOG_LEVEL.DEBUG, "    Calculating the transfer list - " + (int)(nbTimeIntervals+1) + " transfers expected");

        for(int i = (int)nbTimeIntervals; (i >= 0) && !maxVelocityExcessed; i--) // loop from the end to be able to stop if the transfers at the sooner dates happens to be over-energetic
        {
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet, targetAltitude, transferType);

            // calculate a transfer for that arrival time
            double arrival_time = arrival_time_start_window + i * time_interval;
            anaisTransfer.calculateTransfer(start_time, arrival_time);

            if(anaisTransfer.transferOrbit != null)
            {
                // Add it to the list
                LOG(LOG_LEVEL.DEBUG, "      Transfer calculated; arrival time = " + anaisTransfer.transferOrbit.orbitEndTime + "; total ΔV = " + anaisTransfer.total_dv);
                anaisTransferList.Add(anaisTransfer);

                // evaluate velocity - to make it stop to search transfers if we are already too energetic!
                double velocity = anaisTransfer.transferOrbit.GetLocation(anaisTransfer.transferOrbit.orbitStartTime).velocity.magnitude;
                if(velocity > maxTransferVelocity)
                {
                    maxVelocityExcessed = true;
                    LOG(LOG_LEVEL.DEBUG, "    Velocity is too high, stop searching from now on.");
                }
            }
            else
            {
                LOG(LOG_LEVEL.DEBUG, "    transfer is null - arrival time = now + " + (arrival_time - start_time));
            }
        }

        // Reverse the list since it's been calculated from the end to the starting date
        anaisTransferList.Reverse();

        // Insert a transfer close to the last bound
        // ------------------------------------------
        if(anaisTransferList.Count > 1)
        {
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet, targetAltitude, transferType);

            // the time interval between the 2 last points
            double lastTimeInterval = anaisTransferList[anaisTransferList.Count - 1].arrivalTime - anaisTransferList[anaisTransferList.Count - 2].arrivalTime;

            // calculate a transfer for a time value right before the last one
            double arrival_time = anaisTransferList[anaisTransferList.Count - 1].arrivalTime - lastTimeInterval / 8.0;
            anaisTransfer.calculateTransfer(start_time, arrival_time);

            if (anaisTransfer.transferOrbit != null)
            {
                // Add it to the list before last position
                LOG(LOG_LEVEL.DEBUG, "      Additional transfer calculated; arrival time = " + anaisTransfer.transferOrbit.orbitEndTime + "; total ΔV = " + anaisTransfer.total_dv);
                anaisTransferList.Insert(anaisTransferList.Count - 1, anaisTransfer);
            }
        }

        // Insert a transfer close to the first bound
        // ------------------------------------------
        if(anaisTransferList.Count > 1)
        {
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet, targetAltitude, transferType);

            // the time interval between the 2 first points
            double firstTimeInterval = anaisTransferList[1].arrivalTime - anaisTransferList[0].arrivalTime;

            double arrival_time = anaisTransferList[0].arrivalTime + firstTimeInterval / 8.0;
            anaisTransfer.calculateTransfer(start_time, arrival_time);

            if (anaisTransfer.transferOrbit != null)
            {
                // Add it to the list in second position
                LOG(LOG_LEVEL.DEBUG, "      Additional transfer calculated; arrival time = " + anaisTransfer.transferOrbit.orbitEndTime + "; total ΔV = " + anaisTransfer.total_dv);
                anaisTransferList.Insert(1, anaisTransfer);
            }
        }

        return anaisTransferList;
    }


    // CALCULATES THE OPTIMAL TRANSFER WITHIN THE SPECIFIED TIME RANGE
    // ---------------------------------------------------------------
    public static AnaisTransfer calculateAnaisTransferForArrivalTimeRange(Orbit playerOrbit, Orbit targetOrbit, Planet targetPlanet, double targetAltitude, double departureTime, double start_arrivalTime, double end_arrivalTime, ANAIS_Settings.E_TRANSFER_TYPE transferType, double timeWarpFactor)
    {
        // Controls
        // --------
        if (!(start_arrivalTime < end_arrivalTime)) return null; // time range must be defined consistently

        if (!(departureTime < start_arrivalTime)) return null; // arrival must be past departure no matter what happens


        // Calculate several transfers with arrival times evenly dispatched on the time range specified
        // --------------------------------------------------------------------------------------------
        LOG(LOG_LEVEL.DEBUG, "  Building transfer list");

        List<AnaisTransfer> anaisTransferList = calculateTransferList(playerOrbit, targetOrbit, targetPlanet, targetAltitude, departureTime, start_arrivalTime, end_arrivalTime, transferType);

        if (anaisTransferList.Count < 3)
        {
            LOG(LOG_LEVEL.DEBUG, "  Less than 3 transfers were calculated - give up!");
            return null;
        }

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
        double dvMini = double.PositiveInfinity;
        int i_dvMini = 0;

        for (int i = 0; i < anaisTransferList.Count; i++)
        {
            LOG(LOG_LEVEL.DEBUG, "AnaisTransfer[" + i + "] : t = " + anaisTransferList[i].arrivalTime + "; dv = " + anaisTransferList[i].GetDeltaV_valueToMinimize(transferType));

            if (anaisTransferList[i].GetDeltaV_valueToMinimize(transferType) < dvMini)
            {
                dvMini = anaisTransferList[i].GetDeltaV_valueToMinimize(transferType);
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
            LOG(LOG_LEVEL.DEBUG, "  Minimum found at position " + i_dvMini + "; target ΔV = " + anaisTransferList[i_dvMini].GetDeltaV_valueToMinimize(transferType));
        }

        // In the case of an interplanetary transfer, evaluate accurately the serious candidates
        // -------------------------------------------------------------------------------------
        bool needToRefineTransfer = playerOrbit.Planet.parentBody == targetOrbit.Planet; // transfer needs to be refined in case of interplanetary transfer only

        // if player orbit is periodic but the period is too low compared to the time-warp factor, then only the main transfer is took into account
        if ((playerOrbit.pathType == PathType.Eternal) && (playerOrbit.period > 0.0) && (playerOrbit.period / timeWarpFactor < C_MINIMUM_PERIOD_FOR_ACCURATE_CALCULATION)) needToRefineTransfer = false;

        if (needToRefineTransfer)
        {
            // First refine the 3 best transfers
            // transfer 1
            bool transferSuccessfullyRefined = anaisTransferList[i_dvMini - 1].refineTransferCalculation();

            if(!transferSuccessfullyRefined)
            {
                LOG(LOG_LEVEL.DEBUG, "  Failed to refine interplanetary transfer calculation - Give up");
                return null;
            }

            double deltaV_transfer1 = anaisTransferList[i_dvMini - 1].GetDeltaV_valueToMinimize(transferType);

            // transfer 2
            transferSuccessfullyRefined = anaisTransferList[i_dvMini].refineTransferCalculation();

            if (!transferSuccessfullyRefined)
            {
                LOG(LOG_LEVEL.DEBUG, "  Failed to refine interplanetary transfer calculation - Give up");
                return null;
            }

            double deltaV_transfer2 = anaisTransferList[i_dvMini].GetDeltaV_valueToMinimize(transferType);

            // transfer 3
            transferSuccessfullyRefined = anaisTransferList[i_dvMini + 1].refineTransferCalculation();

            if (!transferSuccessfullyRefined)
            {
                LOG(LOG_LEVEL.DEBUG, "  Failed to refine interplanetary transfer calculation - Give up");
                return null;
            }

            double deltaV_transfer3 = anaisTransferList[i_dvMini + 1].GetDeltaV_valueToMinimize(transferType);

            while ((i_dvMini > 0) && (i_dvMini < anaisTransferList.Count - 1) &&                        // i_dvMini is in range (must have a predecessor and a successor)...
                   ((deltaV_transfer1 < deltaV_transfer2) || (deltaV_transfer3 < deltaV_transfer2)))   // ... and deltaV_transfer2 is not the minimum value
            {
                // case when the third transfer is the cheapest
                if (deltaV_transfer3 < deltaV_transfer2)
                {
                    // transfer 3 is the cheapest - evaluate the following one to bound the minimum
                    i_dvMini++;

                    if (i_dvMini < anaisTransferList.Count - 1)
                    {
                        // shift all values
                        deltaV_transfer1 = deltaV_transfer2;
                        deltaV_transfer2 = deltaV_transfer3;

                        // evaluate accurately the new "transfer 3" 
                        transferSuccessfullyRefined = anaisTransferList[i_dvMini + 1].refineTransferCalculation();

                        if (!transferSuccessfullyRefined)
                        {
                            LOG(LOG_LEVEL.DEBUG, "  Failed to refine interplanetary transfer calculation - Give up");
                            return null;
                        }

                        deltaV_transfer3 = anaisTransferList[i_dvMini + 1].GetDeltaV_valueToMinimize(transferType);
                    }
                }

                // case when the first transfer is the cheapest
                else if (deltaV_transfer1 < deltaV_transfer2)
                {
                    // transfer 1 is the cheapest - evaluate the previous one to bound the minimum
                    i_dvMini--;

                    if(i_dvMini > 0)
                    {
                        // shift all values
                        deltaV_transfer3 = deltaV_transfer2;
                        deltaV_transfer2 = deltaV_transfer1;

                        // evaluate accurately the new "transfer 1" 
                        transferSuccessfullyRefined = anaisTransferList[i_dvMini - 1].refineTransferCalculation();

                        if (!transferSuccessfullyRefined)
                        {
                            LOG(LOG_LEVEL.DEBUG, "  Failed to refine interplanetary transfer calculation - Give up");
                            return null;
                        }

                        deltaV_transfer1 = anaisTransferList[i_dvMini - 1].GetDeltaV_valueToMinimize(transferType);
                    }
                }
            } // end while i_dvMini in range and transfer 2 is not the minimal transfer
        }

        // Search for a more accurate minimum around the located minimum
        // -------------------------------------------------------------
        NumericalMinimumCalculator minCalculator = NumericalMinimumCalculator.createInstance(anaisTransferList[i_dvMini - 1].arrivalTime, anaisTransferList[i_dvMini - 1].GetDeltaV_valueToMinimize(transferType),
                                                                                             anaisTransferList[i_dvMini].arrivalTime    , anaisTransferList[i_dvMini].GetDeltaV_valueToMinimize(transferType),
                                                                                             anaisTransferList[i_dvMini + 1].arrivalTime, anaisTransferList[i_dvMini + 1].GetDeltaV_valueToMinimize(transferType));

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
            AnaisTransfer anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, targetPlanet, targetAltitude, transferType);
            anaisTransfer.calculateTransfer(departureTime, newArrivalTime);

            if (anaisTransfer.transferOrbit == null)
            {
                LOG(LOG_LEVEL.DEBUG, "  Error - Anais transfer calculation failed for t = " + newArrivalTime + "; Give-up") ;
                return null;
            }

            if(needToRefineTransfer)
            {
                anaisTransfer.refineTransferCalculation();
            }
            
            // memorize the best transfer
            if (anaisTransfer.GetDeltaV_valueToMinimize(transferType) < bestAnaisTransfer.GetDeltaV_valueToMinimize(transferType)) bestAnaisTransfer = anaisTransfer;

            // Insert appropriately the new transfer in the calculator data
            bool isOK = minCalculator.insertNewPoint(anaisTransfer.arrivalTime, anaisTransfer.GetDeltaV_valueToMinimize(transferType));

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
    public static AnaisTransfer calculateTransferToTarget(Orbit playerOrbit, Orbit targetOrbit, Planet targetPlanet, double targetAltitude, double timeNow, ANAIS_Settings.E_TRANSFER_TYPE transferType, double timeWarpFactor)
    {
        LOG(LOG_LEVEL.DEBUG, "-------------------------------------------");
        LOG(LOG_LEVEL.DEBUG, "---   NEW ANAIS TRANSFER CALCULATION    ---");
        LOG(LOG_LEVEL.DEBUG, "-------------------------------------------");

        // Deal with the return to mother planet case first
        if ((targetPlanet != null) && ((playerOrbit.Planet == targetPlanet) || (playerOrbit.Planet.parentBody == targetPlanet)))
        {
            // Return to mother planet
            AnaisTransfer anaisReturnTransfer = new AnaisTransfer(playerOrbit, null, targetPlanet, targetAltitude, transferType); // No destination orbit in this case

            // if player orbit is periodic but the period is too low compared to the time-warp factor, then only the main transfer is took into account
            bool refineTransfer = true;
            if ((playerOrbit.pathType == PathType.Eternal) && (playerOrbit.period > 0.0) && (playerOrbit.period / timeWarpFactor < C_MINIMUM_PERIOD_FOR_ACCURATE_CALCULATION)) refineTransfer = false;
            
            anaisReturnTransfer.calculateReturnToPlanet(timeNow, refineTransfer);

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
        AnaisTransfer newAnaisTransfer = calculateAnaisTransferForArrivalTimeRange(playerOrbit, targetOrbit, targetPlanet, targetAltitude, timeNow, startTime, endTime, transferType, timeWarpFactor);
        return newAnaisTransfer;
    }


    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.ANAIS_TRANSFER_CALCULATOR, level, message);
    }
}
