using System.Runtime.CompilerServices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using SFS.World;
using SFS.WorldBase;

class AnaisTransfer
{
    // the initial and the target orbit
    public Orbit originOrbit;
    public Orbit destinationOrbit;
    public Planet targetPlanet;

    // all calculated data
    public Orbit transferOrbit;

    public double departureTime;
    public double arrivalTime;

    private Double2 deltaV_start;
    private Double2 deltaV_end;

    public double dv_start;
    public double dv_end;
    public double total_dv;

    public double transfer_efficiency;

    // data for when the transfer starts from a child planet
    bool needsToCalculateExitTrajectory = false;
    public Orbit childTransferOrbit = null;

    public AnaisTransfer(Orbit _originOrbit, Orbit _destinationOrbit, Planet _targetPlanet)
    {
        originOrbit = _originOrbit;
        destinationOrbit = _destinationOrbit;
        targetPlanet = _targetPlanet;

        // init all the shit
        transferOrbit = null;
        departureTime = arrivalTime = 0.0;
        deltaV_start.x = 0.0;
        deltaV_start.y = 0.0;
        deltaV_end.x = 0.0;
        deltaV_end.y = 0.0;
        dv_start = dv_end = total_dv = 0.0;
        transfer_efficiency = 0.0;

        // data for a possible exit trajectory
        needsToCalculateExitTrajectory = false;
        childTransferOrbit = null;
    }

    private bool checkDates(double startTime, double endTime)
    {
        if ((startTime < originOrbit.orbitStartTime) || (endTime > destinationOrbit.orbitEndTime))
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    private bool checkPlanetaryConfiguration()
    {
        if((targetPlanet != null) && ((originOrbit.Planet == targetPlanet) || (originOrbit.Planet.parentBody == targetPlanet)))
        {
            // return to mother planet configuration -> not the good method to call
            return false;
        }
        else if (originOrbit.Planet == destinationOrbit.Planet)
        {
            // origin and destination are in the same SOI
            return true;
        }
        else
        {
            Planet originParentPlanet = originOrbit.Planet.parentBody;

            if( (originParentPlanet != null) && (originParentPlanet == destinationOrbit.Planet) )
            {
                // We are in the SOI of the child planet, we will need to calculate an exit trajectory
                needsToCalculateExitTrajectory = true;
                return true;
            }
        }

        return false;
    }

    private void calculateDeltaVandTransferEfficiency(Double2 startVelocity, Double2 endVelocity, out Double2 DeltaV, out double dv, out double transferEfficiency)
    {
        if (transferOrbit != null)
        {
            DeltaV = endVelocity - startVelocity;
            dv = DeltaV.magnitude;

            if (dv > 1.0)
            {
                // usual calculation (deltaV must not be null)
                transferEfficiency = Math.Abs(endVelocity.magnitude - startVelocity.magnitude) / dv;
            }
            else
            {
                // not significant
                transferEfficiency = 1.0;
            }
        }
        else
        {
            // No tranfer exists -> return default values
            DeltaV = new Double2(0.0, 0.0);
            dv = 0.0;
            transferEfficiency = 0.0;
        }
    }

    private double calculateInsertionDeltaV(double dv_mainTransfer)
    {
        if(targetPlanet == null)
        {
            return dv_mainTransfer;
        }
        else
        {
            double periapsis = targetPlanet.TimewarpRadius_Descend; // Will be the supposed insertion radius

            // Calculate the arrival speed at periapsis using the specific energy formula
            double arrivalDV = dv_mainTransfer * dv_mainTransfer + 2.0 * targetPlanet.mass * (1.0 / periapsis - 1.0 / targetPlanet.SOI);
            arrivalDV = Math.Sqrt(Math.Max(0.0, arrivalDV));

            // calculate orbital speed at periapsis
            double orbitalSpeed = Math.Sqrt(targetPlanet.mass / periapsis);

            // Calculate insertion DV as the difference of the two (Absolute value should be useless, but in case of a REALLY weird situation it will avoid breaking everything)
            double insertionDV = Math.Abs(arrivalDV - orbitalSpeed);

            LOG(LOG_LEVEL.DEBUG, "calculateInsertionDeltaV: periapsis = " + periapsis + "; dv_mainTransfer = " + dv_mainTransfer + "; arrivalDV = " + arrivalDV + "; orbitalSpeed = " + orbitalSpeed + "; insertionDV = " + insertionDV);

            return insertionDV;
        }
    }

    public void calculateTransfer(double startTime, double endTime)
    {
        // CHECK DATA COMPATIBILITY
        // ------------------------
        if ((originOrbit == null) || (destinationOrbit == null))
        {
            LOG(LOG_LEVEL.ERROR, "AnaisTransfer::calculateTransfer: originOrbit or destinationOrbit is null");
            return;
        }

        if (!checkPlanetaryConfiguration())
        {
            LOG(LOG_LEVEL.WARNING, "AnaisTransfer::calculateTransfer: planetary configuration is not good");
            return;
        }

        // Check dates compatibility
        if (checkDates(startTime, endTime) == false)
        {
            LOG(LOG_LEVEL.WARNING, "AnaisTransfer::calculateTransfer: dates are not consistent");
            return;
        }

        departureTime = startTime;
        arrivalTime = endTime;

        if(!needsToCalculateExitTrajectory)
        {
            LOG(LOG_LEVEL.DEBUG, "AnaisTransfer::calculateTransfer: Calculate local transfer");
            calculateTransferInSameSOI();
        }
        else
        {
            LOG(LOG_LEVEL.DEBUG, "AnaisTransfer::calculateTransfer: Calculate interplanetary transfer");
            calculateTransferFromChildPlanet();
        }
    }

    private void calculateTransferInSameSOI()
    {
        // Get starting and ending locations
        Location playerLocation = originOrbit.GetLocation(departureTime);
        Location targetLocation = destinationOrbit.GetLocation(arrivalTime);

        // Calculate the transfer
        transferOrbit = LambertSolver.CalculateTrajectory(originOrbit.Planet, playerLocation.position, targetLocation.position, departureTime, arrivalTime, originOrbit.direction);

        if (transferOrbit != null)
        {
            // Get starting and ending locations on the transfer orbits
            Location startTransferLocation = transferOrbit.GetLocation(departureTime);
            Location endTransferLocation = transferOrbit.GetLocation(arrivalTime);

            // Calculate deltaV and transfer efficiency at start
            calculateDeltaVandTransferEfficiency(playerLocation.velocity, startTransferLocation.velocity, out deltaV_start, out dv_start, out double start_transfer_eff);

            // Calculate deltaV and transfer efficiency at arrival
            calculateDeltaVandTransferEfficiency(endTransferLocation.velocity, targetLocation.velocity, out Double2 deltaV_end_main, out double dv_end_main, out double end_main_transfer_eff);


            // Calculate deltaV and transfer efficiency at arrival
            double end_transfer_eff = 0.0;
            bool noArrivalPlanet = false;
            if ((targetPlanet == null) || noArrivalPlanet)
            {
                // Target is punctual, the last data is what we need
                deltaV_end = deltaV_end_main;
                dv_end = dv_end_main;
            }
            else
            {
                // Target is a planet, we must take into account its SOI
                deltaV_end = deltaV_end_main; // wrong value, but it has no sense and is not used from now on anyway
                dv_end = calculateInsertionDeltaV(dv_end_main); // Calculate the insertion speed after factoring in the planet's gravity
            }

            end_transfer_eff = end_main_transfer_eff;

            // global deltaV and efficiency
            total_dv = dv_start + dv_end;

            if (total_dv > 1.0)
            {
                transfer_efficiency = (start_transfer_eff * dv_start + end_transfer_eff * dv_end) / total_dv;
            }
            else
            {
                // not significant
                transfer_efficiency = 1.0;
            }

            LOG(LOG_LEVEL.DEBUG, "  calculateTransferInSameSOI: Transfer calculated with success");
            LOG(LOG_LEVEL.DEBUG, "    ΔVstart = " + dv_start + "; ΔVend = " + dv_end);
        }
        else
        {
            LOG(LOG_LEVEL.DEBUG, "  calculateTransferInSameSOI: Lambert calculation failed - OUCH!");
            LOG(LOG_LEVEL.DEBUG, "    startTransferTime = " + departureTime + "; endTransferTime = " + arrivalTime);
        }
    }

    private void calculateTransferFromChildPlanet()
    {
        // Sont initialisés: originOrbit, destinationOrbit, departureTime, arrivalTime

        // The planets
        Planet parentPlanet = destinationOrbit.Planet;
        Planet childPlanet = originOrbit.Planet;

        // The orbit of the child body
        Orbit childPlanetOrbit = childPlanet.orbit;

        // the date at which the ship will exit the SOI of the child object (initialized with departureTime for the first iteration)
        double childPlanetExitTime = departureTime;

        // Location of ship at start time
        Location originLocation = originOrbit.GetLocation(departureTime);
        double startRadius = originLocation.position.magnitude;
        double startArgument = originLocation.position.AngleRadians;

        // location of target at arrival time
        Location targetLocation = destinationOrbit.GetLocation(arrivalTime);

        // Planet location when the ship exists the SOI
        Location childPlanetLocationAtExitTime = childPlanetOrbit.GetLocation(childPlanetExitTime);

        // Calculate the transfer without taking into account the SOI, taking the planet position as the starting position
        transferOrbit = LambertSolver.CalculateTrajectory(parentPlanet, childPlanetLocationAtExitTime.position, targetLocation.position, childPlanetExitTime, arrivalTime, childPlanetOrbit.direction);

        // Stop calculation immediately if it failed
        if (transferOrbit == null)
        {
            childTransferOrbit = null;
            LOG(LOG_LEVEL.DEBUG, "ERROR: initial transfer calculation failed");
            return;
        }

        // Calculate velocity at the exit of the SOI
        Location locationAtMainTransferStart = transferOrbit.GetLocation(childPlanetExitTime);
        Double2 exitVelocity = locationAtMainTransferStart.velocity - childPlanetLocationAtExitTime.velocity;

        Double2 oldExitVelocity;
        double oldChildPlanetExitTime;

        bool success = false;

        uint nbIterations = 0;

        do
        {
            nbIterations++;

            // Calculate the exit trajectory that satisfies all those parameters: starting position, exit velocity
            childTransferOrbit = EjectionTrajectoryCalculator.calculateEjectionTrajectories(childPlanet, startRadius, startArgument, departureTime, exitVelocity.magnitude, exitVelocity.AngleRadians, originOrbit.direction, exit: true);

            if (childTransferOrbit == null)
            {
                transferOrbit = null;
                LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate Ejection trajectory (iteration " + nbIterations + ")");
                return;
            }

            // Recalculate exit time and exit position
            oldChildPlanetExitTime = childPlanetExitTime;
            childPlanetExitTime = childTransferOrbit.orbitEndTime;
            Location exitLocation_ChildPlanetRef = childTransferOrbit.GetLocation(childPlanetExitTime);
            childPlanetLocationAtExitTime = childPlanetOrbit.GetLocation(childPlanetExitTime);

            LOG(LOG_LEVEL.DEBUG, " Iteration " + nbIterations + " Original exit velocity: V = " + exitVelocity.magnitude + "; heading = " + exitVelocity.AngleRadians);
            LOG(LOG_LEVEL.DEBUG, " Iteration " + nbIterations + "      New exit velocity: V = " + exitLocation_ChildPlanetRef.velocity.magnitude + "; heading = " + exitLocation_ChildPlanetRef.velocity.AngleRadians);

            // New starting position for the main trajectory calculation
            Double2 startingPosition = childPlanetLocationAtExitTime.position + exitLocation_ChildPlanetRef.position;

            // Calculate the trajectory from that new exit point and starting date
            transferOrbit = LambertSolver.CalculateTrajectory(parentPlanet, startingPosition, targetLocation.position, childPlanetExitTime, arrivalTime, childPlanetOrbit.direction);

            // Stop calculation immediately if it failed
            if (transferOrbit == null)
            {
                childTransferOrbit = null;
                LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate transfer trajectory (iteration " + nbIterations + ")");
                return;
            }

            // Calculate velocity at the exit of the SOI
            locationAtMainTransferStart = transferOrbit.GetLocation(childPlanetExitTime);
            oldExitVelocity = exitVelocity;
            exitVelocity = locationAtMainTransferStart.velocity - childPlanetLocationAtExitTime.velocity;

            // Fixer un critère de sortie: variation négligeable de la date de sortie devant la durée du transfert, variation négligeable de la vitesse d'éjection
            double exitTimeVariation = Math.Abs(childPlanetExitTime - oldChildPlanetExitTime) / (arrivalTime - departureTime);
            double exitSpeedVariation = (exitVelocity - oldExitVelocity).magnitude;

            LOG(LOG_LEVEL.DEBUG, "Iteration " + nbIterations + ": ");
            LOG(LOG_LEVEL.DEBUG, "  Child planet exit time = " + (childPlanetExitTime - departureTime) + " (old = " + (oldChildPlanetExitTime - departureTime) + "); exitTimeVariation = " + exitTimeVariation);
            LOG(LOG_LEVEL.DEBUG, "  Exit speed: V = " + exitVelocity.magnitude + ", heading = " + exitVelocity.AngleRadians + " (old : V = " + oldExitVelocity.magnitude + ", heading = " + oldExitVelocity.AngleRadians + "); exitSpeedVariation = " + exitSpeedVariation);

            // exit condition: exit time and speed must not vary too much
            if (/*(exitTimeVariation < 0.0001) &&*/ (exitSpeedVariation < 0.1))
            {
                success = true;
            }
            else
            {
                // TEST: loop again by taking as exit velocity the median between the 2 calculated exit velocities
                // This is an attempt to avoid looping alternatively between 2 values
                Double2 exitSpeedDelta = exitVelocity - oldExitVelocity;
                exitVelocity = 0.5 * (exitVelocity + oldExitVelocity);



                // Use a weighted median between both velocities
                //exitVelocity = ((childPlanetExitTime - departureTime) * oldExitVelocity + (arrivalTime - childPlanetExitTime) * exitVelocity) / (arrivalTime - departureTime);

                LOG(LOG_LEVEL.DEBUG, "  New exit velocity: V = " + exitVelocity.magnitude + ", heading = " + exitVelocity.AngleRadians + "; exit time = " + (childPlanetExitTime - departureTime) + "; total time = " + (arrivalTime - departureTime));
                LOG(LOG_LEVEL.DEBUG, "  ΔVx = " + exitSpeedDelta.x + "; ΔVy = " + exitSpeedDelta.y);

                
            }

        } while ((nbIterations < 20) && (success == false));


        if(success == false)
        {
            // Something went wrong, make sure we have invalid data
            transferOrbit = null;
            childTransferOrbit = null;
            LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate interplanetary transfer: too many iterations");
            return;
        }

        // Calculate all data and efficiency now
        // -------------------------------------
        Location locationAfterStartBurn = childTransferOrbit.GetLocation(departureTime);
        Location locationBeforeEndBurn = transferOrbit.GetLocation(arrivalTime);

        // Calculate deltaV and transfer efficiency at start
        calculateDeltaVandTransferEfficiency(originLocation.velocity, locationAfterStartBurn.velocity, out deltaV_start, out dv_start, out double start_transfer_eff);

        // Calculate deltaV and transfer efficiency at main transfer start
        calculateDeltaVandTransferEfficiency(childPlanetLocationAtExitTime.velocity, locationAtMainTransferStart.velocity, out Double2 deltaV_start_main, out double dv_start_main, out double start_main_transfer_eff);

        // Calculate deltaV and transfer efficiency at main transfer arrival
        calculateDeltaVandTransferEfficiency(locationBeforeEndBurn.velocity, targetLocation.velocity, out Double2 deltaV_end_main, out double dv_end_main, out double end_main_transfer_eff);

        // Calculate deltaV and transfer efficiency at arrival
        double end_transfer_eff = 0.0;
        bool noArrivalPlanet = false;
        if ((targetPlanet == null) || noArrivalPlanet)
        {
            // Target is punctual, the last data is what we need
            deltaV_end = deltaV_end_main;
            dv_end = dv_end_main;
        }
        else
        {
            // Target is a planet, we must take into account its SOI
            deltaV_end = deltaV_end_main; // wrong value, but it has no sense ans is not used from now on anyway
            dv_end = calculateInsertionDeltaV(dv_end_main); // Calculate the insertion speed after factoring in the planet's gravity
        }

        // global deltaV and efficiency
        total_dv = dv_start + dv_end;
        start_transfer_eff = start_transfer_eff * start_main_transfer_eff;
        end_transfer_eff = end_main_transfer_eff; // the planet insertion efficiency is supposed to be 1 if target is a planet, so no multiplication to make

        if (total_dv > 1.0)
        {
            transfer_efficiency = (start_transfer_eff * dv_start + end_transfer_eff * dv_end) / total_dv;
        }
        else
        {
            // not significant
            transfer_efficiency = 1.0;
        }
    }


    public void calculateReturnToPlanet(double departureTime_)
    {
        if (originOrbit == null) return;

        Orbit mainOrbit = null;
        departureTime = departureTime_;
        double startEfficiency = 1.0; // default value

        // Check planetary configuration
        if (targetPlanet == null)
        {
            return;
        }
        else if(originOrbit.Planet == targetPlanet)
        {
            mainOrbit = originOrbit;
        }
        else if(originOrbit.Planet.parentBody == targetPlanet)
        {
            mainOrbit = originOrbit.Planet.orbit;
        }
        else
        {
            // invalid configuration for a return to mother planet
            return;
        }

        // The radius we will target for the destination
        double targetRadius = targetPlanet.TimewarpRadius_Descend;
        double targetOrbitalSpeed = Math.Sqrt(targetPlanet.mass / targetRadius);

        // The target orbit: a perfectly circular orbit 
        Orbit targetOrbit = Orbit_Utils.CreateOrbit(targetRadius, 0.0, 0.0, mainOrbit.direction, mainOrbit.Planet, PathType.Eternal, null);

        if(targetOrbit == null)
        {
            return; // the risk is purely theoretical actually...
        }

        Location startMainLocation = mainOrbit.GetLocation(departureTime);

        // Calculate the Hohmann transfer as the return trajectory
        Orbit_Utils.CalculateHohmannTransfer(mainOrbit, targetOrbit, startMainLocation, out transferOrbit, out double _);

        if(transferOrbit == null)
        {
            return;
        }


        // Part for transfer from a satellite
        if(originOrbit.Planet.parentBody == targetPlanet)
        {
            // We are in the situation in which we return to the mother planet from one of its satellites
            int nbIterations = 0;
            double diffExitVelocity = 1.0;

            Orbit childOrbit = originOrbit;
            Location startShipLocation = childOrbit.GetLocation(departureTime);

            Location ShipLocationAtExitTime = transferOrbit.GetLocation(departureTime);
            Double2 exitVelocity = ShipLocationAtExitTime.velocity - startMainLocation.velocity;

            do
            {
                nbIterations++;

                // Calculate the corresponding ejection trajectory
                childTransferOrbit = EjectionTrajectoryCalculator.calculateEjectionTrajectories(childOrbit.Planet, startShipLocation.position.magnitude, startShipLocation.position.AngleRadians, departureTime, exitVelocity.magnitude, exitVelocity.AngleRadians, childOrbit.direction);
                if (childTransferOrbit == null)
                {
                    transferOrbit = null;
                    return;
                }

                // Get ejection date
                double ejectionDate = childTransferOrbit.orbitEndTime;

                LOG(LOG_LEVEL.DEBUG, "Iteration " + nbIterations + ": exit speed = " + exitVelocity.magnitude + " m/s; " + exitVelocity.AngleRadians + " rad; transfer duration = " + (ejectionDate - departureTime) + " s");

                // Get locations of satellite and ship at exit time
                startMainLocation = mainOrbit.GetLocation(ejectionDate);
                Location shipLocationAtExitTime = childTransferOrbit.GetLocation(ejectionDate);

                // Ship location in the main planet's frame of reference when it exits
                Location shipLocationAtExitTimeInMainPlanet = new Location(ejectionDate, targetPlanet, startMainLocation.position + shipLocationAtExitTime.position, startMainLocation.velocity + shipLocationAtExitTime.velocity);

                // Calculate the resulting orbit in the main body, by supposing its velocity matches the planet velocity in direction
                // This is because we want to calculate a Hohmann transfer from an orbit which speed matches the planet velocity in direction, the orbit itself doesn't matter
                Location fictiveShipLocationAtExitTimeInMainPlanet = new Location(ejectionDate, targetPlanet, startMainLocation.position + shipLocationAtExitTime.position, startMainLocation.velocity);

                Orbit fictiveExitTrajectory = Orbit.TryCreateOrbit(fictiveShipLocationAtExitTimeInMainPlanet, true, false, out bool success);
                if (!success)
                {
                    transferOrbit = null;
                    childTransferOrbit = null;
                    return;
                }

                // Calculate the Hohmann transfer from that fictive orbit: our exit trajectory should match that Hohmann transfer
                Orbit_Utils.CalculateHohmannTransfer(fictiveExitTrajectory, targetOrbit, fictiveShipLocationAtExitTimeInMainPlanet, out transferOrbit, out double _);

                if(transferOrbit == null)
                {
                    transferOrbit = null;
                    childTransferOrbit = null;
                    return;
                }

                // Calculate the difference in exit velocity
                ShipLocationAtExitTime = transferOrbit.GetLocation(transferOrbit.orbitStartTime);
                exitVelocity = ShipLocationAtExitTime.velocity - startMainLocation.velocity;
                diffExitVelocity = (ShipLocationAtExitTime.velocity - shipLocationAtExitTimeInMainPlanet.velocity).magnitude;
                LOG(LOG_LEVEL.DEBUG, "Iteration " + nbIterations + ": diffExitVelocity = " + diffExitVelocity + "; ShipLocationAtExitTime.velocity = " + ShipLocationAtExitTime.velocity + "; ship velocity = " + shipLocationAtExitTimeInMainPlanet.velocity.magnitude);

            } while ((nbIterations < 20) && (diffExitVelocity > 0.1));

            // Exit if research fails
            if(nbIterations == 20)
            {
                LOG(LOG_LEVEL.DEBUG, "ERROR: AnaisTransfer calculation failed by excessive number of iterations");
                transferOrbit = null;
                childTransferOrbit = null;
                return;
            }

            // Calculate DV and efficiency at start
            Location startShipLocationAfterBurn = childTransferOrbit.GetLocation(departureTime);
            calculateDeltaVandTransferEfficiency(startShipLocation.velocity, startShipLocationAfterBurn.velocity, out _, out dv_start, out startEfficiency);
        }
        else
        {
            // dv_start and efficiency calculation for a local transfer
            Location startLocation_mainTransfer = transferOrbit.GetLocation(departureTime);
            dv_start = Math.Abs(startMainLocation.velocity.magnitude - startLocation_mainTransfer.velocity.magnitude);
            startEfficiency = 1.0;
        }


        arrivalTime = transferOrbit.orbitEndTime;
        Location endLocation_mainTransfer = transferOrbit.GetLocation(arrivalTime);
        
        dv_end = Math.Abs(endLocation_mainTransfer.velocity.magnitude - targetOrbitalSpeed);

        total_dv = dv_start + dv_end;

        if(total_dv > 1.0)
        {
            transfer_efficiency = (startEfficiency * dv_start + dv_end) / total_dv;
        }
        else
        {
            transfer_efficiency = 1.0; // not significant
        }
        
    }

    public bool isValid()
    {
        return (transferOrbit != null);
    }

    // Local log function
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void LOG(LOG_LEVEL level, string message)
    {
#if ACTIVE_LOGS
        AnaisLogger.Log(LOG_CATEGORY.ANAIS_TRANSFER, level, message);
#endif
    }
}
