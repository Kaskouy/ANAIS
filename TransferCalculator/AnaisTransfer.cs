using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

using SFS.World;
using SFS.WorldBase;
using UnityEngine;

class AnaisTransfer
{
    // C_DELTA_V_SIGNIFICANCY_THRESHOLD: Below this threshold, a delta-V value isn't considered significant for efficiency calculations
    private const double C_DELTA_V_SIGNIFICANCY_THRESHOLD = 0.1;
    
    // the initial and the target orbit
    public Orbit originOrbit;
    public Orbit destinationOrbit;
    public Planet targetPlanet;
    public double targetAltitude;

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

    public ANAIS_Settings.E_TRANSFER_TYPE transferType;

    // data for when the transfer starts from a child planet
    bool needsToCalculateExitTrajectory = false;
    public Orbit childTransferOrbit = null;


    // --------------------------------------------------------------------------------------------
    //                                    AnaisTransfer
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   The constructor. Initializes all data for the future calculations. Note that data
    //   consistency is not checked at that point, the caller has to make sure the data is
    //   appropriate.
    // 
    // RESULT:
    //   Nothing, it's a constructor...
    // --------------------------------------------------------------------------------------------
    public AnaisTransfer(Orbit _originOrbit, Orbit _destinationOrbit, Planet _targetPlanet, double _targetAltitude, ANAIS_Settings.E_TRANSFER_TYPE transferType)
    {
        originOrbit = _originOrbit;
        destinationOrbit = _destinationOrbit;
        targetPlanet = _targetPlanet;
        targetAltitude = _targetAltitude;

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
        this.transferType = transferType;
    }

    public Double2 GetDeltaV_start()
    {
        return deltaV_start;
    }

    public double GetDeltaV_valueToMinimize(ANAIS_Settings.E_TRANSFER_TYPE transferType)
    {
        if (transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS) return total_dv;
        else return dv_start;
    }


    // --------------------------------------------------------------------------------------------
    //                               areDeltaV_valuesSignificant
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function tells, once the transfers have been calculated, if the deltaV values are
    //   significant. In particular, for an interplanetary transfer that hasn't been refined, this
    //   will not be the case. This is to warn that the delta-V value calculated is not accurate
    //   and may not be displayed.
    // 
    // RESULT:
    //   Returns true if deltaV values are significant, false otherwise.
    // --------------------------------------------------------------------------------------------
    public bool areDeltaV_valuesSignificant()
    {
        if(transferOrbit == null)
        {
            return false;
        }
        else if(needsToCalculateExitTrajectory && (childTransferOrbit == null))
        {
            return false;
        }
        else
        {
            return true;
        }
    }


    // --------------------------------------------------------------------------------------------
    //                               checkPlanetaryConfiguration
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   Checks if the planet configuration is correct in the case of a transfer (not a return to
    //   planet trajectory). Also checks if the ejection trajectory has to be calculated.
    // 
    // RESULT:
    //   true if planetary configuration is good, false otherwise.
    // --------------------------------------------------------------------------------------------
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


    // --------------------------------------------------------------------------------------------
    //                             calculateDeltaVandTransferEfficiency
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function calculates the delta-V associated to a burn and its efficiency, based on the
    //   start and end velocity.
    // 
    // RESULT:
    //   Returns the delta-V and efficiency associated to the burn.
    // --------------------------------------------------------------------------------------------
    private void calculateDeltaVandTransferEfficiency(Double2 startVelocity, Double2 endVelocity, out Double2 DeltaV, out double dv, out double transferEfficiency)
    {
        DeltaV = endVelocity - startVelocity;
        dv = DeltaV.magnitude;

        if (dv > C_DELTA_V_SIGNIFICANCY_THRESHOLD)
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


    // --------------------------------------------------------------------------------------------
    //                                  calculateInsertionDeltaV
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function calculates the delta-V that will be needed for an insertion into low orbit.
    //   It's simply calculated from the speed at which the ship enters the SOI and using the
    //   specific energy formula.
    // 
    // RESULT:
    //   Returns the delta-V required for insertion.
    // --------------------------------------------------------------------------------------------
    private double calculateInsertionDeltaV(double dv_mainTransfer)
    {
        if(targetPlanet == null)
        {
            return dv_mainTransfer;
        }
        else
        {
            double periapsis = targetPlanet.Radius + targetAltitude; // Will be the supposed insertion radius

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


    // --------------------------------------------------------------------------------------------
    //                                    calculateTransfer
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function is used to calculate a transfer trajectory for a set departure time and
    //   arrival time. It detects automatically the configuration (local/interplanetary transfer)
    //   and calls the appropriate function. In the case of an interplanetary transfer, only the
    //   main transfer will be calculated though, refineTransferCalculation has to be called to
    //   calculate the ejection trajectory. This is done in 2 steps for optimization purposes:
    //   calculating the ejection trajectory is a costly process, so in practice it's done only if
    //   the main transfer is promising enough.
    //   The AnaisTransfer has to be initialized through the constructor.
    // 
    // RESULT:
    //   In case of success, the trajectory is calculated in transferOrbit. The transfer (deltaV)
    //   and efficiency data are calculated aswell.
    //   If the calculation failed, transferOrbit is set to null and the function returns.
    // --------------------------------------------------------------------------------------------
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

        departureTime = startTime;
        arrivalTime = endTime;

        if(!needsToCalculateExitTrajectory)
        {
            //LOG(LOG_LEVEL.DEBUG, "AnaisTransfer::calculateTransfer: Calculate local transfer");
            calculateTransferInSameSOI();
        }
        else
        {
            //LOG(LOG_LEVEL.DEBUG, "AnaisTransfer::calculateTransfer: Calculate interplanetary transfer");
            calculateTransferFromChildPlanet();
        }
    }


    // --------------------------------------------------------------------------------------------
    //                             calculateTransferInSameSOI
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function is used to calculate a local transfer.
    //   The AnaisTransfer has to be initialized through a constructor, and a departure and arrival
    //   date have to be provided.
    // 
    // RESULT:
    //   In case of success, the trajectory is calculated in transferOrbit. The transfer (deltaV)
    //   and efficiency data are calculated aswell.
    //   If the calculation failed, transferOrbit is set to null and the function returns.
    // --------------------------------------------------------------------------------------------
    private void calculateTransferInSameSOI()
    {
        // Get starting and ending locations
        Location playerLocation = originOrbit.GetLocation(departureTime);
        Location targetLocation = destinationOrbit.GetLocation(arrivalTime);

        // Calculate the transfer
        transferOrbit = LambertSolver.CalculateTrajectory(originOrbit.Planet, playerLocation.position, targetLocation.position, departureTime, arrivalTime, originOrbit.direction);

        if (transferOrbit != null)
        {
            LOG(LOG_LEVEL.DEBUG, "transferOrbit: slr = " + transferOrbit.slr + "; ecc = " + transferOrbit.ecc + "; sma = " + transferOrbit.sma);
            LOG(LOG_LEVEL.DEBUG, "transferOrbit: phi = " + transferOrbit.arg + "; Tp = " + transferOrbit.periapsisPassageTime);

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

            if(transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS)
            {
                if (total_dv > C_DELTA_V_SIGNIFICANCY_THRESHOLD)
                {
                    transfer_efficiency = (start_transfer_eff * dv_start + end_transfer_eff * dv_end) / total_dv;
                }
                else
                {
                    // not significant
                    transfer_efficiency = 1.0;
                }
            }
            else
            {
                transfer_efficiency = start_transfer_eff;
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


    // --------------------------------------------------------------------------------------------
    //                             calculateTransferFromChildPlanet
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function is used to calculate a first approximation of an interplanetary transfer.
    //   The interplanetary transfer is treated as a local transfer from the orbited planet position
    //   as a first approximation.
    //   The AnaisTransfer has to be initialized through a constructor, and a departure and arrival
    //   date have to be provided.
    //   In practice, the result is considered sufficient on high time-warp ratios. Otherwise, and
    //   provided the calculation is a success, refineTransferCalculation has to be called after to
    //   take into account the ejection trajectory.
    // 
    // RESULT:
    //   In case of success, the trajectory is calculated in transferOrbit. The transfer (deltaV)
    //   and efficiency data are calculated aswell.
    //   If the calculation failed, transferOrbit is set to null and the function returns.
    // --------------------------------------------------------------------------------------------
    private void calculateTransferFromChildPlanet()
    {
        // The planets
        Planet parentPlanet = destinationOrbit.Planet;
        Planet childPlanet = originOrbit.Planet;

        // The orbit of the child body
        Orbit childPlanetOrbit = childPlanet.orbit;

        // location of target at arrival time
        Location targetLocation = destinationOrbit.GetLocation(arrivalTime);

        // Planet location when the ship exists the SOI
        Location childPlanetLocation = childPlanetOrbit.GetLocation(departureTime);

        // Calculate the transfer without taking into account the SOI, taking the planet position as the starting position
        transferOrbit = LambertSolver.CalculateTrajectory(parentPlanet, childPlanetLocation.position, targetLocation.position, departureTime, arrivalTime, childPlanetOrbit.direction);

        if(transferOrbit != null) // Calculate transfer efficiency like if it was a transfer performed in the main body's frame from the position of the child planet
        {
            // Get starting and ending locations on the transfer orbits
            Location startTransferLocation = transferOrbit.GetLocation(departureTime);
            Location endTransferLocation = transferOrbit.GetLocation(arrivalTime);

            // Calculate deltaV and transfer efficiency at start
            calculateDeltaVandTransferEfficiency(childPlanetLocation.velocity, startTransferLocation.velocity, out deltaV_start, out dv_start, out double start_transfer_eff);

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

            if (transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS)
            {
                if (total_dv > C_DELTA_V_SIGNIFICANCY_THRESHOLD)
                {
                    transfer_efficiency = (start_transfer_eff * dv_start + end_transfer_eff * dv_end) / total_dv;
                }
                else
                {
                    // not significant
                    transfer_efficiency = 1.0;
                }
            }
            else
            {
                transfer_efficiency = start_transfer_eff;
            }

            //LOG(LOG_LEVEL.DEBUG, "  calculateTransferFromChildPlanet: Main transfer calculated with success");
            //LOG(LOG_LEVEL.DEBUG, "    ΔVstart = " + dv_start + "; ΔVend = " + dv_end);

            /*bool success = refineTransferCalculation();

            if(!success)
            {
                LOG(LOG_LEVEL.DEBUG, "    Transfer from child planet: failed to refine transfer calculation - Abort");
            }*/
        }
        else
        {
            LOG(LOG_LEVEL.DEBUG, "    Transfer from child planet: failed to calculate main transfer orbit - Abort");
        }
    }


    // --------------------------------------------------------------------------------------------
    //                                  refineTransferCalculation
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function is used to refine the calculation in the case of an interplanetary transfer.
    //   This function is to be called after calculateTransferFromChildPlanet, provided that one
    //   has succeeded. In this context, transferOrbit is initialized.
    //   The purpose is to take into account the ejection trajectory. The
    //   calculateTransferFromChildPlanet function calculates a transfer from the orbited planet as
    //   if it was a local transfer as a first approximation. This function starts from this
    //   calculation as a basis.
    // 
    // RESULT:
    //   In case of success, the ejection trajectory is calculated in childTransferOrbit. The
    //   transfer (deltaV) and efficiency data are calculated aswell.
    //   If the calculation failed, transferOrbit and childTransferOrbit are both set to null and
    //   the function returns.
    // --------------------------------------------------------------------------------------------
    public bool refineTransferCalculation()
    {
        // The planets
        Planet parentPlanet = destinationOrbit.Planet;
        Planet childPlanet = originOrbit.Planet;

        // The orbit of the child body
        Orbit childPlanetOrbit = childPlanet.orbit;

        // location of target at arrival time
        Location targetLocation = destinationOrbit.GetLocation(arrivalTime);

        // Location of ship at start time
        Location startLocation = originOrbit.GetLocation(departureTime);

        // time at which the ship exits the SOI (first initialized at departure time)
        double exitTime = departureTime;

        // Planet location when the ship exists the SOI
        Location childPlanetLocation = childPlanetOrbit.GetLocation(exitTime);

        // exitVelocity_LambertArc: the velocity at which the ship exits the SOI (in the parent body's frame of reference) - evaluated from the transfer evaluated by the Lambert solver
        Double2 exitVelocity_LambertArc = transferOrbit.GetLocation(exitTime).velocity - childPlanetLocation.velocity;

        // exitVelocity_Ejection: same thing, but evaluated from the ejection trajectory from the child body - initialised with the velocity calculated above
        Double2 exitVelocity_Ejection = exitVelocity_LambertArc;

        uint nbIterations = 0;
        bool success = false;

        Double2 deltaV = new Double2(0.0, 0.0);
        Vector2D_FixedPointSolver theSolver = new Vector2D_FixedPointSolver();

        // Iterate until we find a compatible transfer and ejection trajectory
        // -------------------------------------------------------------------
        do
        {
            nbIterations++;

            // Calculate the exit trajectory that satisfies all those parameters: starting position, exit velocity
            //LOG(LOG_LEVEL.DEBUG, "    refineTransferCalculation: calculate ejection trajectory");
            childTransferOrbit = EjectionTrajectoryCalculator.calculateEjectionTrajectories(childPlanet, startLocation, departureTime, exitVelocity_Ejection, originOrbit.direction, exit: true);

            if (childTransferOrbit == null)
            {
                transferOrbit = null;
                LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate Ejection trajectory (iteration " + nbIterations + ")");
                return false;
            }

            // calculate exit time
            exitTime = childTransferOrbit.orbitEndTime;

            // calculate exit position (which will be the starting position for the Lambert arc that will be calculated right after)
            childPlanetLocation = childPlanetOrbit.GetLocation(exitTime);
            Double2 startingPosition = childPlanetLocation.position + childTransferOrbit.GetLocation(exitTime).position;

            // Calculate the transfer in the main body's frame of reference from the newly calculated starting position
            //LOG(LOG_LEVEL.DEBUG, "    refineTransferCalculation: calculate main Lambert arc");
            transferOrbit = LambertSolver.CalculateTrajectory(parentPlanet, startingPosition, targetLocation.position, exitTime, arrivalTime, childPlanetOrbit.direction);

            if(transferOrbit == null)
            {
                LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate transfer trajectory (iteration " + nbIterations + ")");
                childTransferOrbit = null;
                return false;
            }

            // evaluate new exit velocity
            exitVelocity_LambertArc = transferOrbit.GetLocation(exitTime).velocity - childPlanetLocation.velocity;
            deltaV = exitVelocity_LambertArc - exitVelocity_Ejection;

            // evaluate success criteria : the exit velocity that results from ejection has to match the exit velocity deduced from the Lambert arc calculation
            if (deltaV.magnitude < 0.1)
            {
                LOG(LOG_LEVEL.DEBUG, "   CALCULATION SUCCESSFUL after " + nbIterations + " iterations (delta = " + deltaV.magnitude + " m/s)");
                LOG(LOG_LEVEL.DEBUG, "     Final value: Vx = " + exitVelocity_Ejection.x + "; Vy = " + exitVelocity_Ejection.y);
                success = true;
            }

            // evaluate new exit velocity vector for next iteration
            if(!success)
            {
                theSolver.addPoint(exitVelocity_Ejection, exitVelocity_LambertArc);

                exitVelocity_Ejection = theSolver.GetEstimatedSolution();
            }

        } while ((nbIterations < 20) && (success == false));

        // Exit if the transfer couldn't be calculated
        // -------------------------------------------
        if (!success)
        {
            LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate interplanetary transfer - too many iterations");
            transferOrbit = null;
            childTransferOrbit = null;
            return false;
        }

        
        // Calculate all data and efficiency now
        // -------------------------------------
        // Calculate deltaV and transfer efficiency at start, from ship orbit to child transfer orbit
        calculateDeltaVandTransferEfficiency(startLocation.velocity, childTransferOrbit.GetLocation(departureTime).velocity, out deltaV_start, out dv_start, out double start_transfer_eff);

        // Calculate deltaV and transfer efficiency for start at main transfer: from child planet velocity to main transfer velocity
        calculateDeltaVandTransferEfficiency(childPlanetLocation.velocity, transferOrbit.GetLocation(exitTime).velocity, out Double2 deltaV_start_main, out double dv_start_main, out double start_main_transfer_eff);

        // Calculate deltaV and transfer efficiency at arrival
        calculateDeltaVandTransferEfficiency(transferOrbit.GetLocation(arrivalTime).velocity, targetLocation.velocity, out Double2 deltaV_end_main, out double dv_end_main, out double end_main_transfer_eff);

        // Calculate deltaV and transfer efficiency at arrival
        double end_transfer_eff = 0.0;
        if (targetPlanet == null)
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

        if (transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS)
        {
            if (total_dv > C_DELTA_V_SIGNIFICANCY_THRESHOLD)
            {
                transfer_efficiency = (start_transfer_eff * dv_start + end_transfer_eff * dv_end) / total_dv;
            }
            else
            {
                // not significant
                transfer_efficiency = 1.0;
            }
        }
        else
        {
            transfer_efficiency = start_transfer_eff;
        }

        return true;
    }


    // --------------------------------------------------------------------------------------------
    //                                    calculateReturnToPlanet
    // --------------------------------------------------------------------------------------------
    // DESCRIPTION:
    //   This function is used to calculate a return-to-planet trajectory. It works for both a 
    //   local transfer and for a transfer from a moon.
    //   All parameters of the AnaisTransfer has to be initialized through the constructor, except
    //   the destinationOrbit data that is useless in this case.
    // 
    // RESULT:
    //   In case of success, the main transfer is calculated in transferOrbit, and the ejection
    //   trajectory (if transfer from a moon) is calculated in childTransferOrbit. In both cases,
    //   the transfer (deltaV) and efficiency data are calculated aswell.
    //   If the configuration is not a return to planet configuration, the function returns without
    //   doing anything.
    //   If the calculation failed, transferOrbit and childTransferOrbit are both set to null and
    //   the function returns.
    // --------------------------------------------------------------------------------------------
    public void calculateReturnToPlanet(double departureTime_, bool refineTransfer)
    {
        LOG(LOG_LEVEL.DEBUG, "calculateReturnToPlanet called");
        if (originOrbit == null) return;

        departureTime = departureTime_;
        Orbit mainOrbit = null;
        bool suggestOrbitInsertion;

        // Check planetary configuration
        if (targetPlanet == null)
        {
            return;
        }
        else if(originOrbit.Planet == targetPlanet)
        {
            LOG(LOG_LEVEL.DEBUG, "calculateReturnToPlanet: Return to mother planet needed");
            needsToCalculateExitTrajectory = false;
            mainOrbit = originOrbit;
            suggestOrbitInsertion = true;
        }
        else if(originOrbit.Planet.parentBody == targetPlanet)
        {
            LOG(LOG_LEVEL.DEBUG, "calculateReturnToPlanet: Return to mother planet needed with escape trajectory");
            needsToCalculateExitTrajectory = true;
            mainOrbit = originOrbit.Planet.orbit;
            suggestOrbitInsertion = false;
        }
        else
        {
            // invalid configuration for a return to mother planet
            return;
        }

        // Make a copy of this because we are going to mess with this...
        Location tmpLocation = mainOrbit.GetLocation(departureTime);
        Location startMainLocation = new Location(tmpLocation.time, tmpLocation.planet, tmpLocation.position, tmpLocation.velocity);


        LOG(LOG_LEVEL.DEBUG, "calculateReturnToPlanet: calculating data");
        double startEfficiency = 1.0; // default value
        double targetRadius = targetPlanet.Radius + targetAltitude;
        double targetOrbitalSpeed = Math.Sqrt(targetPlanet.mass / targetRadius);
        bool includeInsertionCost = (transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS);


        if ((originOrbit.Planet == targetPlanet) || ((originOrbit.Planet.parentBody == targetPlanet) && !refineTransfer))
        {
            // LOCAL TRANSFER (or interplanetary unrefined tranfer)
            // --------------
            bool success = ReturnToPlanetCalculator.calculateTrajectory(targetPlanet, startMainLocation, targetRadius, includeInsertionCost, suggestOrbitInsertion, out transferOrbit);

            if (!success)
            {
                LOG(LOG_LEVEL.WARNING, "Return to planet trajectory could not be calculated - abort calculation");
                return;
            }

            // dv_start and efficiency calculation for a local transfer
            Location startLocation_originOrbit = mainOrbit.GetLocation(departureTime);
            Location startLocation_mainTransfer = transferOrbit.GetLocation(departureTime);
            
            calculateDeltaVandTransferEfficiency(startLocation_originOrbit.velocity, startLocation_mainTransfer.velocity, out deltaV_start, out dv_start, out startEfficiency);

            total_dv = dv_start;

            // calculate dv_end and associated efficiency if relevant
            if (transferOrbit.PathType == PathType.Encounter)
            {
                arrivalTime = transferOrbit.orbitEndTime;
                Location endLocation_mainTransfer_ = transferOrbit.GetLocation(arrivalTime);

                dv_end = Math.Abs(endLocation_mainTransfer_.velocity.magnitude - targetOrbitalSpeed);

                total_dv += dv_end;
            }
            else
            {
                dv_end = 0.0;
            }

            // calculate efficiency
            if (total_dv > C_DELTA_V_SIGNIFICANCY_THRESHOLD)
            {
                transfer_efficiency = (startEfficiency * dv_start + dv_end) / total_dv;
            }
            else
            {
                transfer_efficiency = 1.0; // not significant
            }

            return;
        }
        else 
        {
            // INTERPLANETARY TRANSFER
            // -----------------------
            bool success = ReturnToPlanetCalculator.calculateTrajectory(targetPlanet, startMainLocation, targetRadius, includeInsertionCost, suggestOrbitInsertion: false, out transferOrbit);

            if (!success)
            {
                LOG(LOG_LEVEL.WARNING, "Return to planet trajectory could not be calculated - abort calculation");
                return;
            }

            // The orbit of the child body
            Orbit childPlanetOrbit = originOrbit.Planet.orbit;

            // Location of ship at start time
            Location startLocation = originOrbit.GetLocation(departureTime);

            // time at which the ship exits the SOI (first initialized at departure time)
            double exitTime = departureTime;

            // Planet location when the ship exists the SOI
            Location childPlanetLocation = childPlanetOrbit.GetLocation(exitTime);

            // exitVelocity_ReturnTrajectory: the velocity at which the ship exits the SOI (in the parent body's frame of reference) - evaluated from the transfer evaluated by ReturnToPlanetCalculator
            Double2 exitVelocity_ReturnTrajectory = transferOrbit.GetLocation(exitTime).velocity - childPlanetLocation.velocity;

            // exitVelocity_Ejection: same thing, but evaluated from the ejection trajectory from the child body - initialised with the velocity calculated above
            Double2 exitVelocity_Ejection = exitVelocity_ReturnTrajectory;

            uint nbIterations = 0;
            success = false;

            Double2 deltaV = new Double2(0.0, 0.0);

            Vector2D_FixedPointSolver theSolver = new Vector2D_FixedPointSolver();

            // Iterate until we find a compatible transfer and ejection trajectory
            // -------------------------------------------------------------------
            do
            {
                nbIterations++;
                LOG(LOG_LEVEL.DEBUG, " Calculate return to planet trajectory: iteration " + nbIterations + ": exitVelocity_Ejection = (" + exitVelocity_Ejection.x + ", " + exitVelocity_Ejection.y + ")");

                // Calculate the exit trajectory that satisfies all those parameters: starting position, exit velocity
                childTransferOrbit = EjectionTrajectoryCalculator.calculateEjectionTrajectories(originOrbit.Planet, startLocation, departureTime, exitVelocity_Ejection, originOrbit.direction, exit: true);

                if (childTransferOrbit == null)
                {
                    transferOrbit = null;
                    LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate Ejection trajectory (iteration " + nbIterations + ")");
                    return;
                }

                // calculate exit time
                exitTime = childTransferOrbit.orbitEndTime;

                // calculate location of child planet at exit time
                childPlanetLocation = childPlanetOrbit.GetLocation(exitTime);

                // calculate exit location in the child planet reference frame
                Location exitLocation = childTransferOrbit.GetLocation(exitTime);
                
                // calculate the starting position and velocity for the return to mother planet trajectory
                startMainLocation.position = childPlanetLocation.position + exitLocation.position;
                startMainLocation.velocity = childPlanetLocation.velocity; // We consider the exiting ship has velocity of the planet; the calculator is expected to give a transfer which matches ejection speed
                startMainLocation.time = exitTime;

                // calculate opimal return to planet trajectory from the ejection point
                bool calculationOK = ReturnToPlanetCalculator.calculateTrajectory(targetPlanet, startMainLocation, targetRadius, includeInsertionCost, suggestOrbitInsertion: false, out transferOrbit);

                if (!calculationOK)
                {
                    LOG(LOG_LEVEL.WARNING, "Return to planet trajectory could not be calculated - abort calculation");
                    transferOrbit = null;
                    childTransferOrbit = null;
                    return;
                }

                // evaluate new exit velocity
                exitVelocity_ReturnTrajectory = transferOrbit.GetLocation(exitTime).velocity - childPlanetLocation.velocity;
                deltaV = exitVelocity_ReturnTrajectory - exitVelocity_Ejection;

                // evaluate success criteria : the exit velocity that results from ejection has to match the exit velocity deduced from the return trajectory calculation
                if (deltaV.magnitude < 0.1)
                {
                    LOG(LOG_LEVEL.DEBUG, "   CALCULATION SUCCESSFUL after " + nbIterations + " iterations (delta = " + deltaV.magnitude + " m/s)");
                    LOG(LOG_LEVEL.DEBUG, "     Final value: Vx = " + exitVelocity_Ejection.x + "; Vy = " + exitVelocity_Ejection.y);
                    success = true;
                }

                // evaluate new exit velocity vector for next iteration
                if (!success)
                {
                    theSolver.addPoint(exitVelocity_Ejection, exitVelocity_ReturnTrajectory);
                    exitVelocity_Ejection = theSolver.GetEstimatedSolution();
                }

            } while ((nbIterations < 20) && (success == false));


            // Exit if the transfer couldn't be calculated
            // -------------------------------------------
            if (!success)
            {
                LOG(LOG_LEVEL.DEBUG, "ERROR: Failed to calculate interplanetary transfer - too many iterations");
                transferOrbit = null;
                childTransferOrbit = null;
                return;
            }

            // Calculate all data and efficiency now
            // -------------------------------------
            // Calculate deltaV and transfer efficiency at start, from ship orbit to child transfer orbit
            calculateDeltaVandTransferEfficiency(startLocation.velocity, childTransferOrbit.GetLocation(departureTime).velocity, out deltaV_start, out dv_start, out startEfficiency);

            // Calculate deltaV and transfer efficiency for start at main transfer: from child planet velocity to main transfer velocity
            calculateDeltaVandTransferEfficiency(childPlanetLocation.velocity, transferOrbit.GetLocation(exitTime).velocity, out Double2 deltaV_start_main, out double dv_start_main, out double startMainEfficiency);

            total_dv = dv_start;

            // calculate dv_end and associated efficiency if relevant
            if (transferOrbit.PathType == PathType.Encounter)
            {
                arrivalTime = transferOrbit.orbitEndTime;
                Location endLocation_mainTransfer_ = transferOrbit.GetLocation(arrivalTime);

                dv_end = Math.Abs(endLocation_mainTransfer_.velocity.magnitude - targetOrbitalSpeed);

                total_dv += dv_end;
            }
            else
            {
                dv_end = 0.0;
            }


            // global efficiency - end efficiency is 1 since it's a tangent insertion
            startEfficiency *= startMainEfficiency;

            // calculate efficiency
            if (total_dv > C_DELTA_V_SIGNIFICANCY_THRESHOLD)
            {
                transfer_efficiency = (startEfficiency * dv_start + dv_end) / total_dv;
            }
            else
            {
                transfer_efficiency = 1.0; // not significant
            }

            return;
        }
    }

    public bool isValid()
    {
        return (transferOrbit != null);
    }

    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.ANAIS_TRANSFER, level, message);
    }
}
