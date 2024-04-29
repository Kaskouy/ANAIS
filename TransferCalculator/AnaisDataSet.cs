using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using SFS.World.Maps;
using SFS.World;
using SFS.WorldBase;

using UnityEngine;


public class AnaisDataSet
{
    // ----------------------------------------------------------------
    // --------                                                --------
    // --------                 ATTRIBUTES                     --------
    // --------                                                --------
    // ----------------------------------------------------------------

    // Input data
    // ----------
    private bool _hasInputData = false;
    private double _timeNow = double.NegativeInfinity;
    private List<Orbit> _listPlayerOrbits = new List<Orbit>();
    private List<Orbit> _listTarget = new List<Orbit>();
    private Planet _targetPlanet = null;

    private ANAIS_Settings _settings = new ANAIS_Settings(); // the player settings from the ANAIS control panel
    private NavigationVariables _NavVariables = new NavigationVariables(); // Used to store some temp variables transmitted to the calculation algorithms

    private double timeWarpFactor = 1.0;

    // Output data
    // -----------
    public bool _finalApproachMode = false;
    ClosestApproachCalculator.T_ApproachData _approachData_finalApproach = new ClosestApproachCalculator.T_ApproachData(dummy: true);
    private AnaisTransfer _anaisTransfer;
    ClosestApproachCalculator.T_ApproachData _approachData_singleTurn = new ClosestApproachCalculator.T_ApproachData(dummy: true);
    ClosestApproachCalculator.T_ApproachData _approachData_multiTurn = new ClosestApproachCalculator.T_ApproachData(dummy: true);
    uint _nbTurns = 0;
    bool _entryInTargetSOIdetected = false;

    // Internal constants
    // ------------------
    private static Color C_LightBlue_Color = new Color(0.627f, 0.784f, 1.0f, 0.8f);
    private static Color C_LightGreen_Color = new Color(0.706f, 1.0f, 0.902f, 0.8f);

    // ----------------------------------------------------------------
    // --------                                                --------
    // --------                   METHODS                      --------
    // --------                                                --------
    // ----------------------------------------------------------------


    // ----------------------------------------------------------------
    // ------                  SetInputData                      ------
    // ----------------------------------------------------------------
    // This method allows to extract and store all the data needed in
    // prevision for the transfer calculations.
    // ----------------------------------------------------------------
    public void SetInputData(MapPlayer mapPlayer, SelectableObject target)
    {
        LOG(LOG_LEVEL.DEBUG, "SetInputData called");

        // Both player and target must be defined
        if ((mapPlayer == null) || (target == null))
        {
            return;
        }

        // Start with a reset to cleanse some previously memorized data
        Reset();

        // Set current time
        _timeNow = WorldTime.main.worldTime;

        // If target is a planet, memorize it
        // ----------------------------------
        MapPlanet targetMapPlanet = target as MapPlanet;
        if (targetMapPlanet != null) _targetPlanet = targetMapPlanet.planet;

        // Retrieve orbits lists
        // ---------------------
        _listPlayerOrbits = Orbit_Utils.GetListOrbits(mapPlayer.Trajectory);
        _listTarget = Orbit_Utils.GetListOrbits(target.Trajectory);

        // Make a copy of the player settings from the ANAIS control panel
        // ----------------------------------
        _settings.copySettings(ANAIS_Panel._settings);

        // Make a copy of some temp variables for ANAIS, and reset them
        // ----------------------------------
        _NavVariables.Copy(ANAIS_Panel._NavVariables);
        ANAIS_Panel._NavVariables.Reset();

        // time warp index
        timeWarpFactor = SFS.World.WorldTime.main.timewarpSpeed;

        _hasInputData = true;
    }


    // ----------------------------------------------------------------
    // ------               CalculateTransfer                    ------
    // ----------------------------------------------------------------
    // This method runs all the transfers and approach lines
    // calculations if input data is available.
    // Returns true if the calculation occured, false otherwise.
    // WARNING: This method is to be called by the separate ANAIS task
    //     --> No global unprotected data must be used in it!
    // ----------------------------------------------------------------
    public bool CalculateTransfer(ref bool isCurrentlyInFinalApproachMode, ref bool iscurrentlyOnEncounterTrajectory, ref double encounterDate, ref double preferredTimeOfArrivalAtNode1, ref double preferredTimeOfArrivalAtNode2, ref bool ANAIScalculationAllowed)
    {
        // If no input data, leave.
        if (!_hasInputData) return false;

        // make sure there's a player orbit
        if (_listPlayerOrbits.Count == 0)
        {
            return false;
        }

        // Reset preferred times of arrival at nodes if that was asked (if player changed the number of turns through the control panel)
        if(_NavVariables._approachLinesPreferredDatesNeedReset)
        {
            preferredTimeOfArrivalAtNode1 = double.NegativeInfinity;
            preferredTimeOfArrivalAtNode2 = double.NegativeInfinity;
        }

        double targetAltitude = _NavVariables._targetAltitude;

        // RETRIEVE PLAYER/TARGET ORBIT + Evaluate transfer configuration
        // ----------------------------
        Orbit playerOrbit = _listPlayerOrbits[0]; // the current player orbit
        Orbit targetOrbit = null; // the target orbit - will always exist EXCEPT if the target is the main body (the Sun)

        bool isLocalTransfer = false; // This will be true if player and target orbit the same body, or if target is the orbited body itself (the transfer occurs in the current SOI)
        bool isInterplanetaryTransfer = false; // This will be true if the body orbited by the player and the target orbits the same body, or if target is the parent body of the orbited body (the transfer makes the player exit the current SOI and get an encounter in the parent SOI)
        bool isReturnToMotherPlanet = false; // This will be true if the player aims for the orbited body or its parent (we aim for low orbit of such body)
        bool isRendezVousTrajectory = false; // This will be true if the trajectory calculated will lead to an encounter
        
        if ((_targetPlanet != null) && (playerOrbit.Planet == _targetPlanet))
        {
            // Return to mother planet through local transfer
            isReturnToMotherPlanet = true;
            isLocalTransfer = true;

            // Retrieve target orbit if exists (it won't exist if the targeted body is the main body - aka the Sun)
            if (_listTarget.Count > 0) targetOrbit = _listTarget[0];

            LOG(LOG_LEVEL.DEBUG, "  Transfer configuration: return to parent planet (local transfer)");
        }
        else if ((_targetPlanet != null) && (playerOrbit.Planet.parentBody == _targetPlanet))
        {
            // Return to mother planet through interplanetary transfer
            isReturnToMotherPlanet = true;
            isInterplanetaryTransfer = true;

            // Retrieve target orbit if exists (it won't exist if the targeted body is the main body - aka the Sun)
            if (_listTarget.Count > 0) targetOrbit = _listTarget[0];

            LOG(LOG_LEVEL.DEBUG, "  Transfer configuration: return to parent planet (interplanetary transfer)");
        }
        else if (_listTarget.Count > 0)
        {
            // Not a return to mother planet configuration --> Check if it's a rendez-vous trajectory 
            targetOrbit = _listTarget[0];
            
            if (playerOrbit.Planet == targetOrbit.Planet) isLocalTransfer = true;
            else if (playerOrbit.Planet.parentBody == targetOrbit.Planet) isInterplanetaryTransfer = true;
            
            isRendezVousTrajectory = isLocalTransfer || isInterplanetaryTransfer;

            if (isLocalTransfer || isInterplanetaryTransfer)
            {
                LOG(LOG_LEVEL.DEBUG, "  Transfer configuration: Rendez-vous; Type: " + (isLocalTransfer ? "Local" : "Interplanetary") );
            }
            else
            {
                LOG(LOG_LEVEL.DEBUG, "  Transfer configuration: Degraded (ANAIS can't be used)");
            }
        }
        else
        {
            // Not managed situation
            LOG(LOG_LEVEL.DEBUG, "  Transfer configuration: Not managed");
            return false;
        }
        

        // CHECK IF WE ARE IN FINAL APPROACH MODE
        // --------------------------------------
        _finalApproachMode = false;

        if (isLocalTransfer && isRendezVousTrajectory) // This is only for a local transfer (the targeted object is in the same SOI as player)
        {
            // Search if we have a minimal approach in the next 60 seconds
            double endTime = _timeNow + AnaisManager.C_TIME_THRESHOLD_FOR_FINAL_ENCOUNTER_MODE;
            _approachData_finalApproach = ClosestApproachCalculator.CalculateClosestApproachOnShortPeriod(playerOrbit, targetOrbit, _timeNow, endTime);

            if ((_approachData_finalApproach.validity == true) && (_approachData_finalApproach.dist < AnaisManager.C_DISTANCE_THRESHOLD_FOR_FINAL_ENCOUNTER_MODE))
            {
                // There's a closest approach at less than 5000 meters in the next 60 seconds => We are in final approach mode
                _finalApproachMode = true;
            }
            else if (isCurrentlyInFinalApproachMode)
            {
                // There's no longer such a closest approach, if we were previously in final approach mode we stay in this mode provided we are still close to the target
                Location playerLocation = playerOrbit.GetLocation(_timeNow);
                Location targetLocation = targetOrbit.GetLocation(_timeNow);
                double distanceToTarget = (targetLocation.position - playerLocation.position).magnitude;

                if (distanceToTarget < AnaisManager.C_DISTANCE_THRESHOLD_FOR_FINAL_ENCOUNTER_MODE)
                {
                    _finalApproachMode = true;
                }
            }
            //else we are not in final approach mode

            if (!_finalApproachMode)
            {
                _approachData_finalApproach.validity = false; // To avoid using accidentally inappropriate data (if a closest approach was found but at a too large distance)
            }
        }

        if (_finalApproachMode != isCurrentlyInFinalApproachMode)
        {
            LOG(LOG_LEVEL.INFO, "Final approach mode switched to " + _finalApproachMode);
            isCurrentlyInFinalApproachMode = _finalApproachMode;
        }


        // CHECK ENCOUNTER EVENTS
        // ----------------------
        // This is to check if an encounter is found and if a gravity assist is used
        if ((!_finalApproachMode) && (targetOrbit != null) && isRendezVousTrajectory )
        {
            // Search the first orbit that can potentially lead to an encounter with target (this orbit and the target's one must have the same parent body)
            // We search it that way in case some unexpected bodies went in our way, so that we don't consider them (example: we target Venus, we're on a transfer from Earth to Venus, but the Moon comes on our path - this allows to ignore the Moon)
            int i_firstCompatibleTransfer = 0;
            bool gravityAssistDetected = false;

            while (i_firstCompatibleTransfer < _listPlayerOrbits.Count)
            {
                if (_listPlayerOrbits[i_firstCompatibleTransfer].Planet == targetOrbit.Planet) break;
                else i_firstCompatibleTransfer++;
            }

            // Loop through all orbits from the first compatible one to review all encountered bodies
            for(int i_playerOrbit = i_firstCompatibleTransfer; i_playerOrbit < _listPlayerOrbits.Count; i_playerOrbit++)
            {
                if(_listPlayerOrbits[i_playerOrbit].Planet == targetOrbit.Planet) // Orbit is compatible with an encounter with target or a "sister" planet...
                {
                    // ... but doesn't lead to an encounter; We either exit the SOI or turn eternally around the parent body, no more encounters can be found
                    if (_listPlayerOrbits[i_playerOrbit].PathType != PathType.Encounter) break;
                }
                else
                {
                    // We are presumably in the target's SOI or one of its "sister"'s SOI
                    if(_listPlayerOrbits[i_playerOrbit].Planet == _targetPlanet)
                    {
                        // An entry into target SOI is detected
                        _entryInTargetSOIdetected = true;
                        break; // we have all we wanted to know
                    }
                    else
                    {
                        // An entry into a SOI that's not the target's one has been detected
                        gravityAssistDetected = true; // Especially, no previous encounter with the target occured - we would have breaked otherwise
                    }
                }
            }

            // If ANAIS calculation was restricted, we only maintain the restriction if another body is on our way
            // We suppose it's the player's intention to use that body as a gravity assist opportunity, so calculating the ANAIS transfer would disturb him.
            //if(!ANAIScalculationAllowed && !gravityAssistDetected) ANAIScalculationAllowed = true;
            ANAIScalculationAllowed = true; // Forced to true since the player now has the control over this through the control panel
        }
        else
        {
            if (!_finalApproachMode) ANAIScalculationAllowed = true; // ANAIS is allowed to run in any other situation (except in final approach mode)
        }


        // CHECK IF WE ARE ON ENCOUNTER MODE (Player is practically on an encounter trajectory)
        // ---------------------------------
        if (!_finalApproachMode && isRendezVousTrajectory && iscurrentlyOnEncounterTrajectory && (encounterDate < targetOrbit.orbitEndTime) && ANAIScalculationAllowed && _settings._showTransfer)
        {
            // calculate the transfer for the memorized encounter date
            _anaisTransfer = new AnaisTransfer(playerOrbit, targetOrbit, _targetPlanet, targetAltitude, _settings._transferType);
            _anaisTransfer.calculateTransfer(_timeNow, encounterDate);

            // check if we are still on an encounter trajectory - if not, we exit that mode
            if ((_anaisTransfer == null) || !_anaisTransfer.isValid() || (_anaisTransfer.dv_start > AnaisManager.C_DELTAV_THRESHOLD_EXIT_ENCOUNTER_MODE))
            {
                LOG(LOG_LEVEL.INFO, "Encounter trajectory mode switched to false");
                iscurrentlyOnEncounterTrajectory = false;
                encounterDate = -double.NegativeInfinity;
                _anaisTransfer = null;
            }
        }
        else
        {
            iscurrentlyOnEncounterTrajectory = false; // If we are in final approach mode or not on a rendez-vous trajectory type, we are not on encounter mode...
            encounterDate = -double.NegativeInfinity;
        }

        // CALCULATE ANAIS TRANSFER (if not in final approach / on encounter trajectory mode, and if the configurations is compatible with ANAIS)
        // ------------------------
        if (!_finalApproachMode && !iscurrentlyOnEncounterTrajectory && (isRendezVousTrajectory || isReturnToMotherPlanet) && ANAIScalculationAllowed && _settings._showTransfer)
        {
            // Calculate the ANAIS transfer
            _anaisTransfer = AnaisTransferCalculator.calculateTransferToTarget(playerOrbit, targetOrbit, _targetPlanet, targetAltitude, _timeNow, _settings._transferType, timeWarpFactor);

            // If the calculation is successful, check if we are on encounter mode
            if ((_anaisTransfer != null) && _anaisTransfer.isValid() && isRendezVousTrajectory && (_anaisTransfer.dv_start < AnaisManager.C_DELTAV_THRESHOLD_ENTER_ENCOUNTER_MODE))
            {
                LOG(LOG_LEVEL.INFO, "Encounter trajectory mode switched to true");
                iscurrentlyOnEncounterTrajectory = true;
                encounterDate = _anaisTransfer.arrivalTime;
            }
        }


        // MANAGE CLOSEST APPROACH LINES
        // -----------------------------
        if (!_finalApproachMode && (targetOrbit != null))
        {
            int i_firstCompatibleTransfer = 0;

            // Search for the first compatible orbit
            while (i_firstCompatibleTransfer < _listPlayerOrbits.Count)
            {
                if (_listPlayerOrbits[i_firstCompatibleTransfer].Planet == targetOrbit.Planet) break;
                else i_firstCompatibleTransfer++;
            }

            // Browse the following orbits, and calculate the approach lines on the compatible ones
            for (int i_playerOrbit = i_firstCompatibleTransfer; i_playerOrbit < _listPlayerOrbits.Count; i_playerOrbit++)
            {
                // If we find an encounter with the target planet (if it's a planet...), reset the single turn approach data (we want the best after the last known encounter)
                if (_listPlayerOrbits[i_playerOrbit].Planet == _targetPlanet)
                {
                    _approachData_singleTurn.validity = false;
                }

                // if this orbit is around the same body as target, calculate the approach lines
                if (_listPlayerOrbits[i_playerOrbit].Planet == targetOrbit.Planet)
                {
                    // Check if that trajectory is handled by ANAIS
                    bool isHandledByANAIS = false;
                    if ((_anaisTransfer != null) && _anaisTransfer.isValid())
                    {
                        if (isLocalTransfer && (i_playerOrbit == 0)) isHandledByANAIS = true;
                        else if (isInterplanetaryTransfer && (i_playerOrbit == 1)) isHandledByANAIS = true;
                    }

                    // Calculate the closest approach - if this trajectory is not already handled by ANAIS
                    // ------------------------------
                    if(!isHandledByANAIS && (_settings._nbMaxTurns > 0))
                    {
                        ClosestApproachCalculator.T_ApproachData currentApproachData = ClosestApproachCalculator.CalculateClosestApproach(_listPlayerOrbits[i_playerOrbit], targetOrbit, Math.Max(_timeNow, _listPlayerOrbits[i_playerOrbit].orbitStartTime));

                        // If it's the best found until now, memorize it
                        if (currentApproachData.validity)
                        {
                            if(!_approachData_singleTurn.validity || (_approachData_singleTurn.dist > currentApproachData.dist))
                            {
                                _approachData_singleTurn = currentApproachData;
                            }
                        }
                    }

                    // If this orbit escapes the parent body or turns eternally around it, this is the last one we consider -> we calculate the green lines on this one
                    if (_listPlayerOrbits[i_playerOrbit].PathType != PathType.Encounter)
                    {
                        // Calculate the closest approach over several turns (skip if this is the trajectory handled by ANAIS and we are close to getting an encounter)
                        if (!isHandledByANAIS || !iscurrentlyOnEncounterTrajectory)
                        {
                            ClosestApproachCalculator.CalculateClosestApproach_MultiTurn(_listPlayerOrbits[i_playerOrbit], targetOrbit, _approachData_singleTurn, Math.Max(_timeNow, _listPlayerOrbits[i_playerOrbit].orbitStartTime), out _approachData_multiTurn, out _nbTurns, ref preferredTimeOfArrivalAtNode1, ref preferredTimeOfArrivalAtNode2, _settings._nbMaxTurns, _settings._preferredNode);
                        }

                        break;
                    }
                }
            }
        }

        // Reset the closest approach lines data if nothing was calculated
        if (_approachData_multiTurn.validity == false)
        {
            preferredTimeOfArrivalAtNode1 = double.NegativeInfinity;
            _nbTurns = 0;
        }

        // This is to return if the data has been calculated or not
        return true;
    }


    // ----------------------------------------------------------------
    // ------                  DrawTransfer                      ------
    // ----------------------------------------------------------------
    // This method allows to display the transfers and the approach
    // lines if they were calculated.
    // ----------------------------------------------------------------
    public void DrawTransfer()
    {
        if (_finalApproachMode)
        {
            if (_approachData_finalApproach.validity)
            {
                // We are in final approach mode and an encounter is detected in the next 60 seconds --> show the encounter text with the ΔV and the remaining time
                // ---------------------------------------------------------------------------------
                Color textColor = Drawing_Utils.GetEfficiencyTransferColor(1.0);
                string speedText = Units.ToVelocityString(ClosestApproachCalculator.GetApproachSpeed(_approachData_finalApproach), true);

                string timerText;
                double durationBeforeApproach = _approachData_finalApproach.date - WorldTime.main.worldTime;

                if (durationBeforeApproach < 1.0)
                {
                    timerText = " (T-0s)";
                }
                else
                {
                    timerText = " (T-" + Units.ToTimestampString(durationBeforeApproach, true, false) + ")";
                }

                string endLabel;
                endLabel = "Encounter: ΔV = " + speedText + timerText;

                MapDrawer.DrawPointWithText(15, textColor, endLabel, 40, textColor, MapDrawer.GetPosition(_approachData_finalApproach.locPlayer), _approachData_finalApproach.locPlayer.position.normalized, 4, 4);
            }
            else
            {
                // We are in final approach mode but no encounter has been detected --> show the distance to the target
                // ----------------------------------------------------------------
                Location playerLocation = _listPlayerOrbits[0].GetLocation(WorldTime.main.worldTime);
                Location targetLocation = _listTarget[0].GetLocation(WorldTime.main.worldTime);
                Double2 position = targetLocation.position - playerLocation.position;
                double distance = position.magnitude;
                string label = Units.ToDistanceString(distance);
                label = "Distance: " + label;

                Color textColor = Drawing_Utils.GetEfficiencyTransferColor(1.0);

                MapDrawer.DrawPointWithText(15, textColor, label, 40, textColor, MapDrawer.GetPosition(targetLocation), position.normalized, 4, 4);
            }
        }
        else if ((_anaisTransfer != null) && _anaisTransfer.isValid())
        {
            // There's an ANAIS transfer planned --> show the transfer
            // ---------------------------------
            Orbit transferOrbit = _anaisTransfer.transferOrbit;
            Orbit ejectionOrbit = _anaisTransfer.childTransferOrbit;

            Location locationStart = _anaisTransfer.originOrbit.GetLocation(_anaisTransfer.departureTime);

            // start label
            string speedText = Units.ToVelocityString(_anaisTransfer.dv_start, true);
            string startLabel;

            if ((_anaisTransfer.transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS) && (_anaisTransfer.targetPlanet != null) && (transferOrbit.PathType == PathType.Eternal)) 
            {
                // Special label for an insertion in orbit (only if we have an eternal trajectory in rendez-vous mode)
                startLabel = "Insertion orbit: ΔV = " + speedText; 
            } 
            else { startLabel = "Transfer: ΔV = " + speedText; }

            // end label
            speedText = Units.ToVelocityString(_anaisTransfer.dv_end, true);
            string endLabel;

            if (_anaisTransfer.targetPlanet == null)
            {
                endLabel = "Encounter: ΔV = " + speedText;
            }
            else
            {
                endLabel = "Insertion orbit: ΔV = " + speedText;
            }

            // transfer color
            Color textColor = Drawing_Utils.GetEfficiencyTransferColor(_anaisTransfer.transfer_efficiency);
            Color lineColor = textColor;
            lineColor.a = 0.9f;

            if ((_anaisTransfer.dv_start > 0.1) && (_anaisTransfer.areDeltaV_valuesSignificant()))
            {
                MapDrawer.DrawPointWithText(15, textColor, startLabel, 40, textColor, MapDrawer.GetPosition(locationStart), locationStart.position.normalized, 4, 4);
            }
            else
            {
                MapDrawer.DrawPoint(15, textColor, MapDrawer.GetPosition(locationStart), 4, true, 4);
            }

            if(transferOrbit.PathType == PathType.Encounter)
            {
                Location locationEnd = transferOrbit.GetLocation(transferOrbit.orbitEndTime);

                if ((_settings._transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS) && (_anaisTransfer.areDeltaV_valuesSignificant()))
                {
                    // Rendez-vous mode: Show the arrival point and the arrival delta-V
                    MapDrawer.DrawPointWithText(15, textColor, endLabel, 40, textColor, MapDrawer.GetPosition(locationEnd), locationEnd.position.normalized, 4, 4);
                }
                else
                {
                    // Fly-by mode (or non significant deltaV values due to high time-warp mode): only show the arrival point
                    MapDrawer.DrawPoint(15, textColor, MapDrawer.GetPosition(locationEnd), 4, true, 4);
                }
            }
            
            transferOrbit.DrawDashed(drawStats: false, drawStartText: false, drawEndText: false, lineColor);

            if (ejectionOrbit != null)
            {
                ejectionOrbit.DrawDashed(drawStats: false, drawStartText: false, drawEndText: false, lineColor);
            }

            // If a planet is targeted in rendez-vous-mode, render the targeted orbit - supposing we are not already at the targeted altitude (in this case, the transfer orbit calculated is the final orbit, hence the third condition)
            if((_anaisTransfer.transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS) && (_anaisTransfer.targetPlanet != null) && (transferOrbit.PathType == PathType.Encounter))
            {
                double finalOrbitRadius = _anaisTransfer.targetPlanet.Radius + _anaisTransfer.targetAltitude;
                Orbit finalOrbit = Orbit_Utils.CreateOrbit(finalOrbitRadius, 0.0, 0.0, -1, _anaisTransfer.targetPlanet, PathType.Eternal, null);
                finalOrbit.periapsisPassageTime = SFS.World.WorldTime.main.worldTime;
                finalOrbit.orbitStartTime = finalOrbit.periapsisPassageTime;

                finalOrbit.DrawDashed(drawStats: false, drawStartText: false, drawEndText: false, lineColor);
            }
        }

        // CLOSEST APPROACH LINES
        // ----------------------
        if (!_finalApproachMode)
        {
            string closestApproachText;

            // Approach line on a single turn (blue line)
            if (_approachData_singleTurn.validity)
            {
                if ((_targetPlanet != null) && (_approachData_singleTurn.dist < _targetPlanet.SOI))
                {
                    closestApproachText = "Encounter (next pass)";
                }
                else
                {
                    closestApproachText = "Best approach (next pass)";
                }

                Drawing_Utils.DrawDashedLine(_approachData_singleTurn.locPlayer.planet, _approachData_singleTurn.locPlayer, _approachData_singleTurn.locTarget, C_LightBlue_Color, null, closestApproachText);
            }

            // Approach line on multi-turns
            if (_approachData_multiTurn.validity)
            {
                Color lineColor;

                if ((_targetPlanet != null) && (_approachData_multiTurn.dist < _targetPlanet.SOI))
                {
                    closestApproachText = "Encounter";
                }
                else
                {
                    closestApproachText = "Best approach";
                }

                if(_nbTurns > 1)
                {
                    closestApproachText += " (" + _nbTurns + " turns)";
                    lineColor = C_LightGreen_Color;
                }
                else
                {
                    closestApproachText += " (next pass)";
                    lineColor = C_LightBlue_Color;
                }

                Drawing_Utils.DrawDashedLine(_approachData_multiTurn.locPlayer.planet, _approachData_multiTurn.locPlayer, _approachData_multiTurn.locTarget, lineColor, null, closestApproachText);
            }
        }
    }


    // ----------------------------------------------------------------
    // ------              getVelocityArrowData                  ------
    // ----------------------------------------------------------------
    // This method allows to retrieve the needed data to display the
    // appropriate velocity arrows in ship mode.
    // ----------------------------------------------------------------
    // Returns true if there's data available, false otherwise.
    // - isInFinalApproachMode will always indicate the correct value
    // - approachData will only be filled in in final approach mode (it
    //   can still be invalid if the target is close but has no
    //   encounter found over the next 60 seconds)
    // - startingVelocity will be the burn vector if a transfer is
    //   planned (only if isInFinalApproachMode is false)
    // - entrySOIdetected indicates if the ship is on a trajectory that
    //   makes it enter the targetted SOI - always valid; will be false
    //   if the target is another ship (because no SOI...)
    public bool getVelocityArrowData(out bool isInFinalApproachMode, out ClosestApproachCalculator.T_ApproachData approachData, out Double2 startingVelocity, out bool entrySOIdetected)
    {
        startingVelocity = new Double2();
        approachData = new ClosestApproachCalculator.T_ApproachData();
        approachData.validity = false;
        isInFinalApproachMode = _finalApproachMode;
        entrySOIdetected = _entryInTargetSOIdetected;

        if (_finalApproachMode)
        {
            if (_approachData_finalApproach.validity)
            {
                approachData.validity = true;
                approachData.date = _approachData_finalApproach.date;
                approachData.dist = _approachData_finalApproach.dist;
                approachData.sepSpeed = _approachData_finalApproach.sepSpeed;
                approachData.locPlayer = new Location(_approachData_finalApproach.locPlayer.time, _approachData_finalApproach.locPlayer.planet, _approachData_finalApproach.locPlayer.position, _approachData_finalApproach.locPlayer.velocity);
                approachData.locTarget = new Location(_approachData_finalApproach.locTarget.time, _approachData_finalApproach.locTarget.planet, _approachData_finalApproach.locTarget.position, _approachData_finalApproach.locTarget.velocity);
            }

            return true;
        }
        else
        {
            if ((_anaisTransfer != null) && _anaisTransfer.isValid())
            {
                startingVelocity.x = _anaisTransfer.GetDeltaV_start().x;
                startingVelocity.y = _anaisTransfer.GetDeltaV_start().y;
                return true;
            }
        }

        return false;
    }


    // ----------------------------------------------------------------
    // ------                      Reset                         ------
    // ----------------------------------------------------------------
    // Completely reinitializes all data from the object.
    // ----------------------------------------------------------------
    public void Reset()
    {
        _hasInputData = false;

        // Input data
        _timeNow = -double.NegativeInfinity;
        _listPlayerOrbits.Clear();
        _listTarget.Clear();
        _targetPlanet = null;
        _settings.SetDefaultValues();
        _NavVariables.Reset();
        timeWarpFactor = 1.0;

        // Output data
        _finalApproachMode = false;
        _approachData_finalApproach.validity = false;
        _anaisTransfer = null;
        _approachData_singleTurn.validity = false;
        _approachData_multiTurn.validity = false;
        _nbTurns = 0;
        _entryInTargetSOIdetected = false;
    }


    // ----------------------------------------------------------------
    // ------                      Swap                          ------
    // ----------------------------------------------------------------
    // Static method that swaps 2 instances of AnaisDataSet
    // ----------------------------------------------------------------
    public static void Swap(ref AnaisDataSet data1, ref AnaisDataSet data2)
    {
        AnaisDataSet tempData = data1;
        data1 = data2;
        data2 = tempData;
    }


    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.ANAIS_DATA_SET, level, message);
    }
}

