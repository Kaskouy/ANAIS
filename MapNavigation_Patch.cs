using System.Runtime.CompilerServices;

using HarmonyLib;
using SFS.World;
using SFS.WorldBase;
using SFS.World.Maps;
using UnityEngine;
using UnityEngine.UI;


// The 2 following patches allow to reset the latest approach data when the player changes the target
// --------------------------------------------------------------------------------------------------
[HarmonyPatch(typeof(MapNavigation), nameof(MapNavigation.SetTarget))]
public class MapNavigation_SetTarget_Patch
{
	[HarmonyPostfix]
	public static void SetTarget_postfix()
    {
		ClosestApproachCalculator.resetLastApproachDataValidity();
		ClosestApproachCalculator.resetLastApproachDataValidity_MultiTurn();
		ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode1();
		ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode2();

		AnaisTransferCalculator.resetMemorizedTimeOfArrival();

#if ACTIVE_LOGS
		AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, LOG_LEVEL.INFO, "Detected target set to an object");
#endif
	}
}

[HarmonyPatch(typeof(MapNavigation), nameof(MapNavigation.ToggleTarget))]
public class MapNavigation_ToggleTarget_Patch
{
	[HarmonyPostfix]
	public static void ToggleTarget_postfix()
	{
		ClosestApproachCalculator.resetLastApproachDataValidity();
		ClosestApproachCalculator.resetLastApproachDataValidity_MultiTurn();
		ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode1();
		ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode2();

		AnaisTransferCalculator.resetMemorizedTimeOfArrival();

#if ACTIVE_LOGS
		AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, LOG_LEVEL.INFO, "Detected target set to another object");
#endif
	}
}

// This patch allows to reset the latest approach data when the player switches to another ship
// --------------------------------------------------------------------------------------------
[HarmonyPatch(typeof(MapNavigation), "OnPlayerChange")]
public class MapNavigation_OnPlayerChange_Patch
{
	[HarmonyPostfix]
	public static void OnPlayerChange_postfix(Player newPlayer)
	{
		ClosestApproachCalculator.setPlayer(newPlayer);
		ClosestApproachCalculator.resetLastApproachDataValidity();
		ClosestApproachCalculator.resetLastApproachDataValidity_MultiTurn();
		ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode1();
		ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode2();

		AnaisTransferCalculator.resetMemorizedTimeOfArrival();

#if ACTIVE_LOGS
		AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, LOG_LEVEL.INFO, "Detected control set to another ship");
#endif
	}
}


// Patch for the DrawNavigation method: calculates and display the closest approach line
// -------------------------------------------------------------------------------------
[HarmonyPatch(typeof(MapNavigation), nameof(MapNavigation.DrawNavigation))]
public class MapNavigation_Patch
{
	[HarmonyPrefix]
	public static bool DrawNavigation_prefix(MapNavigation __instance, SelectableObject ___target)
    {
		Player value = PlayerController.main.player.Value;
		if (value == null || value.mapPlayer == null)
		{
			return false; // nothing to do anyway
		}

		if (___target == null)
		{
			return false; // no target
		}

		if ( !GetOrbits(value.mapPlayer.Trajectory, out var list_PlayerOrbits) || !GetOrbits(___target.Trajectory, out var list_targetOrbits) )
		{
			LOG(LOG_LEVEL.DEBUG, "No Orbit could be retrieved for either the player or the target");
			return false; // no target
		}

		Orbit targetOrbit = list_targetOrbits[0];
		Orbit playerOrbit = list_PlayerOrbits[0];
		bool imminentEncounterPlanned = false;

		// If target is a planet, identify the concerned planet
		MapPlanet targetMapPlanet = ___target as MapPlanet;
		Planet targetPlanet = null;

		if(targetMapPlanet != null)
        {
			targetPlanet = targetMapPlanet.planet;
		}

		// Check if we can use ANAIS
		bool useANAIS = AnaisTransferCalculator.isAnaisTransferCalculationDoable(playerOrbit, targetOrbit, targetPlanet);

		bool AnaisUsedForLocalTransfer = false;
		bool AnaisUsedForInterplanetaryTransfer = false;


		ClosestApproachCalculator.T_ApproachData approachData = new ClosestApproachCalculator.T_ApproachData();
		approachData.validity = false;

		// Try to calculate a transfer
		// ---------------------------
		if (useANAIS)
        {
			// Calculate approach over the next 60 seconds
			ClosestApproachCalculator.T_ApproachData approachData_shortPeriod = new ClosestApproachCalculator.T_ApproachData();
			approachData_shortPeriod.validity = false;

			if (playerOrbit.Planet == targetOrbit.Planet)
            {
				approachData_shortPeriod = ClosestApproachCalculator.CalculateClosestApproachOnShortPeriod(playerOrbit, targetOrbit, WorldTime.main.worldTime, WorldTime.main.worldTime + 60.0);
			}

			if(approachData_shortPeriod.validity && (approachData_shortPeriod.dist < 10000.0))
            {
				// Player and target are about to encounter very shortly
				// -----------------------------------------------------
				imminentEncounterPlanned = true;

				Color textColor = Drawing_Utils.GetEfficiencyTransferColor(1.0);
				Color lineColor = textColor;
				lineColor.a = 0.9f;

				string speedText = Units.ToVelocityString(ClosestApproachCalculator.GetApproachSpeed(approachData_shortPeriod), true);
				double durationBeforeApproach = approachData_shortPeriod.date - WorldTime.main.worldTime;
				string timerText;

				if (durationBeforeApproach < 1.0)
				{
					timerText = " (T-0s)";
				}
				else
				{
					timerText = " (T-" + Units.ToTimestampString(durationBeforeApproach, true, false) + ")";
				}

				string encounterLabel = "Encounter: ΔV = " + speedText + timerText;

				MapDrawer.DrawPointWithText(15, textColor, encounterLabel, 40, textColor, MapDrawer.GetPosition(approachData_shortPeriod.locTarget), approachData_shortPeriod.locTarget.position.normalized, 4, 4);
			}
			else
            {
				// Time for ANAIS to shine!
				// ------------------------
				AnaisTransfer anaisTransfer = AnaisTransferCalculator.calculateTransferToTarget(playerOrbit, targetOrbit, targetPlanet);

				if ((anaisTransfer != null) && anaisTransfer.isValid())
				{
					Orbit transferOrbit = anaisTransfer.transferOrbit;
					Orbit ejectionOrbit = anaisTransfer.childTransferOrbit;

					Location locationStart = playerOrbit.GetLocation(WorldTime.main.worldTime);
					Location locationEnd = transferOrbit.GetLocation(transferOrbit.orbitEndTime);

					string speedText = Units.ToVelocityString(anaisTransfer.dv_start, true);
					string startLabel = "Transfer: ΔV = " + speedText;
					speedText = Units.ToVelocityString(anaisTransfer.dv_end, true);
					string endLabel;

					if(targetPlanet == null)
                    {
						endLabel = "Encounter: ΔV = " + speedText;
					}
					else
                    {
						endLabel = "Insertion low orbit: ΔV = " + speedText;
					}

					Color textColor = Drawing_Utils.GetEfficiencyTransferColor(anaisTransfer.transfer_efficiency);
					Color lineColor = textColor;
					lineColor.a = 0.9f;

					if(anaisTransfer.dv_start > 0.1)
                    {
						MapDrawer.DrawPointWithText(15, textColor, startLabel, 40, textColor, MapDrawer.GetPosition(locationStart), locationStart.position.normalized, 4, 4);
					}
					
					MapDrawer.DrawPointWithText(15, textColor, endLabel, 40, textColor, MapDrawer.GetPosition(locationEnd), locationEnd.position.normalized, 4, 4);
					transferOrbit.DrawDashed(drawStats: false, drawStartText: false, drawEndText: false, lineColor);

					if(ejectionOrbit != null)
                    {
						ejectionOrbit.DrawDashed(drawStats: false, drawStartText: false, drawEndText: false, lineColor);
						AnaisUsedForInterplanetaryTransfer = true;
					}
					else
                    {
						AnaisUsedForLocalTransfer = true;
					}
				}
				else
				{
					LOG(LOG_LEVEL.INFO, "Anais transfer calculation failed - The closest approach line will be used as a backup");
					useANAIS = false; // to allow the use of the closest approach line if ANAIS failed
				}
			}
		}



		// Search a trajectory that makes an encounter possible with the target (for a potential closest approach calculation)
		// --------------------------------------------------------------------
		targetOrbit = list_targetOrbits[0]; // redefine those
		playerOrbit = null;
		bool skipBlueLine = false;
		for (int i_orbit = 0; i_orbit < list_PlayerOrbits.Length; i_orbit++)
		{
			// search for the first trajectory that orbits the same body as target
			if (list_PlayerOrbits[i_orbit].Planet == targetOrbit.Planet)
			{
				// check if the player object doesn't already have an encounter planned with the targeted object (if a planet is targeted...)
				if (!((i_orbit < list_PlayerOrbits.Length - 1) && (list_PlayerOrbits[i_orbit + 1].Planet.mapPlanet == ___target)))
				{
					playerOrbit = list_PlayerOrbits[i_orbit];
					skipBlueLine = imminentEncounterPlanned || (AnaisUsedForLocalTransfer && (i_orbit == 0)) || (AnaisUsedForInterplanetaryTransfer && (i_orbit == 1));
					break;
				}
			}
		}

		if (playerOrbit == null)
		{
			return false; // Nothing we can do
		}

		

		if (!skipBlueLine) // If Anais was already used to deal with the situation, don't show the blue line
        {
			// Calculate the closest approach
			// ------------------------------
			approachData = ClosestApproachCalculator.CalculateClosestApproach(playerOrbit, targetOrbit);

			if (approachData.validity)
			{
				// If a valid approach is found, display it
				// ----------------------------------------
				Color lineColor = Drawing_Utils.GetClosestApproachLineColor(approachData.dist);
				string closestApproachText = Drawing_Utils.GetClosestApproachText(approachData.dist, ClosestApproachCalculator.GetApproachSpeed(approachData), approachData.date);

				Drawing_Utils.DrawDashedLine(playerOrbit, approachData.locPlayer, approachData.locTarget, lineColor, null, closestApproachText);
			}
		}


		// Calculate closest approach on several turns
		// -------------------------------------------
		if (!(imminentEncounterPlanned || ((targetPlanet == null) && AnaisTransferCalculator.isEncounterPlanned())))
		{
			uint nbTurns1 = 0, nbTurns2 = 0;
			ClosestApproachCalculator.T_ApproachData approachData_multi1, approachData_multi2;
			ClosestApproachCalculator.CalculateClosestApproach_MultiTurn(playerOrbit, targetOrbit, approachData, out approachData_multi1, out nbTurns1, out approachData_multi2, out nbTurns2);

			Color lightGreen = new Color(0.706f, 1.0f, 0.902f, 0.8f); // light green

			// Show approach line for node 1
			if (approachData_multi1.validity && nbTurns1 > 1) // don't show it if nbTurns is 1, would be redundant with classic closest approach line
			{
				string closestApproachText = "Best approach: " + approachData_multi1.dist.ToDistanceString() + " (" + nbTurns1 + " turns)";

				Drawing_Utils.DrawDashedLine(playerOrbit, approachData_multi1.locPlayer, approachData_multi1.locTarget, lightGreen, null, closestApproachText);
			}

			// Show approach line for node 2
			if (approachData_multi2.validity && nbTurns2 > 1) // don't show it if nbTurns is 1, would be redundant with classic closest approach line
			{
				string closestApproachText = "Best approach: " + approachData_multi2.dist.ToDistanceString() + " (" + nbTurns2 + " turns)";

				Drawing_Utils.DrawDashedLine(playerOrbit, approachData_multi2.locPlayer, approachData_multi2.locTarget, lightGreen, null, closestApproachText);
			}
		}
		else
        {
			// Reset memorized number of turns
			ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode1();
			ClosestApproachCalculator.resetNbOfMemorizedTurnsOnNode2();
		}

		return false;
	}


	// copy of the GetOrbits method in the original file
	static bool GetOrbits(Trajectory a, out Orbit[] orbits)
	{
		orbits = new Orbit[a.paths.Count];
		for (int i = 0; i < a.paths.Count; i++)
		{
			if (!(a.paths[i] is Orbit orbit3))
			{
				return false;
			}
			orbits[i] = orbit3;
		}
		return orbits.Length != 0;
	}

	// Local log function
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	private static void LOG(LOG_LEVEL level, string message)
	{
#if ACTIVE_LOGS
		AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, level, message);
#endif
	}
}
