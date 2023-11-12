using System.Runtime.CompilerServices;

using HarmonyLib;
using SFS.World;
using SFS.WorldBase;
using SFS.World.Maps;
using UnityEngine;
using UnityEngine.UI;
using System.Diagnostics;



// The 2 following patches allow to reset the latest approach data when the player changes the target
// --------------------------------------------------------------------------------------------------
[HarmonyPatch(typeof(MapNavigation), nameof(MapNavigation.SetTarget))]
public class MapNavigation_SetTarget_Patch
{
	[HarmonyPostfix]
	public static void SetTarget_postfix(MapNavigation __instance)
    {
		AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, LOG_LEVEL.INFO, "SetTarget: Detected target set to a new object - Target exists: " + ((__instance.target != null) ? "yes":"no"));
		AnaisManager.notifyNewTarget(__instance.target);
	}
}


// This patch allows to reset the latest approach data when the player switches to another ship
// --------------------------------------------------------------------------------------------
[HarmonyPatch(typeof(MapNavigation), "OnPlayerChange")]
public class MapNavigation_OnPlayerChange_Patch
{
	[HarmonyPostfix]
	public static void OnPlayerChange_postfix(MapNavigation __instance, Player newPlayer)
	{
        AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, LOG_LEVEL.INFO, "OnPlayerChange: Detected control set to another ship - Target exists: " + ((__instance.target != null) ? "yes" : "no"));
        AnaisManager.notifyNewTarget(__instance.target);
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

		// Send input data to the ANAIS manager for calculations
		AnaisManager.SetInputData(value.mapPlayer, ___target);

        // Draw the available transfers after the calculations made by ANAIS Manager
		// (Note that it won't be the transfer from the data transmitted above as there's not enough time for the calculation to run between those instructions)
        AnaisManager.DrawAnaisTransfer();

		// Don't run the original method
		return false;
	}


	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.MAP_NAVIGATION, level, message);
	}
}


