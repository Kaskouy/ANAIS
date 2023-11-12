using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using HarmonyLib;
using SFS.World;

[HarmonyPatch(typeof(WorldView), "Update")]
public class WorldView_Update_Patch
{
	[HarmonyPostfix]
	public static void Update_postfix()
	{
		AnaisManager.setIsGameRunning();

		ANAIS_Panel.ANAIS_PanelEventListener.CheckEvent();
    }
}
