using HarmonyLib;
using ModLoader;
using SFS.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using UITools;



namespace ANAIS
{
    public class Main : Mod , IUpdatable
    {
        const string C_STR_MOD_ID = "ANAIS";
        const string C_STR_MOD_NAME = "ANAIS";
        const string C_STR_AUTHOR = "Altaïr";
        const string C_STR_GAME_VERSION = "1.6.00.16";
        const string C_STR_MOD_VERSION = "v1.4.3";
        const string C_STR_MOD_DESCRIPTION = "Advanced NAvigation Innovative System\nReplaces the original navigation system with a more elaborated one.";

        private const string C_STR_CLOSEST_APPROACH_LINE_MOD_ID = "CLOSEST_APPROACH_LINE";

        public override string ModNameID => C_STR_MOD_ID;

        public override string DisplayName => C_STR_MOD_NAME;

        public override string Author => C_STR_AUTHOR;

        public override string MinimumGameVersionNecessary => C_STR_GAME_VERSION;

        public override string ModVersion => C_STR_MOD_VERSION;

        public override string Description => C_STR_MOD_DESCRIPTION;

        public override string IconLink => "https://github.com/Kaskouy/ANAIS/blob/main/Resources/LogoANAIS.png?raw=true"; // link to the logo

        // Set the dependencies
        public override Dictionary<string, string> Dependencies { get; } = new Dictionary<string, string> { { "UITools", "1.0" } };

        public Dictionary<string, FilePath> UpdatableFiles => new Dictionary<string, FilePath> { { "https://github.com/Kaskouy/ANAIS/releases/latest/download/ANAIS.dll", new FolderPath(ModFolder).ExtendToFile("ANAIS.dll") } };
        //public Dictionary<string, FilePath> UpdatableFiles => new Dictionary<string, FilePath> (); // Used to remove auto-update for dev and tests purposes


        public Main() : base()
        {
            
        }

        // This initializes the patcher. This is required if you use any Harmony patches
        public static Harmony patcher;


        public override void Load()
        {
            //UnityEngine.Debug.Log("Load called for ANAIS");
            // Tells the loader what to run when your mod is loaded

            // If closest approach line is active, this one or ANAIS must be disabled!
            if (ModsSettings.main.settings.modsActive.TryGetValue(C_STR_CLOSEST_APPROACH_LINE_MOD_ID, out bool closestApproachActive) && closestApproachActive)
            {
                SFS.UI.MenuGenerator.OpenConfirmation(SFS.Input.CloseMode.Current, () => TextLabel(), () => TextConfirm(), DisableClosestApproach, () => TextCancel(), DisableANAIS);
            }

            // Initialize ANAIS settings
            ANAIS_Config.Init(new FolderPath(ModFolder).ExtendToFile("Config.txt"));

            // Initialize logs
            // NOTE: To enable the logs, set "debug" to true, and compile in Debug configuration. Logs are disabled in release.
            // (logs activation is controlled by the ACTIVE_LOGS conditional compilation symbol - right-click on ANAIS, Properties, then "Build" panel)
            // NOTE2: To customize the logs, see the AnaisLogger.cs file
            AnaisLogger.Init(debug: true, new FolderPath(ModFolder).ExtendToFile("Logs_ANAIS.txt"));

            // Register scene changer events
            ModLoader.Helpers.SceneHelper.OnWorldSceneLoaded += new Action(AnaisManager.setWorldSceneActive);
            ModLoader.Helpers.SceneHelper.OnWorldSceneLoaded += new Action(ANAIS_Panel.ShowGUI);
            ModLoader.Helpers.SceneHelper.OnWorldSceneUnloaded += new Action(AnaisManager.setWorldSceneInactive);


            // Start ANAIS thread
            AnaisManager.StartTask();


            string TextLabel()
            {
                return "ANAIS and Closest approach line are not compatible. Which one do you want to keep?";
            }

            string TextConfirm()
            {
                return "ANAIS";
            }

            string TextCancel()
            {
                return "Closest approach";
            }

            void DisableClosestApproach()
            {
                ModsSettings.main.settings.modsActive[C_STR_CLOSEST_APPROACH_LINE_MOD_ID] = false;
                ModsSettings.main.SaveAll();
                ApplicationUtility.Relaunch();
            }

            void DisableANAIS()
            {
                ModsSettings.main.settings.modsActive[C_STR_MOD_ID] = false;
                ModsSettings.main.SaveAll();
                ApplicationUtility.Relaunch();
            }
        }

        public override void Early_Load()
        {
            // This method runs before anything from the game is loaded. This is where you should apply your patches, as shown below.
            ApplyPatch();
        }

        void ApplyPatch()
        {
            //UnityEngine.Debug.Log("Patching ANAIS");

            // The patcher uses an ID formatted like a web domain
            Main.patcher = new Harmony($"{C_STR_MOD_ID}.{C_STR_MOD_NAME}.{C_STR_AUTHOR}");

            // This pulls your Harmony patches from everywhere in the namespace and applies them.
            Main.patcher.PatchAll();
        }
    }
}

