using System;
using System.Collections.Generic;

using SFS.IO;
using ModLoader;
using UITools;
using HarmonyLib;


namespace ANAIS
{
    public class Main : Mod , IUpdatable
    {
        const string C_STR_MOD_ID = "ANAIS";
        const string C_STR_MOD_NAME = "ANAIS";
        const string C_STR_AUTHOR = "Altaïr";
        const string C_STR_GAME_VERSION = "1.5.10.2";
        const string C_STR_MOD_VERSION = "v1.3.2";
        const string C_STR_MOD_DESCRIPTION = "Advanced NAvigation Innovative System\nReplaces the original navigation system with a more elaborated one.";

        private const string C_STR_CLOSEST_APPROACH_LINE_MOD_ID = "CLOSEST_APPROACH_LINE";

        public override string ModNameID => C_STR_MOD_ID;

        public override string DisplayName => C_STR_MOD_NAME;

        public override string Author => C_STR_AUTHOR;

        public override string MinimumGameVersionNecessary => C_STR_GAME_VERSION;

        public override string ModVersion => C_STR_MOD_VERSION;

        public override string Description => C_STR_MOD_DESCRIPTION;

        public override string IconLink => "https://i.imgur.com/JDBeEJD.png"; // link to the logo

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
            // NOTE: To enable the logs, under Visual right-click on the project ("ANAIS"), select Properties, then the "build" category,
            // then add ACTIVE_LOGS as a conditional compilation symbol. Remove it to disable all logs.
            // NOTE2: To customize the logs, see the AnaisLogger.cs file
            //AnaisLogger.Init(debug: false, "C:\\Users\\JB\\Desktop\\Jeux\\SFS PC\\ANAIS\\Logs_ANAIS.txt");
            System.DateTime dateNow = DateTime.Now;
            AnaisLogger.Init(debug: false, new FolderPath(ModFolder).ExtendToFile("Logs_ANAIS_" + dateNow.Day + dateNow.Month + dateNow.Year + "_" + dateNow.Hour + dateNow.Minute + dateNow.Second + ".txt"));

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

