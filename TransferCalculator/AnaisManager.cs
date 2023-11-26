using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using SFS.Variables;
using SFS.World;
using SFS.World.Maps;
using SFS.World.PlanetModules;
using SFS.WorldBase;
using UnityEngine;

class AnaisManager
{
    
    // This class allows to keep track of the game state and know if the thread dedicated to ANAIS should run or not.
    // This class is completely thread-safe
    private class AnaisManagerState
    {
        // -------------------
        // ---  VARIABLES  ---
        // -------------------

        // Serves as a lock to protect the class from concurrent access
        private readonly object AnaisStateLock = new object();

        // those variables allow to know the game state.
        // ANAIS must calculate transfers if we are in world scene, if the game is running (not paused), and a target is set.
        private bool isWorldSceneActive = false;
        private bool isGameRunning = false;
        private bool isNavActive = false;

        // When the target is reseted (changed/removed), this flags memorizes the reset request so that the ANAIS task 
        // can handle the reset on the working data set (anaisWorkingDataSet) in time (only the ANAIS task is allowed to touch that dataset!)
        private bool needsReset = false;

        // Allows to know if the ANAIS task is busy calculating. In case the player exits the game, this allows to delay the exit
        // until all calculations are over (because otherwise Unity will remove some crucial resources while we are still using them!)
        private bool isAnaisCalculating = false;


        // -------------------
        // ---   METHODS   ---
        // -------------------

        // Those functions are used to change the game state as it's known by AnaisManager.
        // They are meant to be called from the main thread. An internal lock allows to avoid any concurrent access
        public void setWorldSceneActive(bool active)
        {
            lock (AnaisStateLock)
            {
                LOG(LOG_LEVEL.INFO, "World scene: " + (active?"Active": "Inactive"));
                isWorldSceneActive = active;

                if(!active)
                {
                    isNavActive = false; // Because we can't really pretend the navigation is still active once we exited the game
                    isGameRunning = false; // same reason about game being running
                }
            }
        }

        public void setIsGameRunning(bool running)
        {
            lock (AnaisStateLock)
            {
                if(isGameRunning != running)
                {
                    LOG(LOG_LEVEL.INFO, "Game state: " + (running ? "Running": "Paused"));
                }
                isGameRunning = running;
            }
        }

        public void setNavigationActive(bool active)
        {
            lock (AnaisStateLock)
            {
                if(isNavActive != active)
                {
                    LOG(LOG_LEVEL.INFO, "Navigation: " + (active ? "Active" : "Inactive"));
                }
                isNavActive = active;
            }
        }

        // Call this to notify the ANAIS task that a reset is needed on all calculated data
        public void setIsResetNeeded(bool needed)
        {
            lock (AnaisStateLock)
            {
                needsReset = needed;
            }
        }

        // This function allows to tell if ANAIS should run given the current state (game paused/not paused, navigation active/inactive...)
        // If yes, it memorizes that ANAIS is now running and should not be interrupted (if game is exited while ANAIS is calculating...)
        public bool ActivateAnais()
        {
            lock (AnaisStateLock)
            {
                bool active = isWorldSceneActive && isGameRunning && isNavActive;
                //LOG(LOG_LEVEL.DEBUG, "Run ANAIS: " + (active ? "Yes" : "No") + " - World scene active: " + (isWorldSceneActive ? "Yes" : "No") + "; Game running: " + (isGameRunning ? "Yes" : "No") + "; Navigation active: " + (isNavActive ? "Yes" : "No"));
                if (active)
                {
                    isAnaisCalculating = true;
                }
                return active;
            }
        }

        // This function reports that ANAIS now goes back to a sleeping state, and that it's now safe to exit the game if needed
        public void EndAnaisActivity()
        {
            lock (AnaisStateLock)
            {
                isAnaisCalculating = false;
            }
        }

        // Tells if ANAIS is currently running
        public bool isAnaisActive()
        {
            lock (AnaisStateLock)
            {
                return isAnaisCalculating;
            }
        }

        // Accessor on the reset needed flag
        public bool NeedsReset()
        {
            lock (AnaisStateLock)
            {
                return needsReset;
            }
        }
    }


    // The ANAIS calculation periodicity in milliseconds
    private const int C_ANAIS_TASK_PERIOD = 33;
    private const int C_ANAIS_TASK_MIN_COOLDOWN = 10;

    private static Thread AnaisThread = null;

    // Allows to track the game state (in world scene, game paused, navigation active) to know if ANAIS should run or not
    private static AnaisManagerState anaisManagerState = new AnaisManagerState();

    // Internal data for ANAIS, so that it only works with a copy of the game data
    private static AnaisDataSet anaisInputDataSet   = new AnaisDataSet(); // This set is dedicated to receive the input data for the algorithm
    private static AnaisDataSet anaisWorkingDataSet = new AnaisDataSet(); // This set is used for the calculations - ONLY the ANAIS task is allowed to use that set!
    private static AnaisDataSet anaisOutputDataSet  = new AnaisDataSet(); // This set contains the result of the calculations to display the transfers

    // Lock objects to protect the input and output data sets that will be potentially accessed by 2 threads at the same time!
    private static readonly object AnaisInputDataLock  = new object();
    private static readonly object AnaisOutputDataLock = new object();


    // INTERNAL DATA (this is meant to be used by the ANAIS task only)
    // This is used to keep track of some sort of global variables implied in the calculations
    // -------------
    private static bool finalApproachMode = false; // indicates if the ship is close to the target, so the final approach informations should be displayed

    private static bool isOnEncounterTrajectory = false; // Allows to know if we are close to an encounter trajectory or not
    private static double encounterDate = -double.NegativeInfinity; // The encounter date if we are in encounter mode

    // Allows to remember the preferred number of turns for encounters on several turns
    private static double preferredTimeOfArrivalAtNode1;
    private static double preferredTimeOfArrivalAtNode2;

    // This is to know if the ANAIS transfer calculations should be inhibited (for situations in which the approach lines are more appropriate)
    private static bool allowANAIStransferCalculation = false;

    // GLOBAL CONSTANTS FOR THE CALCULATIONS
    // -------------------------------------
    public const double C_DELTAV_THRESHOLD_ENTER_ENCOUNTER_MODE = 2.5;
    public const double C_DELTAV_THRESHOLD_EXIT_ENCOUNTER_MODE = 3.0;
    public const double C_DISTANCE_THRESHOLD_FOR_FINAL_ENCOUNTER_MODE = 5000.0;
    public const double C_TIME_THRESHOLD_FOR_FINAL_ENCOUNTER_MODE = 60.0;

    
    // STATE MODIFIERS
    // ---------------
    // Those functions are used to change the game state as it's known by AnaisManager.
    // They are meant to be called from the main thread. An internal lock allows to avoid any concurrent access
    public static void setWorldSceneActive()
    {
        anaisManagerState.setWorldSceneActive(true);
    }

    public static void setWorldSceneInactive()
    {
        notifyNewTarget(null); // to force to reset the target and all data
        anaisManagerState.setWorldSceneActive(false);

        const int MAX_WAIT_TIME = 1000; // 1 second; In case a problem happened, the game won't be blocked because of this
        int time = 0;
        
        // Wait until ANAIS has stopped its calculations (to force the game to wait before unloading the world scene)
        while((time < MAX_WAIT_TIME) && anaisManagerState.isAnaisActive())
        {
            LOG(LOG_LEVEL.INFO, "ANAIS is running; waiting...");
            Thread.Sleep(C_ANAIS_TASK_MIN_COOLDOWN);
            time += C_ANAIS_TASK_MIN_COOLDOWN;
        }

        LOG(LOG_LEVEL.INFO, "ANAIS activity has ended; ready to exit game");
    }

    public static void setIsGameRunning()
    {
        anaisManagerState.setIsGameRunning(UnityEngine.Time.timeScale > 0.0f);
    }

    public static void notifyNewTarget(SelectableObject target)
    {
        LOG(LOG_LEVEL.INFO, "Asking RESET output data");
        anaisManagerState.setNavigationActive(target != null);

        // Reset the input and the output data sets
        lock (AnaisInputDataLock)
        {
            anaisInputDataSet.Reset();
        }

        lock (AnaisOutputDataLock)
        {
            anaisOutputDataSet.Reset();
        }

        // For the working set, the ANAIS task is notified through a flag, it will reset it itself.
        anaisManagerState.setIsResetNeeded(true);
    }

    // ANAIS INPUT/OUTPUT DATA MANAGEMENT
    // ----------------------------------
    public static void SetInputData(MapPlayer mapPlayer, SelectableObject target)
    {
        lock (AnaisInputDataLock)
        {
            anaisInputDataSet.SetInputData(mapPlayer, target);
        }
    }

    public static void DrawAnaisTransfer()
    {
        lock (AnaisOutputDataLock)
        {
            anaisOutputDataSet.DrawTransfer();
        }
    }

    public static bool getVelocityArrowData(out bool isInFinalApproachMode, out ClosestApproachCalculator.T_ApproachData approachData, out Double2 startingVelocity, out bool entrySOIdetected)
    {
        lock (AnaisOutputDataLock)
        {
            return anaisOutputDataSet.getVelocityArrowData(out isInFinalApproachMode, out approachData, out startingVelocity, out entrySOIdetected);
        }
    }

    // ANAIS TASK MANAGEMENT
    // ---------------------

    // -----------------------------------------------
    // -----         isAnaisThreadRunning         ----
    // -----------------------------------------------
    // A simple accessor that allows any function to know if it's
    // called by ANAIS or by the main thread (to prevent multithread issues...)
    // -----------------------------------------------
    public static bool isAnaisThreadRunning()
    {
        return (Thread.CurrentThread == AnaisThread);
    }

    // -----------------------------------------------
    // -----              StartTask               ----
    // -----------------------------------------------
    // This function allows to start the ANAIS calculation dedicated thread.
    // It's meant to be called once at the beginning of the program.
    // -----------------------------------------------
    public static void StartTask()
    {
        LOG(LOG_LEVEL.INFO, "Activate ANAIS thread");

        try
        {
            AnaisThread = new Thread(new ThreadStart(AnaisTask));
            AnaisThread.Name = "SFS.AnaisManager";
            AnaisThread.IsBackground = true;

            AnaisThread.Start();
        }
        catch (Exception ex)
        {
            LOG(LOG_LEVEL.ERROR, "Failed to start thread: " + ex.Message);
        }
    }


    // -----------------------------------------------
    // -----              AnaisTask               ----
    // -----------------------------------------------
    // The main function run by the dedicated ANAIS thread. Once it started, it loops endlessly
    // and calculates in real time all transfers. It evaluates by itself the game state to know
    // if the calculations should be done or not, and works on a fully independant set of data
    // to avoid any clash due to multithreading.
    // -----------------------------------------------
    private static void AnaisTask()
    {
        LOG(LOG_LEVEL.INFO, "Starting ANAIS process");

        while (Thread.CurrentThread.IsAlive)
        {
            Stopwatch stopwatch = null;

            try
            {
                stopwatch = Stopwatch.StartNew();

                // Reset output data if needed (if target has been changed or game has been exited)
                // ---------------------------
                if (anaisManagerState.NeedsReset())
                {
                    ApplyReset();
                }
                
                // Perform ANAIS calculations if needed
                // ------------------------------------
                if (anaisManagerState.ActivateAnais())
                {
                    // Swap data sets: take the input data set as our working data 
                    lock (AnaisInputDataLock)
                    {
                        AnaisDataSet.Swap(ref anaisInputDataSet, ref anaisWorkingDataSet);
                    }

                    // All good, GO ANAIS!
                    //Stopwatch stopwatch_calculation = Stopwatch.StartNew();
                    bool transferCalculated = anaisWorkingDataSet.CalculateTransfer(ref finalApproachMode, ref isOnEncounterTrajectory, ref encounterDate, ref preferredTimeOfArrivalAtNode1, ref preferredTimeOfArrivalAtNode2, ref allowANAIStransferCalculation);
                    //stopwatch_calculation.Stop();
                    

                    if (anaisManagerState.NeedsReset())
                    {
                        // If a reset order arrived while we were calculating, delete everything (OUCH!)
                        ApplyReset();
                    }
                    else if (transferCalculated)
                    {
                        LOG(LOG_LEVEL.DEBUG, "ANAIS calculation performed successfully!!!!!!!");
                        //LOG(LOG_LEVEL.INFO, "  ANAIS calculation processed in " + stopwatch_calculation.ElapsedTicks / 10.0 + " microseconds");

                        // Swap the working data set with the output to make it available
                        lock (AnaisOutputDataLock)
                        {
                            AnaisDataSet.Swap(ref anaisOutputDataSet, ref anaisWorkingDataSet); // The data is ready for display, swap with output data set
                        }
                        
                        // Initialize the new working data set
                        anaisWorkingDataSet.Reset();
                    }
                }
            }
            catch (System.Threading.ThreadAbortException ex)
            {
                LOG(LOG_LEVEL.INFO, "Stopping ANAIS thread");
            }
            catch (Exception ex)
            {
                LOG(LOG_LEVEL.ERROR, "EXCEPTION - ANAIS thread encountered an error: " + ex.Message);
                UnityEngine.Debug.Log("EXCEPTION - ANAIS thread encountered an error: " + ex.Message);
            }
            finally
            {
                // Report that Anais now goes back to sleeping state
                anaisManagerState.EndAnaisActivity();

                // Leave the thread sleep for the remaining time
                // ---------------------------------------------
                stopwatch.Stop();

                int remainingTime = C_ANAIS_TASK_PERIOD - (int)(stopwatch.ElapsedMilliseconds);
                if (remainingTime < C_ANAIS_TASK_MIN_COOLDOWN) remainingTime = C_ANAIS_TASK_MIN_COOLDOWN;

                Thread.Sleep(remainingTime);
            }
        }

        // Resets all calculated data and all parameters from the calculation algorithm - follows a toggle target order
        void ApplyReset()
        {
            LOG(LOG_LEVEL.INFO, "Reset ANAIS work data");
            anaisManagerState.setIsResetNeeded(false);
            anaisWorkingDataSet.Reset();

            // Reset Anais Manager internal data
            finalApproachMode = false;
            isOnEncounterTrajectory = false;
            encounterDate = double.NegativeInfinity;
            preferredTimeOfArrivalAtNode1 = double.NegativeInfinity;
            preferredTimeOfArrivalAtNode2 = double.NegativeInfinity;
            allowANAIStransferCalculation = false;
        }
    }
    


    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.ANAIS_MANAGER, level, message);
    }
}
