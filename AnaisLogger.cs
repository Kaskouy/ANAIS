using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

using HarmonyLib;


public enum LOG_LEVEL
{
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    NONE
}

public enum LOG_CATEGORY
{
    LAMBERT_SOLVER,
    HOHMANN_TRANSFER,
    ORBIT,
    VELOCITY_ARROW,
    MAP_NAVIGATION,
    CLOSEST_APPROACH,
    EJECTION_TRAJECTORY,
    RETURN_TRAJECTORY,
    NUMERICAL_MINI_CALCULATOR,
    ANAIS_TRANSFER,
    ANAIS_TRANSFER_CALCULATOR,
    ANAIS_MANAGER,
    ANAIS_DATA_SET,
    VECTOR2D_FIXED_POINT_SOLVER
}

class AnaisLogger
{
    private static readonly string[] STR_LOG_LEVEL = { "[DEBUG]", "[INFO ]", "[WARN ]", "[ERROR]", "[NONE ]" };

    private static readonly string[] STR_LOG_CATEGORY = { "LAMBERT_SOLVER", 
                                                          "HOHMANN_TRANSFER", 
                                                          "ORBIT", 
                                                          "VELOCITY_ARROW",
                                                          "MAP_NAVIGATION",
                                                          "CLOSEST_APPROACH",
                                                          "EJECTION_TRAJECTORY",
                                                          "RETURN_TRAJECTORY",
                                                          "NUMERICAL_MINI_CALCULATOR",
                                                          "ANAIS_TRANSFER",
                                                          "ANAIS_TRANSFER_CALCULATOR",
                                                          "ANAIS_MANAGER",
                                                          "ANAIS_DATA_SET",
                                                          "VECTOR2D_FIXED_POINT_SOLVER" };

    private static Dictionary<LOG_CATEGORY, LOG_LEVEL> ListLogLevels = new Dictionary<LOG_CATEGORY, LOG_LEVEL>();

    // Serves as a lock so that only one thread can log data at a given time
    private static readonly object AnaisLoggerLock = new object();

    public static void Init(bool debug, string fileLogPath)
    {
#if ACTIVE_LOGS

        // To specify the log file path...
        Environment.SetEnvironmentVariable("HARMONY_LOG_FILE", fileLogPath);

        HarmonyLib.Tools.HarmonyFileLog.FileWriterPath = fileLogPath;
        HarmonyLib.Tools.HarmonyFileLog.Enabled = false; // instead of debug passed as parameter; setting it to true messes everything...


        // To make sure an entry exists for each category - No log by default
        foreach (LOG_CATEGORY logCategory in Enum.GetValues(typeof(LOG_CATEGORY)))
        {
            ListLogLevels.Add(logCategory, LOG_LEVEL.NONE);
        }

        // Define the desired level of log for each category
        ListLogLevels[LOG_CATEGORY.LAMBERT_SOLVER]              = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.HOHMANN_TRANSFER]            = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.ORBIT]                       = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.VELOCITY_ARROW]              = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.MAP_NAVIGATION]              = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.CLOSEST_APPROACH]            = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.EJECTION_TRAJECTORY]         = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.RETURN_TRAJECTORY]           = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.NUMERICAL_MINI_CALCULATOR]   = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.ANAIS_TRANSFER]              = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.ANAIS_TRANSFER_CALCULATOR]   = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.ANAIS_MANAGER]               = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.ANAIS_DATA_SET]              = LOG_LEVEL.WARNING;
        ListLogLevels[LOG_CATEGORY.VECTOR2D_FIXED_POINT_SOLVER] = LOG_LEVEL.WARNING;

#else
        Environment.SetEnvironmentVariable("HARMONY_NO_LOG", "1");

        foreach(LOG_CATEGORY logCategory in Enum.GetValues(typeof(LOG_CATEGORY)))
        {
            ListLogLevels.Add(logCategory, LOG_LEVEL.NONE);
        }
#endif
    }

    [Conditional("ACTIVE_LOGS")]
    public static void Log(LOG_CATEGORY category, LOG_LEVEL level, string message)
    {
        lock(AnaisLoggerLock)
        {
            if ((level != LOG_LEVEL.NONE) && (level >= ListLogLevels[category]))
            {
                FileLog.Log(STR_LOG_LEVEL[(int)level] + " " + STR_LOG_CATEGORY[(int)category] + " : " + message);
            }
        }
    }
}
