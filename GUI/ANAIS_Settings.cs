using JetBrains.Annotations;
using SFS.IO;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using UITools;

[Serializable]
public class ANAIS_Settings
{
    public enum E_TRANSFER_TYPE
    {
        RENDEZ_VOUS,
        FLY_BY
    };

    public enum E_PREFERRED_NODE
    {
        ANY,
        ASCENDING,
        DESCENDING
    }

    public const int C_NB_TURNS_MIN = 0;
    public const int C_NB_TURNS_MAX = 24;

    public bool _minimized = false;

    public bool _showTransfer = true;

    public E_TRANSFER_TYPE _transferType = E_TRANSFER_TYPE.RENDEZ_VOUS;

    public uint _nbMaxTurns = 12;

    public E_PREFERRED_NODE _preferredNode = E_PREFERRED_NODE.ANY;

    public bool _showTips = true;

    public ANAIS_Settings()
    {
        SetDefaultValues();
    }

    public void copySettings(ANAIS_Settings other)
    {
        if(other != null)
        {
            _transferType = other._transferType;
            _showTransfer = other._showTransfer;
            _nbMaxTurns = other._nbMaxTurns;
            _preferredNode = other._preferredNode;
        }
    }

    public void SetDefaultValues()
    {
        _transferType = E_TRANSFER_TYPE.RENDEZ_VOUS;
        _showTransfer = true;
        _nbMaxTurns = 12;
        _preferredNode = E_PREFERRED_NODE.ANY;
    }
}

public class NavigationVariables
{
    public bool _approachLinesPreferredDatesNeedReset = false;

    // target altitude related variables
    public string _currentPlanetName = "";
    public double _defaultTargetAltitude = 0.0;
    public double _maxTargetAltitude = 0.0;
    public double _targetAltitude = 0.0;

    public NavigationVariables()
    {
        Reset();
    }

    public void Reset()
    {
        _approachLinesPreferredDatesNeedReset = false;
        // _targetAltitude is not to be reset
    }

    public void Copy(NavigationVariables other)
    {
        _approachLinesPreferredDatesNeedReset = other._approachLinesPreferredDatesNeedReset;
        _targetAltitude = other._targetAltitude;
    }
}

public class ANAIS_Config: ModSettings<ANAIS_Settings>
{
    private static ANAIS_Config _instance = null;

    private static FilePath _configPath = null;

    private static Action saveAction = null;

    public static void Init(FilePath filePath)
    {
        _instance = new ANAIS_Config();
        _configPath = filePath;

        _instance.Initialize();

        // set the reference from the panel to the settings instance
        ANAIS_Panel._settings = settings;
    }

    public static void Save()
    {
        saveAction?.Invoke();
    }

    protected override FilePath SettingsFile { get { return _configPath; } }

    protected override void RegisterOnVariableChange(Action onChange)
    {
        saveAction = onChange;
    }
}
