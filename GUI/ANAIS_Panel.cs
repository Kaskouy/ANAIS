using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Permissions;
using System.Text;

using HarmonyLib;
using SFS.UI;
using SFS.UI.ModGUI;
using TMPro;
using UITools;
using UnityEngine;
using static SFS.Parsers.Ini.IniDataEnv;
using Type = SFS.UI.ModGUI.Type;
using static ANAIS_Panel_Constants;


class ANAIS_Panel
{
    // GameObject to attach the window.
    private static GameObject _windowHolder = null;

    // The set of parameters associated to the panel
    public static ANAIS_Settings _settings = null;

    // Set of variables that will be used by ANAIS' calculations
    public static NavigationVariables _NavVariables = new NavigationVariables();

    // ANAIS Control Panel identifier
    private const string C_ANAIS_PANEL_IDENTIFIER = "ANAIS.Panel";

    // Widget list
    private static ClosableWindow _mainWindow = null;

    private static Toggle _showTransfer_toggle = null;
    private static SFS.UI.ModGUI.Button _transferModeButton = null;
    private static Label _label_text_transfer_mode = null;
    private static Label _label_transfer_mode = null;
    private static TextInput _textInput_target_altitude = null;
    private static TextInput _textInput_nb_max_turns = null;
    private static SFS.UI.ModGUI.Button _decreaseNbMaxTurnsButton = null;
    private static SFS.UI.ModGUI.Button _increaseNbMaxTurnsButton = null;
    private static SFS.UI.ModGUI.Button _changeNodeButton = null;
    private static Label _label_node = null;
    private static SFS.UI.ModGUI.Button _defaultSettingsButton = null;

    private enum E_ANAIS_PANEL_BUTTON
    {
        TRANSFER_MODE_BUTTON,
        DECREASE_NB_MAX_TURNS_BUTTON,
        INCREASE_NB_MAX_TURNS_BUTTON,
        CHANGE_NODE_BUTTON,
        DEFAULT_SETTINGS_BUTTON,
        NONE
    }

    // Tooltip panel
    private static GameObject _tooltipHolder;
    private static Window _tooltipWindow;
    private static Label _tooltipLabel;

    // Creates the panel
    public static void ShowGUI()
    {
        // CREATE THE MAIN WINDOW
        // ----------------------

        // Create the window holder, attach it to the currently active scene so it's removed when the scene changes.
        _windowHolder = Builder.CreateHolder(Builder.SceneToAttach.CurrentScene, "ANAIS panel holder");

        // Instantiate the main window - Important: savePosition must be set to false
        _mainWindow = UITools.UIToolsBuilder.CreateClosableWindow(_windowHolder.transform, Builder.GetRandomID(), C_MAIN_WINDOW_WIDTH, C_MAIN_WINDOW_HEIGHT, posX: 0, posY: 0, draggable: true, savePosition: false, opacity: 0.95f, C_MAIN_WINDOW_TITLE, false /*_settings._minimized*/);

        // Bypassing bug from UITools: creating the window in minimized state doesn't set the correct size for it, so it's created in non-minimized state, then manually minimized if needed.
        _mainWindow.Minimized = _settings._minimized;

        // Bind event to be able to keep track of the minimized state
        _mainWindow.OnMinimizedChangedEvent += onMinimizeClicked;

        // Allows to set the window coordinates to the previously saved ones - Also allows to save the window coordinates in real time when the player moves it
        _mainWindow.RegisterPermanentSaving(C_ANAIS_PANEL_IDENTIFIER);



        // -----   BUILD THE CONTENT OF THE ANAIS CONTROL PANEL   -----
        // ------------------------------------------------------------
        // posX = 0: median position; negative value: moves to the left; positive value: moves to the right
        // posY = 0: high position; negative value: moves downwards; positive value: moves upwards

        // First separator
        Builder.CreateSeparator(_mainWindow, C_SEPARATOR_WIDTH, 0, C_FIRST_SEPARATOR_Y_POS);

        // TRANSFER MODE
        // -------------

        // Label for transfer mode (text)
        _label_text_transfer_mode = Builder.CreateLabel(_mainWindow, C_LABEL_TEXT_TRANSFER_MODE_WIDTH, C_LABEL_TEXT_TRANSFER_MODE_HEIGHT, C_LABEL_TEXT_TRANSFER_MODE_POS_X, C_LABEL_TEXT_TRANSFER_MODE_POS_Y, C_LABEL_TEXT_TRANSFER_MODE_TEXT);
        _label_text_transfer_mode.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // Label for transfer mode (mode itself)
        _label_transfer_mode = Builder.CreateLabel(_mainWindow, C_LABEL_TRANSFER_MODE_WIDTH, C_LABEL_TRANSFER_MODE_HEIGHT, C_LABEL_TRANSFER_MODE_POS_X, C_LABEL_TRANSFER_MODE_POS_Y, getTransferTypeText());
        _label_transfer_mode.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // Toggle transfer mode button
        _transferModeButton = Builder.CreateButton(_mainWindow, C_BUTTON_TOGGLE_TRANSFER_MODE_WIDTH, C_BUTTON_TOGGLE_TRANSFER_MODE_HEIGHT, C_BUTTON_TOGGLE_TRANSFER_MODE_POS_X, C_BUTTON_TOGGLE_TRANSFER_MODE_POS_Y, toggleTransferMode, C_BUTTON_TOGGLE_TRANSFER_MODE_TEXT);
        _transferModeButton.SetAlignment(TMPro.TextAlignmentOptions.Midline);

        // Label "Show transfer"
        Label label_show_transfer = Builder.CreateLabel(_mainWindow, C_LABEL_SHOW_TRANSFER_WIDTH, C_LABEL_SHOW_TRANSFER_HEIGHT, C_LABEL_SHOW_TRANSFER_POS_X, C_LABEL_SHOW_TRANSFER_POS_Y, C_LABEL_SHOW_TRANSFER_TEXT);
        label_show_transfer.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // "Show transfer" toggle
        _showTransfer_toggle = Builder.CreateToggle(_mainWindow, getShowTransferValue, C_TOGGLE_SHOW_TRANSFER_POS_X, C_TOGGLE_SHOW_TRANSFER_POS_Y, toggleShowTransfer);
        _showTransfer_toggle.Size = new Vector2(C_TOGGLE_SHOW_TRANSFER_SIZE_X, C_TOGGLE_SHOW_TRANSFER_SIZE_Y);

        // Label "Target altitude"
        Label label_targetAltitude = Builder.CreateLabel(_mainWindow, C_LABEL_TEXT_TARGET_ALTITUDE_WIDTH, C_LABEL_TEXT_TARGET_ALTITUDE_HEIGHT, C_LABEL_TEXT_TARGET_ALTITUDE_POS_X, C_LABEL_TEXT_TARGET_ALTITUDE_POS_Y, C_LABEL_TEXT_TARGET_ALTITUDE_TEXT);
        label_targetAltitude.TextAlignment = TMPro.TextAlignmentOptions.Left;

        _textInput_target_altitude = Builder.CreateTextInput(_mainWindow, C_TEXT_TARGET_ALTITUDE_WIDTH, C_TEXT_TARGET_ALTITUDE_HEIGHT, C_TEXT_TARGET_ALTITUDE_POS_X, C_TEXT_TARGET_ALTITUDE_POS_Y, AltitudeFromGameUnitToPanelString(_NavVariables._targetAltitude), ANAIS_Panel.onTargetAltitudeChanged);
        _textInput_target_altitude.field.characterValidation = TMPro.TMP_InputField.CharacterValidation.Decimal; // Only allows digits to be entered
        _textInput_target_altitude.field.characterLimit = 11; // Max 10 digits + 1 comma
        _textInput_target_altitude.field.SetGlobalPointSize(24.0f);

        // Label "Target altitude unit"
        Label label_targetAltitudeUnit = Builder.CreateLabel(_mainWindow, C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_WIDTH, C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_HEIGHT, C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_POS_X, C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_POS_Y, C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_TEXT);
        label_targetAltitudeUnit.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // Second separator
        Builder.CreateSeparator(_mainWindow, C_SEPARATOR_WIDTH, posX: 0, C_SECOND_SEPARATOR_Y_POS);


        // APPROACH LINES SECTION
        // ----------------------

        // Label for approach lines ("title")
        Label label_approach_lines = Builder.CreateLabel(_mainWindow, C_LABEL_APPROACH_LINES_WIDTH, C_LABEL_APPROACH_LINES_HEIGHT, C_LABEL_APPROACH_LINES_POS_X, C_LABEL_APPROACH_LINES_POS_Y, C_LABEL_APPROACH_LINES_TEXT);
        label_approach_lines.TextAlignment = TMPro.TextAlignmentOptions.Center;

        // Label "Nb max turns"
        Label label_nb_max_turns = Builder.CreateLabel(_mainWindow, C_LABEL_NB_MAX_TURNS_WIDTH, C_LABEL_NB_MAX_TURNS_HEIGHT, C_LABEL_NB_MAX_TURNS_POS_X, C_LABEL_NB_MAX_TURNS_POS_Y, C_LABEL_NB_MAX_TURNS_TEXT);
        label_nb_max_turns.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // Nb turns text input
        _textInput_nb_max_turns = Builder.CreateTextInput(_mainWindow, C_TEXT_INPUT_NB_MAX_TURNS_WIDTH, C_TEXT_INPUT_NB_MAX_TURNS_HEIGHT, C_TEXT_INPUT_NB_MAX_TURNS_POS_X, C_TEXT_INPUT_NB_MAX_TURNS_POS_Y, _settings._nbMaxTurns.ToString(), ANAIS_Panel.onNbMaxTurnsChanged);
        _textInput_nb_max_turns.field.characterValidation = TMPro.TMP_InputField.CharacterValidation.Digit; // Only allows digits to be entered
        _textInput_nb_max_turns.field.characterLimit = 2; // Max 2 characters
        _textInput_nb_max_turns.field.SetGlobalPointSize(24.0f);

        // "Minus 1" button
        _decreaseNbMaxTurnsButton = Builder.CreateButton(_mainWindow, C_BUTTON_MINUS_ONE_WIDTH, C_BUTTON_MINUS_ONE_HEIGHT, C_BUTTON_MINUS_ONE_POS_X, C_BUTTON_MINUS_ONE_POS_Y, decreaseNbMaxTurns, C_BUTTON_MINUS_ONE_TEXT);
        _decreaseNbMaxTurnsButton.SetAlignment(TMPro.TextAlignmentOptions.Center);

        // "Plus 1" button
        _increaseNbMaxTurnsButton = Builder.CreateButton(_mainWindow, C_BUTTON_PLUS_ONE_WIDTH, C_BUTTON_PLUS_ONE_HEIGHT, C_BUTTON_PLUS_ONE_POS_X, C_BUTTON_PLUS_ONE_POS_Y, increaseNbMaxTurns, C_BUTTON_PLUS_ONE_TEXT);
        _increaseNbMaxTurnsButton.SetAlignment(TMPro.TextAlignmentOptions.Center);

        // Label "Nb max turns"
        Label label_preferred_node = Builder.CreateLabel(_mainWindow, C_LABEL_PREFERRED_NODE_WIDTH, C_LABEL_PREFERRED_NODE_HEIGHT, C_LABEL_PREFERRED_NODE_POS_X, C_LABEL_PREFERRED_NODE_POS_Y, C_LABEL_PREFERRED_NODE_TEXT);
        label_preferred_node.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // Label "Node"
        _label_node = Builder.CreateLabel(_mainWindow, C_LABEL_NODE_WIDTH, C_LABEL_NODE_HEIGHT, C_LABEL_NODE_POS_X, C_LABEL_NODE_POS_Y, GetPreferredNode_Str(_settings._preferredNode));
        _label_node.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // "Change node" button
        _changeNodeButton = Builder.CreateButton(_mainWindow, C_BUTTON_CHANGE_NODE_WIDTH, C_BUTTON_CHANGE_NODE_HEIGHT, C_BUTTON_CHANGE_NODE_POS_X, C_BUTTON_CHANGE_NODE_POS_Y, onChangeNode, C_BUTTON_CHANGE_NODE_TEXT);
        _changeNodeButton.SetAlignment(TMPro.TextAlignmentOptions.Midline);

        // Third separator
        Builder.CreateSeparator(_mainWindow, C_SEPARATOR_WIDTH, posX: 0, C_THIRD_SEPARATOR_Y_POS);


        // SETTINGS SECTION
        // ----------------

        // Label "Show tips"
        Label label_tips = Builder.CreateLabel(_mainWindow, C_LABEL_TIPS_WIDTH, C_LABEL_TIPS_HEIGHT, C_LABEL_TIPS_POS_X, C_LABEL_TIPS_POS_Y, C_LABEL_TIPS_TEXT);
        label_tips.TextAlignment = TMPro.TextAlignmentOptions.Left;

        // Tips toggle
        Toggle showTips_toggle = Builder.CreateToggle(_mainWindow, getShowTipsValue, C_TOGGLE_TIPS_POS_X, C_TOGGLE_TIPS_POS_Y, toggleTips);
        showTips_toggle.Size = new Vector2(C_TOGGLE_TIPS_SIZE_X, C_TOGGLE_TIPS_SIZE_Y);

        // "Default settings" button
        _defaultSettingsButton = Builder.CreateButton(_mainWindow, C_BUTTON_DEFAULT_SETTINGS_WIDTH, C_BUTTON_DEFAULT_SETTINGS_HEIGHT, C_BUTTON_DEFAULT_SETTINGS_POS_X, C_BUTTON_DEFAULT_SETTINGS_POS_Y, applyDefaultSettings, C_BUTTON_DEFAULT_SETTINGS_TEXT);
        _defaultSettingsButton.SetAlignment(TMPro.TextAlignmentOptions.Midline);


        // TOOLTIP
        // -------
        _tooltipHolder = Builder.CreateHolder(Builder.SceneToAttach.CurrentScene, "Tooltip holder");

        _tooltipWindow = Builder.CreateWindow(_tooltipHolder.transform, Builder.GetRandomID(), C_LABEL_WINDOW_WIDTH, C_LABEL_WINDOW_HEIGHT, posX: 0, posY: 0, draggable: false, savePosition: false, opacity: 0.5f, "Help for big noobs");

        Box box = Builder.CreateBox(_tooltipWindow, C_LABEL_WINDOW_WIDTH - 6, C_LABEL_WINDOW_HEIGHT - 60, posX: 0, posY: -120);

        _tooltipLabel = Builder.CreateLabel(box, C_LABEL_TOOLTIP_WIDTH, C_LABEL_TOOLTIP_HEIGHT, 0, -40, "");
        _tooltipLabel.AutoFontResize = false;
        _tooltipLabel.FontSize = 24.0f;
        _tooltipLabel.TextAlignment = TMPro.TextAlignmentOptions.TopLeft;

        // Hide the tooltip for now
        _tooltipWindow.Active = false;
    }


    // EVENT FUNCTIONS
    // ---------------

    // Minimize button
    private static void onMinimizeClicked()
    {
        _settings._minimized = !_settings._minimized;
        ANAIS_Config.Save();
    }
    
    // Show transfer button
    private static bool getShowTransferValue()
    {
        return _settings._showTransfer;
    }

    private static void toggleShowTransfer()
    {
        _settings._showTransfer = !_settings._showTransfer;
        ANAIS_Config.Save();
    }


    // Toggle transfer mode button
    private static string getTransferTypeText()
    {
        switch(_settings._transferType)
        {
            case ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS: return C_LABEL_TRANSFER_MODE_RENDEZ_VOUS_TEXT;
            case ANAIS_Settings.E_TRANSFER_TYPE.FLY_BY: return C_LABEL_TRANSFER_MODE_FLY_BY_TEXT;
            default: return "None";
        }
    }

    private static void toggleTransferMode()
    {
        if(_settings._transferType == ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS)
        {
            _settings._transferType = ANAIS_Settings.E_TRANSFER_TYPE.FLY_BY;
        }
        else
        {
            _settings._transferType = ANAIS_Settings.E_TRANSFER_TYPE.RENDEZ_VOUS;
        }

        _label_transfer_mode.Text = getTransferTypeText();

        ANAIS_Config.Save();
    }

    private static double AltitudeFromPanelToGameUnit(double alt)
    {
        return alt * 1000.0;
    }

    private static double AltitudeFromGameToPanelUnit(double alt)
    {
        //double panelAlt = alt / 1000.0;
        //panelAlt.Round(3);
        return alt / 1000.0;
    }

    private static string AltitudeFromGameUnitToPanelString(double alt)
    {
        double panelAlt = Math.Truncate(alt) / 1000.0;
        //return panelAlt.ToString(3, false);
        return panelAlt.ToString();
    }

    // Target altitude
    private static void onTargetAltitudeChanged(string str)
    {
        bool ok = double.TryParse(str, out double altitude);

        if (ok)
        {
            altitude = AltitudeFromPanelToGameUnit(altitude);

            if (altitude < 0.0)
            {
                // entry is below acceptable minimum; force to lowest value
                altitude = 0.0;
                _textInput_target_altitude.Text = altitude.ToString(); // Warning: calls onTargetAltitudeChanged again
                _NavVariables._targetAltitude = altitude;
            }
            else if (altitude > _NavVariables._maxTargetAltitude)
            {
                // entry is beyond acceptable maximum; force to highest value
                altitude = _NavVariables._maxTargetAltitude;
                _textInput_target_altitude.Text = AltitudeFromGameUnitToPanelString(altitude); // Warning: calls onTargetAltitudeChanged again
                _NavVariables._targetAltitude = altitude;
            }
            else
            {
                // Correct entry: update internal variable
                _NavVariables._targetAltitude = altitude;
            }
        }
        else
        {
            _textInput_target_altitude.Text = "0.0"; // Warning: calls onTargetAltitudeChanged again
            _NavVariables._targetAltitude = 0.0;
        }
    }


    public static void updateTargetAltitude(SFS.World.SelectableObject target)
    {
        // Retrieve planet (null if not a planet)
        SFS.WorldBase.Planet planet = null;

        SFS.World.Maps.MapPlanet targetMapPlanet = target as SFS.World.Maps.MapPlanet;
        if (targetMapPlanet != null) planet = targetMapPlanet.planet;

        // Update navigation data
        if (planet == null)
        {
            // No planet selected - default and max are reset, but the planet name and the current altitude are intentionally kept
            // so that they are available again if the player selects the same planet again later.
            _NavVariables._defaultTargetAltitude = 0.0;
            _NavVariables._maxTargetAltitude = 0.0;

            if(_windowHolder != null) // avoids breaking the game when this function is called at start-up
            {
                double tmp_alt = _NavVariables._targetAltitude; // because setting _textInput_target_altitude.Text to "0.0" resets that variable...
                _textInput_target_altitude.Text = "0.0";
                _NavVariables._targetAltitude = tmp_alt;
            }
        }
        else
        {
            // A planet is selected
            _NavVariables._defaultTargetAltitude = planet.TimewarpRadius_Descend - planet.Radius;
            _NavVariables._maxTargetAltitude = planet.SOI - planet.Radius;

            if (planet.name != _NavVariables._currentPlanetName)
            {
                // The planet selected is different from before --> setting planet name and default altitude
                _NavVariables._targetAltitude = _NavVariables._defaultTargetAltitude;
                _NavVariables._currentPlanetName = planet.name;
            }

            if (_windowHolder != null) // avoids breaking the game when this function is called at start-up
            {
                _textInput_target_altitude.Text = AltitudeFromGameUnitToPanelString(_NavVariables._targetAltitude);
            }
        }
    }


    // Nb. max turns input
    private static void onNbMaxTurnsChanged(string str)
    {
        bool ok = uint.TryParse(str, out uint nbTurns);

        if (ok)
        {
            if (nbTurns < ANAIS_Settings.C_NB_TURNS_MIN)
            {
                // entry is below acceptable minimum; force to lowest value
                nbTurns = ANAIS_Settings.C_NB_TURNS_MIN;
                _textInput_nb_max_turns.Text = nbTurns.ToString(); // Warning: calls onNbMaxTurnsChanged again
            }
            else if (nbTurns > ANAIS_Settings.C_NB_TURNS_MAX)
            {
                // entry is beyond acceptable maximum; force to highest value
                nbTurns = ANAIS_Settings.C_NB_TURNS_MAX;
                _textInput_nb_max_turns.Text = nbTurns.ToString(); // Warning: calls onNbMaxTurnsChanged again
            }

            _settings._nbMaxTurns = nbTurns;
        }
        else
        {
            // Invalid entry; use default value
            _settings._nbMaxTurns = 12;
            _textInput_nb_max_turns.Text = _settings._nbMaxTurns.ToString();
        }

        // Tell the navigation to reset the approach lines preferred dates 
        _NavVariables._approachLinesPreferredDatesNeedReset = true;

        ANAIS_Config.Save();
    }

    // increase/decrease Nb. max turns buttons
    private static void decreaseNbMaxTurns()
    {
        if(_settings._nbMaxTurns > 0)
        {
            _textInput_nb_max_turns.Text = ((int)(_settings._nbMaxTurns - 1)).ToString();
        }
        else
        {
            _textInput_nb_max_turns.Text = "0";
        }
    }

    private static void increaseNbMaxTurns()
    {
        _textInput_nb_max_turns.Text = ((int)(_settings._nbMaxTurns + 1)).ToString();
    }


    // Change node button
    private static void onChangeNode()
    {
        switch (_settings._preferredNode)
        {
            case ANAIS_Settings.E_PREFERRED_NODE.ANY:
                {
                    _settings._preferredNode = ANAIS_Settings.E_PREFERRED_NODE.ASCENDING;
                    break;
                }

            case ANAIS_Settings.E_PREFERRED_NODE.ASCENDING:
                {
                    _settings._preferredNode = ANAIS_Settings.E_PREFERRED_NODE.DESCENDING;
                    break;
                }

            case ANAIS_Settings.E_PREFERRED_NODE.DESCENDING:
                {
                    _settings._preferredNode = ANAIS_Settings.E_PREFERRED_NODE.ANY;
                    break;
                }

            default:
                _settings._preferredNode = ANAIS_Settings.E_PREFERRED_NODE.ANY; break;
        }

        _label_node.Text = GetPreferredNode_Str(_settings._preferredNode);
        ANAIS_Config.Save();
    }

    private static string GetPreferredNode_Str(ANAIS_Settings.E_PREFERRED_NODE preferredNode)
    {
        if (preferredNode == ANAIS_Settings.E_PREFERRED_NODE.ASCENDING)
        {
            return C_LABEL_PREFERRED_NODE_ASCENDING;
        }
        else if (preferredNode == ANAIS_Settings.E_PREFERRED_NODE.DESCENDING)
        {
            return C_LABEL_PREFERRED_NODE_DESCENDING;
        }
        else
        {
            return C_LABEL_PREFERRED_NODE_ANY;
        }
    }

    // Show tips toggle
    private static bool getShowTipsValue()
    {
        return _settings._showTips;
    }

    private static void toggleTips()
    {
        _settings._showTips = !_settings._showTips;
        ANAIS_Config.Save();
    }

    // Default settings button
    private static void applyDefaultSettings()
    {
        _settings.SetDefaultValues();

        // Update "Show transfer" toggle
        _showTransfer_toggle.toggleButton.UpdateUI(true);

        // Update transfer mode
        _label_transfer_mode.Text = getTransferTypeText();

        // Update nb max turns
        _textInput_nb_max_turns.Text = _settings._nbMaxTurns.ToString();

        // update preferred node
        _label_node.Text = GetPreferredNode_Str(_settings._preferredNode);

        // Save config
        ANAIS_Config.Save();
    }


    // ANAIS_PanelEventListener: Allows to track onMouseEnter and onMouseExit events for buttons, to allow showing the tooltips
    // ------------------------
    public class ANAIS_PanelEventListener
    {
        private static E_ANAIS_PANEL_BUTTON _buttonCurrentlyOver = E_ANAIS_PANEL_BUTTON.NONE;

        public static void CheckEvent()
        {
            E_ANAIS_PANEL_BUTTON newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.NONE;

            // Determine which button the mouse pointer is over
            if (_transferModeButton.isMouseOver()) newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.TRANSFER_MODE_BUTTON;
            else if(_decreaseNbMaxTurnsButton.isMouseOver()) newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.DECREASE_NB_MAX_TURNS_BUTTON;
            else if (_increaseNbMaxTurnsButton.isMouseOver()) newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.INCREASE_NB_MAX_TURNS_BUTTON;
            else if (_changeNodeButton.isMouseOver()) newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.CHANGE_NODE_BUTTON;
            else if (_defaultSettingsButton.isMouseOver()) newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.DEFAULT_SETTINGS_BUTTON;
            else newButtonCurrentlyOver = E_ANAIS_PANEL_BUTTON.NONE;

            // If the situation has changed, do what is needed
            if(newButtonCurrentlyOver != _buttonCurrentlyOver)
            {
                _buttonCurrentlyOver = newButtonCurrentlyOver;
                updateTooltip(_buttonCurrentlyOver);
            }
        }
    }


    // Display tooltips functions
    // --------------------------
    private static void updateTooltip(E_ANAIS_PANEL_BUTTON e_button)
    {
        // No tooltip: Hide the window
        if((e_button == E_ANAIS_PANEL_BUTTON.NONE) || (_settings._showTips == false))
        {
            _tooltipWindow.Active = false;
            return;
        }

        // Other values: show tooltip with the appropriate text
        float posY = C_TOOLTIP_WINDOW_RELATIVE_POS_Y;
        string toolTipTitle = "";
        string toolTipText = "";

        switch(e_button)
        {
            case E_ANAIS_PANEL_BUTTON.TRANSFER_MODE_BUTTON:
                posY += (float)C_BUTTON_TOGGLE_TRANSFER_MODE_POS_Y;
                toolTipTitle = C_TOOLTIP_TRANSFER_MODE_TITLE;
                toolTipText = C_TOOLTIP_TRANSFER_MODE_TEXT;
                break;

            case E_ANAIS_PANEL_BUTTON.DECREASE_NB_MAX_TURNS_BUTTON:
                posY += (float)C_BUTTON_MINUS_ONE_POS_Y;
                toolTipTitle = C_TOOLTIP_DECREASE_NB_TURNS_TITLE;
                toolTipText = C_TOOLTIP_DECREASE_NB_TURNS_TEXT;
                break;

            case E_ANAIS_PANEL_BUTTON.INCREASE_NB_MAX_TURNS_BUTTON:
                posY += (float)C_BUTTON_PLUS_ONE_POS_Y;
                toolTipTitle = C_TOOLTIP_INCREASE_NB_TURNS_TITLE;
                toolTipText = C_TOOLTIP_INCREASE_NB_TURNS_TEXT;
                break;

            case E_ANAIS_PANEL_BUTTON.CHANGE_NODE_BUTTON:
                posY += (float)C_BUTTON_CHANGE_NODE_POS_Y;
                toolTipTitle = C_TOOLTIP_CHOOSE_NODE_TITLE;
                toolTipText = C_TOOLTIP_CHOOSE_NODE_TEXT;
                break;

            case E_ANAIS_PANEL_BUTTON.DEFAULT_SETTINGS_BUTTON:
                posY += (float)C_BUTTON_DEFAULT_SETTINGS_POS_Y;
                toolTipTitle = C_TOOLTIP_DEFAULT_SETTINGS_TITLE;
                toolTipText = C_TOOLTIP_DEFAULT_SETTINGS_TEXT;
                break;

            default: break;
        }

        // Show the tooltip with the appropriate parameters
        _tooltipWindow.Position = _mainWindow.Position + new Vector2(C_TOOLTIP_WINDOW_RELATIVE_POS_X, posY);
        _tooltipWindow.Title = toolTipTitle;
        _tooltipLabel.Text = toolTipText;
        _tooltipWindow.Active = true;
    }
}

// UTILITIES
// ---------
public static class GUI_Helper
{
    public static void SetAlignment(this SFS.UI.ModGUI.Button button, TMPro.TextAlignmentOptions alignment)
    {
        Traverse.Create(button).Field("_textAdapter").Field("TMProText").GetValue<TextMeshProUGUI>().alignment = alignment;
    }

    public static bool isMouseOver(this SFS.UI.ModGUI.Button button)
    {
        return Traverse.Create(button).Field("_button").Field("over").GetValue<bool>();
    }
}

