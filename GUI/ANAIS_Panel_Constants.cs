using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

class ANAIS_Panel_Constants
{
    // PANEL DEFINITION: DIMENSIONS + LABELS
    // -------------------------------------
    public const int C_HORIZONTAL_MARGIN = 10;   // How far widgets should be placed from the left and right borders

    // Separators
    public const int C_SEPARATOR_WIDTH = C_MAIN_WINDOW_WIDTH - 2 * C_HORIZONTAL_MARGIN; // Separators cover full window width, with a margin left on the left and on the right

    // Labels
    public const int C_LABEL_HEIGHT = 30;

    // Main window
    public const int C_MAIN_WINDOW_WIDTH = 440;
    public const int C_MAIN_WINDOW_HEIGHT = -(C_BUTTON_DEFAULT_SETTINGS_POS_Y - C_LABEL_HEIGHT - 60); // Adjusted from the size of the most downwards widget
    public const string C_MAIN_WINDOW_TITLE = "ANAIS Panel";

    // First separator position
    public const int C_FIRST_SEPARATOR_Y_POS = -10;

    // "Show transfer" label
    public const int C_LABEL_SHOW_TRANSFER_WIDTH = 200;
    public const int C_LABEL_SHOW_TRANSFER_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_SHOW_TRANSFER_POS_X = C_LABEL_SHOW_TRANSFER_WIDTH / 2 - C_MAIN_WINDOW_WIDTH / 2 + C_HORIZONTAL_MARGIN; // To the left, with a margin from the left side of the main window
    public const int C_LABEL_SHOW_TRANSFER_POS_Y = C_FIRST_SEPARATOR_Y_POS - 30;

    public const string C_LABEL_SHOW_TRANSFER_TEXT = "Show transfer:";

    // "Show transfer" toggle
    //public const int C_TOGGLE_SHOW_TRANSFER_POS_X = C_LABEL_SHOW_TRANSFER_POS_X + C_LABEL_SHOW_TRANSFER_WIDTH / 2 + (int)C_TOGGLE_SHOW_TRANSFER_SIZE_X / 2 + 40; // 40 pixels to the left of the "Show tranfer" label
    public const int C_TOGGLE_SHOW_TRANSFER_POS_X = C_MAIN_WINDOW_WIDTH / 2 - (int)C_TOGGLE_SHOW_TRANSFER_SIZE_X / 2 - C_HORIZONTAL_MARGIN; // 10 pixels from the right side
    public const int C_TOGGLE_SHOW_TRANSFER_POS_Y = C_LABEL_SHOW_TRANSFER_POS_Y;
    public const float C_TOGGLE_SHOW_TRANSFER_SIZE_X = 80.0f;
    public const float C_TOGGLE_SHOW_TRANSFER_SIZE_Y = 30.0f;

    // Transfer mode text label
    public const int C_LABEL_TEXT_TRANSFER_MODE_WIDTH = 180;
    public const int C_LABEL_TEXT_TRANSFER_MODE_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_TEXT_TRANSFER_MODE_POS_X = C_LABEL_TEXT_TRANSFER_MODE_WIDTH / 2 - C_MAIN_WINDOW_WIDTH / 2 + C_HORIZONTAL_MARGIN; // Placed on the left side of the main window, with a 10 pixels margin from the border
    public const int C_LABEL_TEXT_TRANSFER_MODE_POS_Y = C_LABEL_SHOW_TRANSFER_POS_Y - 45; // 45 pixels below the separator above

    public const string C_LABEL_TEXT_TRANSFER_MODE_TEXT = "Transfer mode:";

    // Transfer mode label
    public const int C_LABEL_TRANSFER_MODE_WIDTH = 160;
    public const int C_LABEL_TRANSFER_MODE_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_TRANSFER_MODE_POS_X = C_LABEL_TEXT_TRANSFER_MODE_POS_X + C_LABEL_TEXT_TRANSFER_MODE_WIDTH / 2 + C_LABEL_TRANSFER_MODE_WIDTH / 2 + 10; // To the right of previous label, with 10 pixels of margin
    public const int C_LABEL_TRANSFER_MODE_POS_Y = C_LABEL_TEXT_TRANSFER_MODE_POS_Y; // Same level as the label defined above

    public const string C_LABEL_TRANSFER_MODE_RENDEZ_VOUS_TEXT = "Rendez-vous";
    public const string C_LABEL_TRANSFER_MODE_FLY_BY_TEXT = "Fly-by";

    // Toggle transfer mode button
    public const int C_BUTTON_TOGGLE_TRANSFER_MODE_WIDTH = 60;
    public const int C_BUTTON_TOGGLE_TRANSFER_MODE_HEIGHT = 36;
    public const int C_BUTTON_TOGGLE_TRANSFER_MODE_POS_X = C_MAIN_WINDOW_WIDTH / 2 - C_BUTTON_TOGGLE_TRANSFER_MODE_WIDTH / 2 - C_HORIZONTAL_MARGIN; // 10 pixels from the right of the main window
    public const int C_BUTTON_TOGGLE_TRANSFER_MODE_POS_Y = C_LABEL_TEXT_TRANSFER_MODE_POS_Y - 3; // Same level as the labels defined above

    public const string C_BUTTON_TOGGLE_TRANSFER_MODE_TEXT = "◄►";


    // Target altitude text label
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_WIDTH = 160;
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_POS_X = C_LABEL_TEXT_TARGET_ALTITUDE_WIDTH / 2 - C_MAIN_WINDOW_WIDTH / 2 + C_HORIZONTAL_MARGIN; // Placed on the left side of the main window, with a 10 pixels margin from the border
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_POS_Y = C_LABEL_TEXT_TRANSFER_MODE_POS_Y - 45; // 45 pixels below the separator above

    public const string C_LABEL_TEXT_TARGET_ALTITUDE_TEXT = "Target altitude:";

    // Target altitude
    public const int C_TEXT_TARGET_ALTITUDE_WIDTH = 160;
    public const int C_TEXT_TARGET_ALTITUDE_HEIGHT = 35;
    public const int C_TEXT_TARGET_ALTITUDE_POS_X = C_LABEL_TEXT_TARGET_ALTITUDE_POS_X + C_LABEL_TEXT_TARGET_ALTITUDE_WIDTH / 2 + C_TEXT_TARGET_ALTITUDE_WIDTH / 2 + 30; // 20 pixels from the left of previous label
    public const int C_TEXT_TARGET_ALTITUDE_POS_Y = C_LABEL_TEXT_TARGET_ALTITUDE_POS_Y - 3;

    // Target altitude unit text label
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_WIDTH = 30;
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_POS_X = C_TEXT_TARGET_ALTITUDE_POS_X + C_TEXT_TARGET_ALTITUDE_WIDTH / 2 + C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_WIDTH / 2 + 20; // Placed on the left side of the main window, with a 10 pixels margin from the border
    public const int C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_POS_Y = C_LABEL_TEXT_TARGET_ALTITUDE_POS_Y; // 45 pixels below the separator above

    public const string C_LABEL_TEXT_TARGET_ALTITUDE_UNIT_TEXT = "Km";

    // Second separator
    public const int C_SECOND_SEPARATOR_Y_POS = C_LABEL_TEXT_TARGET_ALTITUDE_POS_Y - 40;

    // "Approach lines" label
    public const int C_LABEL_APPROACH_LINES_WIDTH = 200;
    public const int C_LABEL_APPROACH_LINES_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_APPROACH_LINES_POS_X = 0;
    public const int C_LABEL_APPROACH_LINES_POS_Y = C_SECOND_SEPARATOR_Y_POS - 25;

    public const string C_LABEL_APPROACH_LINES_TEXT = "Approach lines";

    // "Nb max turns" label
    public const int C_LABEL_NB_MAX_TURNS_WIDTH = 240;
    public const int C_LABEL_NB_MAX_TURNS_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_NB_MAX_TURNS_POS_X = C_LABEL_NB_MAX_TURNS_WIDTH / 2 - C_MAIN_WINDOW_WIDTH / 2 + C_HORIZONTAL_MARGIN; // 10 pixels from the left of the main window
    public const int C_LABEL_NB_MAX_TURNS_POS_Y = C_LABEL_APPROACH_LINES_POS_Y - 40;

    public const string C_LABEL_NB_MAX_TURNS_TEXT = "Nb. max turns:";

    // Nb max turns text input
    public const int C_TEXT_INPUT_NB_MAX_TURNS_WIDTH = 60;
    public const int C_TEXT_INPUT_NB_MAX_TURNS_HEIGHT = 35;
    public const int C_TEXT_INPUT_NB_MAX_TURNS_POS_X = C_LABEL_NB_MAX_TURNS_POS_X + C_LABEL_NB_MAX_TURNS_WIDTH / 2 + C_TEXT_INPUT_NB_MAX_TURNS_WIDTH / 2 + 20; // 20 pixels from the left of previous label
    public const int C_TEXT_INPUT_NB_MAX_TURNS_POS_Y = C_LABEL_NB_MAX_TURNS_POS_Y - 3;

    // "Minus 1" button
    public const int C_BUTTON_MINUS_ONE_WIDTH = 36;
    public const int C_BUTTON_MINUS_ONE_HEIGHT = 36;
    //public const int C_BUTTON_MINUS_ONE_POS_X = C_TEXT_INPUT_NB_MAX_TURNS_POS_X + C_TEXT_INPUT_NB_MAX_TURNS_WIDTH / 2 + C_BUTTON_MINUS_ONE_WIDTH / 2 + 20; // 20 pixels from the left of previous text input
    public const int C_BUTTON_MINUS_ONE_POS_X = C_BUTTON_PLUS_ONE_POS_X - C_BUTTON_MINUS_ONE_WIDTH / 2 - C_BUTTON_PLUS_ONE_WIDTH / 2 - 4; // 4 pixels from the left to "+" button
    public const int C_BUTTON_MINUS_ONE_POS_Y = C_TEXT_INPUT_NB_MAX_TURNS_POS_Y + 1;

    public const string C_BUTTON_MINUS_ONE_TEXT = "-";

    // "Plus 1" button
    public const int C_BUTTON_PLUS_ONE_WIDTH = 36;
    public const int C_BUTTON_PLUS_ONE_HEIGHT = 36;
    //public const int C_BUTTON_PLUS_ONE_POS_X = C_BUTTON_MINUS_ONE_POS_X + C_BUTTON_MINUS_ONE_WIDTH / 2 + C_BUTTON_PLUS_ONE_WIDTH / 2 + 4; // 4 pixels from the left of previous button
    public const int C_BUTTON_PLUS_ONE_POS_X = C_MAIN_WINDOW_WIDTH / 2 - C_BUTTON_PLUS_ONE_WIDTH / 2 - C_HORIZONTAL_MARGIN; // To the right side of the window
    public const int C_BUTTON_PLUS_ONE_POS_Y = C_TEXT_INPUT_NB_MAX_TURNS_POS_Y + 1;

    public const string C_BUTTON_PLUS_ONE_TEXT = "+";

    // "Preferred node" label
    public const int C_LABEL_PREFERRED_NODE_WIDTH = 180;
    public const int C_LABEL_PREFERRED_NODE_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_PREFERRED_NODE_POS_X = C_LABEL_PREFERRED_NODE_WIDTH / 2 - C_MAIN_WINDOW_WIDTH / 2 + C_HORIZONTAL_MARGIN; // 10 pixels from the left of the main window
    public const int C_LABEL_PREFERRED_NODE_POS_Y = C_LABEL_NB_MAX_TURNS_POS_Y - 45;

    public const string C_LABEL_PREFERRED_NODE_TEXT = "Preferred node:";

    public const string C_LABEL_PREFERRED_NODE_ANY = "Any";
    public const string C_LABEL_PREFERRED_NODE_ASCENDING = "Ascending";
    public const string C_LABEL_PREFERRED_NODE_DESCENDING = "Descending";

    // "Node" label
    public const int C_LABEL_NODE_WIDTH = 160;
    public const int C_LABEL_NODE_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_NODE_POS_X = C_LABEL_PREFERRED_NODE_POS_X + C_LABEL_PREFERRED_NODE_WIDTH / 2 + C_LABEL_NODE_WIDTH / 2 + 10; // 10 pixels from the left of the previous label
    public const int C_LABEL_NODE_POS_Y = C_LABEL_PREFERRED_NODE_POS_Y;

    // "Change node" button
    public const int C_BUTTON_CHANGE_NODE_WIDTH = 60;
    public const int C_BUTTON_CHANGE_NODE_HEIGHT = 36;
    public const int C_BUTTON_CHANGE_NODE_POS_X = C_MAIN_WINDOW_WIDTH / 2 - C_BUTTON_CHANGE_NODE_WIDTH / 2 - 10; // 10 pixels from the right side of the main window
    public const int C_BUTTON_CHANGE_NODE_POS_Y = C_LABEL_PREFERRED_NODE_POS_Y - 2;

    public const string C_BUTTON_CHANGE_NODE_TEXT = "◄►";

    // third separator
    public const int C_THIRD_SEPARATOR_Y_POS = C_LABEL_PREFERRED_NODE_POS_Y - 40;

    // "Show tips" label
    public const int C_LABEL_TIPS_WIDTH = 150;
    public const int C_LABEL_TIPS_HEIGHT = C_LABEL_HEIGHT;
    public const int C_LABEL_TIPS_POS_X = C_LABEL_TIPS_WIDTH / 2 - C_MAIN_WINDOW_WIDTH / 2 + C_HORIZONTAL_MARGIN; // To the left, with a margin from the left side of the main window
    public const int C_LABEL_TIPS_POS_Y = C_THIRD_SEPARATOR_Y_POS - 30;

    public const string C_LABEL_TIPS_TEXT = "Show tips:";

    // "Show tips" toggle
    public const int C_TOGGLE_TIPS_POS_X = C_LABEL_TIPS_POS_X + C_LABEL_TIPS_WIDTH / 2 + (int)C_TOGGLE_TIPS_SIZE_X / 2 + 10; // 10 pixels to the left of the "Show tips" label
    public const int C_TOGGLE_TIPS_POS_Y = C_LABEL_TIPS_POS_Y;
    public const float C_TOGGLE_TIPS_SIZE_X = 80.0f;
    public const float C_TOGGLE_TIPS_SIZE_Y = 30.0f;

    // Reset to default button
    public const int C_BUTTON_DEFAULT_SETTINGS_WIDTH = 140;
    public const int C_BUTTON_DEFAULT_SETTINGS_HEIGHT = 32;
    public const int C_BUTTON_DEFAULT_SETTINGS_POS_X = C_MAIN_WINDOW_WIDTH / 2 - C_BUTTON_DEFAULT_SETTINGS_WIDTH / 2 - C_HORIZONTAL_MARGIN; // 10 pixels from the right of the main window
    public const int C_BUTTON_DEFAULT_SETTINGS_POS_Y = C_LABEL_TIPS_POS_Y; // Same level as the labels defined above

    public const string C_BUTTON_DEFAULT_SETTINGS_TEXT = "Reset";

    // TOOLTIP WINDOW
    // --------------
    public const int C_LABEL_WINDOW_WIDTH = 800;
    public const int C_LABEL_WINDOW_HEIGHT = 300;
    public const int C_LABEL_TOOLTIP_WIDTH = C_LABEL_WINDOW_WIDTH - 20;
    public const int C_LABEL_TOOLTIP_HEIGHT = C_LABEL_WINDOW_HEIGHT - 20;

    public const float C_TOOLTIP_WINDOW_RELATIVE_POS_X = (float)(C_MAIN_WINDOW_WIDTH / 2 + C_LABEL_WINDOW_WIDTH / 2 + 10);
    public const float C_TOOLTIP_WINDOW_RELATIVE_POS_Y = (float)(C_LABEL_WINDOW_HEIGHT / 2 - 50);


    // TOOLTIPS
    // --------
    // For information: the tooltip window is wide for that line to fit precisely in it: "-----------------------------------------------------------------------------"
    public const string C_TOOLTIP_TRANSFER_MODE_TITLE = "Toggle transfer mode";
    public const string C_TOOLTIP_TRANSFER_MODE_TEXT = "- \"Rendez-vous\": ANAIS optimizes both your departure and arrival ΔV.\n" +
                                                       "      Recommended for most common situations; use this if you're not sure.\n" +
                                                       "- \"Fly-by\": ANAIS optimizes your departure ΔV only.\n" +
                                                       "      Use this if you don't intend to spend the ΔV upon arrival.\n" +
                                                       "      Recommended for: fly-bys, gravity assists, impact missions, if you\n" +
                                                       "      aim for a planet with an atmosphere and you intend to aerobrake.";

    public const string C_TOOLTIP_DECREASE_NB_TURNS_TITLE = "Decrease Nb max turns";
    public const string C_TOOLTIP_DECREASE_NB_TURNS_TEXT = "Decreases how many turns into the future ANAIS will calculate encounters.\n\n" +
                                                           "You will get encounters sooner, but you may miss more efficient ones farther in the future.\n\n" +
                                                           "Tip: set this to 0 to remove all approach lines";

    public const string C_TOOLTIP_INCREASE_NB_TURNS_TITLE = "Increase Nb max turns";
    public const string C_TOOLTIP_INCREASE_NB_TURNS_TEXT = "Increases how many turns into the future ANAIS will calculate encounters.\n\n" +
                                                           "A better encounter may be found farther in the future, but you will have to wait longer.";

    public const string C_TOOLTIP_CHOOSE_NODE_TITLE = "Choose node";
    public const string C_TOOLTIP_CHOOSE_NODE_TEXT = "This is for situations when your orbit crosses the target's orbit. ANAIS will calculate the best approach opportunity over several turns. This allows you to choose which node you want an encounter at.\n" +
                                                     "- Any: selects the best opportunity over both nodes; most common use.\n" +
                                                     "- Ascending: Best for gravity assists when you intend to lower your orbit.\n" +
                                                     "- Descending: Best for gravity assists when you intend to raise your orbit.";

    public const string C_TOOLTIP_DEFAULT_SETTINGS_TITLE = "Reset to default settings";
    public const string C_TOOLTIP_DEFAULT_SETTINGS_TEXT = "Reset all settings to default values.\n\n" +
                                                          "Default settings are fine for most cases, do not worry about changing settings if you do not understand them.";

}
