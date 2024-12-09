using HarmonyLib;
using SFS.World;
using SFS.WorldBase;
using SFS.World.Maps;
using UnityEngine;
using UnityEngine.UI;
using SFS.UI;
using System;
using System.Diagnostics;



[HarmonyPatch(typeof(VelocityArrowDrawer), "OnLocationChange")]
class VelocityArrowDrawer_OnLocationChange_Patch
{
	// Those are set to the 2 arrows instances from the VelocityArrowDrawer class; call InstantiateArrows() to instantiate
	private static VelocityArrowDrawer.Arrow arrow_X = null;
	private static VelocityArrowDrawer.Arrow arrow_Y = null;

	// NAVIGATION MODE
	// ---------------
	// This is a variable to tell how arrows shall be represented depending on the situation. Here is what each value means:
	// - DEFAULT                : The ingame original velocity arrows shall be used.
	// - CLOSE_TO_TARGET        : The target is close to the ship (< 10 km), we shall display the DV arrow and the distance arrow.
	// - FINAL_APPROACH         : A closest approach has been detected in the next 60 seconds, we shall display the DV arrow and the closest approach arrow.
	// - IMMINENT_ENCOUNTER     : The ship heads on the target at low speed. Display the DV arrow and the "Encounter" arrow.
	// - IMMINENT_IMPACT        : Same but at higher speed. Display the DV arrow and the "Impact" arrow.
	// - ANAIS_TRANSFER_PLANNED : ANAIS suggests a transfer. Display the specific transfer arrow.
	private enum E_NAV_MODE
	{
		DEFAULT,
		CLOSE_TO_TARGET,
		FINAL_APPROACH,
		IMMINENT_ENCOUNTER,
		IMMINENT_IMPACT,
		ANAIS_TRANSFER_PLANNED
	}

	private static E_NAV_MODE _navState = E_NAV_MODE.DEFAULT;

	// internal variables
	private static Double2 _relativePosition; // valid in all navigation modes but DEFAULT
	private static Double2 _relativeVelocity; // valid in all navigation modes but DEFAULT
	private static Double2 _closestApproachPosition;   // only valid in the following navigation modes: FINAL_APPROACH, IMMINENT_ENCOUNTER, IMMINENT_IMPACT
	private static double  _timeBeforeClosestApproach; // only valid in the following navigation modes: FINAL_APPROACH, IMMINENT_ENCOUNTER, IMMINENT_IMPACT
	private static bool _entrySOIdetected; // valid in all navigation modes but default

	private static double _lastEvaluationTime = 0.0;
	//private static DataSmoother _dataSmoother = new DataSmoother(25, 1.0, 4.0, 0.8);
    private static DataSmoother _dataSmoother = new DataSmoother(8, 1.0, 4.0, 0.8);

    // colors
    private static Color _defaultArrowColor = new Color(1.0f, 1.0f, 1.0f, 0.9019608f);

	private static DynamicColor _arrowX_color = new DynamicColor(_defaultArrowColor, _defaultArrowColor, DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
	private static DynamicColor _arrowY_color = new DynamicColor(_defaultArrowColor, _defaultArrowColor, DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
	private static DynamicColor _textX_color  = new DynamicColor(_defaultArrowColor, _defaultArrowColor, DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
	private static DynamicColor _textY_color  = new DynamicColor(_defaultArrowColor, _defaultArrowColor, DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);

	// constants
	private const double C_IMPACT_DISTANCE_THRESHOLD = 1.0;       // meters  - a target is considered on an impact trajectory if the closest approach is below that value
	private const double C_IMPACT_VELOCITY_THRESHOLD = 5.0;       // m/s     - below that speed value, we display an encounter instead of an impact
	private const double C_MIN_DISTANCE = 20.0;                   // meters  - in CLOSE_TO_TARGET mode, the distance arrow won't be showed below this value
	private const double C_MIN_VELOCITY = 0.1;                    // m/s     - the minimum velocity above which the velocity arrow (and the impact/encounter arrow) will be shown


	// Used to retrieve the game arrows
	public static void InstantiateArrows()
	{
		arrow_X = arrow_Y = null; // default value

		VelocityArrowDrawer theVelocityArrowDrawer = UnityEngine.Object.FindObjectOfType(typeof(VelocityArrowDrawer)) as VelocityArrowDrawer;
		if (theVelocityArrowDrawer != null)
		{
			arrow_X = theVelocityArrowDrawer.velocity_X;
			arrow_Y = theVelocityArrowDrawer.velocity_Y;
			LOG(LOG_LEVEL.INFO, "Velocity arrows correctly initialized");
		}
		else
        {
			LOG(LOG_LEVEL.ERROR, "Velocity arrows could not be instantiated");
		}
	}

	public static void notifyNewTarget(SelectableObject target)
	{
		_dataSmoother.Clear();
    }

	// FUNCTIONS TO EVALUATE THE NAVIGATION MODE
	// -----------------------------------------
	private static void evaluateNavigationState()
    {
		E_NAV_MODE previousNavState = _navState;

		// Initialize default value
		_navState = E_NAV_MODE.DEFAULT;

        if ((PlayerController.main.player.Value != null && PlayerController.main.player.Value.mapPlayer != null) &&    // If the player's rocket exists...
            (Map.navigation.target != null))                                                                           // ...and a target is defined
		{
            if(!Map.manager.mapMode.Value) 
			{
                // Set input data for ANAIS - Do this if view is in ship mode only, because in map mode it's already done in MapNavigation
                AnaisManager.SetInputData(PlayerController.main.player.Value.mapPlayer, Map.navigation.target);
            }

			// NOTE: we intentionally run the rest even if we are in map mode, because SmartSAS uses the data from there when used with ANAIS
			// ----

            // Retrieve output data
            bool hasData = AnaisManager.getVelocityArrowData(out bool finalApproach, out ClosestApproachCalculator.T_ApproachData approachData, out _relativeVelocity, out double evaluationTime, out _entrySOIdetected);

            if (hasData)
            {
                if (finalApproach)
                {
                    // FINAL APPROACH MODE
                    // -------------------
                    _relativePosition = Map.navigation.target.Location.position - PlayerController.main.player.Value.mapPlayer.Location.position;
                    _relativeVelocity = Map.navigation.target.Location.velocity - PlayerController.main.player.Value.mapPlayer.Location.velocity;

                    if (approachData.validity == false)
                    {
                        // We don't have a closest approach in the next 60 seconds anymore but we are still close to the target
						_navState = E_NAV_MODE.CLOSE_TO_TARGET;
                    }
					else
					{
                        // We have an imminent encounter
						_closestApproachPosition = approachData.locTarget.position - approachData.locPlayer.position;
                        _timeBeforeClosestApproach = approachData.date - WorldTime.main.worldTime;

                        if (approachData.dist < C_IMPACT_DISTANCE_THRESHOLD) // approach distance is very low: Impact trajectory
                        {
                            if (_relativeVelocity.magnitude < C_IMPACT_VELOCITY_THRESHOLD)
                            {
                                _navState = E_NAV_MODE.IMMINENT_ENCOUNTER;
                            }
                            else
                            {
                                _navState = E_NAV_MODE.IMMINENT_IMPACT;
                            }
                        }
                        else
                        {
                            _navState = E_NAV_MODE.FINAL_APPROACH;
                        }
                    }
                }
                else
                {
                    // We are far from target
					_navState = E_NAV_MODE.ANAIS_TRANSFER_PLANNED;

					if(_lastEvaluationTime != evaluationTime)
					{
                        _dataSmoother.Add(_relativeVelocity);
						_lastEvaluationTime = evaluationTime;
                    }
					
					_relativeVelocity = _dataSmoother.GetSmoothedData();
                }
            }
        }

        // handle navigation state change if needed
        if (previousNavState != _navState) OnNavigationModeChanged(previousNavState, _navState);
	}

	private static void OnNavigationModeChanged(E_NAV_MODE oldNavState, E_NAV_MODE newNavState)
    {
		// Navigation mode unchanged, nothing to do
		if (oldNavState == newNavState) return;

		// If the new mode is DEFAULT, simply set to all objects the default color since our code won't handle it.
		if(newNavState == E_NAV_MODE.DEFAULT)
        {
			setVelocityArrowColor(arrow_X, _defaultArrowColor, _defaultArrowColor);
			setVelocityArrowColor(arrow_Y, _defaultArrowColor, _defaultArrowColor);
			return;
		}

		// OTHER CASES: arrow X
		// --------------------
		// transition from DEFAULT or ANAIS TRANSFER to one of the encounter modes
		if( ((oldNavState == E_NAV_MODE.DEFAULT) || (oldNavState == E_NAV_MODE.ANAIS_TRANSFER_PLANNED)) && (newNavState != E_NAV_MODE.ANAIS_TRANSFER_PLANNED))
        {
            // arrowX is the velocity arrow: displayed in blinking light blue in this case, to tell the player that he should burn in that direction
            // If the previous mode was something else than DEFAULT or ANAIS_TRANSFER_PLANNED we don't do this, because we use the same color and it would restart the blinking process.
            _arrowX_color.changeColor(DynamicColor.E_COLOR.LIGHT_BLUE_2);
			_arrowX_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_SMOOTH, 1.8f);
			_arrowX_color.restartBlinkingProcess();
			_textX_color.changeColor(DynamicColor.E_COLOR.LIGHT_BLUE);
			_textX_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_SMOOTH, 1.8f);
			_textX_color.restartBlinkingProcess();

            _dataSmoother.Clear();
        }

		if(newNavState == E_NAV_MODE.ANAIS_TRANSFER_PLANNED)
        {
			// arrowX is the velocity arrow: displayed in blinking cyan in this case, to tell the player that he should burn in that direction
			// This is applied no matter what was the previous state because this color is specific to this mode
			_arrowX_color.changeColor(DynamicColor.E_COLOR.CYAN_2);
			_arrowX_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_SMOOTH, 1.8f);
			_arrowX_color.restartBlinkingProcess();
			_textX_color.changeColor(DynamicColor.E_COLOR.CYAN);
			_textX_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_SMOOTH, 1.8f);
			_textX_color.restartBlinkingProcess();
		}


		// OTHER CASES: arrow Y
		// --------------------
		// transition to CLOSE_TO_TARGET mode
		if (newNavState == E_NAV_MODE.CLOSE_TO_TARGET)
        {
			// arrowY is the distance arrow: displayed in non-blinking light green in this case
			_arrowY_color.changeColor(DynamicColor.E_COLOR.LIGHT_GREEN_2);
			_arrowY_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
			_textY_color.changeColor(DynamicColor.E_COLOR.LIGHT_GREEN);
			_textY_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
			return;
		}

		// transition to FINAL APPROACH / IMMINENT ENCOUNTER mode
		if((newNavState == E_NAV_MODE.FINAL_APPROACH) || (newNavState == E_NAV_MODE.IMMINENT_ENCOUNTER))
        {
			// arrowY is the closest approach arrow: displayed in non-blinking light red in this case
			_arrowY_color.changeColor(DynamicColor.E_COLOR.LIGHT_RED_2);
			_arrowY_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
			_textY_color.changeColor(DynamicColor.E_COLOR.LIGHT_RED);
			_textY_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_NONE, 0.0f);
			return;
		}

		// transition to IMMINENT IMPACT mode
		if (newNavState == E_NAV_MODE.IMMINENT_IMPACT)
		{
			// arrowY is the closest approach arrow: displayed in non-blinking red in this case
			_arrowY_color.changeColor(DynamicColor.E_COLOR.RED_2);
			_arrowY_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_BINARY, 1.2f);
			_textY_color.changeColor(DynamicColor.E_COLOR.RED);
			_textY_color.changeBlinkingMode(DynamicColor.E_BLINKING_MODE.BLINKING_BINARY, 1.2f);
			return;
		}

        // transition to ANAIS_TRANSFER_PLANNED mode
        if (newNavState == E_NAV_MODE.ANAIS_TRANSFER_PLANNED)
        {
            // arrowY isn't displayed in this case - color resetted to default
            setVelocityArrowColor(arrow_Y, _defaultArrowColor, _defaultArrowColor);
            return;
        }
    }


	// FUNCTIONS THAT RETURN THE STRINGS TO BE DISPLAYED NEXT TO EACH ARROW
	// --------------------------------------------------------------------
	public static string getArrowX_String(Color color)
    {
		string label = "";

		if (_navState != E_NAV_MODE.DEFAULT)
		{
			if(_entrySOIdetected && (_navState == E_NAV_MODE.ANAIS_TRANSFER_PLANNED))
			{
				label = "(Encounter)"; // Show "Encounter" instead of speed if an encounter is planned, to let the player know that he should switch to map mode
			}
			else
			{
                double speed = _relativeVelocity.magnitude;

                bool showDecimals = (speed < 10.0) ? true : false;
                string speedText = Units.ToVelocityString(speed, showDecimals);

                label = "ΔV = " + speedText;
            }
			
		}

		return label;
	}

	public static string getArrowY_String(Color color)
    {
		string label = "";

		if (_navState == E_NAV_MODE.CLOSE_TO_TARGET)
		{
			float dist = (float)_relativePosition.magnitude;

			bool showDecimals = (dist < 100.0) ? true : false;
			string distText = Units.ToDistanceString(dist, showDecimals);

			label = "Distance = " + distText;
		}
		else if((_navState == E_NAV_MODE.FINAL_APPROACH    ) || 
			    (_navState == E_NAV_MODE.IMMINENT_ENCOUNTER) || 
				(_navState == E_NAV_MODE.IMMINENT_IMPACT   )   )
		{
			// Convert remaining time to string
			string str_remainingSeconds;

			if (_timeBeforeClosestApproach > 1.0) str_remainingSeconds = Units.ToTimestampString(_timeBeforeClosestApproach, true, false);
			else str_remainingSeconds = "0s";

			// compute string to display in each situation
			if (_navState == E_NAV_MODE.FINAL_APPROACH)
			{
				float dist = (float)_closestApproachPosition.magnitude;

				bool showDecimals = (dist < 100.0) ? true : false;
				string distText = Units.ToDistanceString(dist, showDecimals);

				label = "Closest approach = " + distText + " (T-" + str_remainingSeconds + ")";
			}
			else if (_navState == E_NAV_MODE.IMMINENT_IMPACT)
			{
				label = "IMPACT!! (T-" + str_remainingSeconds + ")";
			}
			else if (_navState == E_NAV_MODE.IMMINENT_ENCOUNTER)
			{
				label = "Encounter (T-" + str_remainingSeconds + ")";
			}
		}

		return label;
    }


	// SETS THE VELOCITY ARROW COLOR
	// -----------------------------
	private static void setVelocityArrowColor(VelocityArrowDrawer.Arrow arrow, Color arrowColor, Color textColor)
    {
		arrow.line.color = arrowColor;
		arrow.text.color = textColor;

		try
        {
			// That code will work as long as Stef doesn't modify the transforms (hence why it's in a try catch)
			arrow.line.transform.parent.GetChild(0).gameObject.GetComponent<Image>().color = arrowColor; // Corresponds to Arrow X: tip of the arrow
			arrow.line.transform.GetChild(1).gameObject.GetComponent<Image>().color = arrowColor; // Corresponds to Base X: origin of the arrow
		}
		catch (Exception e)
        {
			LOG(LOG_LEVEL.ERROR, "Velocity arrows components could not be retrieved");
		}
	}


	// FUNCTIONS TO DETERMINE THE ARROWS LENGTH
	// ----------------------------------------
	private static float GetGenericArrowLength(float ref_val, float val)
    {
		// ref_val is the value of val around which the length variation will be maximal
		float x = val / ref_val;

		// the value returned is always between 0 and 1, it's a normalized length
		if(x < 1.0)
        {
			return 0.25f * x * (x + 1.0f);
        }
		else
        {
			return (1.0f - 1.0f / (3.0f * x - 1.0f));
        }
    }
	
	private static float GetVelocityArrowLength(double velocity)
    {
		return 0.2f + 0.8f * GetGenericArrowLength(50.0f, (float)velocity);
	}

	private static float GetDistanceArrowLength(double dist)
	{
		return 0.2f + 0.8f * GetGenericArrowLength(200.0f, (float)dist);
	}

	private static float GetArrowLengthFactor()
    {
		double maxRadius = 0.9 * GameCamerasManager.main.world_Camera.camera.pixelHeight / 2.0;
		return (float)(0.4 * maxRadius);
	}

	private static Vector2 getArrowPos(Double2 val, double length, double distanceFromOrigin)
    {
		Double2 center;

		center.x = GameCamerasManager.main.world_Camera.camera.pixelWidth / 2.0;
		center.y = GameCamerasManager.main.world_Camera.camera.pixelHeight / 2.0;

		double maxRadius = 0.9 * GameCamerasManager.main.world_Camera.camera.pixelHeight / 2.0;

		// Rotate the velocity vector so that the calculated position matches the camera orientation
		val = val.Rotate(-GameCamerasManager.main.world_Camera.CameraRotationRadians);

		double magnitude = val.magnitude;

		double distFromOrigin = distanceFromOrigin * maxRadius + 0.5f * length;

		Vector2 hitPos;
		hitPos.x = (float)(center.x + distFromOrigin * val.x / magnitude);
		hitPos.y = (float)(center.y + distFromOrigin * val.y / magnitude);

		return hitPos;
	}


	// COMPUTE AND DISPLAY THE ARROWS (the navigation mode mustn't be DEFAULT!)
	// ------------------------------
	private static void displayDeltaVarrow(Location location)
	{
		// the length desired for arrows, in percentage of the maximal arrow length
		const float C_VELOCITY_ARROW_LENGTH_FACTOR = 0.5f;
		const float C_DISTANCE_ARROW_LENGTH_FACTOR = 0.65f;

		// the distance of the arrow origin from the center of the screen, in percentage of the maximal length
		const float C_VELOCITY_ARROW_DISTANCE_FROM_CENTER = 0.25f;
		const float C_DISTANCE_ARROW_DISTANCE_FROM_CENTER = 0.6f;

		// If player doesn't control or isn't in ship view, don't display the arrows
		if (!(PlayerController.main.player.Value is Rocket) || (bool)Map.manager.mapMode)
		{
			arrow_X.SetActive(active: false);
			arrow_Y.SetActive(active: false);
			return;
		}

		// If zoom level is too low, don't display the arrows
		float sizeRadius = PlayerController.main.player.Value.GetSizeRadius();
		if ((float)WorldView.main.viewDistance > sizeRadius * 50f + 50f)
		{
			arrow_X.SetActive(active: false);
			arrow_Y.SetActive(active: false);
			return;
		}

		// DISPLAY ARROW X (if active)
		// ---------------
		float speed = (float)_relativeVelocity.magnitude;

		if (speed > C_MIN_VELOCITY)
		{
			// Determine the arrow position and length
			float velocityArrowLength = GetVelocityArrowLength(speed) * GetArrowLengthFactor() * C_VELOCITY_ARROW_LENGTH_FACTOR;
			Double2 directionNormal = _relativeVelocity.Rotate(0f - GameCamerasManager.main.world_Camera.CameraRotationRadians) / speed;
			Vector2 arrowPos = getArrowPos(_relativeVelocity, velocityArrowLength, 0.4f);
			
			// display the arrow
			setVelocityArrowColor(arrow_X, _arrowX_color.getColor(), _textX_color.getColor());
			arrow_X.Position(getArrowX_String, velocityArrowLength, arrowPos, directionNormal);
		}
		else
		{
			arrow_X.SetActive(active: false);
		}

		// DISPLAY ARROW Y
		// ---------------

		// Determine the arrow position in each situation
		if (_navState == E_NAV_MODE.CLOSE_TO_TARGET)
		{
			float dist = (float)_relativePosition.magnitude;

			if (dist > C_MIN_DISTANCE)
			{
				// Determine the arrow position and length
				float distanceArrowLength = GetDistanceArrowLength(dist) * GetArrowLengthFactor() * C_DISTANCE_ARROW_LENGTH_FACTOR;
				Double2 directionNormal = _relativePosition.Rotate(0f - GameCamerasManager.main.world_Camera.CameraRotationRadians) / dist;
				Vector2 arrowPos = getArrowPos(_relativePosition, distanceArrowLength, C_DISTANCE_ARROW_DISTANCE_FROM_CENTER);

				// display the arrow
				setVelocityArrowColor(arrow_Y, _arrowY_color.getColor(), _textY_color.getColor());
				arrow_Y.Position(getArrowY_String, distanceArrowLength, arrowPos, directionNormal);
			}
			else
            {
				arrow_Y.SetActive(false);
			}
		}
		else if(_navState == E_NAV_MODE.FINAL_APPROACH)
        {
			float dist = (float)_closestApproachPosition.magnitude;

			if (dist > C_IMPACT_DISTANCE_THRESHOLD)
			{
				// Determine the arrow position and length
				float distanceArrowLength = GetDistanceArrowLength(dist) * GetArrowLengthFactor() * C_DISTANCE_ARROW_LENGTH_FACTOR;
				Double2 directionNormal = _closestApproachPosition.Rotate(0f - GameCamerasManager.main.world_Camera.CameraRotationRadians) / dist;
				Vector2 arrowPos = getArrowPos(_closestApproachPosition, distanceArrowLength, C_DISTANCE_ARROW_DISTANCE_FROM_CENTER);

				// display the arrow
				setVelocityArrowColor(arrow_Y, _arrowY_color.getColor(), _textY_color.getColor());
				arrow_Y.Position(getArrowY_String, distanceArrowLength, arrowPos, directionNormal);
			}
            else 
			{
				arrow_Y.SetActive(false);
			}
		}
		else if((speed > C_MIN_VELOCITY) && ((_navState == E_NAV_MODE.IMMINENT_ENCOUNTER) || (_navState == E_NAV_MODE.IMMINENT_IMPACT)))
        {
			// Determine the arrow position (arrow length is minimal)
			float distanceArrowLength = GetDistanceArrowLength(0.0) * GetArrowLengthFactor() * C_DISTANCE_ARROW_LENGTH_FACTOR;
			Double2 directionNormal = -_relativeVelocity.Rotate(0f - GameCamerasManager.main.world_Camera.CameraRotationRadians) / speed; // As direction, we take the opposite of velocity direction
			Vector2 arrowPos = getArrowPos(directionNormal, distanceArrowLength, C_DISTANCE_ARROW_DISTANCE_FROM_CENTER);

			// display the arrow
			setVelocityArrowColor(arrow_Y, _arrowY_color.getColor(), _textY_color.getColor());
			arrow_Y.Position(getArrowY_String, distanceArrowLength, arrowPos, directionNormal);
		}
        else
        {
			arrow_Y.SetActive(false);
		}
	}

	

	[HarmonyPrefix]
	public static bool OnLocationChange_Prefix(Location _, Location location)
	{
		// Compute and display custom velocity arrows (if applicable)
		// ------------------------------------------
		if ((arrow_X != null) && (arrow_Y != null) /*&& !Map.manager.mapMode.Value*/) // arrows exist and view is in ship mode
		{
			evaluateNavigationState();

			if (!Map.manager.mapMode.Value) // view is in ship mode
            {
				if (_navState != E_NAV_MODE.DEFAULT)
				{
					try
					{
						displayDeltaVarrow(location);
						return false;
					}
					catch (Exception ex)
					{
						LOG(LOG_LEVEL.ERROR, "Exception : " + ex.Message);
						return true;
					}
				}
				else
				{
					// default mode: call usual code
					return true;
				}
			}
		}

		// call the original method
		return true;
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.VELOCITY_ARROW, level, message);
	}
}


[HarmonyPatch(typeof(VelocityArrowDrawer), "Start")]
class VelocityArrowDrawer_Start_Patch
{
	[HarmonyPostfix]
	public static void Start_postfix()
    {
		// Initialize the reference on arrows each time the game is started
		VelocityArrowDrawer_OnLocationChange_Patch.InstantiateArrows();
	}
}

