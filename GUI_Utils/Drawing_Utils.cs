using SFS.World;
using SFS.World.Maps;
using System;
using UnityEngine;

// UTILITY FUNCTIONS (used for diplay purpose)
// -----------------
public class Drawing_Utils
{
	public static string GetClosestApproachText(double approachDistance, double approachSpeed, double approachTime)
	{
		string closestApproachText = "closest approach: " + approachDistance.ToDistanceString();
		double durationBeforeApproach = approachTime - WorldTime.main.worldTime;

		const double C_MAX_DISTANCE_SHOW_TIMER = 10000.0;

		// Display the remaining time 
		if ((approachDistance < C_MAX_DISTANCE_SHOW_TIMER) && (durationBeforeApproach < 60.0) && (durationBeforeApproach > 0.0))
		{
			bool showDecimals = (approachSpeed < 10.0) ? true : false;
			string speedText = Units.ToVelocityString(approachSpeed, showDecimals);
			
			if (durationBeforeApproach < 1.0)
			{
				closestApproachText = closestApproachText + " (ΔV = " + speedText + "; T-0s)";
			}
			else
			{
				closestApproachText = closestApproachText + " (ΔV = " + speedText + "; T-" + Units.ToTimestampString(durationBeforeApproach, true, false) + ")";
			}
		}

		return closestApproachText;
	}

	public static Color GetClosestApproachLineColor(double distance)
	{
		if (distance > 100.0)
		{
			return new Color(0.627f, 0.784f, 1.0f, 0.8f); // light blue
		}
		else if (distance > 20.0)
		{
			// Ensures a smooth transition from orange to red
			Color orange = new Color(1.0f, 0.549f, 0.235f, 0.8f);
			Color red = new Color(1.0f, 0.1f, 0.0f, 0.8f);
			return Color.Lerp(red, orange, (float)((distance - 20.0) / (100.0 - 20.0)));
		}
		else
		{
			return new Color(1.0f, 0.1f, 0.0f, 0.8f); // red
		}
	}


	public static Color GetEfficiencyTransferColor(double efficiency)
    {
		Color blue = new Color(0.4f, 1.0f, 1.0f, 1.0f);
		Color green = new Color(0.6f, 1.0f, 0.8f, 1.0f);
		Color lemon = new Color(0.624f, 0.714f, 0.255f, 1.0f);

		const double C_TRANSITION_TO_BLUE_THRESHOLD = 0.95;
		const double C_TRANSITION_GREEN_TO_BLUE_THRESHOLD = 0.90;
		const double C_TRANSITION_TO_GREEN_THRESHOLD = 0.8;
		const double C_TRANSITION_LEMON_TO_GREEN_THRESHOLD = 0.75;

		if (efficiency > C_TRANSITION_TO_BLUE_THRESHOLD)
        {
			return blue;
        }
		else if(efficiency > C_TRANSITION_GREEN_TO_BLUE_THRESHOLD)
        {
			return Color.Lerp(green, blue, (float)((efficiency - C_TRANSITION_GREEN_TO_BLUE_THRESHOLD) /(C_TRANSITION_TO_BLUE_THRESHOLD - C_TRANSITION_GREEN_TO_BLUE_THRESHOLD)));
		}
		else if(efficiency > C_TRANSITION_TO_GREEN_THRESHOLD)
        {
			return green;
        }
		else if(efficiency > C_TRANSITION_LEMON_TO_GREEN_THRESHOLD)
        {
			return Color.Lerp(lemon, green, (float)((efficiency - C_TRANSITION_LEMON_TO_GREEN_THRESHOLD) / (C_TRANSITION_TO_GREEN_THRESHOLD - C_TRANSITION_LEMON_TO_GREEN_THRESHOLD)));
		}
		else
        {
			return lemon;
        }
	}

	public static void DrawDashedLine(Orbit orbit, Location start, Location end, Color color, string startText, string endText)
	{
		Vector3[] points = new Vector3[2];
		const double scaleMultiplier = 0.001;

		points[0].x = (float)(start.position.x * scaleMultiplier);
		points[0].y = (float)(start.position.y * scaleMultiplier);
		points[0].z = 0.0f;

		points[1].x = (float)(end.position.x * scaleMultiplier);
		points[1].y = (float)(end.position.y * scaleMultiplier);
		points[1].z = 0.0f;

		Map.dashedLine.DrawLine(points, orbit.Planet, color * new Color(1f, 1f, 1f, 0.5f), color * new Color(1f, 1f, 1f, 0.5f));

		Vector2 unitVector = (end.position - start.position).ToVector2.normalized;

		if (startText != null)
		{
			MapDrawer.DrawPointWithText(15, color, startText, 40, color, orbit.Planet.mapHolder.position + points[0], -unitVector, 4, 4);
		}
		if (endText != null)
		{
			MapDrawer.DrawPointWithText(15, color, endText, 40, color, orbit.Planet.mapHolder.position + points[1], unitVector, 4, 4);
		}
	}
}
