using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

using SFS.World;

public class OrbitDrawer
{
	public Orbit orbit;
	
	public OrbitDrawer(Orbit _orbit)
    {
		orbit = _orbit;
	}
	
	public Vector3[] GetPoints(double fromTrueAnomaly, double toTrueAnomaly, int resolution, double scaleMultiplier)
    {
        if(orbit.ecc < 1.0)
        {
			return GetPoints_Ellipse(fromTrueAnomaly, toTrueAnomaly, resolution, scaleMultiplier);
		}
        else
        {
			return GetPoints_Hyperbola(fromTrueAnomaly, toTrueAnomaly, resolution, scaleMultiplier);
        }
    }

	private Vector3[] GetPoints_Ellipse(double fromTrueAnomaly, double toTrueAnomaly, int resolution, double scaleMultiplier)
	{
		double startTrueAnomaly = Kepler.NormalizeAngle(fromTrueAnomaly);
		double endTrueAnomaly = Kepler.NormalizeAngle(toTrueAnomaly);

		bool fullTurn = (fromTrueAnomaly == toTrueAnomaly);

		double radius, trueAnomaly;
		int nbPoints;

		// Get the ellipse array, it lists the reference points from periapsis to apoapsis
		PolarPoint[] ellipseArray = GetEllipseArray(resolution);

		// Then search the first and last points from the ellipse to be added
		// ------------------------------------------------------------------

		// Those indexes are indexes for the ellipse array, but they can also be negative!
		// This is because the ellipse is browsed from apoapsis to periapsis to apoapsis.
		// The array only lists points from periapsis to apoapsis, the other part is deduced by symmetry
		// negative indexes refer to the symmetric part (apo -> peri), the corresponding true anomaly is simply negated
		int i_ellipse_start; // index of the first point following startTrueAnomaly
		int i_ellipse_stop;  // index of the first point following endTrueAnomaly

		if (orbit.direction > 0)
		{
			// search start index : i_ellipse_start
			i_ellipse_start = -(ellipseArray.Length - 1);
			trueAnomaly = -ellipseArray[Math.Abs(i_ellipse_start)].Argument;

			while (trueAnomaly < startTrueAnomaly && (i_ellipse_start < ellipseArray.Length - 1))
			{
				i_ellipse_start++;
				trueAnomaly = Math.Sign(i_ellipse_start) * ellipseArray[Math.Abs(i_ellipse_start)].Argument;
			}

			// search stop index : i_ellipse_stop
			if (fullTurn)
			{
				// full turn: start and stop indexes are equal
				i_ellipse_stop = i_ellipse_start;
				nbPoints = 2 * ellipseArray.Length + 1; // (2l-1) values for the ellipse array, 2 points for startTrueAnomaly and stopTrueAnomaly
			}
			else
			{
				// trajectory cut: search for the stop index
				i_ellipse_stop = -(ellipseArray.Length - 1);
				trueAnomaly = -ellipseArray[Math.Abs(i_ellipse_stop)].Argument;

				while (trueAnomaly < endTrueAnomaly && (i_ellipse_stop < ellipseArray.Length - 1))
				{
					i_ellipse_stop++;
					trueAnomaly = Math.Sign(i_ellipse_stop) * ellipseArray[Math.Abs(i_ellipse_stop)].Argument;
				}

				nbPoints = i_ellipse_stop - i_ellipse_start;
				if (nbPoints < 0) { nbPoints += 2 * ellipseArray.Length - 1; }
				nbPoints += 2; // number of points from the ellipse, plus 2 points for startTrueAnomaly and stopTrueAnomaly
			}
		}
		else // direction < 0
		{
			// search start index : i_ellipse_start
			i_ellipse_start = ellipseArray.Length - 1;
			trueAnomaly = ellipseArray[Math.Abs(i_ellipse_start)].Argument;

			while (trueAnomaly > startTrueAnomaly && (i_ellipse_start > -(ellipseArray.Length - 1)))
			{
				i_ellipse_start--;
				trueAnomaly = Math.Sign(i_ellipse_start) * ellipseArray[Math.Abs(i_ellipse_start)].Argument;
			}

			// search stop index : i_ellipse_stop
			if (fullTurn)
			{
				// full turn: start and stop indexes are equal
				i_ellipse_stop = i_ellipse_start;
				nbPoints = 2 * ellipseArray.Length + 1; // (2l-1) values for the ellipse array, 2 points for startTrueAnomaly and stopTrueAnomaly
			}
			else
			{
				// trajectory cut: search for the stop index
				i_ellipse_stop = ellipseArray.Length - 1;
				trueAnomaly = ellipseArray[Math.Abs(i_ellipse_stop)].Argument;

				while (trueAnomaly > endTrueAnomaly && (i_ellipse_stop > -(ellipseArray.Length - 1)))
				{
					i_ellipse_stop--;
					trueAnomaly = Math.Sign(i_ellipse_stop) * ellipseArray[Math.Abs(i_ellipse_stop)].Argument;
				}

				nbPoints = i_ellipse_start - i_ellipse_stop;
				if (nbPoints < 0) { nbPoints += 2 * ellipseArray.Length - 1; }
				nbPoints += 2; // number of points from the ellipse, plus 2 points for startTrueAnomaly and stopTrueAnomaly
			}
		}

		// ------------------------------
		// --- Now generate the array ---
		// ------------------------------
		Vector3[] array = new Vector3[nbPoints];
		int i_point = 0;

		// add "startTrueAnomaly"
		// ----------------------
		radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, startTrueAnomaly);
		array[i_point] = Kepler.GetPosition(radius, startTrueAnomaly, orbit.arg) * scaleMultiplier;
		i_point++;

		// Add the points from the ellipse array
		// -------------------------------------
		bool forceFirstIteration = fullTurn; // Otherwise the following loop doesn't run even once
		int i_ellipse = i_ellipse_start;

		while ((i_ellipse != i_ellipse_stop) || forceFirstIteration)
		{
			forceFirstIteration = false;

			trueAnomaly = Math.Sign(i_ellipse) * ellipseArray[Math.Abs(i_ellipse)].Argument;
			radius = ellipseArray[Math.Abs(i_ellipse)].Radius;

			array[i_point] = Kepler.GetPosition(radius, trueAnomaly, orbit.arg) * scaleMultiplier;
			i_point++;

			if (orbit.direction > 0)
			{
				i_ellipse++;
				if (i_ellipse == ellipseArray.Length) { i_ellipse = -(ellipseArray.Length - 1); }
			}
			else
			{
				i_ellipse--;
				if (i_ellipse == -ellipseArray.Length) { i_ellipse = ellipseArray.Length - 1; }
			}
		}

		// Finally, add "endTrueAnomaly"
		// -----------------------------
		radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, endTrueAnomaly);
		array[i_point] = Kepler.GetPosition(radius, endTrueAnomaly, orbit.arg) * scaleMultiplier;
		i_point++;

		return array;
	}

	public PolarPoint[] GetEllipseArray(int resolution)
	{
		const double C_MinAngleVariation = 0.000001; // Avoids creating too many extra points when eccentricity is close to 1

		double angularRes = 2.0 * Math.PI / (double)resolution;

		PolarPoint[] array = new PolarPoint[resolution];
		uint nbPoints = 0;
		double radius, previous_radius, previous_trueAnomaly;
		double C_MaxDeltaRadius = 0.05 * orbit.sma;

		// Start from apoapsis
		// -------------------
		double trueAnomaly = Math.PI;
		radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, trueAnomaly);

		// Store point under the form [radius, argument]
		array[nbPoints] = new PolarPoint(radius, trueAnomaly);
		nbPoints++;

		// calculate next true anomaly
		previous_trueAnomaly = trueAnomaly;
		trueAnomaly -= Math.Max((2.0 - radius / orbit.sma), C_MinAngleVariation) * angularRes; // calculate next true anomaly

		// Calculate and add points one by one
		// -----------------------------------
		while ((trueAnomaly > 0.0) && (nbPoints < resolution - 1)) // make sure to leave one free element to add the periapsis (though the limit should never be reached)
		{
			previous_radius = radius;

			radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, trueAnomaly);

			// If delta between previous radius and this one is too high, add intermediate points
			if (previous_radius - radius > C_MaxDeltaRadius)
			{
				double deltaRadius = previous_radius - radius;
				int nbAdditionalPoints = 0;

				// We basically cut the interval in half and add an intermediate point, until the radius difference between 2 points is satisfying
				while (deltaRadius > C_MaxDeltaRadius)
				{
					nbAdditionalPoints = 2 * nbAdditionalPoints + 1; // 1, 3, 7, 15...
					deltaRadius /= 2.0;
				}

				// difference in true anomaly between 2 successive points
				double deltaTrueAnomaly = (previous_trueAnomaly - trueAnomaly) / (nbAdditionalPoints + 1);

				// Complete the list with the number of necessary points
				for (int i_addPoint = 1; i_addPoint <= nbAdditionalPoints; i_addPoint++)
				{
					double temp_trueAnomaly = previous_trueAnomaly - i_addPoint * deltaTrueAnomaly;
					double temp_radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, temp_trueAnomaly);

					array[nbPoints] = new PolarPoint(temp_radius, temp_trueAnomaly);
					nbPoints++;
				}
			}

			// Store point under the form [radius, argument]
			array[nbPoints] = new PolarPoint(radius, trueAnomaly);
			nbPoints++;

			// calculate next true anomaly
			previous_trueAnomaly = trueAnomaly;
			trueAnomaly -= Math.Max((2.0 - radius / orbit.sma), C_MinAngleVariation) * angularRes; // calculate next true anomaly
		}

		// Add periapsis
		// -------------
		array[nbPoints] = new PolarPoint(orbit.periapsis, 0.0);
		nbPoints++;

		// Now create the final array, from periapsis to apoapsis
		// ------------------------------------------------------
		PolarPoint[] finalArray = new PolarPoint[nbPoints];

		for (uint i_point = 0; i_point < nbPoints; i_point++)
		{
			finalArray[i_point] = array[nbPoints - i_point - 1];
		}

		return finalArray;
	}


	// Also works for parabolas despite the name
	public Vector3[] GetPoints_Hyperbola(double fromTrueAnomaly, double toTrueAnomaly, int resolution, double scaleMultiplier)
	{
		double startTrueAnomaly = Kepler.NormalizeAngle(fromTrueAnomaly);
		double endTrueAnomaly = Kepler.NormalizeAngle(toTrueAnomaly);

		double radius, trueAnomaly;
		int nbPoints;

		double trueAnomaly_maxRadius = Math.Max(Math.Abs(startTrueAnomaly), Math.Abs(endTrueAnomaly));
		double radius_max = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, trueAnomaly_maxRadius);

		// Get the hyperbola array, it lists the reference points from periapsis to apoapsis
		PolarPoint[] hyperbolaArray = GetHyperbolaArray(resolution, radius_max);

		// Then search the first and last points from the ellipse to be added
		// ------------------------------------------------------------------

		// Those indexes are indexes for the hyperbola array, but they can also be negative!
		// This is because the hyperbola is browsed from apoapsis to periapsis to apoapsis.
		// The array only lists points from periapsis to apoapsis, the other part is deduced by symmetry
		// negative indexes refer to the symmetric part (apo -> peri), the corresponding true anomaly is simply negated
		int i_start; // index of the first point following startTrueAnomaly
		int i_stop;  // index of the first point following endTrueAnomaly

		if (orbit.direction > 0)
		{
			// search start index : i_start
			i_start = -(hyperbolaArray.Length - 1);
			trueAnomaly = -hyperbolaArray[Math.Abs(i_start)].Argument;

			while (trueAnomaly < startTrueAnomaly && (i_start < hyperbolaArray.Length - 1))
			{
				i_start++;
				trueAnomaly = Math.Sign(i_start) * hyperbolaArray[Math.Abs(i_start)].Argument;
			}

			// search stop index : i_stop
			i_stop = -(hyperbolaArray.Length - 1);
			trueAnomaly = -hyperbolaArray[Math.Abs(i_stop)].Argument;

			while (trueAnomaly < endTrueAnomaly && (i_stop < hyperbolaArray.Length - 1))
			{
				i_stop++;
				trueAnomaly = Math.Sign(i_stop) * hyperbolaArray[Math.Abs(i_stop)].Argument;
			}

			nbPoints = i_stop - i_start;
			nbPoints += 2; // number of points from the hyperbola, plus 2 point for startTrueAnomaly and stopTrueAnomaly.
		}
		else // direction < 0
		{
			// search start index : i_start
			i_start = hyperbolaArray.Length - 1;
			trueAnomaly = hyperbolaArray[Math.Abs(i_start)].Argument;

			while (trueAnomaly > startTrueAnomaly && (i_start > -(hyperbolaArray.Length - 1)))
			{
				i_start--;
				trueAnomaly = Math.Sign(i_start) * hyperbolaArray[Math.Abs(i_start)].Argument;
			}

			// search stop index : i_stop
			i_stop = hyperbolaArray.Length - 1;
			trueAnomaly = hyperbolaArray[Math.Abs(i_stop)].Argument;

			while (trueAnomaly > endTrueAnomaly && (i_stop > -(hyperbolaArray.Length - 1)))
			{
				i_stop--;
				trueAnomaly = Math.Sign(i_stop) * hyperbolaArray[Math.Abs(i_stop)].Argument;
			}

			nbPoints = i_start - i_stop;
			nbPoints += 2; // number of points from the hyperbola, plus 1 point for startTrueAnomaly and stopTrueAnomal.
		}

		// ------------------------------
		// --- Now generate the array ---
		// ------------------------------
		Vector3[] array = new Vector3[nbPoints];
		int i_point = 0;

		// add "startTrueAnomaly"
		// ----------------------
		radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, startTrueAnomaly);
		array[i_point] = Kepler.GetPosition(radius, startTrueAnomaly, orbit.arg) * scaleMultiplier;
		i_point++;

		// Add the points from the hyperbola array
		// -------------------------------------
		int i_hyperbola = i_start;

		while (i_hyperbola != i_stop)
		{
			trueAnomaly = Math.Sign(i_hyperbola) * hyperbolaArray[Math.Abs(i_hyperbola)].Argument;
			radius = hyperbolaArray[Math.Abs(i_hyperbola)].Radius;

			array[i_point] = Kepler.GetPosition(radius, trueAnomaly, orbit.arg) * scaleMultiplier;
			i_point++;

			if (orbit.direction > 0)
			{
				i_hyperbola++;
			}
			else
			{
				i_hyperbola--;
			}
		}

		// add "endTrueAnomaly"
		// --------------------
		radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, endTrueAnomaly);
		array[i_point] = Kepler.GetPosition(radius, endTrueAnomaly, orbit.arg) * scaleMultiplier;
		i_point++;

		return array;
	}


	public PolarPoint[] GetHyperbolaArray(int resolution, double r_max)
	{
		double angularRes = 2.0 * Math.PI / (double)resolution;
		double C_trueAnomalyMax = Math.Acos(-1.0 / orbit.ecc);

		// the angle at which we stop iterating over the angle (because we are getting dangerously close to the escape angle)
		double C_TrueAnomalyMax_ForIterations = C_trueAnomalyMax - 4.0 * angularRes;

		//double C_MaxDeltaRadius = 0.0625;
		double C_MaxDeltaRadius = Math.Sqrt(2.0) - 1.0;

		PolarPoint[] array = new PolarPoint[resolution];
		int nbPoints = 0;

		double radius, trueAnomaly;
		double previous_radius;

		// start at periapsis
		trueAnomaly = 0.0;
		radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, trueAnomaly);

		array[nbPoints] = new PolarPoint(radius, trueAnomaly);
		nbPoints++;

		// iterate
		// - trueAnomaly must not be too close to the max angle (C_trueAnomalyMax)
		// - radius must stay below r_max
		// - we must not overtake the max number of points
		while ((trueAnomaly < C_TrueAnomalyMax_ForIterations) && (radius < r_max) && (nbPoints < resolution))
		{
			previous_radius = radius;

			// calculate new points
			trueAnomaly += 2.0 * angularRes;
			radius = Kepler.GetRadiusAtTrueAnomaly(orbit.slr, orbit.ecc, trueAnomaly);

			// if new point too far from the previous one, add some intermediate points
			if (radius > (1.0 + C_MaxDeltaRadius) * previous_radius)
			{
				double deltaRadius = radius - previous_radius;
				int nbAdditionalPoints = 0;

				// We basically cut the interval in half and add an intermediate point, until the radius difference between 2 points is satisfying
				while (deltaRadius > C_MaxDeltaRadius * previous_radius)
				{
					nbAdditionalPoints = 2 * nbAdditionalPoints + 1; // 1, 3, 7, 15...
					deltaRadius /= 2.0;
				}

				// Complete the list with the number of necessary points
				for (int i_addPoint = 1; (i_addPoint <= nbAdditionalPoints) && (nbPoints < resolution); i_addPoint++)
				{
					double temp_radius = previous_radius + i_addPoint * deltaRadius;
					double temp_trueAnomaly = Kepler.GetTrueAnomalyAtRadius(orbit, temp_radius);

					array[nbPoints] = new PolarPoint(temp_radius, temp_trueAnomaly);
					nbPoints++;
				}
			}

			// Add the point
			if (nbPoints < resolution)
			{
				array[nbPoints] = new PolarPoint(radius, trueAnomaly);
				nbPoints++;
			}
		}

		// Now we either reached the maximum safe angle, or the max radius --> Make sure we stop at the max radius
		/*while((radius < r_max) && (nbPoints < resolution))
        {
			//radius = Math.Min(1.0625 * radius, r_max); // We increase the previous radius by a constant factor until we reach r_max

			radius = 1.0625 * radius;

			if (radius < r_max) // the radius max is not added: it will be managed later as the first or the last point
            {
				trueAnomaly = Kepler.GetTrueAnomalyAtRadius(this, radius);

				array[nbPoints] = new PolarPoint(radius, trueAnomaly);
				nbPoints++;
			}
		}*///(1.0 + C_MaxDeltaRadius) * previous_radius;

		double fraction_Rmax = 1.0 / Math.Pow(1.05, 4.0);
		bool Rmax1_reached = !(radius < fraction_Rmax * r_max);

		while (!Rmax1_reached && (nbPoints < resolution)) // add points until we reach approx. 82% of max radius
		{
			radius = (1.0 + C_MaxDeltaRadius) * radius;

			if (!(radius < fraction_Rmax * r_max))
			{
				radius = fraction_Rmax * r_max;
				Rmax1_reached = true;
			}

			trueAnomaly = Kepler.GetTrueAnomalyAtRadius(orbit, radius);

			array[nbPoints] = new PolarPoint(radius, trueAnomaly);
			nbPoints++;
		}

		bool Rmax2_reached = !(radius < r_max);

		while (!Rmax2_reached && (nbPoints < resolution))
		{
			radius = 1.05 * radius;

			if (!(radius < r_max))
			{
				radius = r_max;
				Rmax2_reached = true;
			}

			trueAnomaly = Kepler.GetTrueAnomalyAtRadius(orbit, radius);

			array[nbPoints] = new PolarPoint(radius, trueAnomaly);
			nbPoints++;
		}


		// Create the final array
		PolarPoint[] final_array = new PolarPoint[nbPoints];
		for (int i = 0; i < nbPoints; i++)
		{
			final_array[i] = array[i];
		}

		return final_array;
	}


	public PolarPoint[] GetHyperbolaInfiniteBranch_Array(double r_start, double r_max)
	{
		PolarPoint[] array = new PolarPoint[100];
		int nbPoints = 0;

		double r_intermediate = 0.8 * r_max;
		uint nbTotalPoints;

		// We allocate 4 points from r_start to r_intermediate (= 0.8*r_max), and 4 more points from r_intermediate to r_max
		if (r_intermediate < r_start)
		{
			nbTotalPoints = 8;
		}
		else
		{
			nbTotalPoints = 4;
		}

		PolarPoint[] array2 = new PolarPoint[nbTotalPoints];

		double radius;
		double trueAnomaly;

		if (r_intermediate < r_start)
		{
			double multiplier = Math.Pow(r_intermediate / r_start, 1.0 / 4.0);
			radius = r_start;

			for (int i = 0; i < 4; i++)
			{
				radius *= multiplier;
				trueAnomaly = Kepler.GetTrueAnomalyAtRadius(orbit, radius);

				array[nbPoints] = new PolarPoint(radius, trueAnomaly);
				nbPoints++;
			}
		}

		radius = r_intermediate;
		double multiplier2 = Math.Pow(r_max / r_intermediate, 1.0 / 4.0);

		for (int i = 0; i < 4; i++)
		{
			radius *= multiplier2;
			trueAnomaly = Kepler.GetTrueAnomalyAtRadius(orbit, radius);

			array[nbPoints] = new PolarPoint(radius, trueAnomaly);
			nbPoints++;
		}

		return array;
	}
}
