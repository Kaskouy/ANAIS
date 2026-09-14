using SFS.World;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using UnityEngine;


// --------------------------------------------------------------------------------------------------------------------
// -                                             CLASS OrbitDrawer                                                    -
// --------------------------------------------------------------------------------------------------------------------
// This class encapsulates the trajectory calculation.
// Use it as follows:
// 1) Instantiate the OrbitDrawer by calling the constructor.
//    The constructor builds an internal array of precalculated reference points that will be used as a basis to render 
//    all future trajectories.
// 2) On each frame, call GetPoints to quickly calculate a trajectory corresponding to the specified starting and
//    ending true anomaly. This largely uses the precalculated points, so the process is fast.
// --------------------------------------------------------------------------------------------------------------------
//                                             INTERNAL ALGORITHM
//
// The OrbitDrawer first creates an array of evenly dispatched points, so that the result looks as smooth as possible.
// This is done by two means:
// - Points are spaced so that the tangent to the curve has turned by a fixed angle between two consecutive points.
//   This ensures that more points are allocated to highly curved areas, so that the result looks smooth.
//   For this, all calculations are based on the tangential angle (usually noted psi)
// - A tessellation algorithm is applied to add intermediate points to the portions of the curve that still show an
//   important loss in terms of quality. This is based on the surface between the conic arc and the segment joining
//   two consecutive points. If this surface appears to be important compared to the whole conic area, an intermediate
//   point is added until the rendering is satisfying, or the limit has been reached in terms of number of points.
//
// Restrictions: This cannot be used for rectilinear trajectories.
// Outside of this, the algorithm used is universal and works with no instability for ellipses, parabolas, hyperbolas.
// --------------------------------------------------------------------------------------------------------------------
public class OrbitDrawer
{
    // ----------------------------------------------------------------------------------------------------------------
    // -                                                                                                              -
    // -                                         INTERNAL VARIABLES                                                   -
    // -                                                                                                              -
    // ----------------------------------------------------------------------------------------------------------------
    public Orbit orbit;

	private uint nbMaxSegments;
    private double scaleMultiplier;
    private double psiMax;
    private double deltaPsi;
    private double focusDistance = 0.0f;

    private uint nbMaxTessellationSteps = 0;
    private double toleranceThresholdForTessellation = 0.0;

    // cached values to avoid recalculating them many times
    private double cosPhi;
    private double sinPhi;
    private bool closedConic = false;

    // The base point list
    private OrbitPoint[] orbitPointList = null;


    // ----------------------------------------------------------------------------------------------------------------
    // -                                                                                                              -
    // -                                              CONSTANTS                                                       -
    // -                                                                                                              -
    // ----------------------------------------------------------------------------------------------------------------

    // Used to force a somewhat reasonable number of points
    private const int C_MIN_NB_SEGMENTS = 4;
    private const int C_MAX_NB_SEGMENTS = 1800;

    // The minimum angle between two consecutive trajectory segments
    // The algorithm won't create additional points if the angle between the new segments would be below that.
    // NOTE : the value is chosen so that C_MAX_NB_SEGMENTS * C_MIN_DELTA_PSI remains below 2*Pi
    private const double C_MIN_DELTA_PSI = Math.PI / (double)C_MAX_NB_SEGMENTS; // A tenth of a degree



    // ----------------------------------------------------------------------------------------------------------------
    // -                                                                                                              -
    // -                                       INTERNAL STRUCTURES                                                    -
    // -                                                                                                              -
    // ----------------------------------------------------------------------------------------------------------------

    // struct OrbitBasePoint
    // ---------------------
    // Represents an orbit point by its cartesian coordinates in double precision and its associated psi angle.
    // Used for intermediate calculations
    private struct OrbitBasePoint
    {
        public double psi;
        public double x;
        public double y;

        public OrbitBasePoint(double psi, double x, double y)
        {
            this.psi = psi;
            this.x = x;
            this.y = y;
        }
    }

    // struct OrbitPoint
    // -----------------
    // Represents an orbit point by its cartesian coordinates in the game format and its associated psi angle.
    // This is to be stored to be ready for final use.
    private struct OrbitPoint
    {
        public double psi;
        public Vector3 point;

        public OrbitPoint(double psi, Vector3 point)
        {
            this.psi = psi;
            this.point = point;
        }

        public void Set(double psi, Vector3 point)
        {
            this.psi = psi;
            this.point = point;
        }
    }


    // ----------------------------------------------------------------------------------------------------------------
    // -                                                                                                              -
    // -                                               METHODS                                                        -
    // -                                                                                                              -
    // ----------------------------------------------------------------------------------------------------------------

    // ----------------------------------------------------------------------------------------------------------------
    // Constructor
    // ----------------------------------------------------------------------------------------------------------------
    // Initializes all internal data, and builds the base array to serve as a basis to speed up future calculations.
    // Parameters are as follow:
    // - orbit: the orbit object for which we should compute a trajectory.
    // - resolution: will be interpreted as the number of segments that compose a conic - will be rounded down to the
    //               closest multiple of 4.
    // - scaleMultiplier: the scale to apply to calculated points.
    // - areaLossToleranceFactor: tells the program how many times the reference area of a conic you are ready to lose.
    //                            If more is lost, the program will use the tesselation process to add more points and
    //                            reduce the area visually lost.
    // - nbMaxTessellationSteps: how many tessellation steps you want to use at most. Tessellation is the process of
    //                           adding intermediate points so that a curve looks smoother. 3 is enough in practice.
    // ----------------------------------------------------------------------------------------------------------------
    public OrbitDrawer(Orbit orbit, int resolution, double scaleMultiplier, double areaLossToleranceFactor = 10.0, uint nbMaxTessellationSteps = 0)
    {
        this.orbit = orbit;
        this.scaleMultiplier = scaleMultiplier;
        this.nbMaxTessellationSteps = nbMaxTessellationSteps;

        // To avoid that somebody breaks my code by entering stupid values...
        if (resolution < 4) resolution = C_MIN_NB_SEGMENTS;
        if (resolution > 1800) resolution = C_MAX_NB_SEGMENTS;

        // The (maximum) number of segments of which the conic will be made of (rounded down to the nearest multiple
        // of 4 so that we can exploit all the symmetries of an ellipse).
        nbMaxSegments = ((uint)resolution / 4) * 4;

        // Get psiMax, the maximum tangential angle (which is either PI in the case of a full ellipse, or the value that corresponds to the edge of a SOI)
        psiMax = getPsiMax();

        closedConic = (psiMax == Math.PI);

        // the difference of tangential angle between two consecutive points
        deltaPsi = 2.0f * Math.PI / nbMaxSegments;

        // Calculate focus distance (for ellipses only - it won't be used in other cases)
        if (orbit.ecc < 1.0f)
        {
            focusDistance = 2.0f * orbit.ecc * orbit.slr / (1.0f - orbit.ecc * orbit.ecc) * scaleMultiplier;
        }

        // precalculate this values
        cosPhi = Math.Cos(orbit.arg);
        sinPhi = Math.Sin(orbit.arg);

        // calculate the tolerance threshold for the tessellation step
        if(nbMaxTessellationSteps > 0)
        {
            toleranceThresholdForTessellation = 2.0 * Math.Pow(areaLossToleranceFactor * Math.PI * Math.PI * getConicArea(), 1.0 / 3.0) / nbMaxSegments;
        }

        LOG(LOG_LEVEL.DEBUG, "Internal data : psiMax = " + psiMax + "; nbMaxSegments = " + nbMaxSegments + "; deltaPsi = " + deltaPsi + "; focusDistance = " + focusDistance);

        // Build the base points - this will serve as a basis for all future trajectories calculations
        orbitPointList = BuildFullBaseArray();
    }


    // ----------------------------------------------------------------------------------------------------------------
    // getRMax
    // ----------------------------------------------------------------------------------------------------------------
    // Returns the maximum radius up to which it's relevant to draw the trajectory.
    // It's the SOI limit in practice, or determined by the camera position if the SOI is infinite.
    // ----------------------------------------------------------------------------------------------------------------
    private double getRMax()
	{
		if(Double.IsInfinity(orbit.Planet.SOI))
		{
            SFS.World.Maps.MapView.View worldView = SFS.World.Maps.Map.view.view;
            return ((worldView.target.Value.Location.GetSolarSystemPosition(WorldTime.main.worldTime) + worldView.position.Value).magnitude + worldView.distance.Value * (double)SFS.Cameras.ActiveCamera.Camera.ViewSizeNormal * 1.0625);
        }
		else
		{
			return orbit.Planet.SOI;
        }
	}


    // ----------------------------------------------------------------------------------------------------------------
    // getPsiMax
    // ----------------------------------------------------------------------------------------------------------------
    // Returns the maximum psi angle for the given trajectory. It takes into account the SOI limit to decide about the
	// psi max angle. It also works if there's no limit, but a small value is substracted so that using that value for
	// calculating a position doesn't give infinite values.
    // ----------------------------------------------------------------------------------------------------------------
    private double getPsiMax()
	{
        double psiMax = 0.0f;
		double alpha = 1.0 - orbit.slr / orbit.Planet.SOI;

		if (orbit.ecc + alpha < 0.0f)
		{
			// The periapsis is higher than max radius - nothing to display
			psiMax = 0.0f;
		}
		else if (orbit.ecc - alpha < 0.0f)
		{
			// the conic is an ellipse and can be fully displayed
			psiMax = Math.PI;
		}
		else
		{
			// the conic is either a parabola or an hyperbola, or an ellipse that's too large to be fully displayed
			double sinPsiMax = (orbit.ecc + alpha) * (orbit.ecc - alpha);
			sinPsiMax = Math.Sqrt(sinPsiMax);

			double cosPsiMax = orbit.ecc * orbit.ecc - alpha;

			psiMax = Math.Atan2(sinPsiMax, cosPsiMax);

            // In case there's no SOI, reduce a little psiMax so that it doesn't crash the trajectory calculation (the last point would be at infinity otherwise)
            if (Double.IsInfinity(orbit.Planet.SOI))
            {
                if (psiMax > 2.0 * C_MIN_DELTA_PSI)
                {
                    psiMax -= C_MIN_DELTA_PSI;
                }
                else
                {
                    psiMax = psiMax / 2.0; // avoids returning a value that could be negative otherwise (if psiMax was REALLY small...)
                    LOG(LOG_LEVEL.WARNING, "psiMax is really small - reducing it to " + psiMax);
                }
            }
        }

        return psiMax;
    }


    // -------------------------------------------------------------------------------------------------------------------
    // getPsiAtTrueAnomaly
    // -------------------------------------------------------------------------------------------------------------------
    // Returns the psi angle corresponding to the specified true anomaly
    // -------------------------------------------------------------------------------------------------------------------
    private double getPsiAtTrueAnomaly(double trueAnomaly)
	{
		return Math.Atan2(Math.Sin(trueAnomaly),  orbit.ecc + Math.Cos(trueAnomaly));
	}


    // -------------------------------------------------------------------------------------------------------------------
    // calculateXYatPsi
    // -------------------------------------------------------------------------------------------------------------------
    // Calculates the x and y coordinates of a point given its tangential angle (psi)
    // The function works for any psi angle, and also applies the scaling factor.
    // The (x, y) coordinates assume no rotation, a rotation has to be applied to take into account argument of periapsis.
    // -------------------------------------------------------------------------------------------------------------------
    private (double x, double y) calculateXYatPsi(double psi)
    {
        // Check if it's a point in the upper part of an ellipse.
        bool reverseX = false;

        if (Math.Abs(psi) > Math.PI / 2.0)
        {
            // psi is not between -PI/2 and PI/2. We take the corresponding value in this interval and will take the symmetrical value for x (y is unchanged)
            psi = Math.Sign(psi) * Math.PI - psi;
            reverseX = true;
        }

        // calculate x and y
        double sinPsi = Math.Sin(psi);
		double cosPsi = Math.Cos(psi);

        double eSinPsi = orbit.ecc * sinPsi;
        double rho = Math.Sqrt(1.0 - eSinPsi * eSinPsi);
        double y = orbit.slr * sinPsi / rho;
        double x = orbit.periapsis - sinPsi * y / (cosPsi + rho);

        // Apply scale multiplier
        x *= scaleMultiplier;
        y *= scaleMultiplier;

        // Apply symmetry if point is on the other side of the minor axis (if it's an ellipse)
        if (reverseX)
        {
            x = -x - focusDistance;
        }

        return (x, y);
    }


    // ----------------------------------------------------------------------------------------------------------------
    // rotatePoint
    // ----------------------------------------------------------------------------------------------------------------
    // This method applies a rotation to the specified coordinates to take into account the trajectory orientation.
    // ----------------------------------------------------------------------------------------------------------------
    private (double x, double y) rotatePoint(double x, double y)
	{
		double new_x = x * cosPhi - y * sinPhi;
		double new_y = x * sinPhi + y * cosPhi;

		return (new_x, new_y);
    }


    // ----------------------------------------------------------------------------------------------------------------
    // CalculatePointAtPsi
    // ----------------------------------------------------------------------------------------------------------------
    // Directly calculates a point that corresponds to the specified psi angle.
    // ----------------------------------------------------------------------------------------------------------------
    private Vector3 CalculatePointAtPsi(double psi)
	{
        // calculate the point
        (double x, double y) = calculateXYatPsi(psi);

		// Apply the rotation to take into account the orientation of the conic
        (x, y) = rotatePoint(x, y);

		return new Vector3((float)x, (float)y);
    }


    // ----------------------------------------------------------------------------------------------------------------
    // getConicArea
    // ----------------------------------------------------------------------------------------------------------------
    // Returns the adimensionned conic area (this is the actual area divided by slr squared).
    // Open conics are bound to psiMax (which means in practice that the area that goes beyond the SOI is not counted).
    // ----------------------------------------------------------------------------------------------------------------
    private double getConicArea()
    {
        double area = 0.0;

        (double xMax, double _) = calculateXYatPsi(psiMax);
        xMax /= scaleMultiplier; // IMPORTANT: scaling must be cancelled in this particular case !

        double X = orbit.periapsis - xMax;
        double delta = (1.0 - orbit.ecc * orbit.ecc) * X / orbit.slr;

        if (Math.Abs(delta) < 0.0001)
        {
            // Quasi-parabolic case
            double V = delta * (2.0 - delta);
            area = X * (2.0 - delta) / orbit.slr;
            area *= Math.Sqrt(area);
            area *= (1680 + (504 + (270 + 175 * V) * V) * V) / 2520;
        }
        else
        {
            double U = 1 - delta;

            if (orbit.ecc < 1.0)
            {
                // elliptic case
                if (U < -1.0) U = -1.0; // to fix imprecisions just in case...
                area = (Math.Acos(U) - U * Math.Sqrt(1.0 - U * U)) / Math.Pow(1.0 - orbit.ecc * orbit.ecc, 1.5);
            }
            else
            {
                // hyperbolic case
                double sqrt_U2_minus_1 = Math.Sqrt(U * U - 1.0);
                area = (U * sqrt_U2_minus_1 - Math.Log(U + sqrt_U2_minus_1)) / Math.Pow(orbit.ecc * orbit.ecc - 1.0, 1.5);
            }
        }

        return area;
    }


    // ----------------------------------------------------------------------------------------------------------------
    // BuildBasePoints
    // ----------------------------------------------------------------------------------------------------------------
    // Build the original set of points that will be used for future operations.
    // Points are regularly dispatched so that they are separated by a constant deltaPsi angle.
    // Orientation is not taken into account at that stage - a rotation will have to be applied later.
    // Only necessary points are registered. It's necessary to make use of symmetries to deduce all points.
    // ----------------------------------------------------------------------------------------------------------------
    private List<OrbitBasePoint> BuildBasePoints()
	{
        List<OrbitBasePoint> listBasePts = new List<OrbitBasePoint>();

        int i_end;

        if(psiMax < Math.PI / 2.0)
        {
            i_end = (int)(psiMax / deltaPsi);
        }
        else
        {
            // This corresponds to psi = Pi/2. It's calculated that way because the formula above can give a lower value due to numerical precision errors.
            i_end = (int)nbMaxSegments / 4;
        }

		// build the initial list of points : all points are separated by a constant deltaPsi value to have a good basis
		// --------------------------------
		for(int i=0; i <= i_end; i++)
		{
			double psi = i * deltaPsi;
			(double x, double y) = calculateXYatPsi(psi);

            listBasePts.Add(new OrbitBasePoint(psi, x, y));
        }

        return listBasePts;
    }


    // ----------------------------------------------------------------------------------------------------------------
    // ApplyTesselation
    // ----------------------------------------------------------------------------------------------------------------
    // Applies a step of tesselation on a pre-built array of points : if the area lost on a particular segment is too
    // important, an intermediate point is calculated and inserted between those.
    // ----------------------------------------------------------------------------------------------------------------
    private uint ApplyTesselation(ref List<OrbitBasePoint> listBasePts)
	{
		uint nbPointsAdded = 0;

		// Evaluate the deltaPsi that will be used for this tesselation step (basically the previous deltaPsi divided by two)
		// ---------------------
		double currentDeltaPsi = 0.0;

		if(listBasePts.Count >= 2)
		{
			double psi1 = listBasePts[listBasePts.Count - 2].psi;
			double psi2 = listBasePts[listBasePts.Count - 1].psi;

			currentDeltaPsi = (psi2 - psi1) / 2.0;
        }
		else
		{
			// can happen with extremely energetic hyperbolas, for which only the periapsis had been added
			currentDeltaPsi = psiMax / 2.0;
        }

		// Don't add more points if we are already at the maximum resolution allowed 
		if(currentDeltaPsi < C_MIN_DELTA_PSI) return nbPointsAdded;

		// maximum admissible psi value (capped at PI/2 because above that value points are supposed to be obtained through symmetry)
		double psiEnd = Math.Min(Math.PI / 2, psiMax);

		// start at the final point
		int i_curPt = listBasePts.Count - 1;

		// Check if a point needs to be added after the last
		// -------------------------------------------------
		{
			double psi = listBasePts[i_curPt].psi;

			if( psi + currentDeltaPsi < psiEnd) // if there's room for another point...
			{
				double newPsi = psi + currentDeltaPsi;
				double eSinPsi = orbit.ecc * Math.Sin(newPsi);

				if(psiEnd - psi > toleranceThresholdForTessellation * (1.0 - eSinPsi * eSinPsi))
				{
                    (double x, double y) = calculateXYatPsi(newPsi);

                    listBasePts.Add(new OrbitBasePoint(newPsi, x, y));
					nbPointsAdded++;
                }
            }
        }

		i_curPt--;

		// Add as many intermediate points as needed
		// -----------------------------------------
		if(i_curPt >= 0)
        {
			double psi1 = listBasePts[i_curPt].psi;
            double psi2 = listBasePts[i_curPt + 1].psi;
            double newPsi = (psi1 + psi2) / 2;

            double eSinPsi = orbit.ecc * Math.Sin(newPsi);

			while((i_curPt >= 0) && (psi2 - psi1 > toleranceThresholdForTessellation * (1.0 - eSinPsi * eSinPsi)))
			{
                // Add an intermediate point
				(double x, double y) = calculateXYatPsi(newPsi);

				listBasePts.Insert(i_curPt + 1, new OrbitBasePoint(newPsi, x, y));
				nbPointsAdded++;

                // Prepare data for next iteration
                i_curPt--;

				if (i_curPt >= 0)
				{
					psi1 = listBasePts[i_curPt].psi;
					psi2 = listBasePts[i_curPt + 1].psi;
					newPsi = (psi1 + psi2) / 2;

					eSinPsi = orbit.ecc * Math.Sin(newPsi);
				}
            }
        }

        return nbPointsAdded;
	}


    // ----------------------------------------------------------------------------------------------------------------
    // BuildFullBaseArray
    // ----------------------------------------------------------------------------------------------------------------
    // Builds an array of points ideally dispatched along the trajectory. All points are associated to their tangential
    // angle (psi), and are sorted by ascending order.
    // It is a good idea to cache this array and then to build the final array by using these points and inserting a few
    // specific points.
    // ----------------------------------------------------------------------------------------------------------------
    private OrbitPoint[] BuildFullBaseArray()
	{
        // Build list of points on [0 - PI/2] interval range
		// -------------------------------------------------
        List<OrbitBasePoint> listBasePts = BuildBasePoints();


		// Apply tessellation
		// ------------------
		for(uint i=1; i<=nbMaxTessellationSteps; i++)
		{
			uint nbAddedPts = ApplyTesselation(ref listBasePts);

			if(nbAddedPts > 0)
			{
				LOG(LOG_LEVEL.INFO, "Conic (eccentricity = " + orbit.ecc + ") : Added " + nbAddedPts + " points at step " + i + " of tesselation");
			}
        }

		// Calculate total number of points
		// --------------------------------
		int totalNbPoints;

		if(closedConic)
		{
			totalNbPoints = 4 * (listBasePts.Count - 1);
        }
		else
		{
            totalNbPoints = listBasePts.Count;

			if(psiMax > Math.PI / 2.0)
			{
                // Points between pi_minus_psiMax and PI/2 have to be reversed wrt the minor axis.
				// Later, a point which tangential angle is psi will be mirrored in the point of tangential angle (PI - psi)
                double pi_minus_psiMax = Math.PI - psiMax;

                // Note: we start at listBasePts.Count - 2 because listBasePts.Count - 1 corresponds to psi = PI/2, which doesn't need to be mirrored
                for (int i = listBasePts.Count - 2; i >= 0; i--) 
                {
                    if (listBasePts[i].psi >= pi_minus_psiMax) totalNbPoints++;
                    else break;
                }
            }

			// number of points after applying symmetry wrt the major axis
			totalNbPoints *= 2;
			totalNbPoints--; // -1 because periapsis isn't reversed
        }

        // Build full array of points
        // --------------------------
        OrbitPoint[] listOrbitPoints = new OrbitPoint[totalNbPoints];

		// index of point corresponding to psi = 0.0 in array (periapsis)
		int i0 = totalNbPoints / 2;

		// will be used to reference elements in the base array
		int index = 0;
		int nbPointsAdded = 0;

		while (i0 + index < listOrbitPoints.Length)
		{
			double psi, x1, y1, x2, y2;

            if (index < listBasePts.Count)
			{
                psi = listBasePts[index].psi;
                x1 = listBasePts[index].x;
                y1 = listBasePts[index].y;
            }
			else
			{
                // After we've reached the end of listBasePts, we browse it in reverse to build symmetrical points wrt the minor axis
                int actualIndex = 2 * (listBasePts.Count - 1) - index;

                psi = Math.PI - listBasePts[actualIndex].psi;
                x1 = - listBasePts[actualIndex].x - focusDistance;
                y1 = listBasePts[actualIndex].y;
            }

            // point 2 corresponds to -psi (it's the symmetrical wrt the major axis)
            x2 = x1;
            y2 = -y1;

			// register point 1
			(x1, y1) = rotatePoint(x1, y1);
            listOrbitPoints[i0 + index].Set(psi, new Vector3((float)x1, (float)y1));
			nbPointsAdded++;

			// register point 2
            if (index != 0)
			{
                // This is the same point if index = 0 (it's the periapsis)...
                (x2, y2) = rotatePoint(x2, y2);
                listOrbitPoints[i0 - index].Set(-psi, new Vector3((float)x2, (float)y2));
				nbPointsAdded++;
            }
			else if(closedConic)
			{
				//... however if we have a full ellipse we have to register the apoapsis in first position
				x2 = -x2 - focusDistance;
				y2 = 0.0; // should already be the case

                (x2, y2) = rotatePoint(x2, y2);
                listOrbitPoints[0].Set(-Math.PI, new Vector3((float)x2, (float)y2));
				nbPointsAdded++;
            }

			index++;
        }

		if(nbPointsAdded != listOrbitPoints.Length)
		{
			LOG(LOG_LEVEL.ERROR, "GetFullBaseArray : " + nbPointsAdded + " points added - " + listOrbitPoints.Length + " expected !");
		}

        return listOrbitPoints;
    }


    // ----------------------------------------------------------------------------------------------------------------
    // GetLowerPsiIndex
    // ----------------------------------------------------------------------------------------------------------------
    // Finds in the orbit point list the index of the largest psi value that is strictly lower than the psi value
    // passed in parameter.
    // Returns -1 if all values are greater.
    // Note : it assumes that the orbit point list is sorted by order of ascending psi values.
    // ----------------------------------------------------------------------------------------------------------------
    private int GetLowerPsiIndex(double psi)
	{
        int iMin = 0;
        int iMax = orbitPointList.Length - 1;

        // manage extreme cases
        if (psi < orbitPointList[iMin].psi) return -1;
		if (psi > orbitPointList[iMax].psi) return iMax;

        // perform binary search
        while (iMax > iMin + 1)
		{
			int i = (iMax + iMin) / 2;

			if (orbitPointList[i].psi < psi) iMin = i;
			else iMax = i;
        }

		return iMin;
	}


    // ----------------------------------------------------------------------------------------------------------------
    // GetUpperPsiIndex
    // ----------------------------------------------------------------------------------------------------------------
    // Finds in the orbit point list the index of the lowest psi value that is strictly higher than the psi value
    // passed in parameter.
    // Returns orbitPointList.Length if all values are lower.
    // Note : it assumes that the orbit point list is sorted by order of ascending psi values.
    // ----------------------------------------------------------------------------------------------------------------
    private int GetUpperPsiIndex(double psi)
    {
        int iMin = 0;
        int iMax = orbitPointList.Length - 1;

        // manage extreme cases
        if (psi < orbitPointList[iMin].psi) return 0;
        if (psi > orbitPointList[iMax].psi) return iMax+1;

        // perform binary search
        while (iMax > iMin + 1)
        {
            int i = (iMax + iMin) / 2;

            if (orbitPointList[i].psi > psi) iMax = i;
            else iMin = i;
        }

        return iMax;
    }


    // ----------------------------------------------------------------------------------------------------------------
    // GetPoints
    // ----------------------------------------------------------------------------------------------------------------
    // Builds an array of Vector3D that represents the full trajectory, starting and ending at the specified true
    // anomalies.
    // The methods builds the points array by adding in this order :
    // - the point that corresponds to fromTrueAnomaly
    // - all necessary intermediate points (from the pre-calculated array)
    // - the point that corresponds to toTrueAnomaly
    // ----------------------------------------------------------------------------------------------------------------
    public Vector3[] GetPoints(double fromTrueAnomaly, double toTrueAnomaly)
	{
		double psiFrom = getPsiAtTrueAnomaly(fromTrueAnomaly);
        double psiTo = getPsiAtTrueAnomaly(toTrueAnomaly);

		bool fullTurn = closedConic && (psiFrom == psiTo);


		// Calculate start/end indexes for base points that have to be included
		// --------------------------------------------------------------------
        int i_from = 0; // index of first point that comes after psiFrom
        int i_to = 0; // index of the last point that comes before psiTo
        int nbPoints = 0;

        if (orbit.direction > 0)
		{
            // get first value greater than psiFrom in points list
			i_from = GetUpperPsiIndex(psiFrom);

            // get last value lower than psiTo in points list
            i_to = GetLowerPsiIndex(psiTo);

			// take into account that it's a circular array if the conic is closed
			if(closedConic)
			{
				if (i_from == orbitPointList.Length) i_from = 0;
				if (i_to == -1) i_to = orbitPointList.Length - 1;
            }
        }
        else // orbit.direction < 0
        {
			// get first value lower than psiFrom in points list
            i_from = GetLowerPsiIndex(psiFrom);

            // get last value greater than psiTo in points list
            i_to = GetUpperPsiIndex(psiTo);

            // take into account that it's a circular array if the conic is closed
            if (closedConic)
            {
                if (i_from == -1) i_from = orbitPointList.Length - 1;
				if (i_to == orbitPointList.Length) i_to = 0;
            }
        }

		// Calculate the total number of points
		// ------------------------------------

		// calculate number of base points included
		if(orbit.direction > 0)
		{
			nbPoints = i_to - i_from + 1;

			// fix nbPoints if needed (the base point list is used as a circular array for a closed conic, so indexes can be in reverse order)
            if((orbit.ecc < 1.0) && ((nbPoints < 0) || ((nbPoints == 0) && (psiFrom >= psiTo))))
			{
				nbPoints += orbitPointList.Length; // to take into account the fact that it's a circular array for ellipses

                // NOTE : if nbPoints == 0 and psiFrom < psiTo, we just have to draw a small chunk of
                //        trajectory on which no base point is included, so it's ok to have 0 base points.
            }
        }
		else
		{
            nbPoints = i_from - i_to + 1;

            // fix nbPoints if needed (the base point list is used as a circular array for a closed conic, so indexes can be in reverse order)
            if (orbit.ecc < 1.0 && ((nbPoints < 0) || ((nbPoints == 0) && (psiFrom <= psiTo))))
			{
				nbPoints += orbitPointList.Length;

                // NOTE : if nbPoints == 0 and psiFrom > psiTo, we just have to draw a small chunk of
                // trajectory on which no base point is included, so it's ok to have 0 base points.
            }
        }

        // sanity check
        if (nbPoints < 0)
        {
            LOG(LOG_LEVEL.WARNING, "Negative number of points counted (" + nbPoints + ") - psiFrom = " + psiFrom + "; psiTo = " + psiTo + "; dir = " + orbit.direction + "; i_from = " + i_from + "; i_to = " + i_to + "; closedConic = " + closedConic + "; psiMax = " + psiMax);
            nbPoints = 0;
        }

        // add 2 points to take into account "from" and "to"
        nbPoints += 2;


        // Build the array
        // ---------------
        Vector3[] trajectory = new Vector3[nbPoints];
        int trajectoryIndex = 0;

		// add "from"
		trajectory[trajectoryIndex] = CalculatePointAtPsi(psiFrom);
		trajectoryIndex++;

        // add base points
        int i_point = Math.Min(orbitPointList.Length - 1, Math.Max(0, i_from)); // make sure point index is within the base array bounds... 
                                                                                // it can be incorrect in case of incorrect input parameters but at least SFS won't crash because of this!

        while (trajectoryIndex < nbPoints - 1)
		{
			trajectory[trajectoryIndex] = orbitPointList[i_point].point;

			trajectoryIndex++;

			if(orbit.direction > 0)
			{
				i_point++;
				if (i_point == orbitPointList.Length) i_point = 0; // loop if needed
            }
			else
			{
				i_point--;
                if (i_point == -1) i_point = orbitPointList.Length - 1; // loop if needed
            }
        }

        // add "to"
        trajectory[trajectoryIndex] = CalculatePointAtPsi(psiTo);
		
        return trajectory;
    }



    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.ORBIT_DRAWER, level, message);
    }
}
