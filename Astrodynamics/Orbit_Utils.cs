using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

using SFS.World;
using SFS.WorldBase;
using SFS.World.PlanetModules;

static class Orbit_Utils
{
	// The minimal value the slr shall have to ensure the orbit is not rectilinear (or practically)
	public const double C_MINIMAL_ORBIT_SLR = 0.001;
	
	public static double NormalizeAngle(double angle)
	{
		while (angle >  Math.PI) angle -= Math.PI * 2.0;
		while (angle < -Math.PI) angle += Math.PI * 2.0;

		return angle;
	}

	// Deprecated: use orbit.getSpecificEnergy() instead
	public static double GetSpecificEnergy(Orbit orbit)
	{
		double specificEnergy = orbit.Planet.mass * (orbit.ecc * orbit.ecc - 1.0) / (2.0 * orbit.slr);
		return specificEnergy;
	}

	// Deprecated: use orbit.getAngularMomentum() instead
	public static double GetAngularMomentum(Orbit orbit)
	{
		double angularMomentum = Math.Sqrt(orbit.Planet.mass * orbit.slr) * orbit.direction;
		return angularMomentum;
	}

	public static Orbit CreateOrbit(double slr, double ecc, double argOfPeriapsis, int direction, Planet planet, PathType pathType, Planet nextPlanet)
    {
		Orbit_OrbitConstructor2_Patch.useHackedConstructor = true;

		if (slr < C_MINIMAL_ORBIT_SLR)
		{
			AnaisLogger.Log(LOG_CATEGORY.ORBIT, LOG_LEVEL.INFO, "Custom orbit creation aborted - it would be rectilinear (slr = " + slr + ")");
			return null; // The orbit would be rectilinear
		}

		return new Orbit(slr, ecc, argOfPeriapsis, direction, planet, pathType, nextPlanet);
	}

	public static Orbit Clone(this Orbit orbit, Planet parentPlanet)
    {
		Orbit_OrbitConstructor2_Patch.orbitToCopy = orbit;

		// Only the 3 last parameters are important
		LOG(LOG_LEVEL.DEBUG, "Clone called");
		return new Orbit(0.0, 0.0, 0.0, 0, parentPlanet, orbit.pathType, null);
    }

	public static Planet Clone(this Planet planet, UnityEngine.GameObject gameObject)
    {
		Planet newPlanet = gameObject.AddComponent<Planet>();

		// Copy the useful fields from original planet
		//newPlanet.codeName = planet.codeName;
		newPlanet.mass = planet.mass;
		newPlanet.SOI = planet.SOI;
		newPlanet.maxTerrainHeight = planet.maxTerrainHeight;

		// To be assigned manually by the caller if needed
		newPlanet.orbit = null;
		newPlanet.parentBody = null;

		// Initialize all the rest to null/0 by safety
		newPlanet.trajectory = null;
		newPlanet.mapHolder = null;
		newPlanet.mapPlanet = null;
		newPlanet.landmarks = null;
		newPlanet.satellites = null;
		newPlanet.planetTexture = null;
		newPlanet.terrainMaterial = null;
		newPlanet.atmosphereMaterial = null;
		newPlanet.orbitalDepth = 0;
		newPlanet.satelliteIndex = 0;

		// Planet data
		newPlanet.data = new PlanetData();

		// Copy the useful fields from original planet data (basics)
		newPlanet.data.basics = new BasicModule();
		newPlanet.data.basics.radius = planet.data.basics.radius;
		newPlanet.data.basics.gravity = planet.data.basics.gravity;
		newPlanet.data.basics.timewarpHeight = planet.data.basics.timewarpHeight;

		// Copy the useful fields from original planet data (atmospherePhysics)
		newPlanet.data.hasAtmospherePhysics = planet.data.hasAtmospherePhysics;
		if (planet.data.hasAtmospherePhysics)
		{
			newPlanet.data.atmospherePhysics = new Atmosphere_Physics();
			newPlanet.data.atmospherePhysics.height = planet.data.atmospherePhysics.height;
			newPlanet.data.atmospherePhysics.density = planet.data.atmospherePhysics.density;
			newPlanet.data.atmospherePhysics.curve = planet.data.atmospherePhysics.curve;
			newPlanet.data.atmospherePhysics.parachuteMultiplier = planet.data.atmospherePhysics.parachuteMultiplier;
			newPlanet.data.atmospherePhysics.upperAtmosphere = planet.data.atmospherePhysics.upperAtmosphere;
		}
        else
        {
			newPlanet.data.atmospherePhysics = null;
		}

		// Initialize the rest to invalid data
		newPlanet.data.hasAtmosphereVisuals = false;
		newPlanet.data.atmosphereVisuals = null;
		newPlanet.data.hasTerrain = false;
		newPlanet.data.terrain = null;
		newPlanet.data.hasPostProcessing = false;
		newPlanet.data.postProcessing = null;
		newPlanet.data.hasOrbit = false;
		newPlanet.data.orbit = null;
		newPlanet.data.landmarks = null;

		return newPlanet;
    }

	public static bool GetIntersectionAngles(Orbit orbit1, Orbit orbit2, out double angleA, out double angleB)
    {
		double deltaPhi = NormalizeAngle(orbit2.arg - orbit1.arg);
		double cosDeltaPhi = Math.Cos(deltaPhi);
		double A = Math.Sqrt(Math.Pow(orbit1.ecc / orbit1.slr, 2.0) + Math.Pow(orbit2.ecc / orbit2.slr, 2.0) - 2.0 * orbit1.ecc * orbit2.ecc * cosDeltaPhi / (orbit1.slr * orbit2.slr));
		double B = 1.0 / orbit1.slr - 1.0 / orbit2.slr;
		
		// calculate the median argument
		double medianArg = Math.Acos((orbit2.ecc * cosDeltaPhi / orbit2.slr - orbit1.ecc / orbit1.slr) / A);
		if (Math.Sin(deltaPhi) < 0.0)
		{
			medianArg = -medianArg;
		}

		medianArg += orbit1.arg;

		// Solve the equation A * cos(angle - medianArg) = B
		if (A > Math.Abs(B))
		{
			// The equation has 2 solutions -> orbits intersect in two distinct points
			double deltaArg = Math.Acos(B / A);
			angleA = NormalizeAngle(medianArg - deltaArg);
			angleB = NormalizeAngle(medianArg + deltaArg);
			return true;
		}
		else
        {
			// The equation has no solution -> return angle where orbits are closer
			if(B < 0.0)
            {
				medianArg += Math.PI;
				NormalizeAngle(medianArg);
			}
			angleA = medianArg;
			angleB = medianArg;
			return false;
		}
	}

	public static void CalculateHohmannTransfer(Orbit fromOrbit, Orbit targetOrbit, Location location_from, out Orbit hohmannTransfer, out double arg_encounter)
	{
		// Initialize output data
		hohmannTransfer = null;
		arg_encounter = 0.0;

		// Check data consistency
		if((fromOrbit == null) || (targetOrbit == null))
        {
			LOG(LOG_LEVEL.ERROR, "One of the passed orbits is null");
			return;
        }
		else if(fromOrbit.Planet != targetOrbit.Planet)
        {
			LOG(LOG_LEVEL.ERROR, "Both orbits don't have the same parent body");
			return;
		}
		else if( (location_from == null) || (location_from.planet != fromOrbit.Planet))
        {
			//FileLog.Log("ERROR: CalculateHohmannTransfer: origin location is wrong (object is null, or the parent body doesn't correspond)");
			LOG(LOG_LEVEL.ERROR, "origin location is wrong (object is null, or the parent body doesn't correspond)");
			return;
		}

		// local variables to alleviate notations
		Planet thePlanet = fromOrbit.Planet;

		double slr1 = fromOrbit.slr;
		double ecc1 = fromOrbit.ecc;
		double phi1 = fromOrbit.arg;

		double slr2 = targetOrbit.slr;
		double ecc2 = targetOrbit.ecc;
		double phi2 = targetOrbit.arg;

		double arg_start = location_from.position.AngleRadians;
		double R1 = fromOrbit.GetRadiusAtAngle(arg_start);
		double R2 = targetOrbit.GetRadiusAtAngle(arg_start);

		// Calculate Xi = crossing indicator (if negative, orbits don't cross themselves, if positive, orbits cross)
		double Xi = (1.0 - ecc1 * ecc2 * Math.Cos(phi2 - phi1)) / (slr1 * slr2) - 0.5 * ((1 - ecc1 * ecc1) / (slr1 * slr1) + (1 - ecc2 * ecc2) / (slr2 * slr2));

		if (!(Xi < 0.0) && (Math.Abs(R1 - R2) / (R1 + R2) < 0.0001))
		{
			// The orbits cross themselves, and we are too close from the crossing point
			LOG(LOG_LEVEL.INFO, "Orbits cross and we are too close from a crossing point. Calculation is too risky, give up.");
			return;
		}

		// Calculate transfer orbit
		double slr = 1.0 / (1.0 / slr1 + Xi * R1 * R2 / (R2 - R1));

		if (slr < C_MINIMAL_ORBIT_SLR)
		{
			// "impossible hyperbola" (if negative) or rectilinear orbit (if positive but practically null)
			LOG(LOG_LEVEL.INFO, "Hohmann transfer doesn't exist or is rectilinear (slr = " + slr + ")");
			return;
		}

		double ecc = Math.Pow(1.0 / slr - 1.0 / R1, 2.0) + Math.Pow(ecc1 / slr1, 2.0) - Math.Pow(1.0 / slr1 - 1.0 / R1, 2.0);
		if (ecc < 0.0) ecc = 0.0; // To fix for some possible numerical imprecision
		ecc = slr * Math.Sqrt(ecc);

		double trueAnomaly1 = arg_start - phi1;

		while (trueAnomaly1 > Math.PI) { trueAnomaly1 -= 2.0 * Math.PI; }
		while (trueAnomaly1 < -Math.PI) { trueAnomaly1 += 2.0 * Math.PI; }

		double argOfPeriapsis = 0.0;
		double theCosinus = (slr / R1 - 1.0) / ecc;
		if (theCosinus > 1.0) { theCosinus = 1.0; }
		if (theCosinus < -1.0) { theCosinus = -1.0; }
		if (1.0 + ecc > 1.0) { argOfPeriapsis = arg_start - Math.Sign(trueAnomaly1) * Math.Acos(theCosinus); }

		double specificEnergy = thePlanet.mass * (ecc * ecc - 1.0) / slr / 2.0;

		hohmannTransfer = CreateOrbit(slr, ecc, argOfPeriapsis, fromOrbit.direction, thePlanet, PathType.Encounter, null);

		if(hohmannTransfer == null)
        {
			LOG(LOG_LEVEL.ERROR, "Something went wrong when creating the orbit!");
			return;
		}

		// Calculate encounter point
		double R_encounter = 2.0 * (1.0 / slr - 1.0 / slr2) / ((1.0 - ecc * ecc) / (slr * slr) - (1.0 - ecc2 * ecc2) / (slr2 * slr2));

		double deltaPhi = argOfPeriapsis - phi2;

		while (deltaPhi > Math.PI) { deltaPhi -= 2.0 * Math.PI; }
		while (deltaPhi < -Math.PI) { deltaPhi += 2.0 * Math.PI; }

		int sign = 1;

		if (deltaPhi < 0.0) sign = -sign;
		if (1.0 / slr2 - 1.0 / slr < 0.0) sign = -sign;

		// Calculate argument of encounter
		if (1.0 + ecc > 1.0)
		{
			double the_cosinus = (slr / R_encounter - 1.0) / ecc;

			if (the_cosinus > 1.0) { the_cosinus = 1.0; }
			if (the_cosinus < -1.0) { the_cosinus = -1.0; }
			arg_encounter = argOfPeriapsis + sign * Math.Acos(the_cosinus);
		}
		else if (1.0 + ecc2 > 1.0)
		{
			double the_cosinus = (slr2 / R_encounter - 1.0) / ecc2;

			if (the_cosinus > 1.0) { the_cosinus = 1.0; }
			if (the_cosinus < -1.0) { the_cosinus = -1.0; }
			arg_encounter = phi2 + sign * Math.Acos(the_cosinus);
		}
		else
        {
			arg_encounter = arg_start;
		}

		// Set periapsis passage time
		double timeFromPeri = 0.0;
		new KeplerSolver(thePlanet.mass, hohmannTransfer.periapsis, specificEnergy).GetTimeAtPosition(R1, arg_start - hohmannTransfer.arg, ref timeFromPeri);
		hohmannTransfer.periapsisPassageTime = location_from.time - timeFromPeri * (double)hohmannTransfer.direction - 10.0 * hohmannTransfer.period;

		// Set start and end time
		hohmannTransfer.orbitStartTime = location_from.time;
		hohmannTransfer.orbitEndTime = hohmannTransfer.GetNextAnglePassTime(location_from.time, arg_encounter);

		if(!(hohmannTransfer.orbitEndTime > location_from.time))
        {
			// Happens when the trajectory is hyperbolic and the encounter point is passed. In this case the ship will never reach that position
			// --> return an invalid transfer because it's of no use in practice.
			hohmannTransfer = null;
		}
	}

	public static List<Orbit> GetListOrbits(Trajectory a)
	{
		List<Orbit> orbitList = new List<Orbit>();

		for (int i = 0; i < a.paths.Count; i++)
		{
			if (a.paths[i] is Orbit orbit)
			{
				orbitList.Add(orbit);
			}
			else
			{
				break;
			}
		}

		return orbitList;
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.HOHMANN_TRANSFER, level, message);
	}
}
