using HarmonyLib;
using SFS.World;
using SFS.WorldBase;
using System;
using UnityEngine;

using System.Diagnostics;
using System.Runtime.CompilerServices;


// Class Orbit_AdditionalData
// --------------------------
// A class that allows to store additional data and associate them to an orbit instance
public class Orbit_AdditionalData
{
	public double specificEnergy;
	public double angularMomentum;
	public OrbitDrawer orbitDrawer;

    public Orbit_AdditionalData(double specificEnergy, double angularMomentum, Orbit orbit)
    {
        this.specificEnergy = specificEnergy;
        this.angularMomentum = angularMomentum;
        this.orbitDrawer = new OrbitDrawer(orbit);
    }
}


// Class Orbit_AdditionalData_Association
// --------------------------------------
// This class maintains a table of association between all orbit instances and their associated data.
// It also defines some methods to retrieve said data as if they were actual members of the orbit class
// through the use of extension methods.
static class Orbit_AdditionalData_Association
{
	private static ConditionalWeakTable<Orbit, Orbit_AdditionalData> orbitTable = new ConditionalWeakTable<Orbit, Orbit_AdditionalData>();

	// Backup function that allows to recreate a set of associated data if an association appeared to be lost
	// This is only used as a backup, it should never happen in practice.
	private static Orbit_AdditionalData createAdditionalData(Orbit orbit)
    {
		LOG(LOG_LEVEL.WARNING, "createAdditionalData: orbit not found in ConditionalWeakTable; added automatically on the go!");
		return new Orbit_AdditionalData(Orbit_Utils.GetSpecificEnergy(orbit), Orbit_Utils.GetAngularMomentum(orbit), orbit);
	}

	// Function to add an entry in the table; only meant to be called by the Orbit constructor!
	public static void addAdditionalData(this Orbit orbit, double specificEnergy, double angularMomentum)
    {
		orbitTable.Add(orbit, new Orbit_AdditionalData(specificEnergy, angularMomentum, orbit));
	}

	// Accessors to the orbit additional data
	// --------------------------------------
	public static double getSpecificEnergy(this Orbit orbit)
    {
		return orbitTable.GetValue(orbit, createAdditionalData).specificEnergy;
	}

	public static double getAngularMomentum(this Orbit orbit)
	{
		return orbitTable.GetValue(orbit, createAdditionalData).angularMomentum;
	}

	public static void getSpecificEnergyAndAngularMomentum(this Orbit orbit, out double specificEnergy, out double angularMomentum)
	{
		Orbit_AdditionalData orbitData = orbitTable.GetValue(orbit, createAdditionalData);
		
		specificEnergy = orbitData.specificEnergy;
		angularMomentum = orbitData.angularMomentum;
	}

	public static OrbitDrawer getOrbitDrawer(this Orbit orbit)
	{
		return orbitTable.GetValue(orbit, createAdditionalData).orbitDrawer;
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
	}
}

[HarmonyPatch(typeof(Orbit), "TryCreateOrbit")]
public class Orbit_TryCreateOrbit_Patch
{
	[HarmonyPrefix]
	public static bool TryCreateOrbit_Prefix()
	{
		// skip the original method
		return false;
	}

	[HarmonyPostfix]
	public static Orbit TryCreateOrbit_Postfix(Orbit __result, Location location, bool calculateTimeParameters, bool calculateEncounters, ref bool success)
	{
		success = false;
		
		double A = Double3.Cross(location.position, location.velocity).z;

		if (A*A < Orbit_Utils.C_MINIMAL_ORBIT_SLR * location.planet.mass)
        {
			LOG(LOG_LEVEL.INFO, "Orbit creation failed - it would be rectilinear (A = " + A);
			return null;
		}

		if (location.position.magnitude > location.planet.SOI * 1.01)
		{
			LOG(LOG_LEVEL.WARNING, "Orbit creation failed - Out of SOI (" + location.planet.name + "); distance = " + location.position.magnitude + "; SOI radius = " + location.planet.SOI);
			return null;
		}
		Orbit orbit = new Orbit(location, calculateTimeParameters, calculateEncounters);
		
		if (double.IsNaN(orbit.periapsisPassageTime))
		{
			LOG(LOG_LEVEL.ERROR, "Orbit creation failed - periapsisPassageTime is NaN");
			return null;
		}
		
		success = true;
		return orbit;
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
	}
}

[HarmonyPatch(typeof(Orbit), MethodType.Constructor, new Type[] { typeof(Location), typeof(bool), typeof(bool) })]
public class Orbit_OrbitPatch
{
	[HarmonyPrefix]
	public static bool Orbit_prefix(Location location)
    {
		return false;
    }

	[HarmonyPostfix]
	public static void Orbit_postfix(Orbit __instance, Location location, bool calculateTimeParameters, bool calculateEncounters)
    {
		__instance.orbitStartTime = location.time;

		Double3 angularMomentum_vector = Double3.Cross(location.position, location.velocity);
		Double2 eccentricity_vector = (Double2)(Double3.Cross((Double3)location.velocity, angularMomentum_vector) / location.planet.mass) - location.position.normalized;

		double specificEnergy = Math.Pow(location.velocity.magnitude, 2.0) / 2.0 - location.planet.mass / location.Radius;
		double angularMomentum = angularMomentum_vector.z;

		__instance.addAdditionalData(specificEnergy, angularMomentum);

		__instance.ecc = eccentricity_vector.magnitude;
		__instance.slr = angularMomentum * angularMomentum / location.planet.mass;
		__instance.periapsis = __instance.slr / (1.0 + __instance.ecc);
		__instance.arg = eccentricity_vector.AngleRadians;
		__instance.direction = Math.Sign(angularMomentum);

		if (specificEnergy < 0.0)
		{
			// General formulas that also work for rectilinear orbits
			__instance.apoapsis = -location.planet.mass * (1.0 + __instance.ecc) / (2.0 * specificEnergy);
			__instance.sma = location.planet.mass / (-2.0 * specificEnergy);
		}
		else
		{
			// parabolic/hyperbolic orbit --> Infinite
			__instance.apoapsis = double.PositiveInfinity;
			__instance.sma = double.PositiveInfinity;
		}

		double trueAnomaly_Out = Orbit_Utils.NormalizeAngle(location.position.AngleRadians - __instance.arg);
		Traverse.Create(__instance).Field("trueAnomaly_Out").SetValue(trueAnomaly_Out);
		Traverse.Create(__instance).Field("location_Out").SetValue(location);

		if (calculateTimeParameters)
		{
			bool isEscaping = __instance.apoapsis >= location.planet.SOI;

			// calculate period
			__instance.period = (isEscaping ? 0.0 : Kepler.GetPeriod(__instance.sma, location.planet.mass));

			// calculate periapsis passage time
			double timeFromPeri = 0.0;
			new KeplerSolver(location.planet.mass, __instance.periapsis, specificEnergy).GetTimeAtPosition(location.Radius, trueAnomaly_Out, ref timeFromPeri);
			if(double.IsNaN(timeFromPeri))
            {
				LOG(LOG_LEVEL.ERROR, "timeFromPeri is NaN - planet mass = " + location.planet.mass + " - peri = " + __instance.periapsis + " - specEnerg = " + specificEnergy + " - radius = " + location.Radius + " - trueAnomaly = " + trueAnomaly_Out);
			}
			if(double.IsNaN(location.time))
            {
				LOG(LOG_LEVEL.ERROR, "location.time is NaN - planet = " + location.planet.codeName + "; x, y = " + location.position.x + "; " + location.position.y + "; vx, vy = " + location.velocity.x + "; " + location.velocity.y);
			}
			if(double.IsNaN(__instance.period))
            {
				LOG(LOG_LEVEL.ERROR, "period is NaN - sma = " + __instance.sma + "; planet mass = " + location.planet.mass);
			}
			
			__instance.periapsisPassageTime = location.time - timeFromPeri * (double)(__instance.direction) /*- 10.0 * __instance.period*/;


			if (isEscaping)
			{
				// calculate orbit end time
				if (location.planet.SOI == double.PositiveInfinity)
				{
					__instance.orbitEndTime = double.PositiveInfinity;
				}
				else
				{
					double loc_trueAnomaly = Math.Acos((__instance.slr / location.planet.SOI - 1.0) / __instance.ecc);
					new KeplerSolver(location.planet.mass, __instance.periapsis, specificEnergy).GetTimeAtPosition(location.planet.SOI, loc_trueAnomaly, ref timeFromPeri);
					__instance.orbitEndTime = __instance.periapsisPassageTime + timeFromPeri;
				}

				// set path type, next planet
				__instance.pathType = PathType.Escape;
				Traverse.Create(__instance).Field("nextPlanet").SetValue(location.planet.parentBody);
			}
			else
			{
				// set orbit end time, path type, next planet
				__instance.pathType = PathType.Eternal;
				__instance.orbitEndTime = double.PositiveInfinity;
				Traverse.Create(__instance).Field("nextPlanet").SetValue(null);
			}

			// Calculate encounters
			if (calculateEncounters)
			{
				double startWindow = location.time;
				double endWindow = isEscaping ? double.PositiveInfinity : (location.time + __instance.period * 0.99);
				Traverse.Create(__instance).Method("FindEncounters", startWindow, endWindow).GetValue(startWindow, endWindow);
			}
		}
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
	}
}


[HarmonyPatch(typeof(Orbit), MethodType.Constructor, new Type[] { typeof(double), typeof(double), typeof(double), typeof(int), typeof(Planet), typeof(PathType), typeof(Planet) })]
public class Orbit_OrbitConstructor2_Patch
{
	// This is a hack that tells the constructor to interpret the sma parameter as slr instead.
	// In exchange, it allows to create parabolic or hyperbolic orbits, which the original constructor isn't made for
	public static bool useHackedConstructor = false;

	// Another hack: use to act as a copy constructor
	public static Orbit orbitToCopy = null;
	
	[HarmonyPrefix]
	public static bool Orbit_prefix()
	{
		return false;
	}

	[HarmonyPostfix]
	public static void Orbit_postfix(Orbit __instance, double sma, double ecc, double arg, int direction, Planet planet, PathType pathType, Planet nextPlanet)
    {
		LOG(LOG_LEVEL.DEBUG, "Setting location_out");
		Traverse.Create(__instance).Field("location_Out").SetValue(new Location(planet, Double2.zero, Double2.zero));

		__instance.pathType = pathType;
		LOG(LOG_LEVEL.DEBUG, "Setting nextPlanet");
		Traverse.Create(__instance).Field("nextPlanet").SetValue(nextPlanet);

		if (pathType == PathType.Eternal)
		{
			__instance.orbitEndTime = double.PositiveInfinity;
		}

		if(orbitToCopy != null)
        {
			// HACKED CONSTRUCTOR: act as a copy constructor - I hope Stef won't see that mess!
			// ------------------
			LOG(LOG_LEVEL.DEBUG, "Copy constructor called");
			__instance.slr       = orbitToCopy.slr;
			__instance.ecc       = orbitToCopy.ecc;
			__instance.arg       = orbitToCopy.arg;
			__instance.direction = orbitToCopy.direction;
			__instance.periapsis = orbitToCopy.periapsis;
			__instance.apoapsis  = orbitToCopy.apoapsis;
			__instance.sma       = orbitToCopy.sma;
			__instance.period    = orbitToCopy.period;
			__instance.periapsisPassageTime = orbitToCopy.periapsisPassageTime;
			__instance.orbitStartTime = orbitToCopy.orbitStartTime;
			__instance.orbitEndTime = orbitToCopy.orbitEndTime;

			orbitToCopy.getSpecificEnergyAndAngularMomentum(out double specificEnergy, out double angularMomentum);
			__instance.addAdditionalData(specificEnergy, angularMomentum);

			// Reset that variable
			orbitToCopy = null;
		}
		else if (useHackedConstructor)
        {
			// HACKED CONSTRUCTOR: slr has been passed instead of sma (so we can handle the case of parabolas)
			// ------------------
			useHackedConstructor = false; // reset that flag

			__instance.slr = sma;
			__instance.ecc = ecc;
			__instance.arg = arg;
			__instance.direction = direction;
			__instance.periapsis = __instance.slr / (1.0 + ecc);

			if (ecc < 1.0)
            {
				// elliptic orbit: calculate apoapsis, sma, period
				__instance.apoapsis = __instance.slr / (1.0 - ecc);
				__instance.sma = __instance.slr / (1.0 - ecc * ecc);
				__instance.period = 2.0 * Math.PI * Math.Sqrt(__instance.sma * __instance.sma * __instance.sma / planet.mass);
			}
			else
            {
				// parabolic/hyperbolic orbit: data not defined
				__instance.apoapsis = Double.PositiveInfinity;
				__instance.sma = Double.PositiveInfinity;
				__instance.period = 0.0;
			}

			double specificEnergy = planet.mass * (ecc * ecc - 1.0) / (2.0 * __instance.slr);
			double angularMomentum = Math.Sqrt(planet.mass * __instance.slr) * direction;
			__instance.addAdditionalData(specificEnergy, angularMomentum);
		}
		else
        {
			// NORMAL CONSTRUCTOR
			// ------------------
			__instance.sma = sma;
			__instance.ecc = ecc;
			__instance.arg = arg;
			__instance.direction = direction;
			__instance.slr = sma * (1 - ecc * ecc);
			__instance.periapsis = sma * (1 - ecc);
			__instance.period = Kepler.GetPeriod(sma, planet.mass);

			if (ecc < 1)
            {
				__instance.apoapsis = sma * (1 + ecc);
			}
			else
            {
				__instance.apoapsis = Double.PositiveInfinity;
			}

			double specificEnergy = planet.mass * (ecc * ecc - 1.0) / (2.0 * __instance.slr);
			double angularMomentum = Math.Sqrt(planet.mass * __instance.slr) * direction;
			__instance.addAdditionalData(specificEnergy, angularMomentum);
		}
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
	}
}



[HarmonyPatch(typeof(Orbit), "GetVelocityAtTrueAnomaly")]
public class Orbit_GetVelocityAtTrueAnomaly_Patch
{
	[HarmonyPrefix]
	public static bool GetVelocityAtTrueAnomaly_Prefix(double trueAnomaly)
	{
		// skip the original method
		return false;
	}

	[HarmonyPostfix]
	public static void GetVelocityAtTrueAnomaly_Postfix(ref Double2 __result, Orbit __instance, double trueAnomaly)
    {
		__instance.getSpecificEnergyAndAngularMomentum(out double specificEnergy, out double angularMomentum);

		double normalizedTrueAnomaly = Orbit_Utils.NormalizeAngle(trueAnomaly);
		double radiusAtTrueAnomaly = Kepler.GetRadiusAtTrueAnomaly(__instance.slr, __instance.ecc, normalizedTrueAnomaly);
		
		double V_theta = angularMomentum / radiusAtTrueAnomaly;
		double V2 = 2.0 * (specificEnergy + __instance.Planet.mass / radiusAtTrueAnomaly);
		double V_r2 = Math.Max(V2 - V_theta * V_theta, 0.0);
		double V_r = Math.Sqrt(V_r2);
		if (Math.Sign(normalizedTrueAnomaly) != __instance.direction) { V_r = -V_r; }

		Double2 velocity = new Double2(V_r, V_theta);
		velocity = velocity.Rotate(normalizedTrueAnomaly + __instance.arg);
		__result = velocity;
	}
}

[HarmonyPatch(typeof(Orbit), "GetLastTrueAnomalyPassTime")]
public class Orbit_GetLastTrueAnomalyPassTime_Patch
{
	[HarmonyPrefix]
	public static bool GetLastTrueAnomalyPassTime_Prefix(double time, double trueAnomaly)
	{
		// skip the original method
		return false;
	}

	[HarmonyPostfix]
	public static double GetLastTrueAnomalyPassTime_Postfix(double __result, Orbit __instance, double time, double trueAnomaly, Location ___location_Out)
    {
		double timeFromPeri = 0.0;
		double radius = __instance.slr / (1.0 + __instance.ecc * Math.Cos(trueAnomaly));
		new KeplerSolver(___location_Out.planet.mass, __instance.periapsis, __instance.getSpecificEnergy()).GetTimeAtPosition(radius, trueAnomaly, ref timeFromPeri);

		double trueAnomalyPassTime = __instance.periapsisPassageTime + timeFromPeri * (double)(__instance.direction);

		if(__instance.period > 0.0) // period is 0 for hyperbolic or escaping orbits (an object can escape without being on an hyperbolic orbit)
        {
			// the number of periods passed since the memorized periapsis passage time, rounded down
			double nbPeriodsPassed = Math.Floor((time - trueAnomalyPassTime) / __instance.period);
			return trueAnomalyPassTime + nbPeriodsPassed * __instance.period;
		}
		else
        {
			return trueAnomalyPassTime;
		}
	}
}


[HarmonyPatch(typeof(Orbit), "UpdateLocation")]
public class Orbit_UpdateLocation_Patch
{
    [HarmonyPrefix]
    public static bool UpdateLocation_Prefix(double newTime)
    {
        // skip the original method
        return false;
    }

    [HarmonyPostfix]
    public static void UpdateLocation_Postfix(Orbit __instance, double newTime, ref Location ___location_Out, ref double ___trueAnomaly_Out)
    {

		if(double.IsNaN(newTime))
        {
			LOG(LOG_LEVEL.ERROR, "UpdateLocation: newTime is NaN");
			return;
		}
		else if (___location_Out.time != newTime)
		{
			double radius = 0.0;
			double trueAnomaly = 0.0;
			new KeplerSolver(___location_Out.planet.mass, __instance.periapsis, __instance.getSpecificEnergy()).GetPositionAtTime((newTime - __instance.periapsisPassageTime) * (double)(__instance.direction), ref radius, ref trueAnomaly);
			double argument = trueAnomaly + __instance.arg;

			// Calculate position
			Double2 position = new Double2(radius * Math.Cos(argument), radius * Math.Sin(argument));

			// Calculate velocity
			Double2 velocity = __instance.GetVelocityAtTrueAnomaly(trueAnomaly);

			Location location = new Location(newTime, ___location_Out.planet, position, velocity);
			
			if(double.IsNaN(trueAnomaly))
            {
				LOG(LOG_LEVEL.ERROR, "UpdateLocation: trueAnomaly is NaN");
            }
			else if(double.IsNaN(location.position.x) || double.IsNaN(location.position.y))
            {
				LOG(LOG_LEVEL.ERROR, "UpdateLocation: position is wrong, x and/or y is NaN");
			}
			else if(double.IsNaN(location.velocity.x) || double.IsNaN(location.velocity.y))
            {
				LOG(LOG_LEVEL.ERROR, "UpdateLocation: velocity is wrong, x and/or y is NaN");
			}
			else
            {
				___location_Out = location;
				___trueAnomaly_Out = trueAnomaly;
			}
		}
	}

	// Local log function
	[Conditional("ACTIVE_LOGS")]
	private static void LOG(LOG_LEVEL level, string message)
	{
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
	}
}


[HarmonyPatch(typeof(Orbit), "GetPoints")]
public class Orbit_GetPoints_Patch
{
	[HarmonyPrefix]
	public static bool GetPoints_Prefix()
	{
		// skip the original method
		return false;
	}


	[HarmonyPostfix]
	public static Vector3[] GetPoints_Postfix(Vector3[] __result, Orbit __instance, double fromTrueAnomaly, double toTrueAnomaly, int resolution, double scaleMultiplier)
    {
		return __instance.getOrbitDrawer().GetPoints(fromTrueAnomaly, toTrueAnomaly, resolution, scaleMultiplier);
	}
}



