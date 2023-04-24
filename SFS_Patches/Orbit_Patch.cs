using System.Runtime.CompilerServices;

using HarmonyLib;
using SFS.World;
using SFS.WorldBase;
using System;
using UnityEngine;



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
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	private static void LOG(LOG_LEVEL level, string message)
	{
#if ACTIVE_LOGS
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
#endif
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
			KeplerSolver.GetTimeAtPosition(location.planet.mass, __instance.periapsis, specificEnergy, location.Radius, trueAnomaly_Out, ref timeFromPeri);
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
					KeplerSolver.GetTimeAtPosition(location.planet.mass, __instance.periapsis, specificEnergy, location.planet.SOI, loc_trueAnomaly, ref timeFromPeri);
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
}


[HarmonyPatch(typeof(Orbit), MethodType.Constructor, new Type[] { typeof(double), typeof(double), typeof(double), typeof(int), typeof(Planet), typeof(PathType), typeof(Planet) })]
public class Orbit_OrbitConstructor2_Patch
{
	// This is a hack that tells the constructor to interpret the sma parameter as slr instead.
	// In exchange, it allows to create parabolic or hyperbolic orbits, which the original constructor isn't made for
	public static bool useHackedConstructor = false; 
	
	[HarmonyPrefix]
	public static bool Orbit_prefix()
	{
		return false;
	}

	[HarmonyPostfix]
	public static void Orbit_postfix(Orbit __instance, double sma, double ecc, double arg, int direction, Planet planet, PathType pathType, Planet nextPlanet)
    {
		Traverse.Create(__instance).Field("location_Out").SetValue(new Location(planet, Double2.zero, Double2.zero));

		__instance.pathType = pathType;
		Traverse.Create(__instance).Field("nextPlanet").SetValue(nextPlanet);

		if (pathType == PathType.Eternal)
		{
			__instance.orbitEndTime = double.PositiveInfinity;
		}

		if (useHackedConstructor)
        {
			// HACKED CONSTRUCTOR: slr has been passed instead of sma (so we can handle the case of parabolas)
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
		}
		else
        {
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
		}
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
		double angularMomentum = Orbit_Utils.GetAngularMomentum(__instance);
		double specificEnergy = Orbit_Utils.GetSpecificEnergy(__instance);

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
		double specificEnergy = Orbit_Utils.GetSpecificEnergy(__instance);
		double timeFromPeri = 0.0;
		double radius = __instance.slr / (1.0 + __instance.ecc * Math.Cos(trueAnomaly));
		KeplerSolver.GetTimeAtPosition(___location_Out.planet.mass, __instance.periapsis, specificEnergy, radius, trueAnomaly, ref timeFromPeri);

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

		if (___location_Out.time != newTime && !double.IsNaN(newTime))
		{
			double radius = 0.0;
			double trueAnomaly = 0.0;
			double specificEnergy = Orbit_Utils.GetSpecificEnergy(__instance);
			KeplerSolver.GetPositionAtTime(___location_Out.planet.mass, __instance.periapsis, specificEnergy, (newTime - __instance.periapsisPassageTime) * (double)(__instance.direction), ref radius, ref trueAnomaly);
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
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	private static void LOG(LOG_LEVEL level, string message)
	{
#if ACTIVE_LOGS
		AnaisLogger.Log(LOG_CATEGORY.ORBIT, level, message);
#endif
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
		return new OrbitDrawer(__instance).GetPoints(fromTrueAnomaly, toTrueAnomaly, resolution, scaleMultiplier);
	}
}

