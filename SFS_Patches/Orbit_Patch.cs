using HarmonyLib;
using SFS.World;
using SFS.WorldBase;
using System;
using UnityEngine;

using System.Diagnostics;
using System.Runtime.CompilerServices;



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

        // BIG HACK COMING!!!
        // Those 2 variables (specific energy/angular momentum) don't exist in the original Orbit class...
        // Luckily, there are 2 existing variables (semiMinorAxis/meanMotion) in the Orbit class that I don't need myself...
        // ... So yeah, it's ugly, but I store specific energy in semiMinorAxis, and angular momentum in meanMotion
        // I tried using a ConditionalWeakTable to solve that problem properly, but it appeared to occasionally create some breaking bugs.
        Double3 angularMomentum_vector = Double3.Cross(location.position, location.velocity);
		Double2 eccentricity_vector = (Double2)(Double3.Cross((Double3)location.velocity, angularMomentum_vector) / location.planet.mass) - location.position.normalized;

		double specificEnergy = Math.Pow(location.velocity.magnitude, 2.0) / 2.0 - location.planet.mass / location.Radius;
		double angularMomentum = angularMomentum_vector.z;

        __instance.semiMinorAxis = specificEnergy; // HACK
        __instance.meanMotion = angularMomentum;   // HACK

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
			__instance.periapsisPassageTime = location.time - timeFromPeri * (double)(__instance.direction) ; //- 10.0 * __instance.period

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
	public const int C_HACKED_CONSTRUCTOR_VAL = 7777;

    // This is a second hack to tell the constructor to act as a copy constructor
    // Set orbitToCopy to the orbit to duplicate, and pass C_HACK_COPY_CONSTRUCTOR as direction parameter 
    public const int C_HACK_COPY_CONSTRUCTOR = 666;
	public static Orbit orbitToCopy = null;


	[HarmonyPrefix]
	public static bool Orbit_prefix()
	{
		return false;
	}

	[HarmonyPostfix]
	public static void Orbit_postfix(Orbit __instance, double sma, double ecc, double arg, int direction, Planet planet, PathType pathType, Planet nextPlanet)
    {
		if(Math.Abs(direction) == C_HACK_COPY_CONSTRUCTOR)
		{
			// COPY CONSTRUCTOR
			// ----------------
			__instance.sma = orbitToCopy.sma;
            __instance.slr = orbitToCopy.slr;
            __instance.ecc = orbitToCopy.ecc;
            __instance.arg = orbitToCopy.arg;
			__instance.direction = orbitToCopy.direction;
            __instance.period = orbitToCopy.period;
			__instance.periapsis = orbitToCopy.periapsis;
            __instance.apoapsis = orbitToCopy.apoapsis;
            __instance.pathType = orbitToCopy.pathType;
            __instance.orbitStartTime = orbitToCopy.orbitStartTime;
            __instance.orbitEndTime = orbitToCopy.orbitEndTime;
			__instance.periapsisPassageTime = orbitToCopy.periapsisPassageTime;

            Traverse.Create(__instance).Field("location_Out").SetValue(new Location(orbitToCopy.Planet, Double2.zero, Double2.zero));
            Traverse.Create(__instance).Field("nextPlanet").SetValue(orbitToCopy.NextPlanet);

            double specificEnergy = orbitToCopy.Planet.mass * (__instance.ecc * __instance.ecc - 1.0) / (2.0 * __instance.slr);
            double angularMomentum = Math.Sqrt(orbitToCopy.Planet.mass * __instance.slr) * __instance.direction;

			__instance.semiMinorAxis = specificEnergy; // HACK
			__instance.meanMotion = angularMomentum;   // HACK
        }
		else // Usual constructors (either the normal one, or the one that takes slr instead of sma as parameter)
		{ 
			// Determine if we use the nornal or the hacked version of the constructor.
			// The hacked constructor interprets the sma as the slr, so that we can handle the case of parabolas (sma is infinite for those and the calculations are not possible)
			// In this case, the direction parameter (which is usually -1 or 1) is also multiplied by a fixed value. This is how we know that it's the hacked version that we intend to call.
			bool useHackedConstructor = false;

			if(Math.Abs(direction) == C_HACKED_CONSTRUCTOR_VAL)
			{
				useHackedConstructor = true;
			}

			// Use a valid direction parameter in any case
			int dir = (direction > 0)?1:-1;

			//LOG(LOG_LEVEL.DEBUG, "Setting location_out");
			Traverse.Create(__instance).Field("location_Out").SetValue(new Location(planet, Double2.zero, Double2.zero));

			__instance.pathType = pathType;
			//LOG(LOG_LEVEL.DEBUG, "Setting nextPlanet");
			Traverse.Create(__instance).Field("nextPlanet").SetValue(nextPlanet);

			if (pathType == PathType.Eternal)
			{
				__instance.orbitEndTime = double.PositiveInfinity;
			}


			__instance.ecc = ecc;
			__instance.arg = arg;
			__instance.direction = dir;

			if (useHackedConstructor)
			{
				// HACKED CONSTRUCTOR: slr has been passed instead of sma (so we can handle the case of parabolas)
				// ------------------
				__instance.slr = sma;
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
				// NORMAL CONSTRUCTOR
				// ------------------
				__instance.sma = sma;
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

			double specificEnergy = planet.mass * (ecc * ecc - 1.0) / (2.0 * __instance.slr);
			double angularMomentum = Math.Sqrt(planet.mass * __instance.slr) * dir;
            __instance.semiMinorAxis = specificEnergy; // HACK
			__instance.meanMotion = angularMomentum;   // HACK
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
	public static void GetVelocityAtTrueAnomaly_Postfix(ref Double2 __result, Orbit __instance, ref Location ___location_Out, double trueAnomaly)
    {
        double normalizedTrueAnomaly = Orbit_Utils.NormalizeAngle(trueAnomaly);
		double radiusAtTrueAnomaly = Kepler.GetRadiusAtTrueAnomaly(__instance.slr, __instance.ecc, normalizedTrueAnomaly);
		
		double V_theta = __instance.getAngularMomentum() / radiusAtTrueAnomaly;
		double V2 = 2.0 * (__instance.getSpecificEnergy() + ___location_Out.planet.mass / radiusAtTrueAnomaly);
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
	public static double GetLastTrueAnomalyPassTime_Postfix(double __result, Orbit __instance, double time, double trueAnomaly)
    {
        double timeFromPeri = 0.0;
		double radius = __instance.slr / (1.0 + __instance.ecc * Math.Cos(trueAnomaly));
		new KeplerSolver(__instance.Planet.mass, __instance.periapsis, __instance.getSpecificEnergy()).GetTimeAtPosition(radius, trueAnomaly, ref timeFromPeri);

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
        if (___location_Out.time != newTime)
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

            /*if (double.IsNaN(trueAnomaly))
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
            {*/
				___location_Out = location;
				___trueAnomaly_Out = trueAnomaly;
			/*}*/
        }
    }
}


[HarmonyPatch(typeof(Orbit), "GetLocation")]
public class Orbit_GetLocation_Patch
{
	[HarmonyPrefix]
	public static bool GetLocation_Prefix(ref Location __result, Orbit __instance, double time)
	{
		if(AnaisManager.isAnaisThreadRunning() == false)
        {
			// Run normally the original method if this method is called by the main thread
			return true;
        }
		else
        {
			// If the method is called by the ANAIS thread, do the same job as UpdateLocation(), but without updating
			// the internal variables (location_out and trueAnomaly_out) - This is to avoid a clash due to multithreading
			// Not the cleanest way to do this, but probably one of the simplest and the safest (No modification of Stef's code behaviour, no need for mutex)
			double radius = 0.0;
			double trueAnomaly = 0.0;
			new KeplerSolver(__instance.Planet.mass, __instance.periapsis, __instance.getSpecificEnergy()).GetPositionAtTime((time - __instance.periapsisPassageTime) * (double)(__instance.direction), ref radius, ref trueAnomaly);
			double argument = trueAnomaly + __instance.arg;

			// Calculate position
			Double2 position = new Double2(radius * Math.Cos(argument), radius * Math.Sin(argument));

			// Calculate velocity
			Double2 velocity = __instance.GetVelocityAtTrueAnomaly(trueAnomaly);

			// Set the resulting location
			__result = new Location(time, __instance.Planet, position, velocity);

			// Skip the original
			return false;
		}
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



