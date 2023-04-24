using System;
using UnityEngine;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

public struct PolarPoint
{
    private double radius;

    private double argument;

    public double Radius => radius;

    public double Argument => argument;

    public PolarPoint(double rad, double arg)
    {
        radius = rad;
        argument = arg;

        // make sure radius is positive
        if(radius < 0.0) 
        {
            radius = -radius;
            argument += Math.PI;
        }

        // make sure argument is between -Pi and Pi
        NormalizeArgument();
    }

    public PolarPoint(Double2 point)
    {
        radius = point.magnitude;
        argument = point.AngleRadians;
    }

    public void setRadius(double rad)
    {
        radius = rad;

        if(radius < 0.0)
        {
            radius = -radius;
            argument += Math.PI;
            NormalizeArgument();
        }
    }

    public void setArgument(double arg)
    {
        argument = arg;
        NormalizeArgument();
    }

    public void Rotate(double delta_arg)
    {
        argument += delta_arg;
        NormalizeArgument();
    }

    public static double Dot(PolarPoint p1, PolarPoint p2)
    {
        return (p1.radius * p2.radius * Math.Cos(p2.argument - p1.argument));
    }

    public static double CrossProduct(PolarPoint p1, PolarPoint p2)
    {
        return (p1.radius * p2.radius * Math.Sin(p2.argument - p1.argument));
    }

    public Double2 Normalized()
    {
        if (radius > 0.0) { return new Double2(Math.Cos(argument), Math.Sin(argument)); }
        else              { return new Double2(0.0, 0.0); }
    }

    public Vector2 UnitVector()
    {
        if (radius > 0.0)
        {
            return new Vector2((float)Math.Cos(argument), (float)Math.Sin(argument));
        }
        else
        {
            return new Vector2(0.0f, 0.0f);
        }
    }

    public Vector2 ToVector2 => new Vector2((float)(radius * Math.Cos(argument)), (float)(radius * Math.Sin(argument)));

    public Double2 ToDouble2() => new Double2(radius * Math.Cos(argument), radius * Math.Sin(argument));

    public void NormalizeArgument()
    {
        while (argument >  Math.PI) { argument -= 2.0 * Math.PI; }
        while (argument < -Math.PI) { argument += 2.0 * Math.PI; }
    }
}
