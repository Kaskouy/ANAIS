using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

using HarmonyLib;

class NumericalMinimumCalculator<T>
{
    // Internal data
    // -------------
    // The algorithm works with 3 consecutive points. The middle point is the lowest one, so we are sure a minimum exists
    // between point 1 and point 3. 
    private double x1, y1;
    private double x2, y2;
    private double x3, y3;

    // This is additional data associated to the points to allow the caller to store more data
    private T a1;
    private T a2;
    private T a3;

    private bool isLastCalculationAccurate_ = false;

    private uint nbStep = 0;


    // Returns an instance of NumericalMinimumCalculator if data is correct, null otherwise.
    public static NumericalMinimumCalculator<T> createInstance(double x1, double y1, T data1, double x2, double y2, T data2, double x3, double y3, T data3)
    {
        NumericalMinimumCalculator<T> theCalculator = new NumericalMinimumCalculator<T>(x1, y1, data1, x2, y2, data2, x3, y3, data3);

        if (theCalculator.checkDataConsistency())
        {
            LOG(LOG_LEVEL.DEBUG, "NumericalMinimumCalculator: Initialization - List of points is:");
            LOG(LOG_LEVEL.DEBUG, "    x1 = " + x1 + " ; y1 = " + y1);
            LOG(LOG_LEVEL.DEBUG, "    x2 = " + x2 + " ; y2 = " + y2);
            LOG(LOG_LEVEL.DEBUG, "    x3 = " + x3 + " ; y3 = " + y3);

            return theCalculator;
        }
        else
        {
            return null;
        }
    }

    public bool isLastCalculationAccurate()
    {
        return isLastCalculationAccurate_;
    }

    public double getXwidth()
    {
        return x3 - x1;
    }

    public double getMaxYdifference()
    {
        return Math.Max(y1, y3) - y2;
    }

    public void getAdditionalData(out T data1, out T data2, out T data3)
    {
        data1 = a1;
        data2 = a2;
        data3 = a3;
    }


    // Constructor: Initializes the algorithm with 3 points
    private NumericalMinimumCalculator(double x1, double y1, T a1, double x2, double y2, T a2, double x3, double y3, T a3)
    {
        this.x1 = x1;
        this.y1 = y1;
        this.a1 = a1;
        this.x2 = x2;
        this.y2 = y2;
        this.a2 = a2;
        this.x3 = x3;
        this.y3 = y3;
        this.a3 = a3;
    }

    // Checks that the points verify the required conditions for the algorithm to work
    // If this returns false, this most likely means that there was an error from the caller
    private bool checkDataConsistency()
    {
        bool ok = false;

        // Check that x1, x2, x3 are correctly ordered
        if((x1 < x2) && (x2 < x3))
        {
            // Check that point 2 is minimal
            if(!(y2 > y1) && !(y2 > y3))
            {
                ok = true;
            }
        }

        return ok;
    }


    // Inserts a new point (x, y) at the appropriate place. Returns false if the operation failed (point is out of bounds)
    public bool insertNewPoint(double x, double y, T data)
    {
        if((x < x1) || (x > x3))
        {
            // point is out of bounds
            return false;
        }

        double length_orig = x3 - x1;

        // point is within the interval, insert it logically
        if(x < x2)
        {
            if(y < y2)
            {
                // list of points evolves like this: [1, 2, 3] -> [1, new, 2]
                x3 = x2;
                y3 = y2;
                a3 = a2;
                x2 = x;
                y2 = y;
                a2 = data;
            }
            else
            {
                // list of points evolves like this: [1, 2, 3] -> [new, 2, 3]
                x1 = x;
                y1 = y;
                a1 = data;
            }
        }
        else // x > x2
        {
            if (y < y2)
            {
                // list of points evolves like this: [1, 2, 3] -> [2, new, 3]
                x1 = x2;
                y1 = y2;
                a1 = a2;
                x2 = x;
                y2 = y;
                a2 = data;
            }
            else
            {
                // list of points evolves like this: [1, 2, 3] -> [1, 2, new]
                x3 = x;
                y3 = y;
                a3 = data;
            }
        }

        double length_end = x3 - x1;

        nbStep++;

        LOG(LOG_LEVEL.DEBUG, "NumericalMinimumCalculator: Step " + nbStep + " - Length reduction of " + (1.0 - length_end/ length_orig) + " - List of points is:");
        LOG(LOG_LEVEL.DEBUG, "    x1 = " + x1 + " ; y1 = " + y1);
        LOG(LOG_LEVEL.DEBUG, "    x2 = " + x2 + " ; y2 = " + y2);
        LOG(LOG_LEVEL.DEBUG, "    x3 = " + x3 + " ; y3 = " + y3);

        return true;
    }

    // Calculates a new potential x value using a quadratic approximation
    private double calculateNewGuess_QuadraticApproximation()
    {
        double X1 = x2 - x1;
        double X3 = x3 - x2;
        double Y1 = y1 - y2;
        double Y3 = y3 - y2;

        double X1Y3 = X1 * Y3;
        double X3Y1 = X3 * Y1;

        double crossProduct = X3Y1 + X1Y3;

        if(crossProduct < 0.000001 * X1 * X3)
        {
            // Points are practically aligned, we return this by default as we can't do the usual calculation
            LOG(LOG_LEVEL.WARNING, "NumericalMinimumCalculator::calculateNewGuess_QuadraticApproximation : WARNING - points are aligned, return default value");
            return x2;
        }
        else
        {
            double x_min = (X3Y1 * X3 - X1Y3 * X1) / (2.0 * crossProduct) + x2;
            return x_min;
        }
    }

    private double getIntervalReductionFactor(double fractionOfInterval)
    {
        double f = fractionOfInterval;
        
        // to work with the fraction of the greatest part of the interval
        if (f < 0.5) f = 1.0 - f;

        double alpha = Math.Sqrt((1.0 + f) * (1.0 + f) - 4.0 * f * f);
        alpha = (1.0 + f - alpha) / 2.0;
        return alpha;
    }

    public double calculateNewGuess()
    {
        const double C_SMART_CUT_THRESHOLD = 0.2;
        const double C_MAX_REDUCTION_FACTOR = 0.94;

        double length = x3 - x1;
        double ratio = (x2 - x1) / length; // expresses where x2 is located in the [x1, x3] interval (0: x2 = x1; 1: x2 = x3)

        if(ratio < C_SMART_CUT_THRESHOLD)
        {
            // x2 is close to x1 -> lower the right bound
            double newRatio = getIntervalReductionFactor(1.0 - ratio); // calculate the most appropriate reduction factor
            if (newRatio > C_MAX_REDUCTION_FACTOR) newRatio = C_MAX_REDUCTION_FACTOR; // Do not reduce by more than that value (to avoid iterating too much)
            isLastCalculationAccurate_ = false;
            LOG(LOG_LEVEL.DEBUG, "--- New guess : Cut the right part by a factor of " + newRatio);
            return x3 - newRatio * length;
        }
        else if(1.0 - ratio < C_SMART_CUT_THRESHOLD)
        {
            // x2 is close to x3 -> raise the left bound
            double newRatio = getIntervalReductionFactor(ratio); // calculate the most appropriate reduction factor
            if (newRatio > C_MAX_REDUCTION_FACTOR) newRatio = C_MAX_REDUCTION_FACTOR; // Do not reduce by more than that value (to avoid iterating too much)
            isLastCalculationAccurate_ = false;
            LOG(LOG_LEVEL.DEBUG, "--- New guess : Cut the left part by a factor of " + newRatio);
            return x1 + newRatio * length;
        }
        else
        {
            // x2 is in an intermediate position between x1 and x3 -> use the quadratic approximation
            isLastCalculationAccurate_ = true;
            LOG(LOG_LEVEL.DEBUG, "--- New guess : Use quadratic approximation");
            return calculateNewGuess_QuadraticApproximation();
        }
    }

    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.NUMERICAL_MINI_CALCULATOR, level, message);
    }
}
