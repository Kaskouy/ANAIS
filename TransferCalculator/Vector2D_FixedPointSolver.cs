using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

// -------------------------------------------------------------------------------
//                             Vector2D_FixedPointSolver
// -------------------------------------------------------------------------------
// 
// This class offers a high level interface to help to solve equations of type
// [Fx(x, y), Fy(x, y)] = [x, y]
// Basically, it's the same as equations of the form f(x) = x, but in 2 dimensions
// This class is meant to be used when you have to solve such an equation that can
// not be solved analytically, so the only way is to try couples of values (x, y)
// until we find a reasonably close solution.
//
// From a user's point of view, it offers 2 methods of interest:
// - addPoint : allows to specify a new point, i.e. a couple of values (x, y), and
//   the associated values [Fx, Fy]. That point will be registered and used later
//   to propose new estimates
// - GetEstimatedSolution : based on the previously registered points, returns a
//   couple of values that is the most likely to correspond to the solution.
// 
// NOTE: The responsibility to test if a couple of values is close enough to the
//       solution is left to the caller. Some algorithms may fail in case the
//       research is pushed too far.
// -------------------------------------------------------------------------------
class Vector2D_FixedPointSolver
{

    // The representation of a point, and its image under the 2D function we try to solve
    private class Vector2D_ArgAndVal
    {
        public Double2 arg; // The argument of the 2D function
        public Double2 val; // The value of the 2D function for argument arg 

        public Vector2D_ArgAndVal(Double2 arg, Double2 val)
        {
            this.arg = arg;
            this.val = val;
        }
    }

    // An enum used to keep track of the lastly used guessing method
    // Basically, here are the 4 methods we can use:
    // - FIXED_POINT: uses the result of a previous point as the new estimate. Simple, needs only 1 existing point, convergence speed is variable
    // - AVERAGE_DAMPING: takes the mean of a point's argument and its associated value. converges potentially slowerly, but decreases the risk of divergence
    // - LINEAR_APPROXIMATION: guesses a new estimate based on a linear approximation based on 2 points - more accurate than previous methods but needs 2 points
    // - MULTI_SECANT: builds a local linear representation based on 3 non aligned points (similar to the secant method in 1D) - most accurate, but needs 3 non aligned points
    private enum E_GUESSING_METHOD
    {
        NONE = 0,
        FIXED_POINT,
        AVERAGE_DAMPING,
        LINEAR_APPROXIMATION,
        MULTI_SECANT
    }

    // Because we use 3 points that must not be aligned for the multi-secant method, we define the minimal angle that the 3 points must form.
    private const double C_SIN_MIN_ANGLE = 0.04;

    // Sometimes we will intentionally shift a point so that it's not aligned with 2 existing points, to be able to apply the multi-secant method on the next shot.
    // The angle is intentionally higher than the one corresponding to C_SIN_MIN_ANGLE to fix for numerical imprecisions
    private const double C_SIN_SHIFTING_ANGLE = 0.05;
    private const double C_COS_SHIFTING_ANGLE = 0.99874921777190895;

    private List<Vector2D_ArgAndVal> _listPoints;

    private uint numStep;

    private bool isFixedPointMethodSuitable;

    private bool hasLinearApproximationFailed;

    private bool hasMultiSecantFailed;

    private E_GUESSING_METHOD _lastMethodUsed;

    private Vector2D_ArgAndVal _lastPoint;


    public Vector2D_FixedPointSolver()
    {
        _listPoints = new List<Vector2D_ArgAndVal>();
        numStep = 0;
        isFixedPointMethodSuitable = true;
        hasMultiSecantFailed = false;
        _lastMethodUsed = E_GUESSING_METHOD.NONE;
        _lastPoint = null;
    }

    public void addPoint(Double2 arg, Double2 val)
    {
        // Adds this new point to our collection of points
        Vector2D_ArgAndVal point2D = new Vector2D_ArgAndVal(arg, val);

        // update the last point
        _lastPoint = point2D;

        //_listPoints.Add(point2D); // Adds point to the end of the list

        numStep += 1;

        LOG(LOG_LEVEL.DEBUG, "Step " + numStep + " - Added New point : Arg : x = " + point2D.arg.x + "; y = " + point2D.arg.y);
        LOG(LOG_LEVEL.DEBUG, "                           Val : x = " + point2D.val.x + "; y = " + point2D.val.y);
        LOG(LOG_LEVEL.DEBUG, "                           Delta = " + GetDelta(point2D));


        // Insert point in the list so that it's sorted (best point at end of list)
        // --------------------------------------------
        double delta2 = GetDeltaSquared(point2D);
        bool inserted = false;
        int i = _listPoints.Count;

        // Loop from the end of the list - we start from the end because we expect the new points to be better
        while((i > 0) && !inserted)
        {
            // if delta2 improves compared to previous postion, insert it there
            if (delta2 < GetDeltaSquared(_listPoints[i - 1]))
            {
                _listPoints.Insert(i, point2D);
                inserted = true;
            }
            else
            {
                i--; // only decrease i if point was not inserted, so that i is the index at which the point has been inserted
            }
        }

        // insert point to the beginning of the list if it's not inserted yet. In particular for the first point...
        if(!inserted)
        {
            _listPoints.Insert(i, point2D);
        }


        // Lift a warning if the point happens to be a worse estimation than before
        // ------------------------------------------------------------------------
        if(_lastMethodUsed == E_GUESSING_METHOD.FIXED_POINT)
        {
            if(i < _listPoints.Count - 1)
            {
                LOG(LOG_LEVEL.INFO, "  Fixed point method gave a worse result than before - Stop using it from now on.");
                isFixedPointMethodSuitable = false;
            }
        }
        else if(_lastMethodUsed == E_GUESSING_METHOD.AVERAGE_DAMPING)
        {
            if (i < _listPoints.Count - 1)
            {
                LOG(LOG_LEVEL.INFO, "  Average damping method gave a worse result than before - This is going to be hard.");
            }
        }
        else if(_lastMethodUsed == E_GUESSING_METHOD.LINEAR_APPROXIMATION)
        {
            if (i < _listPoints.Count - 2)
            {
                LOG(LOG_LEVEL.INFO, "  Linear based approximation method gave a worse result than before - This is bad!");
                hasLinearApproximationFailed = true;
            }
        }
        else if(_lastMethodUsed == E_GUESSING_METHOD.MULTI_SECANT)
        {
            if (i < _listPoints.Count - 3)
            {
                LOG(LOG_LEVEL.INFO, "  Multi-secant method gave a worse result than the 3 previous points - Now we are doomed!");
                hasMultiSecantFailed = true;
            }
        }


        // Reset "failed" flags if the new point has been inserted in a favorable position
        // --------------------
        if(i >= _listPoints.Count - 3)
        {
            // New point is part of the 3 best: We can attempt again the multi-secant method
            hasMultiSecantFailed = false;

            // New point is part of the 2 best: We can attempt again the linear approximation method
            if (i >= _listPoints.Count - 2) hasLinearApproximationFailed = false;
        }

        //_listPoints.Insert(0, point2D);  // Adds point to the beginning of the list
        //_listPoints.RemoveAt(0); // Removes first point
        //_listPoints.RemoveAt(_listPoints.Count-1); // Removes last point
    }

    public Double2 GetEstimatedSolution()
    {
        Double2 newPoint;

        // Try applying the multi-secant method
        // ------------------------------------
        if (_listPoints.Count > 2) 
        {
            bool success = false;

            if (!hasMultiSecantFailed)
            {
                // Apply the multi-secant method (presumably the best)
                success = ApplyMultiSecantMethod(_listPoints[_listPoints.Count - 3], _listPoints[_listPoints.Count - 2], _listPoints[_listPoints.Count - 1], out newPoint);
            }
            else
            {
                // If the multi-secant method has previously failed, we replace the worse point with the last point - this is generally enough to restore convergence
                success = ApplyMultiSecantMethod(_lastPoint, _listPoints[_listPoints.Count - 2], _listPoints[_listPoints.Count - 1], out newPoint);
            }

            if (success)
            {
                _lastMethodUsed = E_GUESSING_METHOD.MULTI_SECANT;
                return newPoint;
            }
        }


        // If we couldn't, try applying the linear approximation method
        // ------------------------------------------------------------
        if(_listPoints.Count >= 2)
        {
            bool success = false;

            if (!hasLinearApproximationFailed)
            {
                // Apply the linear approximation method method (the best after multi-secant)
                success = ApplyLinearApproximationMethod(_listPoints[_listPoints.Count - 2], _listPoints[_listPoints.Count - 1], out newPoint);
            }
            else
            {
                // If the method previously failed, we replace the worse point with the last point - this is generally enough to restore convergence
                success = ApplyLinearApproximationMethod(_lastPoint, _listPoints[_listPoints.Count - 1], out newPoint);
            }

            if (success)
            {
                _lastMethodUsed = E_GUESSING_METHOD.LINEAR_APPROXIMATION;
                return newPoint;
            }
        }


        // If not applicable, apply a more basic method
        // --------------------------------------------
        if (isFixedPointMethodSuitable)
        {
            _lastMethodUsed = E_GUESSING_METHOD.FIXED_POINT;
            newPoint = ApplyFixedPointMethod(_listPoints[_listPoints.Count - 1]);
        }
        else
        {
            _lastMethodUsed = E_GUESSING_METHOD.AVERAGE_DAMPING;
            newPoint = ApplyAverageDampingMethod(_listPoints[_listPoints.Count - 1]);
        }

        if (_listPoints.Count > 1)
        {
            // If we already have 2 points available, make sure the new point is not aligned with the 2 best points.
            // This is to allow the use of the multi-secant method on the next shot.
            bool shifted = ApplyShifting(_listPoints[_listPoints.Count - 2], _listPoints[_listPoints.Count - 1], ref newPoint);
        }

        return newPoint;
    }

    private double GetDeltaSquared(Vector2D_ArgAndVal point)
    {
        double delta_x = point.val.x - point.arg.x;
        double delta_y = point.val.y - point.arg.y;
        double deltaSqr = delta_x * delta_x + delta_y * delta_y;

        return deltaSqr;
    }

    private double GetDelta(Vector2D_ArgAndVal point)
    {
        double delta_x = point.val.x - point.arg.x;
        double delta_y = point.val.y - point.arg.y;
        double delta = Math.Sqrt(delta_x * delta_x + delta_y * delta_y);

        return delta;
    }

    // Applies the fixed point method: returns the result from the previous point
    // This method is the simplest and the most classic one.
    // Works fine in many cases, but can cause a divergence in tricky situations
    private Double2 ApplyFixedPointMethod(Vector2D_ArgAndVal point) 
    {
        LOG(LOG_LEVEL.INFO, "    Apply fixed point method");
        Double2 estimatedVal = point.val;

        // Return value
        LOG(LOG_LEVEL.DEBUG, "      New point: x = " + estimatedVal.x + "; y = " + estimatedVal.y);
        return estimatedVal;
    }

    // Applies the average damping method: returns the average value between the argument and the resulting value
    // This methods is most likely to converge but may be slower in some cases
    private Double2 ApplyAverageDampingMethod(Vector2D_ArgAndVal point)
    {
        LOG(LOG_LEVEL.INFO, "    Apply average damping method");
        Double2 estimatedVal = new Double2(0.0, 0.0);

        estimatedVal.x = (point.arg.x + point.val.x) / 2.0;
        estimatedVal.y = (point.arg.y + point.val.y) / 2.0;

        // Return value
        LOG(LOG_LEVEL.DEBUG, "      New point: x = " + estimatedVal.x + "; y = " + estimatedVal.y);
        return estimatedVal;
    }

    // Applies the multi-secant method, which consists in building a local linear representation based on 3 non aligned points
    // Gives very good results, but it needs 3 non aligned points
    private bool ApplyMultiSecantMethod(Vector2D_ArgAndVal P1, Vector2D_ArgAndVal P2, Vector2D_ArgAndVal P3, out Double2 result)
    {
        LOG(LOG_LEVEL.INFO, "    Apply multi-secant method");
        result = new Double2(0.0, 0.0);

        double x1 = P1.arg.x;
        double y1 = P1.arg.y;
        double x2 = P2.arg.x;
        double y2 = P2.arg.y;
        double x3 = P3.arg.x;
        double y3 = P3.arg.y;

        double Vx1 = P1.val.x;
        double Vy1 = P1.val.y;
        double Vx2 = P2.val.x;
        double Vy2 = P2.val.y;
        double Vx3 = P3.val.x;
        double Vy3 = P3.val.y;

        double denominator = y1 * (x3 - x2) + y2 * (x1 - x3) + y3 * (x2 - x1);

        // Check denominator is not too close to 0
        // ---------------------------------------
        {
            // Note that those are actually the distances squared, to allow lighter calculations
            double dist1_2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
            double dist1_3 = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1);
            double dist2_3 = (x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2);
            double productDist = 0.0;

            // We calculate the product of the 2 smallest distances
            if((dist1_2 > dist1_3) && (dist1_2 > dist2_3))
            {
                // dist1_2 is the greatest
                productDist = dist1_3 * dist2_3;
            }
            else if((dist1_3 > dist1_2) && (dist1_3 > dist2_3))
            {
                // dist1_3 is the greatest
                productDist = dist1_2 * dist2_3;
            }
            else
            {
                // dist2_3 is the greatest
                productDist = dist1_2 * dist1_3;
            }

            if(denominator * denominator < C_SIN_MIN_ANGLE * C_SIN_MIN_ANGLE * productDist) // means that the angle between 2 vectors is close to 0º or 180º by less than 2.5º (approx.)
            {
                LOG(LOG_LEVEL.INFO, "    Multi-secant method failed - points are closed from alignment");
                return false;
            }
        }

        // calculate the coefficients so that we have:
        // Vx = Axx * x + Axy * y + Cx
        // Vy = Ayx * x + Ayy * y + Cy

        // NOTE: each of those coefficient should be divided by the denominator... But if we pretend
        // that all coefficients are multiplied by that value, then the equation is exactly the same! So we avoid the divisions.
        double Axx = (Vx1 * (y2 - y3) + Vx2 * (y3 - y1) + Vx3 * (y1 - y2)) /*/ denominator*/;
        double Axy = (Vx1 * (x3 - x2) + Vx2 * (x1 - x3) + Vx3 * (x2 - x1)) /*/ denominator*/;
        double Ayx = (Vy1 * (y2 - y3) + Vy2 * (y3 - y1) + Vy3 * (y1 - y2)) /*/ denominator*/;
        double Ayy = (Vy1 * (x3 - x2) + Vy2 * (x1 - x3) + Vy3 * (x2 - x1)) /*/ denominator*/;

        double Cx  = (Vx1 * (x2*y3 - x3*y2) + Vx2 * (x3*y1 - x1*y3) + Vx3 * (x1*y2 - x2*y1)) /*/ denominator*/;
        double Cy  = (Vy1 * (x2*y3 - x3*y2) + Vy2 * (x3*y1 - x1*y3) + Vy3 * (x1*y2 - x2*y1)) /*/ denominator*/;

        // Verification
        //LOG(LOG_LEVEL.DEBUG, "  Vx1 = " + (Axx * x1 + Axy * y1 + Cx) / denominator + " - Expected: " + Vx1);
        //LOG(LOG_LEVEL.DEBUG, "  Vy1 = " + (Ayx * x1 + Ayy * y1 + Cy) / denominator + " - Expected: " + Vy1);
        //LOG(LOG_LEVEL.DEBUG, "  Vx2 = " + (Axx * x2 + Axy * y2 + Cx) / denominator + " - Expected: " + Vx2);
        //LOG(LOG_LEVEL.DEBUG, "  Vy2 = " + (Ayx * x2 + Ayy * y2 + Cy) / denominator + " - Expected: " + Vy2);
        //LOG(LOG_LEVEL.DEBUG, "  Vx3 = " + (Axx * x3 + Axy * y3 + Cx) / denominator + " - Expected: " + Vx3);
        //LOG(LOG_LEVEL.DEBUG, "  Vy3 = " + (Ayx * x3 + Ayy * y3 + Cy) / denominator + " - Expected: " + Vy3);

        // Now we have to solve the matrix equation: [A] * V + C = V, in other terms:
        // Axx * Vx + Axy * Vy + Cx = Vx   <=>   (Axx - 1) * Vx +  Axy      * Vy + Cx = 0
        // Ayx * Vx + Ayy * Vy + Cy = Vy   <=>    Ayx      * Vx + (Ayy - 1) * Vy + Cy = 0

        // calculation of the determinant ("- 1" became "- denominator" because all coefficients were multiplied by denominator to save a few divisions - see optimization move above)
        double det = ((Axx - denominator) * (Ayy - denominator) - Axy * Ayx);

        // Solve for the result - We still have to divide by det, but we will first check that it's not too low.
        result.x = Cy * Axy - Cx * (Ayy - denominator);
        result.y = Cx * Ayx - Cy * (Axx - denominator);

        double xmax = (x2 > x1) ? ((x3 > x2) ? x3 : x2) : ((x3 > x1) ? x3 : x1);
        double ymax = (y2 > y1) ? ((y3 > y2) ? y3 : y2) : ((y3 > y1) ? y3 : y1);

        if( (Math.Abs(result.x) >= 100.0 * Math.Abs(det) * Math.Abs(xmax)) || (Math.Abs(result.y) >= 100.0 * Math.Abs(det) * Math.Abs(ymax)) )
        {
            LOG(LOG_LEVEL.INFO, "    Multi-secant method failed - New point is obviously too far (determinant must be practically null)");
            return false;
        }

        result.x /= det;
        result.y /= det;

        //LOG(LOG_LEVEL.DEBUG, "  Vx = " + result.x + "  Vx(img) = " + (Axx * result.x + Axy * result.y + Cx) / denominator);
        //LOG(LOG_LEVEL.DEBUG, "  Vy = " + result.y + "  Vy(img) = " + (Ayx * result.x + Ayy * result.y + Cy) / denominator);

        LOG(LOG_LEVEL.DEBUG, "      New point: x = " + result.x + "; y = " + result.y);

        return true;
    }

    // An internal method that is used to slightly modify a point to guarantee that it's not aligned with the 2 points
    // given as parameters. The point is moved in the slightiest way possible to break the alignment.
    private bool ApplyShifting(Vector2D_ArgAndVal P1, Vector2D_ArgAndVal P2, ref Double2 result)
    {
        double x1 = P1.arg.x;
        double y1 = P1.arg.y;
        double x2 = P2.arg.x;
        double y2 = P2.arg.y;
        double x3 = result.x;
        double y3 = result.y;

        double xc = (x1 + x2) / 2.0;
        double yc = (y1 + y2) / 2.0;

        double dist_12 = Math.Sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        //double dist_13 = Math.Sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1));
        //double dist_23 = Math.Sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));

        //double dist_c3 = Math.Sqrt((x3 - xc) * (x3 - xc) + (y3 - yc) * (y3 - yc));


        // determine the coordinates of the closest point
        // if the scalar product P1P2.CP3 is negative, P3 is closer to P1, otherwise it's closer to P2.
        double xref, yref;
        double dist_P3_ref;
        
        {
            double scalar_product = (x3 - xc) * (x2 - x1) + (y3 - yc) * (y2 - y1);

            if (scalar_product < 0.0)
            {
                xref = x1;
                yref = y1;
                //LOG(LOG_LEVEL.DEBUG, "   New point is closer to point 1...");
            }
            else
            {
                xref = x2;
                yref = y2;
                //LOG(LOG_LEVEL.DEBUG, "   New point is closer to point 2...");
            }

            dist_P3_ref = Math.Sqrt((x3 - xref) * (x3 - xref) + (y3 - yref) * (y3 - yref));
            //LOG(LOG_LEVEL.DEBUG, "   dist_12 = " + dist_12 + "; dist_P3_ref = " + dist_P3_ref);
        }

        // the point will be shifted so that the angle between P1P2 et CP3 is at least of epsilon
        double cos_eps = C_COS_SHIFTING_ANGLE;
        double sin_eps = C_SIN_SHIFTING_ANGLE;

        double vectorProduct = (x2 - x1) * (y3 - yref) - (x3 - xref) * (y2 - y1);
        //LOG(LOG_LEVEL.DEBUG, "   vectorProduct = " + vectorProduct);

        if (Math.Abs(vectorProduct) > sin_eps * dist_12 * dist_P3_ref)
        {
            // No need of shifting
            //LOG(LOG_LEVEL.DEBUG, "  Shifting not needed: sinus of angle = " + (vectorProduct / (dist_12 * dist_P3_ref)));
            return false;
        }

        LOG(LOG_LEVEL.DEBUG, "  Shifting needed!");

        // Otherwise, shifting is needed!
        double scalarProduct = (x2 - x1) * (x3 - xref) + (y2 - y1) * (y3 - yref);

        

        if(((vectorProduct > 0) && (scalarProduct < 0)) || ((vectorProduct < 0) && (scalarProduct > 0)))
        {
            // the angle formed is slightly below 180º or slightly below 0º
            // --> point will be shifted following anti-trigonometrical sense
            //LOG(LOG_LEVEL.DEBUG, "  Shifting by negative angle");
            sin_eps = -sin_eps;
        }

        // calculate the normalized vector that gives the direction of the line over which will be projected P3
        double ux = ((x2 - x1) * cos_eps - (y2 - y1) * sin_eps) / dist_12;
        double uy = ((x2 - x1) * sin_eps + (y2 - y1) * cos_eps) / dist_12;

        //double cosEpsilon_verif = ((x2 - x1) * ux + (y2 - y1) * uy) / dist_12;
        //double sinEpsilon_verif = ((x2 - x1) * uy - (y2 - y1) * ux) / dist_12;

        //LOG(LOG_LEVEL.DEBUG, "  Verification: cosEps = " + cosEpsilon_verif + "; sinEps = " + sinEpsilon_verif);

        result.x = xref + ((x3 - xref) * ux + (y3 - yref) * uy) * ux;
        result.y = yref + ((x3 - xref) * ux + (y3 - yref) * uy) * uy;

        LOG(LOG_LEVEL.DEBUG, "  Shifting: old x = " + x3 + "; new x = " + result.x);
        LOG(LOG_LEVEL.DEBUG, "  Shifting: old y = " + y3 + "; new y = " + result.y);

        // Verif
        //double new_dist = Math.Sqrt((result.x - xref) * (result.x - xref) + (result.y - yref) * (result.y - yref));

        //double newScalarProduct = (x2 - x1) * (result.x - xref) + (result.y - yref) * (y2 - y1);
        //double newCrossProduct = (x2 - x1) * (result.y - yref) - (result.x - xref) * (y2 - y1);

        //double newSin = newCrossProduct / (dist_12 * new_dist);

        //LOG(LOG_LEVEL.DEBUG, "   Dist from ref to new point = " + new_dist);
        //LOG(LOG_LEVEL.DEBUG, "   newCrossProduct = " + newCrossProduct);
        //LOG(LOG_LEVEL.DEBUG, "   newSin = " + newSin + "; newCos = " + (newScalarProduct / (dist_12 * new_dist)));

        return true;

        /*int rotationDir = 0;

        if(dist_13 < dist_23)
        {
            double vectorProduct = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
            if (vectorProduct > 0) rotationDir = 1;
            else rotationDir = -1;
        }
        else
        {
            double vectorProduct = (x1 - x2) * (y3 - y2) - (x3-x2) * (y1 - y2);
            if (vectorProduct > 0) rotationDir = 1;
            else rotationDir = -1;
        }*/
    }


    // Applies a linear approximation method based on 2 existing points.
    // 3 points are needed to fully model a linear representation in a 2D space, but since we only have 2 we model it
    // along the line joining the 2 points, then we choose the point that would be the closest from the solution on
    // this line, and finally apply the fixed point/average damping method on this point, using the value that would
    // correspond to it following the linear model.
    private bool ApplyLinearApproximationMethod(Vector2D_ArgAndVal P1, Vector2D_ArgAndVal P2, out Double2 result)
    {
        LOG(LOG_LEVEL.INFO, "    Apply linear approximation method");

        result = new Double2(0.0, 0.0);

        double deltaX1 = P1.val.x - P1.arg.x;
        double deltaX2 = P2.val.x - P2.arg.x;
        double deltaY1 = P1.val.y - P1.arg.y;
        double deltaY2 = P2.val.y - P2.arg.y;

        double d_deltaX = deltaX2 - deltaX1;
        double d_deltaY = deltaY2 - deltaY1;

        LOG(LOG_LEVEL.DEBUG, "      Point 1: arg = (" + P1.arg.x + ", " + P1.arg.y + "); val = " + P1.val.x + ", " + P1.val.y + ")");
        LOG(LOG_LEVEL.DEBUG, "      Point 2: arg = (" + P2.arg.x + ", " + P2.arg.y + "); val = " + P2.val.x + ", " + P2.val.y + ")");

        // The behaviour of the function will be simulated over the line joining points 1 and 2
        // lambda is a value that allows to identify the point on this line that's closest to our requirements
        double numerator = -(deltaX1 * d_deltaX + deltaY1 * d_deltaY);
        double denominator = d_deltaX * d_deltaX + d_deltaY * d_deltaY;

        // If the resulting lambda value would be really off (It's ideally expected between 0 and 1, or at least not too far from those bounds), don't try it.
        // (Note that this also ensures that the denominator is not null)
        if(Math.Abs(numerator) >= 10.0 * denominator)
        {
            LOG(LOG_LEVEL.INFO, "      Linear approximation method failed - num = " + numerator + "; den = " + denominator);
            return false;
        }

        double lambda = numerator / denominator;
        LOG(LOG_LEVEL.DEBUG, "      lambda = " + lambda);

        // x, y are the coordinates of the point we identified
        result.x = P1.arg.x + lambda * (P2.arg.x - P1.arg.x);
        result.y = P1.arg.y + lambda * (P2.arg.y - P1.arg.y);

        // val_x, val_y are the values that we expect the function to have at that point
        double val_x = P1.val.x + lambda * (P2.val.x - P1.val.x);
        double val_y = P1.val.y + lambda * (P2.val.y - P1.val.y);

        LOG(LOG_LEVEL.DEBUG, "      New point: arg = (" + result.x + ", " + result.y + "); val = " + val_x + ", " + val_y + ")");

        // We apply either the fixed point method or the average damping method to this simulated point to get a closer approximation
        if(isFixedPointMethodSuitable)
        {
            result.x = val_x;
            result.y = val_y;
        }
        else
        {
            result.x = (result.x + val_x) / 2.0;
            result.y = (result.y + val_y) / 2.0;
        }

        LOG(LOG_LEVEL.DEBUG, "      New guess: x = " + result.x + ", y = " + result.y);

        return true;
    }

    // Local log function
    [Conditional("ACTIVE_LOGS")]
    private static void LOG(LOG_LEVEL level, string message)
    {
        AnaisLogger.Log(LOG_CATEGORY.VECTOR2D_FIXED_POINT_SOLVER, level, message);
    }
}
