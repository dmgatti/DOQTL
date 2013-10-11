/*******************************************************************************
 *  addlog.c
 *  Given two arguments on a log scale, add them on a non-logged scale.
 *  Arguments: x: double, value on a log scale.
 *             y: double, value on a log scale.
 *  Return: double, the value log(exp(x) + exp(y)).
 */
#include <math.h>
double addlog(double x, double y) {
  double retval = 0.0;
  if(!isfinite(x) && x < 0) {
    retval = y;
  } else if(!isfinite(y) && y < 0) { 
    retval =  x;
  } else if (x >= y) {
    retval = x + log1p(exp(y - x));
  } else {
    retval = y + log1p(exp(x - y));
  } /* else */

  return retval;
} /* addlog() */
