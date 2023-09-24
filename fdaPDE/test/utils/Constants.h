#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <limits>

// a set of usefull constants
namespace fdaPDE{
namespace testing{

  // the treshold under which two doubles are considered equal (in the testing framework)
  // ** DO NOT CHANGE THIS TOLERANCE **
  // lower tolerances might not always work for all systems and compilers, higher tolerances are not acceptable.
  // If some test fail then is likely that the implementation is wrong and less likely to have the test to fail 
  // because of some external factors (executing machine, compiler, rounding errors, possible bad dependencies
  // between tests, small numerical instabilities of internal solvers, ...)
  const double DOUBLE_TOLERANCE = std::pow(0.1, 10);
  const double MACHINE_EPSILON  = std::numeric_limits<double>::epsilon();     // approx 2.22*10^-16

  // hardcoded value of pi
  constexpr double pi = 3.14159265358979323846;
}}

#endif // __CONSTANTS_H__
