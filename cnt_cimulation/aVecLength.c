#include "Constants_and_libraries.h"
#include "aVecLength.h"
#include <math.h>

// input: a vector in terms of a1 and a2 components.

// output: the length of the vector in absolute units.
// note: the function assumes she has "A" the length of
// a unit aVec (a1 or a2 length).

double aVecLength(aVec vec)
{
	return A * sqrt( (double)((vec.n * vec.n) + (vec.n * vec.m) + (vec.m * vec.m)) );
}
