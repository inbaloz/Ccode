#include "Constants_and_libraries.h"
#include "CalculateIntersection.h"
#include <math.h>

// input: two atoms.
// output: the intersection of the two circles around them on x-y plane.
double CalculateIntersection(Atom atomTube, Atom atomLattice)
{
	double d = sqrt(pow(atomTube.x - atomLattice.x, 2) + pow(atomTube.y - atomLattice.y, 2));
	double r1 = RCCNT;
	double r2 = RCGRAPHENE;
	double intersection;
	// Too far:
	if (d > r1 + r2)
	{
		return 0;
	}
	// Very close:
	else if (d <= ABS(r1 - r2))
	{
		return MIN(M_PI * pow(r1,2), M_PI * pow(r2,2));
	}
	// In between:
	else
	{
		// Sum of to sectors of circles minus the surface of kite created by them:
		intersection =	pow(r1,2) * acos((pow(d,2) + pow(r1,2) - pow(r2,2)) / (2 * d * r1)) + 
						pow(r2,2) * acos((pow(d,2) + pow(r2,2) - pow(r1,2)) / (2 * d * r2)) - 
						// Heron's formula for triangle surface multiplied by 2:
						0.5 * sqrt( (r1 + r2 - d) * (r1 - r2 + d) * (- r1 + r2 + d) * (r1 + r2 + d) );
		return intersection;
	}

}
