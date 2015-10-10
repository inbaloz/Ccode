#include "Constants_and_libraries.h"
#include "CalculateIntersection.h"
#include "HardSphereIntersection.h"
#include "GaussianIntersection.h"
#include <math.h>

// input: two atoms.
// output: the intersection of the two circles around them on x-y plane.
double CalculateIntersection(Atom atomTube, Atom atomLattice)
{
	double r1, r2; // r1 is atomTube's radius and r2 is atomLattice's radius.
	double radiusToSTDEV1, radiusToSTDEV2;

	double intersection;

	// ----- calculating the distance between the atoms ----
	double d = sqrt(pow(atomTube.x - atomLattice.x, 2) + pow(atomTube.y - atomLattice.y, 2));
	
	// ----- setting the radii according to the atoms types (C/B/N) -----
	if (atomTube.type == C_type)
	{
		r1 = RCCNT;
		radiusToSTDEV1 = RADIUS_TO_STDEV_C;
	}
	else if (atomTube.type == B_type)
	{
		r1 = RBTUBE;
		radiusToSTDEV1 = RADIUS_TO_STDEV_B;
	}
	else
	{
		r1 = RNTUBE;
		radiusToSTDEV1 = RADIUS_TO_STDEV_N;
	}

	if (atomLattice.type == C_type)
	{
		r2 = RCGRAPHENE;
		radiusToSTDEV2 = RADIUS_TO_STDEV_C;
	}
	else if (atomLattice.type == B_type)
	{
		r2 = RBLATTICE;
		radiusToSTDEV2 = RADIUS_TO_STDEV_B;
	}
	else // atomeLattice.type == 'N'
	{
		r2 = RNLATTICE;
		radiusToSTDEV2 = RADIUS_TO_STDEV_N;
	}
	
// ------------ calculate intersection -------------

	if (USE_GAUSSIAN_INTERSECTION) {
		intersection = GaussianIntersection(r1, r2, d, radiusToSTDEV1, radiusToSTDEV2);
	} else {
		intersection = HardSphereIntersection(r1, r2, d);
	}

	return WeightIntersection(intersection);

}

double WeightIntersection(double intersection) {
	return intersection;
}
