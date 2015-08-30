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
	double intersection;

	// ----- calculating the distance between the atoms ----
	double d = sqrt(pow(atomTube.x - atomLattice.x, 2) + pow(atomTube.y - atomLattice.y, 2));
	
	// ----- setting the radii according to the atoms types (C/B/N) -----
	if (atomTube.type == C_type)
	{
		r1 = RCCNT;
	}
	else if (atomTube.type == B_type)
	{
		r1 = RBTUBE;
	}
	else
	{
		r1 = RNTUBE;
	}

	if (atomLattice.type == C_type)
	{
		r2 = RCGRAPHENE;
	}
	else if (atomLattice.type == B_type)
	{
		r2 = RBLATTICE;
	}
	else // atomeLattice.type == 'N'
	{
		r2 = RNLATTICE;
	}
	
// ------------ calculate intersection -------------

	if (USE_GAUSSIAN_INTERSECTION) {
		intersection = GaussianIntersection(r1, r2, d);
	} else {
		intersection = HardSphereIntersection(r1, r2, d);
	}

	return WeightIntersection(intersection);

}

double WeightIntersection(double intersection) {
	return intersection;
}
