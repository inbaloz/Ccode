#include "Constants_and_libraries.h"
#include "CalculateIntersection.h"
#include "HardSphereIntersection.h"
#include "GaussianIntersection.h"
#include "SetRadius.h"
#include "SetAmplitude.h"
#include <math.h>
#include <stdio.h>

// input: two atoms.
// output: the intersection of the two circles around them on x-y plane.
double CalculateIntersection(Atom atomTube, Atom atomLattice)
{
	double r1, r2; // r1 is atomTube's radius and r2 is atomLattice's radius.
	double gaussianAmplitude1, gaussianAmplitude2;

	double intersection;

	// ----- calculating the distance between the atoms ----
	double d = sqrt(pow(atomTube.x - atomLattice.x, 2) + pow(atomTube.y - atomLattice.y, 2));
	
	// ----- setting the radii according to the atoms types (C/B/N) -----
	r1 = SetRadius(atomTube.type);
	r2 = SetRadius(atomLattice.type);
	gaussianAmplitude1 = SetAmplitude(atomTube.type);
	gaussianAmplitude2 = SetAmplitude(atomLattice.type);
	
// ------------ calculate intersection -------------

	if (USE_GAUSSIAN_INTERSECTION) {
		intersection = GaussianIntersection(r1, r2, d, gaussianAmplitude1, gaussianAmplitude2);
	} else {
		intersection = HardSphereIntersection(r1, r2, d);
	}


	if ((atomTube.type == N_type && atomLattice.type == B_type) || (atomTube.type == B_type && atomLattice.type == N_type))
	{
		intersection = - intersection;
	}

	return intersection;
}
