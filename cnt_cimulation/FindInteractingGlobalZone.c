#include "Constants_and_libraries.h"
#include "FindInteractingNearZones.h"
#include "CalculateIntersection.h"
#include <math.h>
#include <stdio.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

// number of atoms in the different zones.

Atom GlobalZone[SIZE_GLOBAL_ZONE];

double FindInteractingGlobalZone(Atom atomMod)
{
	int i;
	double RI = 0;

	for (i=0; i<SIZE_GLOBAL_ZONE; i++) {
		RI += CalculateIntersection(atomMod, GlobalZone[i]);
	}

	return RI;
}