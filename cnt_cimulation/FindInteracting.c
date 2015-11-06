#include "Constants_and_libraries.h"
#include "FindInteracting.h"
#include "CalculateIntersection.h"
#include "FindInteractingNearZones.h"
#include "FindInteractingGlobalZone.h"
#include <math.h>
#include <stdio.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

// number of atoms in the different zones.

double FindInteracting(Atom atom, double xShift, double yShift, int latticeType)
{
	// Ensuring we have a non-negative value in the right range:
	double xMod = remainder(atom.x + xShift, LATTICE_HORIZD);
	double yMod = remainder(atom.y + yShift, LATTICE_HIGHT);
	double RI = 0; // It is the addition to RI before normalizing.
	
	// xMod and yMod need to be not negative:
	if (xMod < 0)
	{
		xMod = xMod + LATTICE_HORIZD;
	}
	if (yMod < 0)
	{
		yMod = yMod + LATTICE_HIGHT;
	}
	Atom atomMod = { .x = xMod, .y = yMod, .z = atom.z , .type = atom.type};

// ------------------------- calculating intersection ----------------------------------

	if (USE_GLOBAL_ZONE) {
		RI = FindInteractingGlobalZone(atomMod);
	} else {
		RI = FindInteractingNearZones(atomMod, xMod, yMod);
	}
	
	return RI;
}
