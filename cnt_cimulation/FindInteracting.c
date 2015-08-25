#include "Constants_and_libraries.h"
#include "FindInteracting.h"
#include "CalculateIntersection.h"
#include <math.h>
#include <stdio.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

// number of atoms in the different zones.

Atom ZONE_1[SIZE_ZONE_1];
Atom ZONE_2[SIZE_ZONE_2];
Atom ZONE_3[SIZE_ZONE_3];
Atom ZONE_4[SIZE_ZONE_4];
Atom ZONE_5[SIZE_ZONE_5];

double FindInteracting(Atom atom, double xShift, double yShift, int latticeType)
{
	int i;


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

	if (xMod >= 0 && xMod < LATTICE_BL) {
		for (i=0; i<SIZE_ZONE_1; i++) {
			RI += CalculateIntersection(atomMod, ZONE_1[i]);
		}	
	}
	else if	(xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)) {
		for (i=0; i<SIZE_ZONE_2; i++) {
			RI += CalculateIntersection(atomMod, ZONE_2[i]);
		}	
	}
	else if (xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL)) {
		if (yMod >= (0.5 * LATTICE_HIGHT)) {
			for (i=0; i<SIZE_ZONE_3; i++) {
				RI += CalculateIntersection(atomMod, ZONE_3[i]);
			}
		}
		else {
			for (i=0; i<SIZE_ZONE_4; i++) {
				RI += CalculateIntersection(atomMod, ZONE_4[i]);
			}
		}	
	}
	else { // xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
		for (i=0; i<SIZE_ZONE_5; i++) {
			RI += CalculateIntersection(atomMod, ZONE_5[i]);
		}	
	}

	return WeightInteracting(RI);

}

double WeightInteracting(double RI) {
	return RI;
}