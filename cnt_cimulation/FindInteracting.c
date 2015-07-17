#include "Constants_and_libraries.h"
#include "FindInteracting.h"
#include "CalculateIntersection.h"
#include <math.h>
#include <stdio.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

// number of atoms in the different zones.
int SIZE_ZONE_1  = 12; 
int SIZE_ZONE_2  = 12;
int SIZE_ZONE_3_4  = 14;
int SIZE_ZONE_5  = 17;

// ----------------- Creating the 5 zones for h-BN ---------------------------------------------

// Zone 1:  xMod >= 0 && xMod < LATTICE_BL 
Atom ZONE_BN_1[12] = {
	{ .x = 0    * BN_LATTICE_BL, .y = 0    * BN_HEIGHT, .z = 0 , .type = B_type},  // 1
	{ .x = -0.5 * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 2
	{ .x = 0    * BN_LATTICE_BL, .y = 1    * BN_HEIGHT, .z = 0 , .type = B_type},  // 3
	{ .x = 1    * BN_LATTICE_BL, .y = 1    * BN_HEIGHT, .z = 0 , .type = N_type}, // 4
	{ .x = 1.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type},  // 5
	{ .x = 1    * BN_LATTICE_BL, .y = 0    * BN_HEIGHT, .z = 0 , .type = N_type}, // 6
	{ .x = -0.5 * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = N_type}, // 7
	{ .x = -1.5 * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type},  // 8
	{ .x = -0.5 * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 9
	{ .x = 1.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = B_type},  // 10
	{ .x = 2.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 11
	{ .x = 1.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = B_type}   // 12
};

// Zone 2: xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)
Atom ZONE_BN_2[12] = {
	{ .x = 0.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 1
	{ .x = 0.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 2
	{ .x = 1.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = B_type}, // 3
	{ .x = 1.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 4
	{ .x = 1.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 5
	{ .x = 1.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 6
	{ .x = 1.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 7
	{ .x = 2.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = N_type}, // 8
	{ .x = 3.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 9
	{ .x = 2.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 10
	{ .x = 3.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 11
	{ .x = 2.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = N_type}  // 12
};

// Zones 3 and 4: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod > or < (0.5 * LATTICE_HIGHT)
Atom ZONE_BN_3[14] = {
	{ .x = 0.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 1
	{ .x = 0.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 2
	{ .x = 1.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = B_type}, // 3
	{ .x = 1.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 4
	{ .x = 1.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 5
	{ .x = 1.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 6
	{ .x = 1.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 7
	{ .x = 2.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = N_type}, // 8
	{ .x = 3.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 9
	{ .x = 2.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 10
	{ .x = 3.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 11
	{ .x = 2.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 12
	{ .x = 4.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 13
	{ .x = 4.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}  // 14
};

Atom ZONE_BN_4[14] = {
	{ .x = 0.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 1
	{ .x = 0.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 2
	{ .x = 1.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = B_type}, // 3
	{ .x = 1.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 4
	{ .x = 1.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 5
	{ .x = 1.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 6
	{ .x = 1.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 7
	{ .x = 2.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = N_type}, // 8
	{ .x = 3.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 9
	{ .x = 2.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 10
	{ .x = 3.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 11
	{ .x = 2.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 12
	{ .x = 4.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 13
	{ .x = 4.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}  // 14
};

// Zones 5: xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
Atom ZONE_BN_5[17] = {
	{ .x = 1.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = B_type}, // 1
	{ .x = 1.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 2
	{ .x = 1.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 3
	{ .x = 1.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 4
	{ .x = 1.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 5 
	{ .x = 3.0  * BN_LATTICE_BL, .y = -1.0 * BN_HEIGHT, .z = 0 , .type = B_type}, // 6
	{ .x = 2.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = N_type}, // 7
	{ .x = 3.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 8
	{ .x = 2.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 9
	{ .x = 3.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 10
	{ .x = 2.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = N_type}, // 11
	{ .x = 3.0  * BN_LATTICE_BL, .y = 2.0  * BN_HEIGHT, .z = 0 , .type = B_type}, // 12
	{ .x = 4.5  * BN_LATTICE_BL, .y = -0.5 * BN_HEIGHT, .z = 0 , .type = B_type}, // 13
	{ .x = 4.0  * BN_LATTICE_BL, .y = 0.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 14
	{ .x = 4.5  * BN_LATTICE_BL, .y = 0.5  * BN_HEIGHT, .z = 0 , .type = B_type}, // 15
	{ .x = 4.0  * BN_LATTICE_BL, .y = 1.0  * BN_HEIGHT, .z = 0 , .type = N_type}, // 16
	{ .x = 4.5  * BN_LATTICE_BL, .y = 1.5  * BN_HEIGHT, .z = 0 , .type = B_type}  // 17
};

// ----------------- Creating the 5 zones for graphene ---------------------------------------------

// Zone 1:  xMod >= 0 && xMod < LATTICE_BL 
Atom ZONE_G_1[12] = {
	{ .x = 0    * GRAPHENE_BL, .y = 0    * G_HEIGHT, .z = 0 , .type = C_type},  // 1
	{ .x = -0.5 * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 2
	{ .x = 0    * GRAPHENE_BL, .y = 1    * G_HEIGHT, .z = 0 , .type = C_type},  // 3
	{ .x = 1    * GRAPHENE_BL, .y = 1    * G_HEIGHT, .z = 0 , .type = C_type}, // 4
	{ .x = 1.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type},  // 5
	{ .x = 1    * GRAPHENE_BL, .y = 0    * G_HEIGHT, .z = 0 , .type = C_type}, // 6
	{ .x = -0.5 * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 7
	{ .x = -1.5 * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type},  // 8
	{ .x = -0.5 * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 9
	{ .x = 1.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type},  // 10
	{ .x = 2.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 11
	{ .x = 1.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}   // 12
};

// Zone 2: xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)
Atom ZONE_G_2[12] = {
	{ .x = 0.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 1
	{ .x = 0.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 2
	{ .x = 1.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 3
	{ .x = 1.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 4
	{ .x = 1.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 5
	{ .x = 1.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 6
    { .x = 1.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 7
	{ .x = 2.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 8
	{ .x = 3.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 9
	{ .x = 2.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 10
	{ .x = 3.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 11
	{ .x = 2.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}  // 12
};

// Zone 3: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod > (0.5 * LATTICE_HIGHT)
Atom ZONE_G_3[14] = {
	{ .x = 0.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 1
	{ .x = 0.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 2
	{ .x = 1.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 3
	{ .x = 1.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 4
	{ .x = 1.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 5
	{ .x = 1.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 6
	{ .x = 1.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 7
	{ .x = 2.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 8
	{ .x = 3.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 9
	{ .x = 2.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 10
	{ .x = 3.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 11
	{ .x = 2.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 12
	{ .x = 4.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 13
	{ .x = 4.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}  // 14
};

// Zone 4: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod < (0.5 * LATTICE_HIGHT)
Atom ZONE_G_4[14] = {
	{ .x = 0.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 1
	{ .x = 0.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 2
	{ .x = 1.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 3
	{ .x = 1.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 4
	{ .x = 1.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 5
	{ .x = 1.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 6
	{ .x = 1.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 7
	{ .x = 2.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 8
	{ .x = 3.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 9
	{ .x = 2.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 10
	{ .x = 3.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 11
	{ .x = 2.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 12
	{ .x = 4.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 13
	{ .x = 4.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}  // 14
};

// Zones 5: xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
Atom ZONE_G_5[17] = {
	{ .x = 1.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 1
	{ .x = 1.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 2
	{ .x = 1.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 3
	{ .x = 1.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 4
	{ .x = 1.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 5 
	{ .x = 3.0  * GRAPHENE_BL, .y = -1.0 * G_HEIGHT, .z = 0 , .type = C_type}, // 6
	{ .x = 2.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 7
	{ .x = 3.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 8
	{ .x = 2.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 9
	{ .x = 3.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 10
	{ .x = 2.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 11
	{ .x = 3.0  * GRAPHENE_BL, .y = 2.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 12
	{ .x = 4.5  * GRAPHENE_BL, .y = -0.5 * G_HEIGHT, .z = 0 , .type = C_type}, // 13
	{ .x = 4.0  * GRAPHENE_BL, .y = 0.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 14
	{ .x = 4.5  * GRAPHENE_BL, .y = 0.5  * G_HEIGHT, .z = 0 , .type = C_type}, // 15
	{ .x = 4.0  * GRAPHENE_BL, .y = 1.0  * G_HEIGHT, .z = 0 , .type = C_type}, // 16
	{ .x = 4.5  * GRAPHENE_BL, .y = 1.5  * G_HEIGHT, .z = 0 , .type = C_type}  // 17
};


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

	if (latticeType == 1) {
		if (xMod >= 0 && xMod < LATTICE_BL) {
			for (i=0; i<SIZE_ZONE_1; i++) {
				RI += CalculateIntersection(atomMod, ZONE_BN_1[i]);
			}	
		}
		else if	(xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL))
		{
			for (i=0; i<SIZE_ZONE_2; i++) {
				RI += CalculateIntersection(atomMod, ZONE_BN_2[i]);
			}	
		}
		else if (xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL)) {
			if (yMod >= (0.5 * LATTICE_HIGHT)) {
				for (i=0; i<SIZE_ZONE_3_4; i++) {
					RI += CalculateIntersection(atomMod, ZONE_BN_3[i]);
				}
			}
			else {
				for (i=0; i<SIZE_ZONE_3_4; i++) {
					RI += CalculateIntersection(atomMod, ZONE_BN_4[i]);
				}
			}	
		}
		else { // xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
			for (i=0; i<SIZE_ZONE_5; i++) {
				RI += CalculateIntersection(atomMod, ZONE_BN_5[i]);
			}	
		}
	}
	else if (latticeType == 0) {
		if (xMod >= 0 && xMod < LATTICE_BL) {
			for (i=0; i<SIZE_ZONE_1; i++) {
				RI += CalculateIntersection(atomMod, ZONE_G_1[i]);
			}	
		}
		else if	(xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL))
		{
			for (i=0; i<SIZE_ZONE_2; i++) {
				RI += CalculateIntersection(atomMod, ZONE_G_2[i]);
			}	
		}
		else if (xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL)) {
			if (yMod >= (0.5 * LATTICE_HIGHT)) {
				for (i=0; i<SIZE_ZONE_3_4; i++) {
					RI += CalculateIntersection(atomMod, ZONE_G_3[i]);
				}
			}
			else {
				for (i=0; i<SIZE_ZONE_3_4; i++) {
					RI += CalculateIntersection(atomMod, ZONE_G_4[i]);
				}
			}	
		}
		else { // xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
			for (i=0; i<SIZE_ZONE_5; i++) {
				RI += CalculateIntersection(atomMod, ZONE_G_5[i]);
			}	
		}

	}

	
	return RI;

}
