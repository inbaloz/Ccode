#include "Constants_and_libraries.h"
#include "CutUnitcell.h"
#include "Rotate.h"
#include "Move.h"
#include "aVecLength.h"
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "ArrayToFile.h"


// input: a lattice that needs to be cut, the size of it, chirality vector,
// translational vector, the angle between x axis and the chirality vector,
// the expected amount of atoms after cutting.
// note: the input lattice remains unchanged.

// output: a new atom array that includes the tube's unitcell centered
// and parallel to y axis.

// How this function works:
// 1.	We rotate the atoms to cause the chirality vector, to become
// parallel to the x-axis and therefor, translational vector parallel
// to the y-axis.
// 2.	Now, we have a comftrable rectangle made of the Ch and T vectors.
// We can check by measuring the distance from the middle which atoms
// are required and which aren't.
// 3.	Moving the cut lattice so the it's center will be (0,0) instead of
// (xMax/2, yMax/2).
// 4.	Making a tube from the cut lattice. Note: the center of the
// tube is (xMax/2, yMax/2) and not (0,0).
// 5.	We rotate everything back to the normal coordiante
// system. Note: we keep our tube centered (we rotate him back without
// moving it).

Atom* CutUnitcell(Atom* lattice, int latticeN, aVec Ch, aVec T, double shiftAngle, int cutLatticeN)
{
	int count = 0;
	int i;					// Counter for loop
	int missingAtomsN;		// Number of missing atoms (part 3)
	double pWorst;			// Potentially worst atom that fits (part 3)
	double xMax, yMax;		// Rectangle borders (part 2)
	double* distMissing;	// Missing atoms distances for the rectangle borders (part 3)
	double epsilon = 1e-5;	// Disntace of equality
	Atom* cutLattice = malloc(cutLatticeN * sizeof(Atom));
	Atom* missingAtoms;		// Missing atoms array (part 3)

//********************** PART 1 ******************************
	Rotate(lattice, latticeN, 3, -shiftAngle);	// Rotating the input lattice.
	// Note: the output atoms haven't been located yet, so we don't need to
	// rotate them.
//********************** PART 2 ******************************
	xMax = aVecLength(Ch);
	yMax = aVecLength(T);

	// Checking for atoms in the rectangle.
	for ( i = 0; i < latticeN; i++ )
	{
			// We want to include only the left border without the right border
			// to prevent duplication of atoms. Same thing goes for bottom atoms
			// without top atoms.
		if ((lattice[i].x >= 0 - epsilon) && (lattice[i].x <= xMax - epsilon) && 
			(lattice[i].y >= 0 - epsilon) && (lattice[i].y <= yMax - epsilon)	)
		{
			if (count < latticeN)
			{
				cutLattice[count].x = lattice[i].x;
				cutLattice[count].y = lattice[i].y;
				cutLattice[count].z = lattice[i].z;
				cutLattice[count].type = lattice[i].type;
				count++;
			}
			else
			{
				printf("An error has occured, too many atoms in the unitcell.\n");
			}
		}
	}



//********************** PART 3 ******************************
	Move(cutLattice, cutLatticeN, -(xMax/2), -(yMax/2), 0);

//********************** PART 4 ******************************
	// Note: we only need to change the x and z values. y remains
	// unchanged.
	// The tube's radius is xMax/(2*M_PI).
	for ( i = 0; i < cutLatticeN; i++ )
	{
		cutLattice[i].z = ILD + (xMax/(2*M_PI)) * ( 1 - cos(2 * M_PI * cutLattice[i].x /xMax) );
		cutLattice[i].x = (xMax/(2*M_PI)) * sin(2 * M_PI * cutLattice[i].x /xMax);
	}

//********************** PART 5 ******************************
	// We don't rotate back the tube because we want y axis to be
	// parallel to the tube's axis.
	Rotate(lattice, latticeN, 3, shiftAngle);		// Rotating back the input lattice.
	return cutLattice;
}
