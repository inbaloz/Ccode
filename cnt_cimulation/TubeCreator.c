#include "Constants_and_libraries.h"
#include "WriteCoordinates.h"
#include "aVecLength.h"
#include "LatticeCreator.h"
#include "CutUnitcell.h"
#include "DuplicateTube.h"
#include "TubeCreator.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>


int TubeCreator(aVec T, aVec Ch, double teta, int tubeType, int tubeUnitN,
				int unitcellN, Atom** tube)
{	
	Atom* lattice;
	int latticeN, tubeN;
	Atom* tubeUnit;			// The tube's unitcell

	// The function needs the future unitcell borders (before rotating and moving)
	// to create a lattice in a suitable size:
	{
	double xMin = - aVecLength(T) * sin( (M_PI / 6) - teta);
	double xMax = aVecLength(Ch) * cos((M_PI / 6) - teta);
	double yMin = 0;
	double yMax = ( aVecLength(Ch) * sin((M_PI / 6) - teta) + aVecLength(T) * cos( (M_PI / 6) - teta));


	// Creating the lattice (to make a tube of):
	latticeN = LatticeCreator(&lattice, xMin, yMin, xMax, yMax, tubeType);
	}
	// Making the tube's unitcell from the lattice:
	tubeUnit = CutUnitcell(lattice, latticeN, Ch, T, (M_PI / 6) - teta, tubeUnitN);

	// Duplicating the unitcell:
	*tube = DuplicateTube(tubeUnit, tubeUnitN, unitcellN, aVecLength(T));
	tubeN = tubeUnitN * unitcellN;

	// cut last part of tube if requested
	//if (percentTruncated != 0.0) {
	//	CutLastPartOfTube(tube, &tubeN, percentTruncated);
	//}

	printf("number of atoms: %d\n",tubeN);

	// Free the unit tube, it isn't needed anymore:
	free(tubeUnit);

	// This lattice was used to create the tube. Now we free it for later reuse:
	free(lattice);
	return tubeN;
}


		