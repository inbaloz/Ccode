#include "Constants_and_libraries.h"
#include "WriteCoordinates.h"
#include "AtomsToFile.h"
#include <stdio.h>
#include <math.h>


void WriteCoordinates(Atom *tube, int tubeN, Atom *surfaceLattice,
					  int surfaceN, double xShift, double yShift, int step,
					  char *prefix)
{
	 if (step > 0) {
	 	return;
	 }

	int append = 0; // 0 - write, 1 - append
	char tubeFile[105];		// Tube file name
	double xMod = remainder(xShift, LATTICE_HORIZD);
	double yMod = remainder(yShift, LATTICE_HIGHT);

	// print the lattice and the tube
	sprintf(tubeFile, "%s - atoms %d", prefix, step);
	AtomsToFile(tube, tubeN, tubeFile, 0.0, 0.0, 0.0, append);
	append = 1;
	AtomsToFile(surfaceLattice, surfaceN, tubeFile, -xMod, -yMod, 0.0, append);
}


		