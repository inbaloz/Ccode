#include "Constants_and_libraries.h"
#include "LatticeCreator.h"
#include "CreateSurface.h"
#include "Move.h"
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This function creates the surface and then moves it so the tube 

int CreateSurface(Atom** surfaceLattice, double length, double radius, int latticeType) 

{

	double xMin = 0.0, yMin = 0.0;
	double xMax, yMax;
	double xMod, yMod;
	int surfaceN;
	double squareEdge;

	xMax = 2 * radius + 2 * (INTERACTION_BUFFER + LATTICE_HORIZD);
	yMax = 1 * length + 2 * (INTERACTION_BUFFER + LATTICE_HIGHT);
	squareEdge = MAX(xMax, yMax);


	surfaceN = LatticeCreator(surfaceLattice, xMin, yMin, squareEdge, squareEdge, latticeType);
	xMod = remainder(-squareEdge/2, LATTICE_HORIZD);
	yMod = remainder(-squareEdge/2, LATTICE_HIGHT);
	Move(*surfaceLattice, surfaceN, -squareEdge/2 - xMod, -squareEdge/2 - yMod, 0);

	return surfaceN;

}
