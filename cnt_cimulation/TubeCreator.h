#ifndef __TubeCreator__
#define __TubeCreator__

#include "Constants_and_libraries.h"
#include "WriteCoordinates.h"
#include "aVecLength.h"
#include "LatticeCreator.h"
#include "CutUnitcell.h"
#include "DuplicateTube.h"


int TubeCreator(aVec T, aVec Ch, double teta, int tubeType, int tubeUnitN,
				 int unitcellN, Atom** tube);

#endif
		