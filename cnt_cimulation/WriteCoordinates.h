#ifndef __WriteCoordinates__
#define __WriteCoordinates__

#include "Constants_and_libraries.h"
void WriteCoordinates(Atom *tube, int tubeN, Atom *surfaceLattice,
					  int surfaceN, double xShift, double yShift, int step,
					  char *prefix);


#endif