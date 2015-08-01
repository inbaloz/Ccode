#ifndef __SPINNINGMOTION__
#define __SPINNINGMOTION__

#include "Constants_and_libraries.h"

void SpinningMotion(double* RI, double spinningStep, double amountOfSteps,
					Atom* tube, int tubeN, Atom* surfaceLattice, int surfaceN,
					double xShift, double yShift, int latticeType, char* prefix);

#endif
