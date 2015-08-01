#ifndef __SLIDINGMOTION__
#define __SLIDINGMOTION__

#include "Constants_and_libraries.h"

void SlidingMotion(double* RI, double xStep, double yStep, double amountOfSteps,
				   Atom* tube, int tubeN, Atom* surfaceLattice, int surfaceN, double xShift, 
				   double yShift, int latticeType,
				   char* prefix);

#endif
