#ifndef __ROTATIONMOTION__
#define __ROTATIONMOTION__

#include "Constants_and_libraries.h"
void RotationMotion(double* RI, double rotationStep, int amountOfSteps,
					Atom* tube, int tubeN, Atom* surfaceLattice, int surfaceN, 
					double radius, double shiftAngle, 
					double xShift, double yShift, int latticeType,
					char* prefix);
#endif
