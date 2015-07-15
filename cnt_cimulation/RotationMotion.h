#ifndef __ROTATIONMOTION__
#define __ROTATIONMOTION__

#include "Constants_and_libraries.h"
void RotationMotion(double* RI, double rotationStep, int amountOfSteps,
					Atom* tube, int tubeN, double radius, double shiftAngle, 
					double xShift, double yShift, int latticeType);
#endif
