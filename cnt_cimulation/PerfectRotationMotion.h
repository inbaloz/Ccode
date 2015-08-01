#ifndef __PERFECTROTATIONMOTION__
#define __PERFECTROTATIONMOTION__

#include "Constants_and_libraries.h"
void PerfectRotationMotion(double* RI, double xStart, double yStart, double xStep, double yStep,
							double rotationStep, double amountOfSteps, 
							Atom* tube, int tubeN, Atom* surfaceLattice, int surfaceN,
							double radius, double shiftAngle, 
							int latticeType, char* prefix);

#endif
