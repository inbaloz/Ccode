#include "Constants_and_libraries.h"
#include "RotationMotion.h"
#include "RotateShift.h"
#include "FindInteracting.h"
#include <math.h>

// input: Array for RI values, rotation step, amount of steps, the tube,
// amount of atoms in the tube, the angle between tube axis and y axis, 
// x,y shifts of the center of the tube from (0,0).

// output: Calculate all RI values and put them in the input array.
// Notes:	1. At the end the tube will be rotated by rotation end angle 
// + one rotation step.
//			2. The function assumes the rotationStart is zero or
// the tube has been already rotated by rotationStart.

void RotationMotion(double* RI, double rotationStep, int amountOfSteps,
					Atom* tube, int tubeN, double radius, double shiftAngle, double xShift, double yShift)
{
	int i, j;
	double effectiveNum;
	// The rotation (RotateShift) is in the end because when i = 0 we want
	// the start angle, so we don't want to rotate by the regular step.
	for (i = 0; i < amountOfSteps; i++)
	{
		RI[i] = 0;
		for (j = 0; j < tubeN; j++)
		{
			effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (RND - ILD) );
			if (effectiveNum > NP)
			{
				RI[i] = RI[i] + effectiveNum * FindInteracting(tube[j], xShift, yShift);
			}
		}
		RotateShift(tube, tubeN, rotationStep, shiftAngle, ILD + radius);
	}
}
