#include "Constants_and_libraries.h"
#include "PerfectRotationMotion.h"
#include "RotateShift.h"
#include "WriteCoordinates.h"
#include "FindInteracting.h"
#include <math.h>


// input: Array for RI values, x starting value, y starting value, x step,
// y step, rotation step, amount of steps, the tube, amount of atoms in the
// tube, the angle between tube axis and y axis.

// output: Calculate all RI values and put them in the input array.
// Notes:	1. At the end the tube will be rotated by rotation end angle 
// + one rotation step.
//			2. Note that xStart and yStart aren't pointers. Therefore,
// their change is local.
//			3. The function assumes the rotationStart is zero or
// the tube has been already rotated by rotationStart.

void PerfectRotationMotion(double* RI, double xStart, double yStart, double xStep, double yStep,
							double rotationStep, double amountOfSteps, 
							Atom* tube, int tubeN, Atom* surfaceLattice, int surfaceN,
							double radius, double shiftAngle, 
							int latticeType, char* prefix)
{
	int i, j;
	double effectiveNum;
	double xShift = xStart;
	double yShift = yStart;
	
	// The rotation (RotateShift) is in the end because when i = 0 we want
	// the start angle, so we don't want to rotate by the regular step.
	// Same goes with x-yStep.
	for (i = 0; i < amountOfSteps; i++)
	{
		RI[i] = 0;
		for (j = 0; j < tubeN; j++)
		{
			effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (RND - ILD) );
			if (tube[j].z < MAX_HEIGHT)
			{
				RI[i] = RI[i] + effectiveNum * FindInteracting(tube[j], xShift, yShift, latticeType);
			}
		}


		WriteCoordinates(tube, tubeN, surfaceLattice, surfaceN, 
						  xShift, yShift, i, prefix); 

		RotateShift(tube, tubeN, rotationStep, shiftAngle, ILD + radius);
		xShift = xShift + xStep;
		yShift = yShift + yStep;


	}
}
