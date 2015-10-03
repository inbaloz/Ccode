#include "Constants_and_libraries.h"
#include "SpinningMotion.h"
#include "Rotate.h"
#include "FindInteracting.h"
#include "WriteCoordinates.h"
#include <math.h>

// input: Array for RI values, spinning step, amount of steps, the tube,
// amount of atoms in the tube, x,y shifts of the center of the tube from (0,0).

// output: Calculate all RI values and put them in the input array.
// Notes:	1. At the end the tube will be spinned by rotSpinEnd angle 
// + one Spinning step.
//			2. The function assumes the spinningStart is zero or
// the tube has been already spinned by rotSpinStart.

void SpinningMotion(double* RI, double spinningStep, double amountOfSteps,
					Atom* tube, int tubeN, Atom* surfaceLattice, int surfaceN,
					double xShift, double yShift, int latticeType, char* prefix)
{
	int i, j;
	double effectiveNum;

	// The spinning (Rotate) is in the end because when i = 0 we want
	// the start angle, so we don't want to spin by the regular step.
	for (i = 0; i < amountOfSteps; i++)
	{
		RI[i] = 0;
		for (j = 0; j < tubeN; j++)
		{
			effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (RND - ILD) );
			if (tube[i].z < MAX_HEIGHT)
			{
				RI[i] = RI[i] + effectiveNum * FindInteracting(tube[j], xShift, yShift, latticeType);
			}
		}
		
		WriteCoordinates(tube, tubeN, surfaceLattice, surfaceN, 
						  xShift, yShift, i, prefix);
		Rotate(tube, tubeN, 3, spinningStep);
	}
}
