#include "Constants_and_libraries.h"
#include "SlidingMotion.h"
#include "FindInteracting.h"
#include "Move.h"
#include <math.h>
// input:

void CalculateRIMinMax(double RIMin, double amountOfSteps, Atom* tube, int tubeN, 
	double xShift, double yShift, int latticeType)
{
	int i, j;
	double effectiveNum;
	double currentInteracting;
	// The sliding (Move) is in the end because when i = 0 we want
	// the start x,y, so we don't want to slide by the regular step.
	for (i = 0; i < amountOfSteps; i++)
	{
		RIMin = 0;
		for (j = 0; j < tubeN; j++)
		{
			effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (MAX_HEIGHT - ILD) ); // find weight of atom
			if (effectiveNum > NP)
			{
				currentInteracting = FindInteracting(tube[j], xShift, yShift, latticeType);
				RIMin = RIMin + effectiveNum * currentInteracting;
			}
		}

	}
}
