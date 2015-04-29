#include "Constants_and_libraries.h"
#include "EffectiveNum.h"
#include <math.h>

// input: array of atoms, the number of atoms, the minimum hight in the array,
// the maximum hight that isn't negligible (1% of effect).

// output: The effective amount of atoms (sum of atom's effectiveness level).

double EffectiveNum(Atom* array, int arrayN, double minHight, double maxHight)
{
	double result = 0;
	int i;		// counter
	double temp;

	for (i = 0; i < arrayN; i++)
	{
		// temp - current atom in tube, normalized by its hight:
		// the lowest will get 1, the closest to maxHight will get 0.
		temp = exp( EXPNORM * (minHight - array[i].z) / (maxHight - minHight) );
		if	(temp > NP)
		{
			result = result + temp;
		}
	}
	return result;
}
