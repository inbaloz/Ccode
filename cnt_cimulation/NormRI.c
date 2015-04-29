#include "Constants_and_libraries.h"
#include "NormRI.h"
#include <math.h>

// input: unnormalized RI data, amount of RI data, RI minimum, RI maximum.
// output: no actual output, RI is normalized.

void NormRI(double* RI, int RIN, double RIMin, double RIMax)
{
	int i;
	for (i = 0; i < RIN; i++)
	{
		RI[i] = (RI[i] - RIMin) / (RIMax - RIMin);
	}
}
