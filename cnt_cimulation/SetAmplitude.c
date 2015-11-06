#include "Constants_and_libraries.h"
#include "SetAmplitude.h"

#include <math.h>

double SetAmplitude(int atomType)
{

	if (atomType == C_type)
	{
		return GAUSSIAN_AMPLITUDE_C;
	}
	else if (atomType == B_type)
	{
		return GAUSSIAN_AMPLITUDE_B;
	}
	else
	{
		return GAUSSIAN_AMPLITUDE_N;
	}
}