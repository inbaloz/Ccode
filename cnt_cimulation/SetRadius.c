#include "Constants_and_libraries.h"
#include "SetRadius.h"

#include <math.h>

double SetRadius(int atomType)
{

	if (atomType == C_type)
	{
		return RCCNT;
	}
	else if (atomType == B_type)
	{
		return RBTUBE;
	}
	else
	{
		return RNTUBE;
	}
}