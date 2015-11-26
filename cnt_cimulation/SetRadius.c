#include "Constants_and_libraries.h"
#include "SetRadius.h"

#include <math.h>

double SetRadius(int atomType)
{

	if (atomType == C_type)
	{
		return RCGRAPHENE_HOMO;
	}
	else if (atomType == B_type)
	{
		return RBLATTICE_HETERO;
	}
	else
	{
		return RNLATTICE_HETERO;
	}


}