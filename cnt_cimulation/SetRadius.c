#include "Constants_and_libraries.h"
#include "SetRadius.h"

#include <math.h>
#include <stdio.h>

double SetRadius(int atomType)
{
	if (TUBETYPE == 0 && LATTICETYPE == 0) {     // 1. CNT on graphene
		return RCGRAPHENE_HOMO;
	} 
	else if (TUBETYPE == 1 && LATTICETYPE == 1) { // 2. BN tube on BN lattice
		if (atomType == B_type) {
			return RBTUBE_HOMO;
		} else {
			return RNTUBE_HOMO;
		} 
	}
	else  {  // 3. Heterojunctions	
		if (atomType == C_type) {
			return RCGRAPHENE_HETERO;
		} else if (atomType == B_type) {
			return RBTUBE_HETERO;
		} else {
			return RNTUBE_HETERO;
		}
	}
}
	