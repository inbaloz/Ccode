#include <math.h>
#include "Constants_and_libraries.h"
#include "SetGlobals.h"
#include "CreateGaussianZone.h"
#include "CreateHardSphereZone.h"
#include <stdio.h>

double RCCNT;
double RCGRAPHENE;
double LATTICE_BL;
double TUBE_BL;
double ILD;

double RBTUBE;
double RBLATTICE;
double RNTUBE;
double RNLATTICE;

void SetGlobals(int tubeType, int latticeType)
{
	if (tubeType == 0 & latticeType == 0) {
		RCCNT = RCCNT_HOMO;
		RCGRAPHENE = RCGRAPHENE_HOMO;
		ILD        = CNT_G_ILD;
		TUBE_BL    = CNT_BL_HOMO;
		LATTICE_BL = GRAPHENE_BL_HOMO;
	}
	else if (tubeType == 0 & latticeType == 1) {
		RCCNT      = RCCNT_HETERO;
		RBLATTICE  = RBLATTICE_HETERO;
		RNLATTICE  = RNLATTICE_HETERO;
		ILD 	   = BNT_G_ILD;
		TUBE_BL    = CNT_BL_HETERO;
		LATTICE_BL = BN_LATTICE_BL_HETERO;
	}
	else if (tubeType == 1 & latticeType == 0) {
		RNTUBE     = RNTUBE_HETERO;
		RBTUBE     = RBTUBE_HETERO;
		RCGRAPHENE = RCCNT_HETERO;
		ILD        = BNT_G_ILD;
		TUBE_BL    = BN_TUBE_BL_HETERO;
		LATTICE_BL = GRAPHENE_BL_HETERO;
	}
	else {
		RNTUBE     = RNTUBE_HOMO;
		RBTUBE     = RBTUBE_HOMO;
		RBLATTICE  = RBLATTICE_HOMO;
		RNLATTICE  = RNLATTICE_HOMO;
		ILD        = BNT_BNL_ILD;
		TUBE_BL    = BN_TUBE_BL_HOMO;
		LATTICE_BL = BN_LATTICE_BL_HOMO;
	}

	if (USE_GAUSSIAN_INTERSECTION) {
		CreateGaussianZone(latticeType);
	}
	else {
		CreateHardSphereZone(latticeType);
	}
}

