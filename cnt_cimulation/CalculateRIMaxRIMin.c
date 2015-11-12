#include "Constants_and_libraries.h"
#include "Rotate.h"
#include "FindInteracting.h"
#include "WriteCoordinates.h"
#include "Move.h"
#include <math.h>
#include <stdio.h>

// This function creates the surface and then moves it so the tube 

void CalculateRIMaxRIMin(double* RIMax, double* RIMin, Atom *surfaceLattice, int surfaceN, 
						 Atom* tube, int tubeN, double teta,
					 	 int tubeType, int latticeType) 
{	
	double tempLatticeBL;
	double xShiftRIMax = 0.0;
    double xShiftRIMin = 0.0;
    double yShiftRIMax = 0.0;
    double yShiftRIMin = 0.0;
    double rotationAngleRIMin = 0.0;
    int	   j = 0;
    double effectiveNum = 0;

    *RIMax = 0.0;
    *RIMin = 0.0;

    tempLatticeBL = LATTICE_BL;
    LATTICE_BL = TUBE_BL;

	// ---- (I) Placing the tube for easy reach of desired configurations ---------

    Rotate(tube, tubeN, 3, (M_PI/6) - teta); // rotation around the z axis (spinning)

    // ---- (II) Setting the correct x and yshifts of desired configurations: ------
    // 1. CNT on graphene: maximal at AA, minimal at AB.
    // 2. BNNT on h-BN	 : maximal at AA. minimal at AA' (All N on B)
    // 3. Heterojunctions: maximal at AA (all eclipsed), minimal when half
    // 	  of the C atoms are atop B atoms and the rest are above hexagon centers
    //	  of the h-BN lattice.
    

    if (tubeType == 0 && latticeType == 0) {      // 1. CNT on graphene
		//RIMax
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

		//RIMin
    	xShiftRIMin = CNT_BL_HOMO;
    	yShiftRIMin = 0.0;
	}
	else if (tubeType == 1 && latticeType == 1) { // 2. BN tube on BN lattice
		//RIMax
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

		//RIMin
    	rotationAngleRIMin = M_PI/3;
    	xShiftRIMin = BN_LATTICE_BL_HOMO;
    	yShiftRIMin = 0.0;
    	printf("BL lattice: %lf\n, BL tube: %lf\n", LATTICE_BL, TUBE_BL);
	}
	else if (tubeType == 0 && latticeType ==1) {  // 3. Heterojunctions: CNT on BN lattice	
		//RIMax
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

		//RIMin
    	xShiftRIMin = -BN_LATTICE_BL_HETERO;
    	yShiftRIMin = 0.0;
	}
	else {  									  // 3. Heterojunctions: BNNT on graphene
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

    	xShiftRIMin = BN_LATTICE_BL_HETERO;
    	yShiftRIMin = 0.0;
	}									  
	

    // ------------ (III) Calculating RIMax and RIMin ----------------------

	// Calculating the RI max
	for (j = 0; j < tubeN; j++)
	{
		effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (RND - ILD) ); // find weight of atom
		if (tube[j].z < MAX_HEIGHT)
		{
			*RIMax = *RIMax + effectiveNum * FindInteracting(tube[j], xShiftRIMax, xShiftRIMax,
														   latticeType);
		}
	}

	WriteCoordinates(tube, tubeN, surfaceLattice, surfaceN, 
						  xShiftRIMax,yShiftRIMax, -2, "normalizedRIMax");
	// Calculating the RI min
	Rotate(tube, tubeN, 3, rotationAngleRIMin);
	WriteCoordinates(tube, tubeN, surfaceLattice, surfaceN, 
						  xShiftRIMin,yShiftRIMin, -1, "normalizedRIMin");

	for (j = 0; j < tubeN; j++)
	{
		effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (RND - ILD) ); // find weight of atom
		if (tube[j].z < MAX_HEIGHT)
		{
			*RIMin = *RIMin + effectiveNum * FindInteracting(tube[j], xShiftRIMin, yShiftRIMin,
														   latticeType);
		}
	}
	printf("RIMin:%lf\n", *RIMin);
	printf("RIMax:%lf\n", *RIMax);
	// ------- (IV) Returning the tube to the original location -----------------
	Rotate(tube, tubeN, 3, -((M_PI/6) - teta + rotationAngleRIMin));

	LATTICE_BL = tempLatticeBL;
}
