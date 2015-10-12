#include "Constants_and_libraries.h"
#include "Rotate.h"
#include "FindInteracting.h"
#include "WriteCoordinates.h"
#include <math.h>
#include <stdio.h>

// This function creates the surface and then moves it so the tube 

void CalculateRIMaxRIMin(double* RIMax, double* RIMin, Atom *surfaceLattice, int surfaceN, 
						 Atom* tube, int tubeN, double teta,
					 	 int tubeType, int latticeType) 
{
	double xShiftRIMax = 0.0;
    double xShiftRIMin = 0.0;
    double yShiftRIMax = 0.0;
    double yShiftRIMin = 0.0;
    double rotationAngleRIMin = 0.0;
    int	   j = 0;
    double effectiveNum = 0;
    double radius_tube_0, radius_tube_1, radius_surface_0, radius_surface_1;
    double halfDifference;

    *RIMax = 0.0;
    *RIMin = 0.0;


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

    	// radii
    	radius_tube_0 = RCCNT_HOMO;
    	radius_tube_1 = RCCNT_HOMO;
    	radius_surface_0 = RCGRAPHENE_HOMO;
    	radius_surface_1 = RCGRAPHENE_HOMO;
	}
	else if (tubeType == 1 && latticeType == 1) { // 2. BN tube on BN lattice
		//RIMax
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

		//RIMin
    	rotationAngleRIMin = M_PI/3;
    	xShiftRIMin = -BN_LATTICE_BL_HOMO;
    	yShiftRIMin = 0.0;

    	// radii
    	radius_tube_0 = RBTUBE_HOMO;
    	radius_tube_1 = RNTUBE_HOMO;
    	radius_surface_0 = RBLATTICE_HOMO;
    	radius_surface_1 = RNLATTICE_HOMO;
	}
	else if (tubeType == 0 && latticeType ==1) {  // 3. Heterojunctions: CNT on BN lattice
		if (CNT_BL_HETERO == BN_LATTICE_BL_HETERO) {
			//RIMax
			xShiftRIMax = 0.0;
			yShiftRIMax = 0.0;

			//RIMin
    		xShiftRIMin = -BN_LATTICE_BL_HETERO;
    		yShiftRIMin = 0.0;
		} else {
			halfDifference = MAX(CNT_BL_HETERO, BN_LATTICE_BL_HETERO)/2;
			//RIMax
			xShiftRIMax = -halfDifference;
			yShiftRIMax = 0.0;

			//RIMin
    		xShiftRIMin = - (BN_LATTICE_BL_HETERO + halfDifference);
    		yShiftRIMin = 0.0;
			
		}

    	// radii
    	radius_tube_0 = RCCNT_HETERO;
    	radius_tube_1 = RCCNT_HETERO;
    	radius_surface_0 = RBTUBE_HETERO;
    	radius_surface_1 = RNTUBE_HETERO;
	}
	else {  									  // 3. Heterojunctions: BNNT on graphene
		
		if (BN_TUBE_BL_HETERO == GRAPHENE_BL_HETERO) {
			
			halfDifference = MAX(BN_TUBE_BL_HETERO, GRAPHENE_BL_HETERO)/2;
			//RIMax
			xShiftRIMax = halfDifference;
			yShiftRIMax = 0.0;

			//RIMin
    		xShiftRIMin = (BN_LATTICE_BL_HETERO + halfDifference);
    		yShiftRIMin = 0.0;
		} else {

		}

    	// radii
    	radius_tube_0 = RBTUBE_HETERO;
    	radius_tube_1 = RNTUBE_HETERO;
    	radius_surface_0 = RBLATTICE_HETERO;
    	radius_surface_1 = RNLATTICE_HETERO;
	}									  
	

    // ------------ (III) Calculating RIMax and RIMin ----------------------

	// Calculating the RI max
	for (j = 0; j < tubeN; j++)
	{
		effectiveNum = exp( EXPNORM * (ILD - tube[j].z) / (RND - ILD) ); // find weight of atom
		if (tube[j].z < MAX_HEIGHT)
		{
			*RIMax = *RIMax + effectiveNum * FindInteracting(tube[j], xShiftRIMax, yShiftRIMax,
														   latticeType);

		}
	}

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
	// ------- (IV) Returning the tube to the original location -----------------
	Rotate(tube, tubeN, 3, -((M_PI/6) - teta) + rotationAngleRIMin);
}
