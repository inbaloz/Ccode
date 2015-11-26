#include "Constants_and_libraries.h"
#include "Rotate.h"
#include "FindInteracting.h"
#include "WriteCoordinates.h"
#include "CreateSurface.h"
#include "TubeCreator.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

// This function creates the surface and then moves it so the tube 

void CalculateRIMaxRIMin(double* RIMax, double* RIMin, Atom *surfaceLattice, int surfaceN, 
						 Atom* tube, int tubeN, double teta,
					 	 int tubeType, int latticeType, aVec T, aVec Ch, int tubeUnitN,
						 int unitcellN, double radius) 
{
	Atom** normaliztionTube; 		// a pointer to the right tube we're going to use for the normaliztion
									// either the hetro tube or the given tube (regular tube)
	int* normalizationTubeN;
	Atom** normalizationSurfaceLattice;	// The lattice that will be used to create the surface.
	int* normalizationSurfaceN;

	Atom* heteroTube;			// tube specificaly for hetero junctions
	int heteroTubeN;
	Atom* heteroSurfaceLattice;	// The lattice that will be used to create the surface.
	int heteroSurfaceN;

	double tempLatticeBL, tempTubeBL;
	double xShiftRIMax = 0.0, xShiftRIMin = 0.0, yShiftRIMax = 0.0, yShiftRIMin = 0.0;
    double rotationAngleRIMin = 0.0;
    int	   j = 0;
    double effectiveNum = 0;
    double totalEffectiveNum = 0; // Variable intended to capture the weight of the entire tube for energy normaliation per atom.

    *RIMax = 0.0;
    *RIMin = 0.0;

	

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

	} else if (tubeType == 1 && latticeType == 1) { // 2. BN tube on BN lattice
		//RIMax
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

		//RIMin
    	rotationAngleRIMin = M_PI/3;
    	xShiftRIMin = BN_LATTICE_BL_HOMO;
    	yShiftRIMin = 0.0;
    	printf("BL lattice: %lf\n, BL tube: %lf\n", LATTICE_BL, TUBE_BL);

	} else if (tubeType == 0 && latticeType ==1) {  // 3. Heterojunctions: CNT on BN lattice	
		//RIMax
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

		//RIMin
    	xShiftRIMin = -BL_NORMALIZATION_HETERO;
    	yShiftRIMin = 0.0;

	} else {  									  // 3. Heterojunctions: BNNT on graphene
		xShiftRIMax = 0.0;
		yShiftRIMax = 0.0;

    	xShiftRIMin = BL_NORMALIZATION_HETERO;
    	yShiftRIMin = 0.0;
	}									  
	
	// creating the hetero tube:

	if (tubeType != latticeType)
	{	
		tempLatticeBL = LATTICE_BL;
		tempTubeBL = TUBE_BL;
		TUBE_BL    = BL_NORMALIZATION_HETERO;
		LATTICE_BL = BL_NORMALIZATION_HETERO;
		heteroTubeN = TubeCreator(T, Ch, teta, tubeType, tubeUnitN, unitcellN, &heteroTube);
		heteroSurfaceN = CreateSurface(&heteroSurfaceLattice, T, unitcellN, radius, latticeType);

		normaliztionTube = &heteroTube;
		normalizationTubeN = &heteroTubeN;
		normalizationSurfaceLattice = &heteroSurfaceLattice;
		normalizationSurfaceN = &heteroSurfaceN;
	} else {
		normaliztionTube = &tube;
		normalizationTubeN = &tubeN;
		normalizationSurfaceLattice = &surfaceLattice;
		normalizationSurfaceN = &surfaceN;
	}
	
	// ---- (I) Placing the tube for easy reach of desired configurations ---------

    Rotate((*normaliztionTube), *normalizationTubeN, 3, (M_PI/6) - teta); // rotation around the z axis (spinning)

    // ------------ (III) Calculating RIMax and RIMin ----------------------

	// Calculating the RI max
	for (j = 0; j < *normalizationTubeN; j++)
	{
		effectiveNum = exp( EXPNORM * (ILD - (*normaliztionTube)[j].z) / (RND - ILD) ); // find weight of atom
		totalEffectiveNum = totalEffectiveNum + effectiveNum;
		if ((*normaliztionTube)[j].z < MAX_HEIGHT)
		{
			*RIMax = *RIMax + effectiveNum * FindInteracting((*normaliztionTube)[j], xShiftRIMax, xShiftRIMax,
														   latticeType);
		}
	}

	WriteCoordinates((*normaliztionTube), *normalizationTubeN, *normalizationSurfaceLattice, *normalizationSurfaceN, 
					 xShiftRIMax, yShiftRIMax, -2, "normalizedRIMax");
	// Calculating the RI min
	Rotate((*normaliztionTube), *normalizationTubeN, 3, rotationAngleRIMin);


	for (j = 0; j < *normalizationTubeN; j++)
	{
		effectiveNum = exp( EXPNORM * (ILD - (*normaliztionTube)[j].z) / (RND - ILD) ); // find weight of atom
		if ((*normaliztionTube)[j].z < MAX_HEIGHT)
		{
			*RIMin = *RIMin + effectiveNum * FindInteracting((*normaliztionTube)[j], xShiftRIMin, yShiftRIMin,
														   latticeType);
		}
	}
	WriteCoordinates((*normaliztionTube), *normalizationTubeN, *normalizationSurfaceLattice, *normalizationSurfaceN, 
					 xShiftRIMin, yShiftRIMin, -1, "normalizedRIMin");

	printf("RIMin:%lf\n", *RIMin);
	printf("RIMax:%lf\n", *RIMax);
	printf("totalEffectiveNum:%lf\n", totalEffectiveNum);

	// ------- (IV) Returning the tube to the original location -----------------
	Rotate((*normaliztionTube), *normalizationTubeN, 3, -((M_PI/6) - teta + rotationAngleRIMin));

	if (tubeType != latticeType) {	
		TUBE_BL    = tempTubeBL;
		LATTICE_BL = tempLatticeBL;
		free(heteroTube);
		free(heteroSurfaceLattice);
	} 

}
