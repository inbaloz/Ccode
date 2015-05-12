#include "Constants_and_libraries.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ArrayToFile.h"
#include "aVecLength.h"
#include "CalculateIntersection.h"
#include "CalculateTranslational.h"
#include "CLI.h"
#include "CutUnitcell.h"
#include "DuplicateTube.h"
#include "EffectiveNum.h"
#include "gcd.h"
#include "LatticeCreator.h"
#include "LoadInPar.h"
#include "NormRI.h"
#include "PerfectRotationMotion.h"
#include "Rad2Deg.h"
#include "RotationMotion.h"
#include "SaveInPar.h"
#include "SlidingMotion.h"
#include "SpinningMotion.h"
#include "TwodDataToFile.h"

#include "FindInteracting.h"

// This program creates a nanotube and lattice to preform
// simulations using the RI (registry index) method.

/*
 *	Step 1: Recieve input data
 *
 *	Step 2: Calculate tube parameters
 *
 *	Step 3: Create the tube
 *
 *	Step 4: Normalize RI
 *
 *	Step 5: Calculate RI
 *
 */

/*
 *	Types of motion:
 *
 *	1. Rotating: start angle, end angle, amount of steps.
 *
 *	2. Spinning: Start angle, end angle, amount of steps.
 *
 *	3. Sliding: Start x-y, end x-y, amount of steps.
 *
 *	4. Rotation without sliding: start angle, end angle,
 *	amount of steps, start x-y.
 *
 */


/*
 *	Notes:
 *
 *	1. Shift angles are implemented by rotating the tube,
 *	however x,y shifts are implemented by moving the lattice.
 *	
 *	2. There is no actual lattice. The lattice is cyclic,
 *	therefore, we calculate (x % HORIZD) and (y % HIGHT).
 *	This way, we can find which lattice atoms would have
 *	interact with the tube atom (if we actually had a lattice)
 *	and create them temporarily.
 *
 *	3. It seems like note 1 and 2 contradict each other.
 *	They don't, the x,y shifts are taken into consideration
 *	by calculating ((x + xShift) % HORIZD) instead of
 *	(x % HORIZD), as well as ((y + yShift) % HIGHT) instead of
 * 	(y % HIGHT).
 *
 *
 */

double LATTICE_BL = 1.44;
double TUBE_BL = 1.44;
double ILD = 3.33;

int main(int argc, char *argv[])
{
	int i;					// Loop counter.
	int runOrGui;			// 1 for no gui, and 0 for gui.
	char prefix[80];		// File prefixes.

//-----------------from configuration ----------------------------

	int motionType;			// Type of tube motion.
	int unitcellN;			// Amount of tube unitcells.
	aVec Ch, T;				// Chirality and translational vectors.
	double shiftAngle = 0;	// The angle between tube axis and y axis (type 1, 3, 4)
	double rotateAngle = 0;	// The angle of rotation of the tube around his axis (type 2, 3)
	double xShift;			// Tube's x axis shift (type 1, 2)
	double yShift;			// Tube's y axis shift (type 1, 2)
	double teta;			// Chiral angle
	double rotSpinStart;	// Starting rotation angle	(type 1, 2, 4)
	double rotSpinEnd;		// Ending rotation angle	(type 1, 2, 4)
	double rotSpinStep;		// Rotation step...			(type 1*, 2*, 4*)
	int amountOfSteps;		// Amount of angle steps	(type 1, 2, 3, 4)
	double* rotSpinValues;	// Rotation angle values	(type 1*, 2*)
	double xStart;			// Starting x.				(type 3, 4)
	double yStart;			// Starting y.				(type 3, 4)
	double xEnd;			// Ending x.				(type 3)
	double yEnd;			// Ending y.				(type 3)
	double percentTruncated;// amount of Truncated tube from the end, in a right-angle triangular fashion.
	int tubeType = 0;		// tube type: 	 0 for CNT, 1 for BN.
	int latticeType = 0;	// lattice type: 0 for graphene, 1 for BN.

//-----------------------------------------------------------------------------

	double xStep;			// x Step...				(type 3*, 4*)
	double yStep;			// y Step...				(type 3*, 4*)
	double* slideValues;	// Slide values				(type 3*, 4*)
	double slideStep;		// Slide step values		(type 3*, 4*)
	double totDist;			// the total distance travelled by the tube (type 4)
	Atom* tubeUnit;			// The tube's unitcell
	int tubeUnitN;			// Number of atoms in the tube's unitcell
	double radius;			// The tube's radius
	double length;			// The tube's length
	Atom* tube;				// The final tube
	int tubeN;				// The number of atoms in the tube
	Atom* lattice;			// The lattice that will be used to create the tube.
	int latticeN;			// The number of atoms in the lattice
	double RIMin;			// The minimum surface
	double RIMax;			// The maximum surface
	// double effTubeN;		// The effective number of atoms (normalized by hight)	
	double* RI;				// The registry index data.
							// type: used to mark what variables are used to
							// each type of motion. If the number is marked
							// by a '*' it means that variable is used but it isn't
							// an input, it is calculated from the other inputs.
	InPar input;			// Input parameters.

//********************** Step 1 - Recieve input data ***************************
	
	// Load data
	switch (argc)
	{
	case 1:
		// Load default data
		strcpy(prefix, "Default");
		break;
	case 2:
		// Load data of the file name:
		strcpy(prefix, argv[1]);
		break;
	default:
		printf("Not a valid 1st argument, please enter -1 or 0 or 1.\n");
		exit(0);
		break;
	}
	input = LoadInPar(prefix);


// The values are duplicated for extra readability:
// (and because there aren't many of them..)
	Ch.n = input.Ch.n;
	Ch.m = input.Ch.m;
	unitcellN = input.unitcellN;
	motionType = input.motionType;
	shiftAngle = input.shiftAngle;
	rotateAngle = input.rotateAngle;
	xShift = input.xShift;
	yShift = input.yShift;
	rotSpinStart = input.rotSpinStart;
	rotSpinEnd = input.rotSpinEnd;
	amountOfSteps = input.amountOfSteps;
	xStart = input.xStart;
	yStart = input.yStart;
	xEnd = input.xEnd;
	yEnd = input.yEnd;
	percentTruncated = input.percentTruncated;
	tubeType = input.tubeType;
	latticeType = input.latticeType;
	


//********************** Step 2 - Calculate tube parameters ********************
	T = CalculateTranslational(Ch);
	tubeUnitN = 4 * ( ((Ch.m) * (Ch.m)) + ((Ch.n) * (Ch.n)) + ((Ch.n) * (Ch.m)) ) / gcd(Ch.n, Ch.m);
	teta = acos( (2 * Ch.n + Ch.m) / (2 * sqrt( (Ch.n * Ch.n) + (Ch.m * Ch.m) + (Ch.n * Ch.m) ) ) );
	radius = aVecLength(Ch) / (2 * M_PI);
	length = aVecLength(T) * unitcellN;

//********************** Step 3 - Create the tube ******************************

	// The function needs the future unitcell borders (before rotating and moving)
	// to create a lattice in a suitable size:
	{
	double xMin = - aVecLength(T) * sin( (M_PI / 6) - teta);
	double xMax = aVecLength(Ch) * cos((M_PI / 6) - teta);
	double yMin = 0;
	double yMax = ( aVecLength(Ch) * sin((M_PI / 6) - teta) + aVecLength(T) * cos( (M_PI / 6) - teta));
	
	// Creating the lattice (to make a tube of):
	latticeN = LatticeCreator(&lattice,	xMin, yMin, xMax, yMax, tubeType);
	}
	// Making the tube's unitcell from the lattice:
	tubeUnit = CutUnitcell(lattice, latticeN, Ch, T, (M_PI / 6) - teta, tubeUnitN);

	// Duplicating the unitcell:
	tube = DuplicateTube(tubeUnit, tubeUnitN, unitcellN, aVecLength(T));
	tubeN = tubeUnitN * unitcellN;
	// Free the unit tube, it isn't needed anymore:
	free(tubeUnit);
	// Initially rotating and spinning the tube as requested (around z axis).
	Rotate(tube, tubeN, 3, shiftAngle);
	RotateShift(tube, tubeN, rotateAngle, shiftAngle, ILD + radius);
	// This lattice was used to create the tube. Now we free it for later reuse:
	free(lattice);

//********************** Step 4 - Normalizie RI ********************************
												
	// effTubeN = EffectiveNum(tube, tubeN, ILD, RND); 
	// RIMin = effTubeN * MIN(M_PI * pow(RCGRAPHENE,2), M_PI * pow(RCCNT,2)) / 2;
	// RIMax = RIMin * 2;

	Rotate(tube, tubeN, 3, -shiftAngle + M_PI/6 - teta); // get to AA
	double AAshiftx = 0;
	double AAshifty = 0;


    //------------------ calculating RImax -----------------------------
	double effectiveNum, currentInteracting;
	RIMax = 0;
	for (i = 0; i < tubeN; i++)
	{
		effectiveNum = exp( EXPNORM * (ILD - tube[i].z) / (RND - ILD) );
		if (effectiveNum > NP)
		{
			currentInteracting = FindInteracting(tube[i], AAshiftx, AAshifty);
			RIMax = RIMax + effectiveNum * currentInteracting;
		}
	}

	printf("RI max: %d\n", RIMax);

	//------------------- calculating RIMin ------------------------------
	RIMin = 0;
	double ABshiftx = AAshiftx - 1.25*LATTICE_BL;
	for (i = 0; i < tubeN; i++)
	{
		effectiveNum = exp( EXPNORM * (ILD - tube[i].z) / (RND - ILD) );
		if (effectiveNum > NP)
		{
			currentInteracting = FindInteracting(tube[i], ABshiftx, AAshifty);
			RIMin = RIMin + effectiveNum * currentInteracting;
		}
	}

	Rotate(tube, tubeN, 3, -(-shiftAngle + M_PI/6 - teta));

	
//********************** Step 5 - Calculate RI *********************************

	switch(motionType)
	{
	case 1: // Rotating: start angle, end angle, amount of steps.
		// Defining data arrays:
		RI = malloc(amountOfSteps * sizeof(double));
		rotSpinValues = malloc(amountOfSteps * sizeof(double));
		// Calculating constants and rotSpinValues:
		rotSpinStep = ( (rotSpinEnd - rotSpinStart) / (amountOfSteps - 1) );
		for (i = 0; i < amountOfSteps; i++)
		{
			rotSpinValues[i] = rotSpinStart + (rotSpinStep * i);
		}
		// Rotating by "rotSpinStart":
		RotateShift(tube, tubeN, rotSpinStart, shiftAngle, ILD + radius);
		// Calculating RI:
		RotationMotion(RI, rotSpinStep, amountOfSteps, tube, tubeN, radius, shiftAngle, xShift, yShift);
		// Printing the RI and rotSpinValues to a file:
		NormRI(RI, amountOfSteps, RIMin, RIMax);
		Rad2Deg(rotSpinValues, amountOfSteps);
		TwodDataToFile(rotSpinValues, RI, amountOfSteps, strcat(prefix, " - Rotation RI Data"));
		free(rotSpinValues);
		free(RI);
		break;
	case 2: // Spinning: Start angle, end angle, amount of steps.
		// Defining data arrays:
		RI = malloc(amountOfSteps * sizeof(double));
		rotSpinValues = malloc(amountOfSteps * sizeof(double));
		// Calculating constants and rotSpinValues:
		rotSpinStep = ( (rotSpinEnd - rotSpinStart) / (amountOfSteps - 1) );
		for (i = 0; i < amountOfSteps; i++)
		{
			rotSpinValues[i] = rotSpinStart + (rotSpinStep * i);
		}
		// Spinning by "rotSpinStart":
		Rotate(tube, tubeN, 3, rotSpinStart);
		// Calculating RI:
		SpinningMotion(RI, rotSpinStep, amountOfSteps, tube, tubeN, xShift, yShift);			
		// Printing the RI and rotSpinValues to a file:
		NormRI(RI, amountOfSteps, RIMin, RIMax);
		Rad2Deg(rotSpinValues, amountOfSteps);
		TwodDataToFile(rotSpinValues, RI, amountOfSteps, strcat(prefix, " - Spinning RI Data"));
		free(rotSpinValues);
		free(RI);
		break;
	case 3: // Sliding: Start x-y, end x-y, amount of steps.
		// Defining data arrays:
		RI = malloc(amountOfSteps * sizeof(double));
		slideValues = malloc(amountOfSteps * sizeof(double));
		// Calculating constants and SlideStep:
		xStep = ( (xEnd - xStart) / (amountOfSteps - 1) );
		yStep = ( (yEnd - yStart) / (amountOfSteps - 1) );
		slideStep = sqrt( (xStep * xStep) + (yStep * yStep) );
		for (i = 0; i < amountOfSteps; i++)
		{
			slideValues[i] = slideStep * i;
		}
		// Calculating RI:
		SlidingMotion(RI, xStep, yStep, amountOfSteps, tube, tubeN, xStart, yStart);			
		// Printing the RI and x-yValues to a file:
		NormRI(RI, amountOfSteps, RIMin, RIMax);
		TwodDataToFile(slideValues, RI, amountOfSteps, strcat(prefix, " - Sliding RI Data"));
		free(slideValues);
		free(RI);
		break;
	case 4: // Rotation without sliding (rolling): start angle, end angle, amount of steps, start x-y.
		// Defining data arrays:
		RI = malloc(amountOfSteps * sizeof(double));
		slideValues = malloc(amountOfSteps * sizeof(double));
		// Calculating constants and rotSpinValues:
		rotSpinStep = ( (rotSpinEnd - rotSpinStart) / (amountOfSteps - 1) );

		//--------------------------------------------
        // I tried adding to the code - hope it's ok

        totDist = (rotSpinEnd - rotSpinStart) * radius;
        xStep = ( (totDist * cos(shiftAngle)) / (amountOfSteps - 1) );
        yStep = ( (totDist * sin(shiftAngle)) / (amountOfSteps - 1) );

        //--------------------------------------------

		slideStep = sqrt( (xStep * xStep) + (yStep * yStep) );

		//xStep = slideStep * radius * cos(shiftAngle); // why not multiply by rotSpinStep?
		//yStep = slideStep * radius * sin(shiftAngle);
		rotSpinStep = ( (rotSpinEnd - rotSpinStart) / (amountOfSteps - 1) ); // why do it twice?

		for (i = 0; i < amountOfSteps; i++)
		{
			slideValues[i] = slideStep * i;
		}
		// Rotating by "rotSpinStart":
		RotateShift(tube, tubeN, rotSpinStart, shiftAngle, ILD + radius);
		// Calculating RI:
		PerfectRotationMotion(RI, xStart, yStart, xStep, yStep, rotSpinStep, amountOfSteps, tube, tubeN, radius, shiftAngle);
		// Printing the RI and x-yValues to a file:
		NormRI(RI, amountOfSteps, RIMin, RIMax);
		TwodDataToFile(slideValues, RI, amountOfSteps, strcat(prefix, " - Perfect Rotation RI Data"));
		free(slideValues);
		free(RI);
		break;
	}
	// Free the tube:
	free(tube);
}
