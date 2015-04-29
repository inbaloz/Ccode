#include "Constants_and_libraries.h"
#include <stdio.h>
#include <string.h>
#include "SaveInPar.h"

// input: input parameters, prefix.

// output: no direct output. However, the function
// writes into all the input data to a file.

// ************* To do: *****************
// In the future if you want to add more properties,
// to prevent syntax errors you should create a common
// header file for "SaveInPar.c" and "LoadInPar.c" with
// constants like this one:
// #define STR_SPIN_START "Rotation or spinning start angle = %e\n"
// To prevent typo errors...

void SaveInPar(InPar input, char prefix[])
{
	char destination[80];
	// Duplicating the prefix to prevent it from changing:
	strcpy(destination, prefix);
	strcat(destination, " - inputData");

	FILE* desFile = fopen(destination, "w");

	// Printing the data:
	fprintf(desFile, "Ch.n = %d\n", input.Ch.n);
	fprintf(desFile, "Ch.m = %d\n", input.Ch.m);
	fprintf(desFile, "Amount of tube unitcells = %d\n", input.unitcellN);
	fprintf(desFile, "Motion type = %d\n", input.motionType);
	fprintf(desFile, "Shift angle = %e\n", input.shiftAngle);
	fprintf(desFile, "Rotate angle = %e\n", input.rotateAngle);
	fprintf(desFile, "x shift = %e\n", input.xShift);
	fprintf(desFile, "y shift = %e\n", input.yShift);
	fprintf(desFile, "Rotation or spinning start angle = %e\n", input.rotSpinStart);
	fprintf(desFile, "Rotation or spinning end angle = %e\n", input.rotSpinEnd);
	fprintf(desFile, "Amount of steps = %d\n", input.amountOfSteps);
	fprintf(desFile, "x start = %e\n", input.xStart);
	fprintf(desFile, "y start = %e\n", input.yStart);
	fprintf(desFile, "x end = %e\n", input.xEnd);
	fprintf(desFile, "y end = %e\n", input.yEnd);

	fclose(desFile);
}
