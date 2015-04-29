#include "Constants_and_libraries.h"
#include "TwodDataToFile.h"
#include <stdio.h>

void TwodDataToFile(double* xValues, double* yValues, int amountOfData, char* destination)
{
	int i;
	FILE* desFile = fopen(destination, WRITE);
	if (PT == 1)
	{
		fprintf(desFile, "x y\n");
	}
	for (i = 0; i < amountOfData; i++)
	{
		fprintf(desFile, "%e %e\n", xValues[i], yValues[i]);
	}
	fclose(desFile);
}
