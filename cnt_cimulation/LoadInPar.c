#include "Constants_and_libraries.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "LoadInPar.h"

// input: prefix.

// output: Input parrameters.

InPar LoadInPar(char prefix[])
{
	char destination[80];
	char line[100], *res;
	int i;			// Simply a counter.
	InPar input;


	// default values for backwards compatability
	input.percentTruncated = 0.0;
	input.tubeType = 0;
	input.latticeType = 0;
	// Duplicating the prefix to prevent it from changing:
	strcpy(destination, prefix);
	strcat(destination, " - inputData");
	FILE *desFile = fopen(destination, "r");
	fread(line, 100, 1, desFile);
	// res = fgets(line, 100, desFile);
	while (res != NULL)
	{
		i = 0;
		while (i != 100 && line[i] != '=')
		{
			i++;
		}

		if (strncmp(line, "Ch.n", 4) == 0)
		{
			input.Ch.n = atoi(&line[i + 2]);
		}
		else if (strncmp(line, "Ch.m", 4) == 0)
		{
			input.Ch.m = atoi(&line[i + 2]);
		}
		else if (strncmp(line, "Amount of tube unitcells", 24) == 0)
		{
			input.unitcellN = atof(&line[i + 2]);
		}
		else if (strncmp(line, "Motion type", 11) == 0)
		{
			input.motionType = atof(&line[i + 2]);
		}
		else if (strncmp(line, "Shift angle", 11) == 0)
		{
			input.shiftAngle = atof(&line[i + 2]);
		}
		else if (strncmp(line, "Rotate angle", 12) == 0)
		{
			input.rotateAngle = atof(&line[i + 2]);
		}
		else if (strncmp(line, "x shift", 7) == 0)
		{
			input.xShift = atof(&line[i + 2]);
		}
		else if (strncmp(line, "y shift", 7) == 0)
		{
			input.yShift = atof(&line[i + 2]);
		}
		else if (strncmp(line, "Rotation or spinning start angle", 32) == 0)
		{
			input.rotSpinStart = atof(&line[i + 2]);
		}
		else if (strncmp(line, "Rotation or spinning end angle", 30) == 0)
		{
			input.rotSpinEnd = atof(&line[i + 2]);
		}
		else if (strncmp(line, "Amount of steps", 15) == 0)
		{
			input.amountOfSteps = atoi(&line[i + 2]);
		}
		else if (strncmp(line, "x start", 7) == 0)
		{
			input.xStart = atof(&line[i + 2]);
		}
		else if (strncmp(line, "y start", 7) == 0)
		{
			input.yStart = atof(&line[i + 2]);
		}
		else if (strncmp(line, "x end", 5) == 0)
		{
			input.xEnd = atof(&line[i + 2]);
		}
		else if (strncmp(line, "y end", 5) == 0)
		{
			input.yEnd = atof(&line[i + 2]);
		} 
		else if (strncmp(line, "percent truncated", 17) == 0)
		{
			input.percentTruncated = atof(&line[i + 2]);
		}
		else if (strncmp(line, "tube type", 9) == 0)
		{
			input.tubeType = atoi(&line[i + 2]);
		}
		else if (strncmp(line, "lattice type", 12) == 0)
		{
			input.latticeType = atoi(&line[i + 2]);
		}
		else
		{
			printf("a messed up line\n");
		}

		res = fgets(line, 100, desFile);
	}
	fclose(desFile);
	return input;
}
