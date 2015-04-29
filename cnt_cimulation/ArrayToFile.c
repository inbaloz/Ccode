#include "Constants_and_libraries.h"
#include "ArrayToFile.h"
#include <stdio.h>

void ArrayToFile(Atom* array, int arrayN, char* destination)
{
	int i;
	FILE* desFile = fopen(destination, "w");
	if (PT == 1)
	{
		fprintf(desFile, "x y z\n");
	}
	for (i = 0; i < arrayN; i++)
	{
		fprintf(desFile, "%e %e %e \n", array[i].x, array[i].y, array[i].z);
	}
	fclose(desFile);
}
