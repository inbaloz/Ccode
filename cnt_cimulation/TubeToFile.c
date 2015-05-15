#include "Constants_and_libraries.h"
#include "TubeToFile.h"
#include <stdio.h>

void TubeToFile(Atom *tube, int tubeN, char *destination)
{
	int i;
	FILE* desFile = fopen(destination, WRITE);
	
	if (PT == 1)
	{
		fprintf(desFile, "x y z\n");
	}

	for (i = 0; i < tubeN; i++)
	{
		fprintf(desFile, "%e %e %e\n", tube[i].x, tube[i].y, tube[i].z);
	}
	fclose(desFile);
}
