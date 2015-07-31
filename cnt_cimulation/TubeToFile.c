#include "Constants_and_libraries.h"
#include "TubeToFile.h"
#include <stdio.h>

// getting types:
// C - 0
// N - 1
// B - 2
char temp;

void TubeToFile(Atom *tube, int tubeN, char *destination)
{
	int i;
	FILE* desFile = fopen(destination, WRITE);
	
	if (PT == 1)
	{
		fprintf(desFile, "x y z type\n");
	}

	for (i = 0; i < tubeN; i++)
	{
		fprintf(desFile, "%lf %lf %lf %d\n", tube[i].x, tube[i].y, tube[i].z, tube[i].type);
	}
			
	fclose(desFile);
}
