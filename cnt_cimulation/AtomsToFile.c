#include "Constants_and_libraries.h"
#include "AtomsToFile.h"
#include <stdio.h>

// getting types:
// C - 0
// N - 1
// B - 2
char temp;

void AtomsToFile(Atom *tube, int tubeN, char *destination, double xShift, double yShift,
				 double zshift, int append)
{
	int i;
	FILE* desFile;

	if (append) {
		desFile = fopen(destination, "a");
	}
	else {
		desFile = fopen(destination, "w");
	}
	
	if (PT == 1) {
		fprintf(desFile, "x y z type\n");
	}

	for (i = 0; i < tubeN; i++) {
		fprintf(desFile, "%lf %lf %lf %d\n", tube[i].x + xShift, tube[i].y 
				+ yShift, tube[i].z + zshift, tube[i].type);
	}
			
	fclose(desFile);
}
