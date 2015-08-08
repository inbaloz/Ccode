#include "Constants_and_libraries.h"
#include "AtomsToFile.h"
#include <stdio.h>

// This funcion creates coordinate filez that will eork with Itai's code.
// They are organized:
// atom_type	layer (0 or 1)	x 	y 	z

// The tube is considered layer 0. It is written first and hence we use "write"
// and not "append". Only atoms with z<MAX_HEIGHT are considered.
// The surface is layer 1.

void AtomsToFile(Atom *tube, int tubeN, char *destination, double xShift, double yShift,
				 double zshift, int append)
{
	int i;
	FILE* desFile;

	if (append) { // append == 1 means we are writing the surface.
		desFile = fopen(destination, "a");

		for (i = 0; i < tubeN; i++) {
			fprintf(desFile, "%d 1 %lf %lf %lf\n", tube[i].type, tube[i].x + xShift, tube[i].y 
				+ yShift, tube[i].z + zshift);
		}
	}
	
	else { // append == 0 means we are writing the tube.
		desFile = fopen(destination, "w");
	
		fprintf(desFile, "0.0 0.0 0.0\n");
		for (i = 0; i < tubeN; i++) {
			if (tube[i].z <= MAX_HEIGHT) {
				fprintf(desFile, "%d 0 %lf %lf %lf\n", tube[i].type, tube[i].x + xShift, 
				tube[i].y + yShift, tube[i].z + zshift);	
			}
		}
	}

			
	fclose(desFile);
}
