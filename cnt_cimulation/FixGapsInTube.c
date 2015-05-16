#include "Constants_and_libraries.h"
#include "FixGapsInTube.h"
#include <stdlib.h>
#include <stdio.h>

/*
 gets a tube with gaps, represented by atoms with dummy type as type
 fixes the tube by shifting the cells left where there are gaps
 algorithm:
  	2 indexes, one for new cells and one for old cells
 	the index for the old cells keeps going right
 	the new cells index is going right until it gets to a real cell (indicated by dummy type)
 	the old cells index keeps going right one step for each iteration
	
*/
void FixGapsInTube(Atom *tube, int *tubeN) {
	int i, j;
	int newTubeN;

	i = 0;
	j = 0;
	newTubeN = *tubeN;
	for (i = j = 0; j < *tubeN; i++, j++) {
		while ((tube[j].type == DUMMY_TYPE) && (j < *tubeN))
		{
			j++;
			newTubeN--;
		}
		if (j < *tubeN)
		{
			// now j is on a real cell (and not a gap cell)
			tube[i] = tube[j];
		}
	}
	*tubeN = newTubeN;
}