#include "Constants_and_libraries.h"
#include "DuplicateTube.h"
#include "Move.h"
#include <malloc.h>

// input: a centered tube unitcell, amount of atoms, amount of duplications,
// unitcell length (aVecLength(T)).
// output: a centered duplicate of the tube with duplicateN unitcells.
Atom * DuplicateTube(Atom* tubeUnit, int tubeUnitN, int duplicateN, double tubeLength)
{
	int i, j;	// counters.
	Atom* tubeDuplicated = malloc(tubeUnitN * duplicateN * sizeof(Atom));
	for (i = 0; i < duplicateN; i++)
	{
		for (j = 0; j < tubeUnitN; j++)
		{
			tubeDuplicated[(i * tubeUnitN) + j].x = tubeUnit[j].x;
			tubeDuplicated[(i * tubeUnitN) + j].z = tubeUnit[j].z;
			// The tube's axis is y axis. Therefore, only y is changed
			// between unitcells.
			tubeDuplicated[(i * tubeUnitN) + j].y = tubeUnit[j].y + (tubeLength * i);
		}
	}
	// The new tube isn't centered, the (0,0,0) point is only in the middle
	// of the first unitcell. Therefore, lets move it:
	Move(tubeDuplicated, tubeUnitN * duplicateN, 0, - ((duplicateN - 1) * tubeLength) / 2.0, 0);
	return tubeDuplicated;
}


