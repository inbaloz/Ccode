#include "Constants_and_libraries.h"
#include "CutLastPartOfTube.h"
#include "FixGapsInTube.h"
#include <string.h>
#include <stdio.h>


/*
 cuts last part of precentTruncated of tube (in a right triangle shape)
 to make the tube look more like the one in the experiments
 assumes the tube aligned with y and x axis (and not rotated yet)

 gets: 
 	tube - the tube we're gonna cut.
 	tubeN - number of atoms in tube
 	precentTruncated - the height precentage we want to cut off (e.g: 0.1)
*/
void CutLastPartOfTube(Atom *tube, int *tubeN, float percentTruncated) {
	int i;
	double a, b; 					// height and width of the right triangle we're gonna cut
	double yMin, yMax, zMin, zMax;	// min and max y value of the tube
	double yCut; 					// the minimal y value that we want to cut from
	double normY, normZ;			// current coordinates of atom
									// normalized by the start of the traingle we want to cut
									// and the height of the tube
	// find ymin and ymax
	yMin = tube[0].y;
	yMax = tube[0].y;
	zMin = tube[0].z;
	zMax = tube[0].z;
	for (i = 0; i < *tubeN; i++)
	{
		if (tube[i].y < yMin)
		{
			yMin = tube[i].y;
		}
		if (tube[i].y > yMax)
		{
			yMax = tube[i].y;
		}
		if (tube[i].z < zMin)
		{
			zMin = tube[i].z;
		}
		if (tube[i].z > zMax)
		{
			zMax = tube[i].z;
		}
	}

	yCut = (1 - percentTruncated) * (yMax - yMin) + yMin;
	a = zMax - zMin;
	b = yMax - yCut;

	for (i = 0; i < *tubeN; i++)
	{
		// we only want to remove atoms that are above yCut
		if (tube[i].y >= yCut)
		{
			normY = tube[i].y - yCut;
			normZ = tube[i].z - zMin;

			if (a*normY + b*normZ - a*b >= 0)
			{
				tube[i].type = DUMMY_TYPE;
			}
		}
	}
	FixGapsInTube(tube, tubeN);
}
