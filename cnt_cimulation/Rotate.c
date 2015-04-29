#include "Constants_and_libraries.h"
#include "Rotate.h"
#include <math.h>

// input: array of atmos, the array's size, the axis that isn't used,
// and the angle of rotation.

// ouput: the array's atoms are rotated in the unspecified axis's by
// the angle "angle".

void Rotate (Atom* array, int arraySize, int unUsedAxis, double angle)
{
	int i;
	double temp;
	if (angle != 0)
	{
		switch (unUsedAxis)
			{
				case 1:
					for (i = 0; i < arraySize; i++)
					{
						temp = array[i].y;
						array[i].y =	(temp * cos(angle)) + (- array[i].z * sin(angle));
						array[i].z =	(temp * sin(angle)) + (array[i].z * cos(angle));
					}
				break;
				case 2:
					for (i = 0; i < arraySize; i++)
					{
						temp = array[i].z;
						array[i].z =	(temp * cos(angle)) + (- array[i].x * sin(angle));
						array[i].x =	(temp * sin(angle)) + (array[i].x * cos(angle));
					}
				break;
				case 3:
					for (i = 0; i < arraySize; i++)
					{
						temp = array[i].x;
						array[i].x =	(temp * cos(angle)) + (- array[i].y * sin(angle));
						array[i].y =	(temp * sin(angle)) + (array[i].y * cos(angle));
					}
				break;
			}
	}
}
