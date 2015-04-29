#include "Constants_and_libraries.h"
#include "Move.h"

void Move(Atom* array, int arraySize, double x, double y, double z)
{
	int i;
	for (i = 0; i < arraySize; i++)
	{
		array[i].x = array[i].x + x;
		array[i].y = array[i].y + y;
		array[i].z = array[i].z + z;
	}
}
