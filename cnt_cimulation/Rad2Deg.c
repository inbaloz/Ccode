#include "Constants_and_libraries.h"
#include "Rad2Deg.h"
#include <math.h>

// input: data in radians, amount of data.
// output: the data is transfered to degrees.

void Rad2Deg(double* data, int dataN)
{
	int i;
	for (i = 0; i < dataN; i++)
	{
		data[i] = 180 * data[i] / M_PI;
	}
}
