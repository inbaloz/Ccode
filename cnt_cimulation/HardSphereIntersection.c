#include "Constants_and_libraries.h"
#include "HardSphereIntersection.h"
#include <math.h>
#include <stdio.h>

double HardSphereIntersection(double radius1, double radius2, double distance)
{
	double intersection;

	// Too far:
	if (distance > radius1 + radius2)
	{
		return 0;
	}
	// Very close:
	else if (distance <= ABS(radius1 - radius2))
	{
		return MIN(M_PI * pow(radius1,2), M_PI * pow(radius2,2));
	}
	
	// In between:
	else
	{
		// Sum of to sectors of circles minus the surface of kite created by them:
		intersection =	pow(radius1,2) * acos((pow(distance,2) + pow(radius1,2) - pow(radius2,2)) / (2 * distance * radius1)) + 
						pow(radius2,2) * acos((pow(distance,2) + pow(radius2,2) - pow(radius1,2)) / (2 * distance * radius2)) - 
						// Heron's formula for triangle surface multiplied by 2:
						0.5 * sqrt( (radius1 + radius2 - distance) * (radius1 - radius2 + distance) * (- radius1 + radius2 + distance) * (radius1 + radius2 + distance) );
		return intersection;
	}

}
