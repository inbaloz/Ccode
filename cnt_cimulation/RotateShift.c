#include "Constants_and_libraries.h"
#include "RotateShift.h"
#include "Rotate.h"
#include "Move.h"
#include <math.h>

// inupt: An atom array, the size of it, angle of rotation, the
// angle of the rotation axis (2d x-y), and the hight of the center.
// Notes:	1. The axis of rotation is on the x-y plane however
// all coordinates can be changed because the rotation is
// around y axis (while shiftAngle = 0).
//			2. The function assumes that the tube is centered by x-y
// but not necessarily in z.

// output: The atom array will be rotated in the "shiftAngle" axis by
// "angle" radians.

// How it works: we lower the array by hight, we rotate the array by
// "-shiftAngle" in x-y plane, then we rotate by "angle" in xz left
// handed plane, then we rotate back by "shiftAngle" in x-y plane,
// then we lift the array by hight.

void RotateShift(Atom* array, int arraySize, double angle, double shiftAngle, double hight)
{
	if (angle != 0)
	{
		Move(array, arraySize, 0, 0, -hight);
		Rotate(array, arraySize, 3, -shiftAngle);
		Rotate(array, arraySize, 2, angle);
		Rotate(array, arraySize, 3, shiftAngle);
		Move(array, arraySize, 0, 0, hight);
	}
}
