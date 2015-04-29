#include "Constants_and_libraries.h"
#include "FindInteracting.h"
#include "CalculateIntersection.h"
#include <math.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

// Note: the function assumes we have a standard lattice. This assumption
// is always true unless you expand the codes functionality by changing
// the way lattices and tubes are located.

double FindInteracting(Atom atom, double xShift, double yShift)
{
	// Ensuring we have a non-negative value in the right range:
	double xMod = remainder(atom.x + xShift, HORIZD);
	double yMod = remainder(atom.y + yShift, HIGHT);
	double RI = 0; // It is the addition to RI before normalizing.
	// xMod and yMod need to be not negative:
	if (xMod < 0)
	{
		xMod = xMod + HORIZD;
	}
	if (yMod < 0)
	{
		yMod = yMod + HIGHT;
	}
	Atom atomMod = { .x = xMod, .y = yMod, .z = atom.z };
	Atom temp = { .x = 0, .y = 0, .z = 0 };

	if (xMod >= 0 && xMod < BL)
	{
		temp.x = 0;		//	**ALREADY THAT VALUE**
		temp.y = 0;		//	**ALREADY THAT VALUE**
		RI += CalculateIntersection(atomMod, temp);
		temp.x = -0.5 * BL;
		temp.y = 0.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 0;
		temp.y = HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = BL;
		temp.y = HIGHT;	//	**ALREADY THAT VALUE**
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 1.5 * BL;
		temp.y = 0.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = BL;
		temp.y = 0;
		RI += CalculateIntersection(atomMod, temp);
	}
	else if	(xMod >= BL && xMod < (1.5 * BL))
	{
		temp.x = 1.5 * BL;
		temp.y = -0.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = BL;
		temp.y = 0;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 1.5 * BL;
		temp.y = 0.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = BL;
		temp.y = HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 1.5 * BL;
		temp.y = 1.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
	}
	else if (xMod >= (1.5 * BL) && xMod < (2.5 * BL))
	{
		if (yMod < (0.5 * HIGHT))
		{
			temp.x = 1.5 * BL;
			temp.y = -0.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = BL;
			temp.y = 0;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 1.5 * BL;
			temp.y = 0.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * BL;
			temp.y = 0.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 3 * BL;
			temp.y = 0;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * BL;
			temp.y = -0.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
		}
		else
		{
			temp.x = 1.5 * BL;
			temp.y = 0.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = BL;
			temp.y = HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 1.5 * BL;
			temp.y = 1.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * BL;
			temp.y = 1.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 3 * BL;
			temp.y = HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * BL;
			temp.y = 0.5 * HIGHT;
			RI += CalculateIntersection(atomMod, temp);		
		}
	}
	else	// xMod >= (2.5 * BL) && xMod < (3 * BL)
	{
		temp.x = 2.5 * BL;
		temp.y = -0.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 3 * BL;
		temp.y = 0;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 2.5 * BL;
		temp.y = 0.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 3 * BL;
		temp.y = HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 2.5 * BL;
		temp.y = 1.5 * HIGHT;
		RI += CalculateIntersection(atomMod, temp);		
	}
	return RI;

}
