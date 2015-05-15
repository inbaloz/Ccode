#include "Constants_and_libraries.h"
#include "FindInteracting.h"
#include "CalculateIntersection.h"
#include <math.h>
#include <stdio.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

// Note: the function assumes we have a standard lattice. This assumption
// is always true unless you expand the codes functionality by changing
// the way lattices and tubes are located.

double FindInteracting(Atom atom, double xShift, double yShift)
{
	// Ensuring we have a non-negative value in the right range:
	double xMod = remainder(atom.x + xShift, LATTICE_HORIZD);
	double yMod = remainder(atom.y + yShift, LATTICE_HIGHT);
	double RI = 0; // It is the addition to RI before normalizing.
	// xMod and yMod need to be not negative:
	if (xMod < 0)
	{
		xMod = xMod + LATTICE_HORIZD;
	}
	if (yMod < 0)
	{
		yMod = yMod + LATTICE_HIGHT;
	}
	Atom atomMod = { .x = xMod, .y = yMod, .z = atom.z , .type = 'C'};
	Atom temp = { .x = 0, .y = 0, .z = 0 , .type = 'C'};

	if (xMod >= 0 && xMod < LATTICE_BL)
	{
		temp.x = 0;		//	**ALREADY THAT VALUE**
		temp.y = 0;		//	**ALREADY THAT VALUE**
		RI += CalculateIntersection(atomMod, temp);
		temp.x = -0.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 0;
		temp.y = LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = LATTICE_BL;
		temp.y = LATTICE_HIGHT;	//	**ALREADY THAT VALUE**
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 1.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = LATTICE_BL;
		temp.y = 0;
		RI += CalculateIntersection(atomMod, temp);
	}
	else if	(xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL))
	{
		temp.x = 1.5 * LATTICE_BL;
		temp.y = -0.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = LATTICE_BL;
		temp.y = 0;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 1.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = LATTICE_BL;
		temp.y = LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 1.5 * LATTICE_BL;
		temp.y = 1.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
	}
	else if (xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL))
	{
		if (yMod < (0.5 * LATTICE_HIGHT))
		{
			temp.x = 1.5 * LATTICE_BL;
			temp.y = -0.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = LATTICE_BL;
			temp.y = 0;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 1.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 3 * LATTICE_BL;
			temp.y = 0;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * LATTICE_BL;
			temp.y = -0.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
		}
		else
		{
			temp.x = 1.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = LATTICE_BL;
			temp.y = LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 1.5 * LATTICE_BL;
			temp.y = 1.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * LATTICE_BL;
			temp.y = 1.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 3 * LATTICE_BL;
			temp.y = LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);
			temp.x = 2.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			RI += CalculateIntersection(atomMod, temp);		
		}
	}
	else	// xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
	{
		temp.x = 2.5 * LATTICE_BL;
		temp.y = -0.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 3 * LATTICE_BL;
		temp.y = 0;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 2.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 3 * LATTICE_BL;
		temp.y = LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
		temp.x = 2.5 * LATTICE_BL;
		temp.y = 1.5 * LATTICE_HIGHT;
		RI += CalculateIntersection(atomMod, temp);
	}
	return RI;

}
