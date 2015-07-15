#include "Constants_and_libraries.h"
#include "FindInteracting.h"
#include "CalculateIntersection.h"
#include <math.h>
#include <stdio.h>

// Input: a tube atom, the x,y requested shift of the tube's center.

// Output: the additiont to RI as a result of the interaction between
// the input atom and the lattice atoms.

double FindInteracting(Atom atom, double xShift, double yShift, int latticeType)
{
	// Types of atoms
	int firstAtom, secondAtom;
	firstAtom = latticeType == 0 ? C_type : B_type; // the BN lattice starts at Boron.
	secondAtom = latticeType == 0 ? C_type : N_type;

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
	Atom atomMod = { .x = xMod, .y = yMod, .z = atom.z , .type = atom.type};
	Atom temp = { .x = 0, .y = 0, .z = 0 , .type = C_type};

	if (xMod >= 0 && xMod < LATTICE_BL)
	{
		// 1
		temp.x = 0;		//	**ALREADY THAT VALUE**
		temp.y = 0;		//	**ALREADY THAT VALUE**
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 2
		temp.x = -0.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 3
		temp.x = 0;
		temp.y = LATTICE_HIGHT;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 4
		temp.x = LATTICE_BL;
		temp.y = LATTICE_HIGHT;	//	**ALREADY THAT VALUE**
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 5
		temp.x = 1.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 6
		temp.x = LATTICE_BL;
		temp.y = 0;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
	}
	else if	(xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL))
	{
		// 1
		temp.x = 1.5 * LATTICE_BL;
		temp.y = -0.5 * LATTICE_HIGHT;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 2
		temp.x = LATTICE_BL;
		temp.y = 0;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 3
		temp.x = 1.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 4
		temp.x = LATTICE_BL;
		temp.y = LATTICE_HIGHT;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 5
		temp.x = 1.5 * LATTICE_BL;
		temp.y = 1.5 * LATTICE_HIGHT;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
	}
	else if (xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL))
	{
		if (yMod < (0.5 * LATTICE_HIGHT))
		{
			// 1
			temp.x = 1.5 * LATTICE_BL;
			temp.y = -0.5 * LATTICE_HIGHT;
			temp.type = firstAtom; // B or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 2
			temp.x = LATTICE_BL;
			temp.y = 0;
			temp.type = secondAtom; // N or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 3
			temp.x = 1.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			temp.type = firstAtom; // B or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 4
			temp.x = 2.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			temp.type = secondAtom; // N or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 5
			temp.x = 3 * LATTICE_BL;
			temp.y = 0;
			temp.type = firstAtom; // B or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 6	
			temp.x = 2.5 * LATTICE_BL;
			temp.y = -0.5 * LATTICE_HIGHT;
			temp.type = secondAtom; // N or C
			RI += CalculateIntersection(atomMod, temp);
		}
		else
		{
			// 1
			temp.x = 1.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			temp.type = firstAtom; // B or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 2
			temp.x = LATTICE_BL;
			temp.y = LATTICE_HIGHT;
			temp.type = secondAtom; // N or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 3
			temp.x = 1.5 * LATTICE_BL;
			temp.y = 1.5 * LATTICE_HIGHT;
			temp.type = firstAtom; // B or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 4
			temp.x = 2.5 * LATTICE_BL;
			temp.y = 1.5 * LATTICE_HIGHT;
			temp.type = secondAtom; // N or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 5
			temp.x = 3 * LATTICE_BL;
			temp.y = LATTICE_HIGHT;
			temp.type = firstAtom; // B or C
			RI += CalculateIntersection(atomMod, temp);
			
			// 6
			temp.x = 2.5 * LATTICE_BL;
			temp.y = 0.5 * LATTICE_HIGHT;
			temp.type = secondAtom; // N or C
			RI += CalculateIntersection(atomMod, temp);		
		}
	}
	else	// xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
	{
		// 1
		temp.x = 2.5 * LATTICE_BL;
		temp.y = -0.5 * LATTICE_HIGHT;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 2
		temp.x = 3 * LATTICE_BL;
		temp.y = 0;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 3
		temp.x = 2.5 * LATTICE_BL;
		temp.y = 0.5 * LATTICE_HIGHT;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 4
		temp.x = 3 * LATTICE_BL;
		temp.y = LATTICE_HIGHT;
		temp.type = firstAtom; // B or C
		RI += CalculateIntersection(atomMod, temp);
		
		// 5
		temp.x = 2.5 * LATTICE_BL;
		temp.y = 1.5 * LATTICE_HIGHT;
		temp.type = secondAtom; // N or C
		RI += CalculateIntersection(atomMod, temp);
	}
	return RI;

}
