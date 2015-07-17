#ifndef __FINDINTERACTING__
#define __FINDINTERACTING__

#include "Constants_and_libraries.h"

#define SQRT3 (1.7320508075688772)
#define G_HEIGHT (GRAPHENE_BL * SQRT3)
#define BN_HEIGHT (BN_LATTICE_BL * SQRT3)

double FindInteracting(Atom atom, double xShift, double yShift, int latticeType);


#endif
