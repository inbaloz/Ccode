#include <math.h>
#include "Constants_and_libraries.h"
#include <stdio.h>

double RCCNT;
double RCGRAPHENE;
double LATTICE_BL;
double TUBE_BL;
double ILD;

double RBTUBE;
double RBLATTICE;
double RNTUBE;
double RNLATTICE;

void SetGlobals(int tubeType, int latticeType)
{
	if (tubeType == 0 & latticeType == 0) {
		RCCNT = RCCNT_HOMO;
		RCGRAPHENE = RCGRAPHENE_HOMO;
		ILD        = CNT_G_ILD;
		TUBE_BL    = CNT_BL_HOMO;
		LATTICE_BL = GRAPHENE_BL_HOMO;

		// Zone 1:  xMod >= 0 && xMod < LATTICE_BL 
		ZONE_1[0] = (Atom){ .x = 0    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_1[1] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_1[2] = (Atom){ .x = 0    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_1[3] = (Atom){ .x = 1    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_1[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_1[5] = (Atom){ .x = 1    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_1[6] = (Atom){ .x = -0.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_1[7] = (Atom){ .x = -1.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_1[8] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_1[9] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_1[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_1[11] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12

		// Zone 2: xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)
		ZONE_2[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_2[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_2[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_2[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_2[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_2[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_2[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_2[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_2[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_2[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_2[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_2[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12

		// Zone 3: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod > (0.5 * LATTICE_HIGHT)
		ZONE_3[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_3[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_3[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_3[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_3[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_3[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_3[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_3[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_3[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_3[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_3[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_3[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12
		ZONE_3[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 13
		ZONE_3[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 14
	

		// Zone 4: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod < (0.5 * LATTICE_HIGHT)
		ZONE_4[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_4[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_4[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_4[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_4[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_4[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_4[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_4[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_4[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_4[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_4[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_4[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12
		ZONE_4[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 13
		ZONE_4[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 14
	

		// Zones 5: xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
		ZONE_5[0] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_5[1] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_5[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_5[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_5[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5 
		ZONE_5[5] = (Atom){ .x = 3.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_5[6] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_5[7] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_5[8] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_5[9] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_5[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_5[11] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 2.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12
		ZONE_5[12] = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 13
		ZONE_5[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 14
		ZONE_5[14] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 15
		ZONE_5[15] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 16
		ZONE_5[16] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 17

	}
	else if (tubeType == 0 & latticeType == 1) {
		RCCNT      = RCCNT_HETERO;
		RBLATTICE  = RBLATTICE_HETERO;
		RNLATTICE  = RNLATTICE_HETERO;
		ILD 	   = BNT_G_ILD;
		TUBE_BL    = CNT_BL_HETERO;
		LATTICE_BL = BN_LATTICE_BL_HETERO;

		// Zone 1:  xMod >= 0 && xMod < LATTICE_BL 
		ZONE_1[0] = (Atom){ .x = 0    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 1
		ZONE_1[1] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 2
		ZONE_1[2] = (Atom){ .x = 0    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 3
		ZONE_1[3] = (Atom){ .x = 1    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 4
		ZONE_1[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 5
		ZONE_1[5] = (Atom){ .x = 1    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 6
		ZONE_1[6] = (Atom){ .x = -0.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 7
		ZONE_1[7] = (Atom){ .x = -1.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 8
		ZONE_1[8] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 9
		ZONE_1[9] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 10
		ZONE_1[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 11
		ZONE_1[11] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 12


		// Zone 2: xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)
		ZONE_2[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_2[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 2
		ZONE_2[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_2[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_2[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5
		ZONE_2[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 6
		ZONE_2[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 7
		ZONE_2[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 8
		ZONE_2[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 9
		ZONE_2[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 10
		ZONE_2[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 11
		ZONE_2[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 12

		// Zones 3 and 4: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod > or < (0.5 * LATTICE_HIGHT)
		ZONE_3[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_3[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 2
		ZONE_3[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_3[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_3[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5
		ZONE_3[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 6
		ZONE_3[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 7
		ZONE_3[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 8
		ZONE_3[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 9
		ZONE_3[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 10
		ZONE_3[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 11
		ZONE_3[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 12
		ZONE_3[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 13
		ZONE_3[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 14

		ZONE_4[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_4[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 2
		ZONE_4[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_4[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_4[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5
		ZONE_4[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 6
		ZONE_4[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 7
		ZONE_4[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 8
		ZONE_4[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 9
		ZONE_4[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 10
		ZONE_4[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 11
		ZONE_4[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 12
		ZONE_4[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 13
		ZONE_4[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 14

		// Zones 5: xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
		ZONE_5[0] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_5[1] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 2
		ZONE_5[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_5[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_5[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5 
		ZONE_5[5] = (Atom){ .x = 3.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 6
		ZONE_5[6] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 7
		ZONE_5[7] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 8
		ZONE_5[8] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 9
		ZONE_5[9] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 10
		ZONE_5[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 11
		ZONE_5[11] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 2.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 12
		ZONE_5[12] = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 13
		ZONE_5[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 14
		ZONE_5[14] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 15
		ZONE_5[15] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 16
		ZONE_5[16] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 17
	}
	else if (tubeType == 1 & latticeType == 0) {
		RNTUBE     = RNTUBE_HETERO;
		RBTUBE     = RBTUBE_HETERO;
		RCGRAPHENE = RCCNT_HETERO;
		ILD        = BNT_G_ILD;
		TUBE_BL    = BN_TUBE_BL_HETERO;
		LATTICE_BL = GRAPHENE_BL_HETERO;
		// Zone 1:  xMod >= 0 && xMod < LATTICE_BL 
		ZONE_1[0] = (Atom){ .x = 0    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_1[1] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_1[2] = (Atom){ .x = 0    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_1[3] = (Atom){ .x = 1    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_1[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_1[5] = (Atom){ .x = 1    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_1[6] = (Atom){ .x = -0.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_1[7] = (Atom){ .x = -1.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_1[8] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_1[9] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_1[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_1[11] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12

		printf("IN GLOBALS: zone 1 atom 1: %lf, %d\n", ZONE_1[1].x, ZONE_1[0].type);

		// Zone 2: xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)
		ZONE_2[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_2[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_2[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_2[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_2[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_2[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_2[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_2[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_2[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_2[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_2[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_2[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12

		// Zone 3: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod > (0.5 * LATTICE_HIGHT)
		ZONE_3[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_3[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_3[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_3[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_3[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_3[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_3[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_3[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_3[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_3[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_3[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_3[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12
		ZONE_3[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 13
		ZONE_3[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 14
	

		// Zone 4: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod < (0.5 * LATTICE_HIGHT)
		ZONE_4[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_4[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_4[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_4[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_4[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5
		ZONE_4[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_4[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_4[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_4[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_4[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_4[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_4[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12
		ZONE_4[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 13
		ZONE_4[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 14
	

		// Zones 5: xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
		ZONE_5[0] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 1
		ZONE_5[1] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 2
		ZONE_5[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 3
		ZONE_5[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 4
		ZONE_5[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 5 
		ZONE_5[5] = (Atom){ .x = 3.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 6
		ZONE_5[6] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 7
		ZONE_5[7] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 8
		ZONE_5[8] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 9
		ZONE_5[9] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 10
		ZONE_5[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 11
		ZONE_5[11] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 2.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 12
		ZONE_5[12] = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 13
		ZONE_5[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 14
		ZONE_5[14] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 15
		ZONE_5[15] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 16
		ZONE_5[16] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = C_type}; // 17
	}
	else {
		RNTUBE     = RNTUBE_HOMO;
		RBTUBE     = RBTUBE_HOMO;
		RBLATTICE  = RBLATTICE_HOMO;
		RNLATTICE  = RNLATTICE_HOMO;
		ILD        = BNT_BNL_ILD;
		TUBE_BL    = BN_TUBE_BL_HOMO;
		LATTICE_BL = BN_LATTICE_BL_HOMO;
		// Zone 1:  xMod >= 0 && xMod < LATTICE_BL 
		ZONE_1[0] = (Atom){ .x = 0    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 1
		ZONE_1[1] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 2
		ZONE_1[2] = (Atom){ .x = 0    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 3
		ZONE_1[3] = (Atom){ .x = 1    * LATTICE_BL, .y = 1    * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 4
		ZONE_1[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 5
		ZONE_1[5] = (Atom){ .x = 1    * LATTICE_BL, .y = 0    * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 6
		ZONE_1[6] = (Atom){ .x = -0.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 7
		ZONE_1[7] = (Atom){ .x = -1.5 * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 8
		ZONE_1[8] = (Atom){ .x = -0.5 * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 9
		ZONE_1[9] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 10
		ZONE_1[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type};  // 11
		ZONE_1[11] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};  // 12


		// Zone 2: xMod >= LATTICE_BL && xMod < (1.5 * LATTICE_BL)
		ZONE_2[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_2[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 2
		ZONE_2[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_2[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_2[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5
		ZONE_2[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 6
		ZONE_2[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 7
		ZONE_2[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 8
		ZONE_2[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 9
		ZONE_2[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 10
		ZONE_2[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 11
		ZONE_2[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 12

		// Zones 3 and 4: xMod >= (1.5 * LATTICE_BL) && xMod < (2.5 * LATTICE_BL), yMod > or < (0.5 * LATTICE_HIGHT)
		ZONE_3[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_3[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 2
		ZONE_3[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_3[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_3[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5
		ZONE_3[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 6
		ZONE_3[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 7
		ZONE_3[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 8
		ZONE_3[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 9
		ZONE_3[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 10
		ZONE_3[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 11
		ZONE_3[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 12
		ZONE_3[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 13
		ZONE_3[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 14

		ZONE_4[0] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_4[1] = (Atom){ .x = 0.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 2
		ZONE_4[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_4[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_4[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5
		ZONE_4[5] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 6
		ZONE_4[6] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 7
		ZONE_4[7] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 8
		ZONE_4[8] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 9
		ZONE_4[9] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 10
		ZONE_4[10] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 11
		ZONE_4[11] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 12
		ZONE_4[12] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 13
		ZONE_4[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 14

		// Zones 5: xMod >= (2.5 * LATTICE_BL) && xMod < (3 * LATTICE_BL)
		ZONE_5[0] = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 1
		ZONE_5[1] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 2
		ZONE_5[2] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 3
		ZONE_5[3] = (Atom){ .x = 1.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 4
		ZONE_5[4] = (Atom){ .x = 1.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 5 
		ZONE_5[5] = (Atom){ .x = 3.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 6
		ZONE_5[6] = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 7
		ZONE_5[7] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 8
		ZONE_5[8] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 9
		ZONE_5[9] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 10
		ZONE_5[10] = (Atom){ .x = 2.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 11
		ZONE_5[11] = (Atom){ .x = 3.0  * LATTICE_BL, .y = 2.0  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 12
		ZONE_5[12] = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 13
		ZONE_5[13] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 14
		ZONE_5[14] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 0.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 15
		ZONE_5[15] = (Atom){ .x = 4.0  * LATTICE_BL, .y = 1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; // 16
		ZONE_5[16] = (Atom){ .x = 4.5  * LATTICE_BL, .y = 1.5  * LATTICE_HIGHT, .z = 0 , .type = B_type}; // 17
	}


}
