#include "Constants_and_libraries.h"
#include "CreateGaussianZone.h"
#include <math.h>
#include <stdio.h>

Atom GaussianZone[SIZE_GAUSSIAN_ZONE];

void CreateGaussianZone(int latticeType)
{
	if (latticeType == 0) {
		//col 1
	GaussianZone[0]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[1]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[2]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[3]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[4]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[5]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 2
	GaussianZone[6]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[7]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[8]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[9]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[10]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[11]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		
		//col 3
	GaussianZone[12]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[13]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[14]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[15]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[16]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[17]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[18]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 4
	GaussianZone[19]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[20]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[21]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[22]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[23]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[24]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[25]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 5
	GaussianZone[26]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[27]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[28]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[29]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[30]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};  
	GaussianZone[31]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 

		//col 6
	GaussianZone[32]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[33]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[34]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[35]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[36]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[37]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 7
	GaussianZone[38]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[39]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[40]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[41]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[42]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[43]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[44]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 8
	GaussianZone[45]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[46]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[47]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[48]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[49]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[50]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[51]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 9
	GaussianZone[52]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[53]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[54]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[55]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[56]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};  
	GaussianZone[57]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 

		//col 10
	GaussianZone[58]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[59]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[60]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[61]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[62]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};  
	GaussianZone[63]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 

		//col 11
	GaussianZone[64]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[65]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[66]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[67]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[68]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[69]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[70]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 12
	GaussianZone[71]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[72]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[73]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[74]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[75]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[76]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[77]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 13
	GaussianZone[78]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[79]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[80]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[81]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[82]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[83]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 14
	GaussianZone[84]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[85]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[86]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[87]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[88]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[89]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 15
	GaussianZone[90]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[91]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[92]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[93]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[94]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[95]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[96]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 16
	GaussianZone[97]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =   3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[98]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =   2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[99]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =   1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[100]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[101]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[102]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[103]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 17
	GaussianZone[104]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[105]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[106]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[107]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[108]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[109]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 18
	GaussianZone[110]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  3 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[111]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  2 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[112]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  1 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[113]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[114]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -1 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[115]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -2 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 19
	GaussianZone[116]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[117]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[118]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[119]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[120]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[121]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[122]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

		//col 20
	GaussianZone[123]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[124]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[125]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[126]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[127]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
	GaussianZone[128]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
	GaussianZone[129]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};	
}

	else {
	//col 1
	GaussianZone[0]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[1]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[2]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[3]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[4]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[5]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 2
	GaussianZone[6]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[7]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[8]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[9]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[10]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[11]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		
		//col 3
	GaussianZone[12]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[13]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[14]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[15]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[16]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[17]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[18]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 4
	GaussianZone[19]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[20]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[21]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[22]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[23]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[24]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[25]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 5
	GaussianZone[26]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[27]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[28]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[29]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[30]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};  
	GaussianZone[31]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 

		//col 6
	GaussianZone[32]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[33]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[34]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[35]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[36]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[37]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 7
	GaussianZone[38]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[39]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[40]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[41]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[42]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[43]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[44]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 8
	GaussianZone[45]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[46] = (Atom){ .x = -0.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[47]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[48]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[49]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[50]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[51]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 9
	GaussianZone[52]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[53]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[54]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[55]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[56]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};  
	GaussianZone[57]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 

		//col 10
	GaussianZone[58]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[59]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[60]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[61]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[62]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};  
	GaussianZone[63]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 

		//col 11
	GaussianZone[64]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[65]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[66]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[67]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[68]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[69]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[70]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 12
	GaussianZone[71]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[72]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[73]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[74]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[75]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[76]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[77]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 13
	GaussianZone[78]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[79]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[80]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[81]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[82]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[83]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 14
	GaussianZone[84]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[85]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[86]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[87]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[88]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[89]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 15
	GaussianZone[90]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[91]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[92]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[93]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[94]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[95]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[96]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 16
	GaussianZone[97]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[98]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[99]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[100]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[101]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[102]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[103]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 17
	GaussianZone[104]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[105]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[106]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[107]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[108]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[109]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 18
	GaussianZone[110]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  3  * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[111]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  2  * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[112]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  1  * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[113]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  0  * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[114]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -1 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[115]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -2 * LATTICE_HIGHT, .z = 0 , .type = N_type};

		//col 19
	GaussianZone[116]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[117]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[118]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[119]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[120]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
	GaussianZone[121]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
	GaussianZone[122]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

		//col 20
	GaussianZone[123]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[124]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[125]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[126]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[127]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
	GaussianZone[128]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	GaussianZone[129]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
	
	}
}	