#include "Constants_and_libraries.h"
#include "CreateGlobalZone.h"
#include <math.h>
#include <stdio.h>

Atom GlobalZone[SIZE_GLOBAL_ZONE];

void CreateGlobalZone(int latticeType)
{
	if (latticeType == 0) {
		//col 1
		GlobalZone[0]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[1]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[2]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[3]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[4]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[5]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 2
		GlobalZone[6]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[7]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[8]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[9]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[10]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[11]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
			
			//col 3
		GlobalZone[12]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[13]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[14]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[15]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[16]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[17]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[18]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 4
		GlobalZone[19]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[20]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[21]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[22]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[23]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[24]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[25]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 5
		GlobalZone[26]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[27]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[28]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[29]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[30]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};  
		GlobalZone[31]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 

			//col 6
		GlobalZone[32]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[33]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[34]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[35]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[36]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[37]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 7
		GlobalZone[38]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[39]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[40]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[41]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[42]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[43]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[44]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 8
		GlobalZone[45]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[46]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[47]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[48]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[49]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[50]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[51]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 9
		GlobalZone[52]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[53]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[54]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[55]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[56]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};  
		GlobalZone[57]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 

			//col 10
		GlobalZone[58]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[59]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[60]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[61]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[62]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};  
		GlobalZone[63]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 

			//col 11
		GlobalZone[64]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[65]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[66]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[67]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[68]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[69]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[70]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 12
		GlobalZone[71]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[72]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[73]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[74]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[75]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[76]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[77]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 13
		GlobalZone[78]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[79]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[80]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[81]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[82]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[83]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 14
		GlobalZone[84]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[85]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[86]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[87]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[88]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[89]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 15
		GlobalZone[90]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[91]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[92]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[93]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[94]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[95]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[96]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 16
		GlobalZone[97]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =   3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[98]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =   2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[99]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =   1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[100]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[101]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[102]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[103]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 17
		GlobalZone[104]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[105]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[106]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[107]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[108]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[109]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 18
		GlobalZone[110]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  3 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[111]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  2 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[112]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  1 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[113]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  0 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[114]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -1 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[115]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -2 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 19
		GlobalZone[116]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[117]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[118]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[119]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[120]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[121]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[122]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};

			//col 20
		GlobalZone[123]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[124]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[125]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[126]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[127]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = C_type}; 
		GlobalZone[128]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};
		GlobalZone[129]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = C_type};	
	} else {
		//col 1
		GlobalZone[0]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[1]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[2]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[3]  = (Atom){ .x = -6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[4]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[5]  = (Atom){ .x = -6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 2
		GlobalZone[6]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[7]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[8]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[9]  = (Atom){ .x = -5.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[10]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[11]  = (Atom){ .x = -5.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
			
			//col 3
		GlobalZone[12]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[13]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[14]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[15]  = (Atom){ .x = -4.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[16]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[17]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[18]  = (Atom){ .x = -4.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 4
		GlobalZone[19]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[20]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[21]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[22]  = (Atom){ .x = -3.5 * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[23]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[24]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[25]  = (Atom){ .x = -3.5 * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 5
		GlobalZone[26]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[27]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[28]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[29]  = (Atom){ .x = -3.0 * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[30]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};  
		GlobalZone[31]  = (Atom){ .x = -3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 

			//col 6
		GlobalZone[32]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[33]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[34]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[35]  = (Atom){ .x = -2.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[36]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[37]  = (Atom){ .x = -2.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 7
		GlobalZone[38]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[39]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[40]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[41]  = (Atom){ .x = -1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[42]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[43]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[44]  = (Atom){ .x = -1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 8
		GlobalZone[45]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[46] = (Atom){ .x = -0.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[47]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[48]  = (Atom){ .x = -0.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[49]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[50]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[51]  = (Atom){ .x = -0.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 9
		GlobalZone[52]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[53]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[54]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[55]  = (Atom){ .x = 0.0  * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[56]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};  
		GlobalZone[57]  = (Atom){ .x = 0.0  * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 

			//col 10
		GlobalZone[58]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  3.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[59]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  2.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[60]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  1.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[61]  = (Atom){ .x = 1.0 * LATTICE_BL, .y =  0.0  * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[62]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};  
		GlobalZone[63]  = (Atom){ .x = 1.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 

			//col 11
		GlobalZone[64]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[65]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[66]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[67]  = (Atom){ .x = 1.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[68]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[69]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[70]  = (Atom){ .x = 1.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 12
		GlobalZone[71]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[72]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[73]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[74]  = (Atom){ .x = 2.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[75]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[76]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[77]  = (Atom){ .x = 2.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 13
		GlobalZone[78]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[79]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[80]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[81]  = (Atom){ .x = 3.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[82]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[83]  = (Atom){ .x = 3.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 14
		GlobalZone[84]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[85]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[86]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[87]  = (Atom){ .x = 4.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[88]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[89]  = (Atom){ .x = 4.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 15
		GlobalZone[90]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[91]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[92]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[93]  = (Atom){ .x = 4.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[94]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[95]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[96]  = (Atom){ .x = 4.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 16
		GlobalZone[97]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[98]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[99]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[100]  = (Atom){ .x = 5.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[101]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[102]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[103]  = (Atom){ .x = 5.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 17
		GlobalZone[104]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  3.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[105]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[106]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[107]  = (Atom){ .x = 6.0 * LATTICE_BL, .y =  0.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[108]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -1.0 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[109]  = (Atom){ .x = 6.0 * LATTICE_BL, .y = -2.0 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 18
		GlobalZone[110]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  3  * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[111]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  2  * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[112]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  1  * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[113]  = (Atom){ .x = 7.0 * LATTICE_BL, .y =  0  * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[114]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -1 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[115]  = (Atom){ .x = 7.0 * LATTICE_BL, .y = -2 * LATTICE_HIGHT, .z = 0 , .type = N_type};

			//col 19
		GlobalZone[116]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[117]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[118]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[119]  = (Atom){ .x = 7.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[120]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = B_type}; 
		GlobalZone[121]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};
		GlobalZone[122]  = (Atom){ .x = 7.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = B_type};

			//col 20
		GlobalZone[123]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  3.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[124]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[125]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[126]  = (Atom){ .x = 8.5  * LATTICE_BL, .y =  0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[127]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -0.5 * LATTICE_HIGHT, .z = 0 , .type = N_type}; 
		GlobalZone[128]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -1.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		GlobalZone[129]  = (Atom){ .x = 8.5  * LATTICE_BL, .y = -2.5 * LATTICE_HIGHT, .z = 0 , .type = N_type};
		
	}
}	