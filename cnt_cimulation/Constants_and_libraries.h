/* The following line prevents this file from being included twice */
#ifndef __CONSTANTS_AND_LIBRARIES__
#define __CONSTANTS_AND_LIBRARIES__

//*************** Flags **********************
#define USE_GAUSSIAN_INTERSECTION 1 // 1 gaussian intersection, 0 for hard sphere.
#define USE_GLOBAL_ZONE 1			// 1 use a global zone, 0 use individual zones
#define WRITE_ENTIRE_TUBE 1 		// 1 write entire tube to coords, 0 write bottom half

//*************** Constants ******************
#define NULL	0

#define CNT_BL_HOMO		     1.42	        // C-C bond length in Angstroms	for a CNT in a homojunction
#define CNT_BL_HETERO	     1.42	        // C-C bond length in Angstroms	for a CNT in a heterojunction
#define GRAPHENE_BL_HOMO     1.42			// C-C bond length in Angstroms	for graphene in a homojunction
#define GRAPHENE_BL_HETERO   1.42			// C-C bond length in Angstroms	for graphene in a heterojunction	
#define BN_TUBE_BL_HOMO	     1.45			// B-N bond length in Angstroms in a tube in a homojunction
#define BN_TUBE_BL_HETERO    1.45			// B-N bond length in Angstroms in a tube in a heterojunction
#define BN_LATTICE_BL_HOMO	 1.45			// B-N bond length in Angstroms in a lattice in a homojunction
#define BN_LATTICE_BL_HETERO 1.45			// B-N bond length in Angstroms in a lattice in a heterojunction
#define BL_NORMALIZATION_HETERO 1.431		// bond length of a hetero system for normalization

#define RCGRAPHENE_HOMO	    (0.5  * GRAPHENE_BL_HOMO)     // Radius of graphene atom (for RI calculations)
#define RCGRAPHENE_HETERO	(0.5  * GRAPHENE_BL_HETERO)   // Radius of graphene atom (for RI calculations)
#define RCCNT_HOMO		    (0.5  * CNT_BL_HOMO)		  // Radius of carbon nanotube atom (for RI Calculations)
#define RCCNT_HETERO		(0.5  * CNT_BL_HETERO)		  // Radius of carbon nanotube atom (for RI Calculations)
#define RBTUBE_HOMO		    (0.15 * BN_TUBE_BL_HOMO)      // Radius of Boron atom in BN tube in a homojunction
#define RBTUBE_HETERO       (0.2  * BN_TUBE_BL_HETERO)    // Radius of Boron atom in BN tube in a heterojunction
#define RBLATTICE_HOMO	    (0.15 * BN_LATTICE_BL_HOMO)	  // Radius of Boron atom in BN hexagonal lattice
#define RBLATTICE_HETERO	(0.2  * BN_LATTICE_BL_HETERO) // Radius of Boron atom in BN hexagonal lattice
#define RNTUBE_HOMO			(0.5  * BN_TUBE_BL_HOMO)	  // Radius of Nitrogen atom in BN tube
#define RNTUBE_HETERO		(0.4  * BN_TUBE_BL_HETERO)	  // Radius of Nitrogen atom in BN tube
#define RNLATTICE_HOMO	    (0.5  * BN_LATTICE_BL_HOMO)	  // Radius of Nitrogen atom in BN hexagonal lattice
#define RNLATTICE_HETERO    (0.4  * BN_LATTICE_BL_HETERO) // Radius of Nitrogen atom in BN hexagonal lattice

#define LATTICE_HIGHT	(LATTICE_BL * sqrt(3))	// The hexagon's hight in the lattice
#define TUBE_HIGHT		(TUBE_BL * sqrt(3))		// The hexagon's hight in the tube
#define LATTICE_WIDTH	(LATTICE_BL * 2)		// The hexagon's width in the lattice (not in use)
#define TUBE_WIDTH		(TUBE_BL * 2)			// The hexagon's width in the tube (not in use)
#define LATTICE_HORIZD	(LATTICE_BL * 3)		// The horizontal symmetry distance of the lattice
#define TUBE_HORIZD		(TUBE_BL * 3)			// The horizontal symmetry distance of the tube
#define LATTICE_HORIZS	(LATTICE_BL * 1.5)		// The horizontal shift between atom layers of the lattice (not in use)
#define TUBE_HORIZS		(TUBE_BL * 1.5)			// The horizontal shift between atom layers of the tube
#define A				(TUBE_BL * sqrt(3))		// The length of a1 and a2 lattice vectors 

#define CNT_G_ILD		3.33			// Interlayer difference between CNT and graphene lattice in Angstroms // TODO Inbal - I found 3.35 for cc and 3.33 for BN
#define BNT_G_ILD		3.33			// Interlayer difference between BN tube and graphene lattice in Angstroms // TODO Inbal - I found 3.35 for cc and 3.33 for BN
#define BNT_BNL_ILD		3.33			// Interlayer difference in Angstroms // TODO Inbal - I found 3.35 for cc and 3.33 for BN

#define NP		(0.01)			// Negligible part;
#define EXPNORM	(4.605170186)	// ABS(log(0.01)) (log = ln)
#define RND		(ILD + 1.53)

#define WRITE	"w"				// Writing type
#define PT		0				// Print title for data files? (1 = yes, 0 = no)
#define DUMMY_TYPE '\xff'		// dummy type for dummy atoms - atoms we want to earase eventually

//*************** atom types ******************
#define C_type 6
#define N_type 7
#define B_type 5

//**************** Gaussian intersection **********
#define RADIUS_TO_STDEV (0.75)	// sigma = radius * RADIUS_TO_STDEV								   	// The larger it is, the wider the gaussian.
//#define RADIUS_TO_STDEV_C (0.75)
//#define RADIUS_TO_STDEV_N (0.75)
//#define RADIUS_TO_STDEV_B (2.15)
#define GAUSSIAN_AMPLITUDE_C (1)
#define GAUSSIAN_AMPLITUDE_N (1)
#define GAUSSIAN_AMPLITUDE_B (1)

//**************** Useful *********************
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define ABS(X)   ((X) >= 0  ? (X) : -(X))

#define	SIZE_ZONE_1   12
#define	SIZE_ZONE_2   12
#define	SIZE_ZONE_3   14
#define	SIZE_ZONE_4   14
#define	SIZE_ZONE_5   17

#define	SIZE_GLOBAL_ZONE   130

//*************** interlayer potential (Itai and Oded's code) *********
#define INTERACTION_BUFFER 4 // in [Angstrom]. Note: The maximal interaction distance 
							 //in Itai's code is 16 [A], I added a buffer.

//*************** Structures *****************

typedef struct
{
	int type; // Carbon, Boron, Nitrogen
	// Coordinates:
	double x;
	double y;
	double z;
} Atom;

typedef struct
{
	int n;	// a1 component
	int m;	// a2 component
} aVec;

typedef struct
{
	aVec Ch;				// Chirality vector
	int	unitcellN;			// Amount of tube unitcells
	int motionType;			// Type of motion
	double shiftAngle;		// The angle between tube axis and y axis (type 1, 3, 4)
	double rotateAngle;		// The angle of rotation of the tube around his axis (type 2, 3)
	double xShift;			// Tube's x axis shift 		(type 1, 2)
	double yShift;			// Tube's y axis shift 		(type 1, 2)
	double rotSpinStart;	// Starting rotation angle	(type 1, 2, 4)
	double rotSpinEnd;		// Ending rotation angle	(type 1, 2, 4)
	int amountOfSteps;		// Amount of angle steps	(type 1, 2, 3, 4)
	double xStart;			// Starting x.				(type 3, 4)
	double yStart;			// Starting y.				(type 3, 4)
	double xEnd;			// Ending x.				(type 3)
	double yEnd;			// Ending y.				(type 3)
	double percentTruncated;// amount of Truncated tube from the end, in a right-angle triangular fashion.
	int tubeType;			// tube type: 	 0 for CNT, 1 for BN.
	int latticeType;		// lattice type: 0 for graphene, 1 for BN.
} InPar;


//**************** Globals *********************
// tube BL, lattice BL
extern double LATTICE_BL;
extern double TUBE_BL;
extern double ILD;
extern double MAX_HEIGHT;
extern double RCGRAPHENE;
extern double RCCNT;
extern double RBTUBE;
extern double RBLATTICE;
extern double RNTUBE;
extern double RNLATTICE;
extern int TUBETYPE;
extern int LATTICETYPE;


extern Atom ZONE_1[SIZE_ZONE_1];
extern Atom ZONE_2[SIZE_ZONE_2];
extern Atom ZONE_3[SIZE_ZONE_3];
extern Atom ZONE_4[SIZE_ZONE_4];
extern Atom ZONE_5[SIZE_ZONE_5];

extern Atom GlobalZone[SIZE_GLOBAL_ZONE];


#endif /* On the #ifndef __CONSTANTS_AND_LIBRARIES__ */