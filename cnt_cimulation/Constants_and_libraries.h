/* The following line prevents this file from being included twice */
#ifndef __CONSTANTS_AND_LIBRARIES__
#define __CONSTANTS_AND_LIBRARIES__

//*************** Constants ******************
#define NULL	0

// TODO: paper 34
#define CC_BL		1.44				// C-C bond length in Angstroms			
#define BN_BL		1.45				// B-N bond length in Angstroms
#define RBTUBE		(0.15 * BN_BL)		// Radius of Boron atom in BN tube
#define RNTUBE		(0.5 * BN_BL)		// Radius of Nitrogen atom in BN tube
#define RBLATTICE	(0.15 * BN_BL)		// Radius of Boron atom in BN hexagonal lattice
#define RNLATTICE	(0.5 * BN_BL)		// Radius of Nitrogen atom in BN hexagonal lattice
#define RCGRAPHENE	(0.5 * CC_BL)		// Radius of graphene atom (for RI calculations)
#define RCCNT		(0.5 * CC_BL)		// Radius of carbon nanotube atom (for RI Calculations)

#define LATTICE_HIGHT	(LATTICE_BL * sqrt(3))	// The hexagon's hight in the lattice
#define TUBE_HIGHT		(TUBE_BL * sqrt(3))		// The hexagon's hight in the tube
#define LATICE_WIDTH	(LATTICE_BL * 2)		// The hexagon's width in the lattice
#define TUBE_WIDTH		(TUBE_BL * 2)			// The hexagon's width in the tube
#define LATTICE_HORIZD	(LATTICE_BL * 3)		// The horizontal symmetry distance of the lattice
#define TUBE_HORIZD		(TUBE_BL * 3)			// The horizontal symmetry distance of the tube
#define LATTICE_HORIZS	(LATTICE_BL * 1.5)		// The horizontal shift between atom layers of the lattice
#define TUBE_HORIZS		(TUBE_BL * 1.5)			// The horizontal shift between atom layers of the tube
#define A				(TUBE_BL * sqrt(3))		// The length of a1 and a2 lattice vectors 

// TODO: paper 28
#define CNT_G_ILD		3.33			// Interlayer difference between CNT and graphene lattice in Angstroms // TODO Inbal - I found 3.35 for cc and 3.33 for BN
#define BNT_G_ILD		3.33			// Interlayer difference between BN tube and graphene lattice in Angstroms // TODO Inbal - I found 3.35 for cc and 3.33 for BN
#define BNT_BNL_ILD		3.33			// Interlayer difference in Angstroms // TODO Inbal - I found 3.35 for cc and 3.33 for BN
#define RND				(ILD + 1.53)	// Maximum hight that isn't negligible

#define NAN		(0/0)			// Nan
#define NP		(0.01)			// Negligible part;
#define EXPNORM	(4.605170186)	// ABS(log(0.01)) (log = ln)

#define WRITE	"w"				// Writing type
#define PT		0				// Print title for data files? (1 = yes, 0 = no)
#define DUMMY_TYPE '\xff'		// dummy type for dummy atoms - atoms we want to earase eventually

//**************** Useful *********************
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define ABS(X)   ((X) >= 0  ? (X) : -(X))

//**************** Globals *********************
// tube BL, lattice BL
extern double LATTICE_BL;
extern double TUBE_BL;
extern double ILD;

//*************** Structures *****************

typedef struct
{
	char type; // 'C'arbon, 'B'oron, 'N'itrogen
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

#endif /* On the #ifndef __CONSTANTS_AND_LIBRARIES__ */
