/* The following line prevents this file from being included twice */
#ifndef __CONSTANTS_AND_LIBRARIES__
#define __CONSTANTS_AND_LIBRARIES__

//*************** Constants ******************
#define NULL	0
#define BL		1.44 			// C-C bond length in Angstroms
#define HIGHT	(BL * sqrt(3))	// The hexagon's hight
#define WIDTH	(BL * 2)		// The hexagon's width
#define HORIZD	(BL * 3)		// The horizontal symmetry distance
#define HORIZS	(BL * 1.5)		// The horizontal shift between atom layers
#define A		(BL * sqrt(3))  // The length of a1 and a2 lattice vectors 
#define NAN		(0/0)			// Nan
#define ILD		3.33			// Interlayer difference
#define RND		(ILD + 1.53)	// Maximum hight that isn't negligible
#define NP		(0.01)			// Negligible part;
#define RCG		(0.5 * BL)		// Radius of graphene atom (for RI calculations)
#define RCC		(0.5 * BL)		// Radius of carbon nanotube atom (for RI Calculations)
#define EXPNORM	(4.605170186)	// ABS(log(0.01)) (log = ln)
#define WRITE	"w"				// Writing type
#define PT		0				// Print title for data files? (1 = yes, 0 = no)
//**************** Useful *********************
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define ABS(X)   ((X) >= 0  ? (X) : -(X))

//*************** Structures *****************

typedef struct
{
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
} InPar;

#endif /* On the #ifndef __CONSTANTS_AND_LIBRARIES__ */
