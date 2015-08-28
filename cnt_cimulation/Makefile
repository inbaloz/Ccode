C_FILES = \
	ArrayToFile.c \
	aVecLength.c \
	CalculateIntersection.c \
	CalculateTranslational.c \
	CLI.c \
	CutUnitcell.c \
	DuplicateTube.c \
	CutLastPartOfTube.c \
	CreateSurface.c \
	FixGapsInTube.c \
	EffectiveNum.c \
	FindInteracting.c \
	HardSphereIntersection.c \
	gcd.c \
	LatticeCreator.c \
	LoadInPar.c \
	GaussianIntersection.c \
	main.c \
	Move.c \
	NormRI.c \
	PerfectRotationMotion.c \
	Rad2Deg.c \
	Rotate.c \
	RotateShift.c \
	RotationMotion.c \
	SaveInPar.c \
	SlidingMotion.c \
	SpinningMotion.c \
	TwodDataToFile.c \
	AtomsToFile.c \
	SetGlobals.c \
	WriteCoordinates.c \

H_FILES = \
	ArrayToFile.h \
	aVecLength.h \
	CalculateIntersection.h \
	CalculateTranslational.h \
	CLI.h \
	Constants_and_libraries.h \
	CutUnitcell.h \
	DuplicateTube.h \
	CutLastPartOfTube.h \
	CreateSurface.h \
	HardSphereIntersection.h \
	FixGapsInTube.h \
	EffectiveNum.h \
	FindInteracting.h \
	gcd.h \
	LatticeCreator.h \
	LoadInPar.h \
	GaussianIntersection.h \
	Move.h \
	NormRI.h \
	PerfectRotationMotion.h \
	Rad2Deg.h \
	Rotate.h \
	RotateShift.h \
	RotationMotion.h \
	SaveInPar.h \
	SlidingMotion.h \
	SpinningMotion.h \
	TwodDataToFile.h \
	AtomsToFile.h \
	SetGlobals.h \
	WriteCoordinates.h \

O_FILES = ${C_FILES:.c=.o} # Take the .c files and replace with .o

all: cnt-simulation

cnt-simulation: $(O_FILES)
	gcc $(O_FILES) -lm -o cnt-simulation

# This rule means that for every file.o, use the file.c
# In the actual rule below, $< is the file.c
.c.o:
	gcc -g -c $<

# Remove everything which is not the original source code
clean:
	rm -f $(O_FILES)
	rm -f cnt-simulation
