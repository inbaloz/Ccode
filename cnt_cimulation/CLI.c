#include "Constants_and_libraries.h"
#include <stdio.h>
#include <math.h>
#include "CLI.h"

// input: Input parameters.

// output: Input parameters after changes made by the user.

InPar CLI(InPar input)
{
	int counter;				// Counter..
	int inputNum, inputValue;	// User's input parameter number and it's value
	double inputValueD;			// User's input parameter value (double)
	for (counter = 0; counter < 100; counter++) // Just to make this loop doesn't run forever...
	{
		printf("The current parameters are:\n");
		printf("1.  Ch.n = %d\n", input.Ch.n);
		printf("2.  Ch.m = %d\n", input.Ch.m);
		printf("3.  Amount of tube unitcells = %d\n", input.unitcellN);
		printf("4.  Motion type = %d\n", input.motionType);
		printf("   (1. Rotating, 2. Spinning (in place), 3. Sliding, 4. Spinning without sliding)\n");
		printf("5.  Shift angle (degrees) = %e\n", 180 * input.shiftAngle / M_PI);
		printf("6.  Rotate angle (degress) = %e\n", 180 * input.rotateAngle / M_PI);
		printf("7.  x shift = %e\n", input.xShift);
		printf("8.  y shift = %e\n", input.yShift);
		printf("9.  Rotation or Spinning start angle (degrees) = %e\n", 180 * input.rotSpinStart / M_PI);
		printf("10. Rotation or Spinning end angle (degrees) = %e\n", 180 * input.rotSpinEnd / M_PI);
		printf("11. Amount of steps = %d\n", input.amountOfSteps);
		printf("12. x start = %e\n", input.xStart);
		printf("13. y start = %e\n", input.yStart);
		printf("14. x end = %e\n", input.xEnd);
		printf("15. y end = %e\n", input.yEnd);
		printf("To change the value of a paramter enter his number and then the requested value.\n");
		printf("To finish changing and continue enter -1.\n");
		scanf("%d", &inputNum);
		if (inputNum == -1)
		{
			return input;
		}
		scanf("%lf", &inputValueD);
		switch(inputNum)
		{
		case 1:
			input.Ch.n = (int)inputValueD;
			break;
		case 2:
			input.Ch.m = (int)inputValueD;
			break;
		case 3:
			input.unitcellN = (int)inputValueD;
			break;
		case 4:
			input.motionType = (int)inputValueD;
			break;
		case 5:
			input.shiftAngle = M_PI * inputValueD / 180;
			break;
		case 6:
			input.rotateAngle = M_PI * inputValueD / 180;
			break;
		case 7:
			input.xShift = inputValueD;
			break;
		case 8:
			input.yShift = inputValueD;
			break;
		case 9:
			input.rotSpinStart = M_PI * inputValueD / 180;
			break;
		case 10:
			input.rotSpinEnd = M_PI * inputValueD / 180;
			break;
		case 11:
			input.amountOfSteps = (int)inputValueD;
			break;
		case 12:
			input.xStart = inputValueD;
			break;
		case 13:
			input.yStart = inputValueD;
			break;
		case 14:
			input.xEnd = inputValueD;
			break;
		case 15:
			input.yEnd = inputValueD;
			break;
		default:
			printf("Not a legal parameter number.\n");
			break;
		}
	}
	return input;
}
