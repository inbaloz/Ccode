#include "Constants_and_libraries.h"
#include "GaussianIntersection.h"
#include <math.h>

/*
 * The intersection between two 2d gaussians is calculated using the method from
 * (1) Grant and Pickup, J. Phys. Chem. 1995, 99, 3503-3510
 *	   http://pubs.acs.org/doi/pdf/10.1021/j100011a016		
 * (2) Grant, Gallardo, and Pickup, J. Comp. Chem. 1996, 17, 1653-1666
 * 
 * The gaussians have the general form: Gaussian = p X exp {-alpha X r^2}
 */
double GaussianIntersection(double radius1, double radius2, double distance,
							double gaussianAmplitude1, double gaussianAmplitude2)
{
	double intersection;

	double stdev1 = RADIUS_TO_STDEV * radius1;
	double stdev2 = RADIUS_TO_STDEV * radius2;
	double p1 = gaussianAmplitude1;
	double p2 = gaussianAmplitude2;
	double alpha1 = 1 / (2 * pow(stdev1,2));
	double alpha2 = 1 / (2 * pow(stdev2,2));
	double delta12 = alpha1 + alpha2;
	double k12 = exp(-(alpha1 * alpha2 * pow(distance, 2)) / delta12);

	// --------- calcualte intersection -----------

	intersection = p1 * p2 * k12 * pow((M_PI / delta12), 1);

	return intersection;
}

	//double p1 = 1 / (2 * M_PI * pow(stdev1,2));
	//double p2 = 1 / (2 * M_PI * pow(stdev2,2));
