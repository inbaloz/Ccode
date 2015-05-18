#include "Constants_and_libraries.h"
#include "CalculateTranslational.h"
#include "gcd.h"

// input: Ch the chirality vector.
// output: T the translational vector.

// The function assumes the Ch is defined correctly (0 <= abs(m) <= n
// and they aren't both 0).
aVec CalculateTranslational(aVec Ch)
{
	/*int dR;
	aVec T;

	dR = gcd(2 * Ch.n + Ch.m, 2 * Ch.m + Ch.n);
	T.n =	( (2 * Ch.m) + Ch.n ) / dR;
	T.m = -	( (2 * Ch.n) + Ch.m ) / dR;
	return T;*/

	int n;
	aVec T;
	n = ABS(Ch.m);
	n = gcd(Ch.n,n);
	T.n =	( (2 * Ch.m) + Ch.n ) / n;
	T.m = -	( (2 * Ch.n) + Ch.m ) / n;
	return T;
}
