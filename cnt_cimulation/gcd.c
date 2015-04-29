#include "Constants_and_libraries.h"
#include "gcd.h"

// input: 2 integers.
// output: their gcd.
// The function assumes both numbers are not negative
// and atleast one of them isn't 0.
int gcd(int n, int m)
{
	int temp;

	if (m > n)
	{
		temp = m;
		m = n;
		n = temp;
	}
	if (m != 0)
	{
		do
		{
			temp = n % m;
			n = m;
			m = temp;
		} while (temp != 0);
	}
	return n;
}
