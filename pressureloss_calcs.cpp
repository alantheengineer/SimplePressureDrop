#include <math.h>
#include "pressureloss_calcs.h"

double pressure_loss::friction_factor(const double * relativeRoughness, const double * reynoldsNumber)
{
	double ff = 0.002;
	double rn = *reynoldsNumber;
	if (rn < 3000.0)
	{
		if (rn < 700.0)
		{
			ff = 64.0/rn;
		}
		else
		{
			ff = 64.0 / rn + 1.008E-05*(rn - 700);
		}
	}
	else
	{
		for (int i = 0; i < 5; i++)
		{
			double rr = *relativeRoughness;
			double x = log10(0.27*rr + 2.52 / rn / sqrt(ff));
			ff = 0.25 / x / x;
		}
	}

	return ff;
}
