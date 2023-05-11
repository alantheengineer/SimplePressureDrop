#pragma once
#include <vector>

namespace pressure_loss
{
	double friction_factor(const double* relativeRoughness, const double* reynoldsNumber);
	std::vector<double> pressure_loss(const double* density, const double* velocity, const double* friction_factor);
}
