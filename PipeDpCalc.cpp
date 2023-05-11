// PipeDpCalc.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <math.h>
#include <string>
#include "enhancedwater.h"
#include "pressureloss_calcs.h"

using namespace std;
using namespace pressure_loss;
using namespace fluid_properties;

int main(int argc, char* argv[])
{
	double dPressureIn = 101.325348872000; // kPa
	double dTemperature = 21.0; // �C
	double dSalinity = 0.0; // ppm
	double dHumidity = 0.0; // %
	double dDiameter = 90.0; // mm
	double dAngleInclination = 0.0; // radians
	double dFlowRate = 10.0; // l/s
	double dRoughness = 1.5E-05; // m
	int iFluid = 1; // 0 - 'water', 1 - 'air'

	if (argc > 1)
	{
		if (string(argv[1]) == "--help")
		{
			cout << "Simple pipe pressure drop calculator:\n"
				<< "Requires at least a rate (-q) as input\n"
				<< "If only rate is entered assumes default values for all else\n"
				<< "Arguments that can be provided:\n"
				<< "\t-f,--fluid: 'water' or 'air', default = 'air'\n"
				<< "\t\t-s,--salinity: input if using water salinity, default = 0 (ppm)\n"
				<< "\t\t-h,--humidity: input if using relative humidity (%) default = 0\n\n"
				<< "\t-p,--pressure: input pressure (kPa), default = 101.325 kPa\n"
				<< "\t-t,--temperature: input temperature (celsius), default = 21 degC\n"
				<< "\t-q,--rate: flow rate @ standard conditions (l/s), default = 10.0 l/s\n"
				<< "\t -d,--diameter: input pipe diameter (mm), default = 10.0 mm\n"
				<< "\t -a,--angle: angle from horizontal (radians), default = 0 rad, min -1.57, max +1.57 (horizontal)\n"
				<< "\t -r,--roughness: pipe roughness (m), default = 1.5e-05 m (pvc pipe)"
				<< "This program returns pipeline pressure gradients for the given flow rates\n"
				<< "assuming single phase flow, steadystate flow as well just to be clear."
				<< "Units: density (kg/m^3), viscosity (mPa.s), compressibility (1/MPa)\n"
				<< endl;
			return 0;
		}

		for (int i = 1; i < argc; ++i)
		{
			string arg = argv[i];
			if ((arg == "-d") || (arg == "--diameter"))
			{
				if (i + 1 < argc) dDiameter = atof(argv[i + 1]);
			}
			if ((arg == "-q") || (arg == "--rate"))
			{
				if (i + 1 < argc) dFlowRate = atof(argv[i + 1]);
			}
			if ((arg == "-a") || (arg == "--angle"))
			{
				if (i + 1 < argc) dAngleInclination = atof(argv[i + 1]);
			}
			if ((arg == "-r") || (arg == "--roughness"))
			{
				if (i + 1 < argc) dRoughness = atof(argv[i + 1]);
			}
			if ((arg == "-h") || (arg == "--humidity"))
			{
				if (i + 1 < argc) dHumidity = atof(argv[i + 1]);
			}
			if ((arg == "-p") || (arg == "--pressure"))
			{
				if (i + 1 < argc) dPressureIn = atof(argv[i + 1]);
			}
			if ((arg == "-t") || (arg == "--temperature"))
			{
				if (i + 1 < argc) dTemperature = atof(argv[i + 1]);
			}
			if ((arg == "-s") || (arg == "--salinity"))
			{
				if (i + 1 < argc) dSalinity = atof(argv[i + 1]);
			}
			if ((arg == "-f") || (arg == "--fluid"))
			{
				std::cout << argv[i + 1];
				if (i + 1 < argc)
				{
					string fluid_arg = argv[i + 1];
					if (fluid_arg == "water") iFluid = 0;
				}
			}
		}
	}

	double dPressureMPa = dPressureIn / 1000.0; //MPa

	double dDensity = 0.0;
	double dStdDensity = 0.0;
	double dPStd = 0.101325348872;
	double dTStd = 21.0;
	double dCompressibility = 0.0;
	double dViscosity = 0.0;
	double dViscosity_metric = 0.0;

	double dDiameterMetric = dDiameter*0.001;

	if (iFluid == 0) printf("Water properties at: %g kPa, %g oC, %g ppm\n", dPressureIn, dTemperature, dSalinity);
	if (iFluid == 1) printf("Air properties at: %g kPa, %g �C\n", dPressureIn, dTemperature);

	if (iFluid == 0)
	{
		// std
		iapws95(dPStd, dTStd, dSalinity, &dStdDensity, &dCompressibility);
		// conditions
		iapws95(dPressureMPa, dTemperature, dSalinity, &dDensity, &dCompressibility);
		dViscosity = kestin_viscosity(dPressureMPa, dTemperature, dSalinity);

		printf("Viscosity: %g mPa.s\n", dViscosity);
		printf("Compressibility: %g 1/MPa\n", dCompressibility);
		dViscosity_metric = dViscosity * 0.001;
	}

	if (iFluid == 1)
	{
		// std
		dStdDensity = air_density(dPStd, dTStd, dHumidity);
		// conditions
		dDensity = air_density(dPressureMPa, dTemperature, dHumidity);
		dCompressibility = air_compressibility_factor(dPressureMPa, dTemperature, dHumidity);
		dViscosity = air_viscosity(dPressureMPa, dTemperature, dHumidity);

		printf("Viscosity: %g uPa.s\n", dViscosity);
		printf("Compressibility (Z): %g \n", dCompressibility);
		dViscosity_metric = dViscosity * 0.000001;
	}

	printf("Density: %g kg/m^3\n", dDensity);

	double dMassRate = dStdDensity * (dFlowRate*0.001);
	double dVelocity = (dMassRate / dDensity) / (3.14159265359 * pow(0.5*dDiameterMetric,2));

	double dRelativeRoughness = dRoughness / dDiameterMetric;
	double dReynoldsNumber = (dDensity*dVelocity*dDiameterMetric) / dViscosity_metric;

	double dFF = friction_factor(&dRelativeRoughness, &dReynoldsNumber);

	double dFrictionGradient = dFF*dDensity*pow(dVelocity, 2) / (2 * dDiameterMetric);
	double dGravityHead = dDensity*sin(dAngleInclination) * 9.81;
	double dAccelerationHead = dDensity*dVelocity*dVelocity / (9.81*1.E06*dPressureMPa);

	if (dAccelerationHead > 0.99) dAccelerationHead = 0.99;

	double dAccelerationGradient = (dFrictionGradient + dGravityHead)*dAccelerationHead;

	double dTotalGradient = dGravityHead + dFrictionGradient + dAccelerationGradient;

	printf("Mass flow rate: %g kg/s\n", dMassRate);
	printf("fluid velocity: %g m/s\n", dVelocity);
	printf("Fluid reynolds number: %g (-)\n", dReynoldsNumber);
	printf("friction factor: %g (-)\n\n", dFF);

	printf("static gradient: %g Pa/m\n", dGravityHead);
	printf("friction gradient: %g Pa/m\n", dFrictionGradient);
	printf("acceleration gradient: %g Pa/m\n", dAccelerationGradient);

	printf("total pressure gradient: %g Pa/m\n", dTotalGradient);


    return 0;
}
