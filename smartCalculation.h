#ifndef SMARTCALCUTATIONH
#define SMARTCALCUTATIONH

#include "commonFunctions.h"
#include "MagneticFieldDependence.h"
#include <deque>

/*
Класс для работы адаптации фильтра
*/



enum parameterType
{
	ELECTRON_MOBILITY=0,
	ELECTRON_DENSITY=1,
	LIGHT_HOLE_MOBILITY=2,
	LIGHT_HOLE_DENSITY=3,
	HEAVY_HOLE_MOBILITY=4,
	HEAVY_HOLE_DENSITY=5
};

struct Results
{
	FilterParams fpHall;
	FilterParams fpRes;
	int polinomPowHall;
	int polinomPowMagnetoresistance;
	long double electronConcentration;
	long double electronMobility;
	long double heavyHoleConcentration;
	long double heavyHoleMobility;
	long double lightHoleConcentration;
	long double lightHoleMobility;
	long double targetFuncValue;
	std::vector <PeaksCriteria> vPC;
    AdditionalData additionalData;
    parameterType pt;
};


class smartCalculation
{
public:
	smartCalculation(MagneticFieldDependence * base);
	~smartCalculation();
	bool processData();
	bool saveResults(std::string filename);

	std::vector<long double> getUpBound();
	std::vector<long double> getLowBound();


private:

	enum vPCType
	{
		HHOLE=0,
		LHOLE=1,
		ELECTRON=2,
		LOGHHOLE=3,
		LOGLHOLE=4,
		LOGELECTRON=5
	};

	size_t numberOfParameters;
	size_t numberOfStoredBestResults;
	long double targetFunction (Results & r, parameterType pt);
	bool isBetterResults(Results & r1, Results & r2); // возвращает true если r1 лучше чем r2
	bool processResults();
    bool isGoodExtremum(Results & r, parameterType pt);

	std::deque < std::deque<Results> > results;

	std::vector < Results > allSpectras;
	MagneticFieldDependence * MFDData;


	std::vector<long double> upBound;
	std::vector<long double> lowBound;


};

#endif
