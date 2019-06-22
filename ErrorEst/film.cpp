//---------------------------------------------------------------------------

#include "film.h"
#include <cmath>

film::film(long double molarCompositionCadmium,
long double Temperature,long double heavyHoleConcentration,long double aFactor,long double kFactor,
long double thickness,long double cbRatio,long double currentIntensity,
long double numberOfCarrierTypes):NumberOfCarrierTypes(numberOfCarrierTypes)
{
	Concentration=new long double[NumberOfCarrierTypes];
	Mobility=new long double[NumberOfCarrierTypes];

	MolarCompositionCadmium=molarCompositionCadmium;
	CurrentTemperature=Temperature;
	Concentration[0]=heavyHoleConcentration;
	AFactor=aFactor;
	KFactor=kFactor;
	Thickness=thickness;
	CBRatio=cbRatio;
    CurrentIntensity=currentIntensity;

	calcutatingCarrierParams();
}

film::~film()
{
	delete[] Concentration;
	Concentration=NULL;
	delete[] Mobility;
	Mobility=NULL;
}

void film::calcutatingCarrierParams()
{
	long double m0=9.10938188E-31; // масса электрона хотя она и не нужна в общем-то:)
	long double mh=0.443*m0; // масса тяжелых дырок
	long double ml=0.001*m0; // масса легких дырок

	Concentration[1]=pow(ml/mh,1.5)*Concentration[0]*50; // концентрация легких дырок домножена на 50
	Concentration[2]=-pow(ownConcentrationOfElectrons(),2)/Concentration[0];

	Mobility[2]=-AFactor*pow(CurrentTemperature/77.0,-KFactor);
	Mobility[1]=-0.1*Mobility[2];
	Mobility[0]=-0.005*Mobility[2];

}

long double film::getNumberOfCarrierTypes()
{
	return NumberOfCarrierTypes;
}

long double film::getConcentration(int index)
{
	if(index>=0 && index<NumberOfCarrierTypes)
		return Concentration[index];
	else
		return 0;
}

long double film::getMobility(int index)
{
	if(index>=0 && index<NumberOfCarrierTypes)
		return Mobility[index];
	else
		return 0;
}

long double film::getThickness()
{
	return Thickness;
}

long double film::getCBRatio()
{
    return CBRatio;
}

long double film::getCurrentIntensity()
{
    return CurrentIntensity;
}

long double film::getAFactor()
{
	return AFactor;

}

long double film::getKFactor()
{
	return KFactor;
}

long double film::getMolarCompositionCadmium()
{
	return MolarCompositionCadmium;
}

long double film::getCurrentTemperature()
{
    return CurrentTemperature;
}

long double film::ownConcentrationOfElectrons()
{
	long double T=CurrentTemperature;
	long double x=MolarCompositionCadmium;
	long double k=1.380648813E-23/1.60217646E-19; // постоянная больцмана в электрон-Вольтах
    long double Eg=-0.302+1.93*x-0.810*x*x+0.832*x*x*x+5.35E-4*(1-2.0*x)*T;
    return  (5.585-3.820*x+1.753E-3*T-1.364E-3*x*T)*1E20*powl(Eg,3/4.)*powl(T,1.5)*exp(-Eg/2./k/T); // собственная концентрация
    // Вопрос, в формуле с журнала 1e14 вместо 1e20
}
