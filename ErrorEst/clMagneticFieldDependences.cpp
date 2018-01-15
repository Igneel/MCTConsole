//---------------------------------------------------------------------------

#include "clMagneticFieldDependences.h"


long double clMagneticFieldDependences::getConcentration(int n)
{
	return carrierParams->getConcentration(n);
}

long double clMagneticFieldDependences::getMobility(int n)
{
	return carrierParams->getMobility(n);
}

void clMagneticFieldDependences::MemoryAllocation()
{
	sxx.resize(NumberOfPoints);
	sxy.resize(NumberOfPoints);
	B.resize(NumberOfPoints);
	Us.resize(NumberOfPoints);
	Uy.resize(NumberOfPoints);
	s_eff.resize(NumberOfPoints);
	Rh_eff.resize(NumberOfPoints);

}

const long double clMagneticFieldDependences::THEALMOSTZERO=0.00000000001;
const long double clMagneticFieldDependences::electronCharge=1.60217646E-19;

clMagneticFieldDependences::clMagneticFieldDependences(int size, long double shag,
long double molarCompositionCadmium,
long double Temperature, long double heavyHoleConcentrerion,
long double AFactor, long double KFactor,
long double thickness, long double cbRatio,
long double currentIntensity, long double numberOfCarrierTypes)
:NumberOfPoints(size)
{
	// инициализация пленки.
    carrierParams = new film(molarCompositionCadmium,Temperature,heavyHoleConcentrerion,
	AFactor,KFactor,thickness,cbRatio,currentIntensity,numberOfCarrierTypes);

	h=shag;
	MemoryAllocation();

	calculateMagneticFieldPoints();
}

clMagneticFieldDependences::clMagneticFieldDependences(int size,long double shag,film * cp)
:NumberOfPoints(size)
{
/*
Вот загрузили мы данные в плёнку, заполнили массив магнитного поля.
*/
	h=shag;
	carrierParams = new film(cp->getMolarCompositionCadmium(),
	cp->getCurrentTemperature(),cp->getConcentration(0),
	cp->getAFactor(),cp->getKFactor(),cp->getThickness(),
	cp->getCBRatio(),cp->getCurrentIntensity(),cp->getNumberOfCarrierTypes());
	MemoryAllocation();

	calculateMagneticFieldPoints();
}

clMagneticFieldDependences::~clMagneticFieldDependences()
{
	delete carrierParams;
}

void clMagneticFieldDependences::calculateMagneticFieldPoints()
{
	//B.push_back(0);
	B[0]=0;
	for(int i=1;i<NumberOfPoints;i++)
	{
		//B.push_back(B[i-1]+h);// это наше магнитное поле
		B[i]=B[i-1]+h;
	}

}

void clMagneticFieldDependences::calculateTenzorFromCarrierParams()
{
	// вызывается автоматически из calculateMagneticFieldDependences()
	for(int i=0;i<NumberOfPoints;i++)
	{
		sxx[i]=sxy[i]=0;
		int NumberOfCarrierTypes=carrierParams->getNumberOfCarrierTypes();

		for(int j=0;j<NumberOfCarrierTypes;j++)
		{
			sxx[i]+=electronCharge*carrierParams->getConcentration(j)*carrierParams->getMobility(j)/(1.0+
			carrierParams->getMobility(j)*carrierParams->getMobility(j)*B[i]*B[i]);

			sxy[i]+=electronCharge*carrierParams->getConcentration(j)*
			carrierParams->getMobility(j)*carrierParams->getMobility(j)*
			B[i]/(1.0+carrierParams->getMobility(j)*
			carrierParams->getMobility(j)*B[i]*B[i]);
		}
	}

}

void clMagneticFieldDependences::calculateEffectiveParamsFromTenzor()
{
	// вызывается автоматически из calculateMagneticFieldDependences()
	for(int i=0;i<NumberOfPoints;i++)
	{
		s_eff[i]=(sxx[i]*sxx[i]+sxy[i]*sxy[i])/sxx[i];
		Rh_eff[i]=sxy[i]/(sxx[i]*sxx[i]+sxy[i]*sxy[i]);
	}

}

void clMagneticFieldDependences::calculateSignalsFromEffectiveParams()
{
	for(int i=0;i<NumberOfPoints;i++)
	{
		Us[i]=carrierParams->getCBRatio()/carrierParams->getThickness()*carrierParams->getCurrentIntensity()/s_eff[i];
		Uy[i]=Rh_eff[i]*carrierParams->getCurrentIntensity()/carrierParams->getThickness();
	}

}
void clMagneticFieldDependences::calculateEffectiveParamsFromSignals()
{
	for (int i = 0; i < NumberOfPoints ; i++)
	{
		if(fabs(Us[i])<THEALMOSTZERO)
			s_eff[i]=0.0;
		else
		{
			s_eff[i]=carrierParams->getCBRatio()/carrierParams->getThickness()*carrierParams->getCurrentIntensity()/Us[i];
		}
		Rh_eff[i]=carrierParams->getThickness()*Uy[i]/carrierParams->getCurrentIntensity();
	}

}

void clMagneticFieldDependences::calculateTenzorFromEffectiveParams()
{
	// вызывается автоматически из calculateMagneticFieldDependences()
	for (int i = 0; i < NumberOfPoints ; i++)
	{
		sxx[i]=s_eff[i]/
			(Rh_eff[i]*Rh_eff[i]*s_eff[i]*s_eff[i]+1.0);
		sxy[i]=s_eff[i]*s_eff[i]*Rh_eff[i]/
		   (Rh_eff[i]*Rh_eff[i]*s_eff[i]*s_eff[i]+1.0);
	}
}

void clMagneticFieldDependences::calculateMagneticFieldDependences()
{
	calculateTenzorFromCarrierParams();
	calculateEffectiveParamsFromTenzor();
	calculateSignalsFromEffectiveParams();
}

TSignal const & clMagneticFieldDependences::getSignalUs()
{
    return Us;
}

TSignal const & clMagneticFieldDependences::getSignalB()
{
	return B;
}

TSignal const & clMagneticFieldDependences::getSignalUy()
{
    return Uy;
}

TSignal const & clMagneticFieldDependences::getSxx()
{
	return sxx;
}

TSignal const & clMagneticFieldDependences::getSxy()
{
    return sxy;
}

int clMagneticFieldDependences::modifySignals
(void (*ShumAdding)(TSignal const &x,TSignal &out,
	TSignal &ret, long double koeff,const int l),
	const TSignal & idealUs,const TSignal & idealUy,
	TSignal &returnData,long double koeff)
{

	ShumAdding(idealUs,Us,returnData,koeff,NumberOfPoints);
	ShumAdding(idealUy,Uy,returnData,koeff,NumberOfPoints);

	calculateEffectiveParamsFromSignals();
	calculateTenzorFromEffectiveParams();

	return returnData.size();
}

int clMagneticFieldDependences::modifySignals(void (*ShumAdding)(TSignal const &x,TSignal &out,
	TSignal &ret, long double koeff,const int l,long double mzr, int amplitude),
	const TSignal & idealB,
	const TSignal & idealUs,
	const TSignal & idealUy,TSignal &returnData,long double koeff,
	long double mzr, int amplitude)

{
	ShumAdding(idealB,B,returnData,koeff,NumberOfPoints,mzr*10,amplitude);
	ShumAdding(idealUs,Us,returnData,koeff,NumberOfPoints,mzr,amplitude);
	ShumAdding(idealUy,Uy,returnData,koeff,NumberOfPoints,mzr,amplitude);


	calculateEffectiveParamsFromSignals();
	calculateTenzorFromEffectiveParams();

	return returnData.size();
}

void clMagneticFieldDependences::setB_Us_Uy(long double *newB, long double * newUs,long double *newUy)
{
	for (int i = 0; i < NumberOfPoints; i++) {
		B[i]=newB[i];
		Us[i]=newUs[i];
		Uy[i]=newUy[i];
	}
	calculateEffectiveParamsFromSignals();
	calculateTenzorFromEffectiveParams();
}

void clMagneticFieldDependences::setB_Us_Uy(TSignal &newB,
	TSignal & newUs,TSignal &newUy)
{
	for (int i = 0; i < NumberOfPoints; i++) {
		B[i]=newB[i];
		Us[i]=newUs[i];
		Uy[i]=newUy[i];
	}
	calculateEffectiveParamsFromSignals();
	calculateTenzorFromEffectiveParams();
}



//---------------------------------------------------------------------------
