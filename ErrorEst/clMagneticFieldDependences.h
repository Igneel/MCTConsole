//---------------------------------------------------------------------------

#ifndef clMagneticFieldDependencesH
#define clMagneticFieldDependencesH

#include "film.h"
#include "NoiseUnit.h"

#include "../commonFunctions.h"

class clMagneticFieldDependences {

	public:
    //static const enum DependenceType {SXX,SXY,US,UY,S_EFF,RH_EFF};
    //static const enum ModifyType {SHUM_ADDING, FILTERING, EXTRAPOLATE};
    
	clMagneticFieldDependences(int size,long double shag,long double molarCompositionCadmium,
    long double Temperature,long double heavyHoleConcentrerion, long double AFactor,long double KFactor,
	long double thickness,long double cbRatio,long double currentIntensity,long double numberOfCarrierTypes);

	clMagneticFieldDependences(int size,long double shag,film * cp);
	

	~clMagneticFieldDependences();

	TSignal const & getSignalB();
	TSignal const & getSignalUs();
	TSignal const & getSignalUy();
	TSignal const & getSxx();
	TSignal const & getSxy();

	long double getConcentration(int n);
	long double getMobility(int n);

	int modifySignals(void (*ShumAdding)(TSignal const &x,TSignal &out,
	TSignal &ret, long double koeff,const int l),const TSignal & idealUs,
	 const TSignal & idealUy,TSignal &returnData,long double koeff);

	int modifySignals(void (*ShumAdding)(TSignal const &x,TSignal &out,
	TSignal &ret, long double koeff,const int l,long double mzr, int amplitude),
	const TSignal & idealB,
	const TSignal & idealUs,
	const TSignal & idealUy,TSignal &returnData,long double koeff,
	long double mzr, int amplitude);

    void setB_Us_Uy(long double *newB, long double * newUs,long double *newUy);

	void setB_Us_Uy(TSignal &newB, TSignal & newUs,TSignal &newUy);

	void calculateEffectiveParamsFromSignals();
	void calculateTenzorFromEffectiveParams();
	void calculateMagneticFieldDependences();

	private:

	film *carrierParams;

	void calculateTenzorFromCarrierParams();
	void calculateEffectiveParamsFromTenzor();
	void calculateSignalsFromEffectiveParams();
	void calculateMagneticFieldPoints(); // расчет точек магнитного поля.
	void MemoryAllocation(); // выделение памяти.

	static const long double THEALMOSTZERO;
	static const long double electronCharge;

	long double h;
	const int NumberOfPoints;

	TSignal sxx;
	TSignal sxy;
	TSignal B;
	TSignal Us; // продольный сигнал
	TSignal Uy; // поперечный сигнал
	TSignal s_eff;
	TSignal Rh_eff;

};
//---------------------------------------------------------------------------
#endif
