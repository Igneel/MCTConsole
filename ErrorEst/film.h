//---------------------------------------------------------------------------

#ifndef filmH
#define filmH

class film {

	public:

	film(long double molarCompositionCadmium,
	long double Temperature,long double heavyHoleConcentration,
	long double AFactor,long double KFactor,
	long double thickness,long double cbRatio,
	long double currentIntensity,long double numberOfCarrierTypes);

	~film();

	long double getNumberOfCarrierTypes();
	long double getConcentration(int index);
	long double getMobility(int index);

	long double getThickness();
	long double getCBRatio();
    long double getCurrentIntensity();
	long double getCurrentTemperature();
	long double getAFactor();
	long double getKFactor();

	long double getMolarCompositionCadmium();

	private:

	void calcutatingCarrierParams();

	long double ownConcentrationOfElectrons();
	long double Thickness;
	long double CBRatio;
	long double MolarCompositionCadmium;
	long double CurrentTemperature;
	long double CurrentIntensity;
	long double AFactor;
	long double KFactor;

	const int NumberOfCarrierTypes;

	long double *Concentration;
	long double *Mobility;
};
//---------------------------------------------------------------------------
#endif
