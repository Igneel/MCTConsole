#ifndef MagneticFieldDependenceH
#define MagneticFieldDependenceH

#include "FilteringUnit.h"

#include "DataSaver.h"
#include "FilterParams.h"
#include "ExtrapolateUnit.h"
#include "commonFunctions.h"
#include "mobilityspectrum.h"
#include "ErrorEst/clMagneticFieldDependences.h"

class MagneticFieldDependence
{
public:

TSignal * getPointerSxx();
//-------------------------------------------------------------------------
TSignal * getPointerSxy(/*DataKind dataKind*/);

	enum FeatType {ODD_FEAT, EVEN_FEAT};

	MagneticFieldDependence(std::string current, std::string temperature, std::string SampleInventoryNumber,
		std::string length, std::string width, std::string Thickness);

    MagneticFieldDependence(MyDataType current, MyDataType temperature, std::string SampleInventoryNumber,
        MyDataType length, MyDataType width, MyDataType Thickness);

	~MagneticFieldDependence();

    TSignal const * getB();

	TSignal const * getBHall();
	TSignal const * getBMagnetoResistance();
    TSignal const * getHallEffect();
    TSignal const * getMagnetoResistance();
  
    TSignal const * getFilteredB();

	TSignal const * getFilteredBHall();
	TSignal const * getFilteredBMagnetoResistance();
    TSignal const * getFilteredHallEffect();
    TSignal const * getFilteredMagnetoResistance();
    TSignal const * getExtrapolatedB();

	TSignal const * getExtrapolatedBHall();
	TSignal const * getExtrapolatedBMagnetoResistance();
    TSignal const * getExtrapolatedHallEffect();
    TSignal const * getExtrapolatedMagnetoResistance();   

    TSignal const * getSxx();
    TSignal const * getSxy();
    TSignal const * getAveragedB();
    TSignal const * getRh_eff();
    TSignal const * getS_eff();

    TSignal const * getMobility();
	TSignal const * getHoleConductivity();
	TSignal const * getElectronConductivity();
	
    TSignal const * getHoleConcentration();
    TSignal const * getHoleMobility();
    TSignal const * getElectronConcentration();
    TSignal const * getElectronMobility();

    void CopyTheorToSignals();
    void CopyTheorToNoisySignals();
	

	//(получение, фильтрация, экстраполяция, увеличение/уменьшение, вырезка и т.п.,
	//    построение графиков, сохранение результатов)

	void loadData(TStringList * tsl);
	void loadSampleDescription(TStringList *Names,TStringList *Values,std::string FileName);

	void setDependence(TSignal::iterator beginB, 
		TSignal::iterator endB, TSignal::iterator beginHall, 
		TSignal::iterator beginResistance);
	// Сохранение данных.
	void setSampleDescription(std::string Temperature, std::string Current, std::string SampleInventoryNumber, std::string length, std::string width, std::string Thickness);
    void setSampleDescription(MyDataType Temperature, MyDataType Current, std::string SampleInventoryNumber, MyDataType length, MyDataType width, MyDataType Thickness);


	void SaveAllData(std::string FileName);

	void shiftCurve(DataKind dataKind,SignalType dependenceType,MyDataType shiftValue,MyDataType leftBound, MyDataType rightBound);

	//-----Расчет тензора проводимости---------------------------------
	
	void calculateTenzor(DataKind dataKind);

	//------Фильтрация результатов-------------------------------------
    bool setFilterParamsHall(std::string samplingFrequecy,std::string bandwidthFrequency,std::string attenuationFrequency, std::string length);
	bool setFilterParamsHall(MyDataType samplingFrequecy,MyDataType bandwidthFrequency,MyDataType attenuationFrequency, int length);
    bool setFilterParamsResistance(string samplingFrequecy, string bandwidthFrequency, string attenuationFrequency, string length);
	bool setFilterParamsResistance(MyDataType samplingFrequecy,MyDataType bandwidthFrequency,MyDataType attenuationFrequency, int length);

	FilterParams const * getFilterParamsHall();
	FilterParams const * getFilterParamsResistance();
	// дискретизации, пропускания, затухания
	void filterData(FilterParams &fPHall, FilterParams &fPResistance);
	void filterData();

	// Экстраполяция результатов

    void setExtrapolateParams(int powPolinowHall,int powPolinomRes);
    void getExtrapolateParams(int & powPolinomHall,int & powPolinomRes);

	//bool extrapolateData(const int polinomPowForMagnetoResistance, const int polinomPowForHallEffect);
	bool extrapolateData( DataKind dataKind, const int polinomPowForMagnetoResistance, const int polinomPowForHallEffect);
	
	void setRoundNeeded(bool needRound);
	void setParamsType(ParamsType pt);
    void multiplySignal(SignalType s, MyDataType x);
    
    // Спектр подвижности
    bool calculateMobilitySpectrum(unsigned int NumberOfDecPlaces);
    bool runSmartCalcutation();
    // Многозонная подгонка
    bool runMultiCarrierFit(long double VesGxx, long double VesGxy);
    bool runSmartMultiCarrierFit(long double VesGxx, long double VesGxy);
    void getMultiCarrierFitResults(InDataSpectr & nMagSpectr, InDataSpectr & nGxxIn, InDataSpectr & GxyIn,
    MyData_spektr & outGxx, MyData_spektr &  outGxy, TStatistic & outValues);

    // Сохранение результатов------------------------------------------
	DataSaver * saver;

	mobilitySpectrum * MobilitySpectrumObj;

    //----Что нужно для консоли?---------------------------------------

    void calculateDependencesFromFilm(size_t NumberOfPoints, MyDataType h, MyDataType eMolarCompositionCadmium, int eTemperature, MyDataType eHeavyHoleConcentration, MyDataType eAFactor, MyDataType eKFactor, MyDataType eSampleThickness, MyDataType eCBRatio, MyDataType eCurrentIntensity, size_t NumberOfCarrierTypes);
    void addNoiseToSignals(long double coef);
    
    long double getTheorMobility(int n);
    long double getTheorConcentration(int n);

    //-----Автоматический генератор данных-----------------------------

private:

	//----Генераторы шума----------------------------------------------

    void QuativativeNoiseGenerator();
    void GaussianNoiseGenerator(long double coef);

    long double noiseCoef;

    std::vector<long double> mzr;

    
	//-----Модель плёнки и исходные сигналы----------------------------

    clMagneticFieldDependences * theoreticalDependences;
    clMagneticFieldDependences * noisyDependences;

    TSignal theorB;
    TSignal theorUs;
    TSignal theorUy;
    TSignal theorSxx;
    TSignal theorSxy;

    //-----------------------------------------------------------------

	ParamsType paramsType;
	TSignal leftBound;
	TSignal rightBound;		

	//------Фильтрация результатов-------------------------------------
	FilterParams *filterParamsHall;
	FilterParams *filterParamsResistance;

    int PowPolinomRes;
    int PowPolinomHall;

    void filterDataHelper(FilterParams &fP,
    SignalType dependenceType,DataKind dataKind);

    //------Загрузка данных-------------------------------------------- 
    void loadDataHelper(TSignal &temp, std::string AnsiS,const std::string delimiter);
	inline void ReplaceDotsToComma(std::string &in, std::string & out);

	//---------Обработка данных----------------------------------------
	void multiplyB(DataKind dataKind);

    void GetEqualNumberOfPoints(TSignal & B,
	TSignal & BHall,TSignal & BRes, TSignal & Hall,
	TSignal & Res);

	void calculateEffectiveParamsFromSignals();
	void calculateTenzorFromEffectiveParams();
	void featData(DataKind dataKind); // усреднение зависимостей, вызывать эту.
    void averageData(TSignal & inY, TSignal &outY, FeatType featType, unsigned long centerOfSimmetry);
	
    void cutData(); // оставляет только положительные значения магнитного поля

	void clearCurrentParams();
	void clearFilteredParams();

	TSignal * getPointerB(DataKind dataKind);
	TSignal * getPointerBHall(DataKind dataKind);
	TSignal * getPointerBResistance(DataKind dataKind);
	TSignal * getPointerHall(DataKind dataKind);
	TSignal * getPointerMagnetoResistance(DataKind dataKind);

	//---------------Много переменных----------------------------------
	// Текущие магнитное поле и эффект Холла/магнитосопротивление,
	// после всяческих преобразований (вырезка, увеличение и т.п.).

    DataSet currentData;

	TSignal B;
	TSignal BHall;
	TSignal BMagnetoResistance;
	TSignal HallEffect;
    TSignal MagnetoResistance;
	// Первоначальные значения, полученные от АЦП сохраняются в текущие.
	// Фильтрованные значения.
	TSignal FilteredB;
	TSignal FilteredBHall;
	TSignal FilteredBMagnetoResistance;
    TSignal FilteredHallEffect;
    TSignal FilteredMagnetoResistance;

	// Экстраполированные значения.
	TSignal ExtrapolatedB;
	TSignal ExtrapolatedBHall;
	TSignal ExtrapolatedBMagnetoResistance;
    TSignal ExtrapolatedHallEffect;
    TSignal ExtrapolatedMagnetoResistance;

    // Усредненные значения.
    TSignal AveragedB;
	TSignal AveragedBHall;
	TSignal AveragedBMagnetoResistance;
	TSignal AveragedHallEffect;
	TSignal AveragedMagnetoResistance;

    // И ещё мы кушаем память для усредненных значений, которые будем
    // кормить спектру подвижности.
    /*
    TSignal AveragedB;
    TSignal sxx;
	TSignal sxy; */

	TSignal s_eff;
	TSignal Rh_eff;
	TSignal sxx;
	TSignal sxy;

	// Результаты работы спектра подвижности
	TSignal mobility;
	TSignal holeConductivity;
	TSignal electronConductivity;
	// В этих хранятся значения, соответствующие пикам для дырок и электронов
    TSignal holeConcentration;
    TSignal holeMobility;
    TSignal electronConcentration;
    TSignal electronMobility;


    // Статистика спектра подвижности, необходимая для расчета многозонной подгонки:
    std::vector<long double> upBound;
    std::vector<long double> lowBound;

    // Сюда сохраняются выходные значения подгонки, смотри функцию getMultiCarrierFitResults
    InDataSpectr nMagSpectr;
    InDataSpectr nGxxIn;
    InDataSpectr nGxyIn;
    MyData_spektr outGxx;
    MyData_spektr outGxy;
    TStatistic outValues;
};


#endif
