#include "smartCalculation.h"
#include <sstream>
#include <fstream>


smartCalculation::smartCalculation(MagneticFieldDependence * base)
{
	/*
	Всё делается на основе данных, которые хранятся и рассчитываются в классе MagneticFieldDependence.
	*/
	numberOfParameters=6;
	numberOfStoredBestResults=25;
	MFDData=base;
	results.resize(numberOfParameters);
}

smartCalculation::~smartCalculation()
{
	;
}

bool isBetterResult(Results & r1, Results & r2)
// возвращает true если r1 лучше чем r2
{
	return fabs(r1.targetFuncValue)<fabs(r2.targetFuncValue);
}

bool smartCalculation::isBetterResults(Results & r1, Results & r2)
// возвращает true если r1 лучше чем r2
{
	return isBetterResult(r1, r2);
}

//bool compareResult (int i,int j) { return (i<j); }

bool smartCalculation::processResults()
{
	// Провести обработку

	Results nr; // наши новые результаты.
	nr.fpHall=*(MFDData->getFilterParamsHall());
	nr.fpRes=*(MFDData->getFilterParamsResistance());
	MFDData->getExtrapolateParams(nr.polinomPowHall,nr.polinomPowMagnetoresistance);

	if(MFDData->getElectronConcentration()->size()>=1)
	{
		nr.electronConcentration=MFDData->getElectronConcentration()->operator[](0);
		nr.electronMobility=MFDData->getElectronMobility()->operator[](0);
	}
	else
	{
		nr.electronConcentration=0;
		nr.electronMobility=0;
	}

	if(MFDData->getHoleConcentration()->size()>=1)
	{
		nr.heavyHoleConcentration=MFDData->getHoleConcentration()->operator[](0);
		nr.heavyHoleMobility=MFDData->getHoleMobility()->operator[](0);
	}
	else
	{
		nr.heavyHoleConcentration=0;
		nr.heavyHoleMobility=0;
	}

	if(MFDData->getHoleConcentration()->size()>=2)
	{
		nr.lightHoleConcentration=MFDData->getHoleConcentration()->operator[](1);
		nr.lightHoleMobility=MFDData->getHoleMobility()->operator[](1);
	}
	else
	{
		nr.lightHoleConcentration=0;
		nr.lightHoleMobility=0;
	}

    nr.vPC=MFDData->MobilitySpectrumObj->getPeaksCriteria();
    nr.additionalData=MFDData->MobilitySpectrumObj->getAdditionalData();

	nr.targetFuncValue=10000;
    allSpectras.push_back(nr);

	// Рассчитать целевые функции для всех видов носителей заряда
    for (size_t i = 0; i < numberOfParameters; ++i)
	{
        parameterType parType;
        switch (i) {
        case 0:
            parType = ELECTRON_MOBILITY;
            break;
        case 1:
            parType = ELECTRON_DENSITY;
            break;
        case 2:
            parType = LIGHT_HOLE_MOBILITY;
            break;
        case 3:
            parType = LIGHT_HOLE_DENSITY;
            break;
        case 4:
            parType = HEAVY_HOLE_MOBILITY;
            break;
        case 5:
            parType = HEAVY_HOLE_DENSITY;
            break;
        }
        nr.targetFuncValue= targetFunction(nr,parType);
        nr.pt=parType;

		if (results[i].size()==0)
		{
			/* (если пусто - добавить так) */
			results[i].push_back(nr);
			sort(results[i].begin(), results[i].end(), isBetterResult);
		}
		else
		{
			// Сравнить с имеющимися данными 
			if (isBetterResults(nr,results[i].back()))
			{
				// Если целевая функция лучше - сохранить
				results[i].push_back(nr);
				sort(results[i].begin(), results[i].end(), isBetterResult);
			}

			// Хранить будем numberOfStoredBestResults лучших значений целевой функции для каждого параметра

			if (results[i].size()>numberOfStoredBestResults)
			{
				results[i].pop_back();
			}
		}
	}
	return true;
}

bool smartCalculation::processData()
{
	
	/*
	Итак, тут обработка данных ведется.
	*/
    unsigned int NumberOfDecPlaces=5;

	// Для каждого значения частоты среза фильтра
	FilterParams fpHall=*(MFDData->getFilterParamsHall());
	FilterParams fpRes=*(MFDData->getFilterParamsResistance());
	int polinomPowHall;
	int polinomPowMagnetoresistance;
	MFDData->getExtrapolateParams(polinomPowHall,polinomPowMagnetoresistance);

	long double step=0.1;
    //Form1->ProgressBar1->Max=20*20;
	size_t iteration=0;
	for (long double BandWidthFHall=0.01, AttenuationFHall=0.1; BandWidthFHall < 10; BandWidthFHall+=step, AttenuationFHall+=step)
	{
		

		if(step==0.1 && BandWidthFHall>1.0)
		{
			step=1.0;
		}
		long double step=0.1;
		for (long double BandWidthFRes=0.01, AttenuationFRes=0.1; BandWidthFRes < 10; BandWidthFRes+=step, AttenuationFRes+=step)
		{
            //Form1->ProgressBar1->Position=iteration++;
            //Application->ProcessMessages();

			if(step==0.1 && BandWidthFRes>1.0)
			{
				step=1.0;
			}

        	fpHall.BandwidthFrequency=BandWidthFHall;
        	fpHall.AttenuationFrequency=AttenuationFHall;
        	fpRes.BandwidthFrequency=BandWidthFRes;
        	fpRes.AttenuationFrequency=AttenuationFRes;

        	MFDData->setFilterParamsResistance(fpRes.SamplingFrequecy, fpRes.BandwidthFrequency, 
        		fpRes.AttenuationFrequency, fpRes.filterLength);
			MFDData->setFilterParamsHall(fpHall.SamplingFrequecy, fpHall.BandwidthFrequency, 
        		fpHall.AttenuationFrequency, fpHall.filterLength);
			MFDData->filterData();
        	MFDData->extrapolateData(FILTERED_DATA,polinomPowMagnetoresistance,polinomPowHall);
    		MFDData->calculateTenzor(FILTERED_DATA);
            MFDData->calculateMobilitySpectrum(NumberOfDecPlaces);
			MFDData->MobilitySpectrumObj->calculatePeaksWeigth();
			if (0/* Нужна подгонка */)
			{
				/* сначала перепиши её в класс MagneticFieldDependence! */
			}
			processResults();
		}
	}
	return true;
}

std::string convertToString(AdditionalData & ad)
{
	std::ostringstream strs;
	strs << "leftPointElectronConductivity " << ad.leftPointElectronConductivity << "\n";
	strs << "leftPointHoleConductivity " << ad.leftPointHoleConductivity << "\n";
	strs << "rightPointElectronConductivity " << ad.rightPointElectronConductivity << "\n";
	strs << "rightPointHoleConductivity " << ad.rightPointHoleConductivity << "\n";
	strs << "leftPointElectronConductivityLog " << ad.leftPointElectronConductivityLog << "\n";
	strs << "leftPointHoleConductivityLog " << ad.leftPointHoleConductivityLog << "\n";
	strs << "rightPointElectronConductivityLog " << ad.rightPointElectronConductivityLog << "\n";
	strs << "rightPointHoleConductivityLog " << ad.rightPointHoleConductivityLog << "\n";
	strs << "leftPointMobility " << ad.leftPointMobility << "\n";
	strs << "rightPointMobility " << ad.rightPointMobility << "\n";

	return strs.str();
}

std::string convertToString(PeaksCriteria & pc)
{
	std::ostringstream strs;
	strs << "peakHeigh " << pc.peakHeigh << "\n";
	strs << "peakWidth " << pc.peakWidth << "\n";
	strs << "peakWidthOrd " << pc.peakWidthOrd << "\n";
	strs << "peakVelocityR " << pc.peakVelocityR << "\n";
	strs << "peakVelocity2R " << pc.peakVelocity2R << "\n";
	strs << "symmetri " << pc.symmetri << "\n";
	strs << "peakVelocityL " << pc.peakVelocityL << "\n";
	strs << "peakVelocity2L " << pc.peakVelocity2L << "\n";
	strs << "peakLeftBorderFirstDerivative " << pc.peakLeftBorderFirstDerivative << "\n";
	strs << "peakFirstDerivative " << pc.peakFirstDerivative << "\n";
	strs << "peakRightBorderFirstDerivative " << pc.peakRightBorderFirstDerivative << "\n";
	strs << "peakLeftMiddleFirstDerivative " << pc.peakLeftMiddleFirstDerivative << "\n";
	strs << "peakRightMiddleFirstDerivative " << pc.peakRightMiddleFirstDerivative << "\n";
	strs << "peakLeftBorderSecondDerivative " << pc.peakLeftBorderSecondDerivative << "\n";
	strs << "peakSecondDerivative " << pc.peakSecondDerivative << "\n";
	strs << "peakRightBorderSecondDerivative " << pc.peakRightBorderSecondDerivative << "\n";
	strs << "peakLeftMiddleSecondDerivative " << pc.peakLeftMiddleSecondDerivative << "\n";
	strs << "peakRightMiddleSecondDerivative " << pc.peakRightMiddleSecondDerivative << "\n";

	return strs.str();
}

bool smartCalculation::saveResults(std::string filename)
{
	std::ostringstream strs;

	for (std::deque < deque <Results> >::iterator res2 = results.begin(); res2 != results.end(); ++res2)
	{
		for (deque <Results> ::iterator res = res2->begin(); res != res2->end(); ++res)
		{
			strs << "fpHall " << res->fpHall.getFilterParams().c_str() << "\n";
			strs << "fpRes " << res->fpRes.getFilterParams().c_str() << "\n";
			strs << "polinomPowHall " << res->polinomPowHall << "\n";
			strs << "polinomPowMagnetoresistance " << res->polinomPowMagnetoresistance << "\n";

			strs << "electronConcentration " << res->electronConcentration << "\n";
			strs << "electronMobility " << res->electronMobility << "\n";
			strs << "heavyHoleConcentration " << res->heavyHoleConcentration << "\n";
			strs << "heavyHoleMobility " << res->heavyHoleMobility << "\n";
			strs << "lightHoleConcentration " << res->lightHoleConcentration << "\n";
			strs << "lightHoleMobility " << res->lightHoleMobility << "\n";

			strs << "targetFuncValue " << res->targetFuncValue << "\n";
			strs << res->pt << "\n";
			strs << "vPC\n";
			for (std::vector<PeaksCriteria>::iterator i = res->vPC.begin(); i != res->vPC.end(); ++i)
			{
				strs << convertToString(*i) << "\n";
			}

			strs << "additionalData " << convertToString(res->additionalData)<<endl;
		}
		strs << "end of first parameterType\n";
	}

	//---------------------------------------------------------
	// Расчет статистики по параметрам из спектра. Из которых мы будем исключать промахи и определять интервалы для многозонной подгонки.
	TSignal x;

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		if (allSpectras[i].electronMobility<90 && allSpectras[i].electronMobility > 0.1)
		{
			x.push_back(allSpectras[i].electronMobility);
		}
	}
	strs << "standard Deviation for electronMobility is " << standardDeviation(x) << "\n";
	strs << "mean Value for electronMobility is " << calculateMeanValue(x) << "\n";

	long double koefForDeviation=1;
	upBound.push_back(calculateMeanValue(x)+koefForDeviation*standardDeviation(x));
	lowBound.push_back(calculateMeanValue(x)-koefForDeviation*standardDeviation(x));

	x.clear();

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		if (allSpectras[i].electronMobility<90 && allSpectras[i].electronMobility > 0.1)
		{
			x.push_back(allSpectras[i].electronConcentration);
		}
	}
	strs << "standard Deviation for electronConcentration is " << standardDeviation(x) << "\n";
	strs << "mean Value for electronConcentration is " << calculateMeanValue(x) << "\n";

	upBound.push_back(calculateMeanValue(x)+koefForDeviation*standardDeviation(x));
	lowBound.push_back(calculateMeanValue(x)-koefForDeviation*standardDeviation(x));
	x.clear();

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		x.push_back(allSpectras[i].lightHoleMobility);
	}
	strs << "standard Deviation for lightHoleMobility is " << standardDeviation(x) << "\n";
	strs << "mean Value for lightHoleMobility is " << calculateMeanValue(x) << "\n";

	koefForDeviation=3;
	upBound.push_back(calculateMeanValue(x)+koefForDeviation*standardDeviation(x));
	lowBound.push_back(calculateMeanValue(x)-koefForDeviation*standardDeviation(x));
	x.clear();

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		x.push_back(allSpectras[i].lightHoleConcentration);
	}
	strs << "standard Deviation for lightHoleConcentration is " << standardDeviation(x) << "\n";
	strs << "mean Value for lightHoleConcentration is " << calculateMeanValue(x) << "\n";

	upBound.push_back(calculateMeanValue(x)+koefForDeviation*standardDeviation(x));
	lowBound.push_back(calculateMeanValue(x)-koefForDeviation*standardDeviation(x));
	x.clear();

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		x.push_back(allSpectras[i].heavyHoleMobility);
	}
	strs << "standard Deviation for heavyHoleMobility is " << standardDeviation(x) << "\n";
	strs << "mean Value for heavyHoleMobility is " << calculateMeanValue(x) << "\n";

	upBound.push_back(calculateMeanValue(x)+koefForDeviation*standardDeviation(x));
	lowBound.push_back(calculateMeanValue(x)-koefForDeviation*standardDeviation(x));
	x.clear();

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		x.push_back(allSpectras[i].heavyHoleConcentration);
	}
	strs << "standard Deviation for heavyHoleConcentration is " << standardDeviation(x) << "\n";
	strs << "mean Value for heavyHoleConcentration is " << calculateMeanValue(x) << "\n";

	upBound.push_back(calculateMeanValue(x)+koefForDeviation*standardDeviation(x));
	lowBound.push_back(calculateMeanValue(x)-koefForDeviation*standardDeviation(x));
	x.clear();

    for (size_t i = 0; i < allSpectras.size(); ++i)
	{
		strs << allSpectras[i].electronMobility << '\t' << allSpectras[i].electronConcentration << '\t' 
		<< allSpectras[i].lightHoleMobility << '\t' << allSpectras[i].lightHoleConcentration << '\t'
		<< allSpectras[i].heavyHoleMobility << '\t' << allSpectras[i].heavyHoleConcentration << '\n';
	}

    std::ofstream fout(filename);


    fout << strs.str();


	// Нужно изменить порядок границ на те, которые используются при подгонке.
	long double upTemp, lowTemp;
	
	std::vector<long double> newUpBound;
	std::vector<long double> newLowBound;

	newUpBound.resize(upBound.size());
	newLowBound.resize(upBound.size());

	 newUpBound[0]=upBound[0];
	newLowBound[0]=lowBound[0];
	 newUpBound[1]=upBound[2];
	newLowBound[1]=lowBound[2];
	 newUpBound[2]=upBound[4];
	newLowBound[2]=lowBound[4];
	 newUpBound[3]=upBound[1];
	newLowBound[3]=lowBound[1];
	 newUpBound[4]=upBound[3];
	newLowBound[4]=lowBound[3];
	 newUpBound[5]=upBound[5];
	newLowBound[5]=lowBound[5];

	upBound=newUpBound;
	lowBound=newLowBound;

	return true;
}

bool smartCalculation::isGoodExtremum(Results & r, parameterType pt)
{
	return true;
}

long double smartCalculation::targetFunction (Results & r, parameterType pt)
{
	/*
	Надо определять тип носителя.
	*/

	

	/*
	hhole
	lhole
	electron
	loghhole
	loglhole
	logelectron

	r.vPC[i].peakHeigh;
	r.vPC[i].peakWidth;
	r.vPC[i].peakWidthOrd;
	r.vPC[i].peakVelocityR;
	r.vPC[i].peakVelocity2R;
	r.vPC[i].symmetri;
	r.vPC[i].peakVelocityL;
	r.vPC[i].peakVelocity2L;
	r.vPC[i].peakLeftBorderFirstDerivative;
	r.vPC[i].peakFirstDerivative;
	r.vPC[i].peakRightBorderFirstDerivative;
	r.vPC[i].peakLeftMiddleFirstDerivative;
	r.vPC[i].peakRightMiddleFirstDerivative;
	r.vPC[i].peakLeftBorderSecondDerivative;
	r.vPC[i].peakSecondDerivative;
	r.vPC[i].peakRightBorderSecondDerivative;
	r.vPC[i].peakLeftMiddleSecondDerivative;
	r.vPC[i].peakRightMiddleSecondDerivative;

	r.additionalData.leftPointElectronConductivity;
	r.additionalData.leftPointHoleConductivity;
	r.additionalData.rightPointElectronConductivity;
	r.additionalData.rightPointHoleConductivity;
	r.additionalData.leftPointElectronConductivityLog;
	r.additionalData.leftPointHoleConductivityLog;
	r.additionalData.rightPointElectronConductivityLog;
	r.additionalData.rightPointHoleConductivityLog;

	*/
	const int VERYBIGTARGETFUNC=10000;
	switch (pt)
	{
		case ELECTRON_MOBILITY:
			if (fabs(r.additionalData.rightPointMobility-r.electronMobility)<1.0 || fabs(r.additionalData.leftPointMobility-r.electronMobility)<1.0)
			{
				return VERYBIGTARGETFUNC;
			}
			if (isGoodExtremum(r,pt))
			{
				//return r.vPC[ELECTRON].peakHeigh*r.additionalData.rightPointElectronConductivity/r.additionalData.rightPointHoleConductivity;
				return r.vPC[HHOLE].peakHeigh*r.vPC[HHOLE].peakWidth+2.94999999812644*r.vPC[ELECTRON].peakWidth*r.vPC[ELECTRON].peakWidth;
			}
			else
			{
				return r.vPC[ELECTRON].peakWidth- r.vPC[ELECTRON].peakLeftMiddleFirstDerivative;
			}
		/*
		ярко выраженный экстремум:
		высота экстремума электронов * величины проводимости электронов на границе вычисляемого спектра подвижности 100 м^2/(В∙с) / величины проводимости дырок на границе вычисляемого спектра подвижности 100 м^2/(В∙с)

		слабо выраженный:
		ширина пика электронов - величина первой производной для середины боковой стороны экстремума, соответствующего электронам
		*/
		case ELECTRON_DENSITY:
			if (fabs(r.additionalData.rightPointMobility-r.electronMobility)<1.0 || fabs(r.additionalData.leftPointMobility-r.electronMobility)<1.0)
			{
				return VERYBIGTARGETFUNC;
			}
			if (isGoodExtremum(r,pt))
			{
				return 1;
			}
			else
			{
				return r.vPC[ELECTRON].peakHeigh*r.additionalData.rightPointElectronConductivity/r.additionalData.rightPointHoleConductivity/r.vPC[ELECTRON].peakWidth;
			}
		/*
		слабо выраженный:
		высота экстремума электронов * величины проводимости электронов на границе вычисляемого спектра подвижности 100 м^2/(В∙с) / величины проводимости дырок на границе вычисляемого спектра подвижности 100 м^2/(В∙с) /ширина экстремума электронов
		*/
		case LIGHT_HOLE_MOBILITY:
			if (isGoodExtremum(r,pt))
			{
				r.vPC[LHOLE].peakWidth*r.vPC[HHOLE].peakWidth/r.vPC[ELECTRON].peakWidth/r.vPC[ELECTRON].peakWidth;
			}
			else
			{
				return 1;
			}
		/*
		ярко выраженный экстремум:
		ширина пика легких дырок* ширина пика тяжелых дырок / квадрат ширины пика электронов
		*/
		case LIGHT_HOLE_DENSITY:
			if (isGoodExtremum(r,pt))
			{
				return r.additionalData.rightPointElectronConductivity/r.additionalData.rightPointHoleConductivity;
			}
			else
			{
				return 1;
			}
		/*
		ярко выраженный экстремум:
		величины проводимости электронов на границе вычисляемого спектра подвижности 100 м^2/(В∙с) / величины проводимости дырок на границе вычисляемого спектра подвижности 100 м^2/(В∙с)
		*/
		case HEAVY_HOLE_MOBILITY:
			if (isGoodExtremum(r,pt))
			{
				return r.vPC[LHOLE].peakWidth/r.vPC[ELECTRON].peakWidth;
			}
			else
			{
				return r.additionalData.rightPointElectronConductivity;
			}
		/*
		ярко выраженный экстремум:
		ширина пика легких дырок / ширина пика электронов

		слабо выраженный:
		величины проводимости электронов на границе вычисляемого спектра подвижности 100 м^2/(В∙с)
		*/
		case HEAVY_HOLE_DENSITY:
			if (isGoodExtremum(r,pt))
			{
				return r.vPC[LHOLE].peakWidth/r.vPC[ELECTRON].peakWidth;
			}
			else
			{
				return 1;
			}
		/*
		ярко выраженный экстремум:
		ширина пика легких дырок / ширина пика электронов
		*/
	}

	return 0;
}


std::vector<long double> smartCalculation::getUpBound()
{
	return upBound;
}

std::vector<long double> smartCalculation::getLowBound()
{
	return lowBound;
}
