#include <iostream>
#include "main.h"
using namespace std;


void getAndDisplayMultiCarrierFitResults(MagneticFieldDependence * p);

int main(int argc, char *argv[])
{
    std::ios::sync_with_stdio(false);

    //cout << __LDBL_EPSILON__ << "\n";

    /*
     * если kNoise=0 - отдавать чистые данные и шум даже не считать.
     * ./main T kNoise current sampleLength sampleWidth sampleThickness eHeavyHoleConcentration molarCadmiunValue AFactor KFactor
     * ./main 77 1 10e-3 30e-3 10e-3 1e-5 1e22 0.21 5 1.3
     *
     * eSamplingFRes,eBandWidthFRes,eAttenuationFRes,eLengthFilterRes eSamplingFHall,eBandWidthFHall,eAttenuationFHall,eLengthFilterHall;
    */

    int Temperature; // начальная температура
    long double coef;// коэффициент шума

    const MyDataType shagB=0.001L;

    MyDataType eMolarCompositionCadmium;
    size_t NumberOfCarrierTypes = 3;
    std::string inventoryNumber="007";

    MyDataType current;

    MyDataType sampleLength,
            sampleWidth,
            sampleThickness;

    MyDataType eHeavyHoleConcentration;

    MyDataType eAFactor = 5; //(5-8)
    MyDataType eKFactor = 1.3L; //(1.3-1.5)


    std::string eSamplingFRes,eBandWidthFRes,eAttenuationFRes,eLengthFilterRes;
    std::string eSamplingFHall,eBandWidthFHall,eAttenuationFHall,eLengthFilterHall;


    unsigned int NumberOfDecPlaces;

    std::ofstream fout("out.txt");

    if (argc<3)
    {
        cout << "Ops. Please specify agruments\n";
        return -1;
    }
    else if (argc == 3)
    {

        std::ifstream fin("data.txt");


        fin >> Temperature;
        fin >> coef;
        fin >> current;
        fin >> sampleLength >> sampleWidth >> sampleThickness;
        fin >> eHeavyHoleConcentration;
        fin >> eMolarCompositionCadmium;
        fin >> eAFactor;
        fin >> eKFactor;

        fin >>eSamplingFRes;
        fin >>eBandWidthFRes;
        fin >>eAttenuationFRes;
        fin >>eLengthFilterRes;
        fin >>eSamplingFHall;
        fin >>eBandWidthFHall;
        fin >>eAttenuationFHall;
        fin >>eLengthFilterHall;

        fin >>NumberOfDecPlaces;

        }
    else
    {
        Temperature = stoi(argv[1]); // Температура
        coef = stold(argv[2]); // Коэффициент шума
        current = stold(argv[3]); // Величина силы тока
        sampleLength = stold(argv[4]); // Длина образца
        sampleWidth = stold(argv[5]); // Ширина образца
        sampleThickness = stold(argv[6]); // Толщина образца
        eHeavyHoleConcentration = stold(argv[7]); // Концентрация тяжелых дырок
        eMolarCompositionCadmium = stold(argv[8]); // Молярный состав кадмия
        eAFactor = stold(argv[9]); // Влияет на подвижность электронов
        eKFactor = stold(argv[10]); // Влияет на подвижность электронов

        // Параметры фильтрации
        eSamplingFRes = argv[11];
        eBandWidthFRes = argv[12];
        eAttenuationFRes = argv[13];
        eLengthFilterRes = argv[14];
        eSamplingFHall = argv[15];
        eBandWidthFHall = argv[16];
        eAttenuationFHall = argv[17];
        eLengthFilterHall = argv[18];

        NumberOfDecPlaces = stoi(argv[19]);
    }


    MagneticFieldDependence * p = new MagneticFieldDependence(current,Temperature,inventoryNumber,
                                                              sampleLength,sampleWidth,sampleThickness);

    size_t NumberOfPoints =static_cast<size_t>(2 / shagB + 1);



    MyDataType CBRatio = sampleLength / sampleWidth;




    p->calculateDependencesFromFilm(NumberOfPoints,shagB,
                                    eMolarCompositionCadmium,
                                    Temperature,eHeavyHoleConcentration,
                                    eAFactor,eKFactor,
                                    sampleThickness,CBRatio,
                                    current,NumberOfCarrierTypes);
    p->CopyTheorToSignals();

    p->addNoiseToSignals(coef);
    /*
    Предупреждение Initializing Params Type with const int - не актуально,
    т.к. значения как раз совпадают с перечислением, но аккуратней в этом месте.
    */
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //p->saver->setParamsType(ResCurveIndex->ItemIndex); // значения их совпадают.
    //p->setParamsType(ResCurveIndex->ItemIndex);


    int PowPolinomHall=1, PowPolinomRes=2;

    p->setFilterParamsResistance(eSamplingFRes, eBandWidthFRes,
    eAttenuationFRes, eLengthFilterRes);
    p->setFilterParamsHall(eSamplingFHall, eBandWidthFHall,
    eAttenuationFHall, eLengthFilterHall);
    p->setExtrapolateParams(PowPolinomHall,PowPolinomRes);



    // filter & extrapolation
    p->filterData();
    //cout << "after filter\n";
    p->extrapolateData(FILTERED_DATA,PowPolinomRes,PowPolinomHall);


    // calculateTenzor
    p->setSampleDescription(to_string(Temperature),to_string(current),
    inventoryNumber,to_string(sampleLength),to_string(sampleWidth),to_string(sampleThickness));
    //p->setParamsType(ResCurveIndex->ItemIndex);

    DataKind dataflag=FILTERED_DATA; //(coef==0?CURRENT_DATA:FILTERED_DATA);
    //cout << "CalculateTenzor\n";
    p->calculateTenzor(dataflag);

    // save Tenzor
    /*
    DataSaver * tenzorSaver=new DataSaver(to_string(Temperature),
    to_string(current), inventoryNumber,to_string(sampleLength),to_string(sampleWidth),to_string(sampleThickness));

    std::string fileName="tempSaveFileName";
    tenzorSaver->SaveData(CURRENT_DATA,p->getAveragedB(),
    p->getSxy(), p->getSxx(), ALL_POINTS,fileName);

    // 11 точек пока не нужны.
    //tenzorSaver->SaveData(CURRENT_DATA,p->getAveragedB(),
    //p->getSxy(), p->getSxx(), POINTS_11,fileName);

    delete tenzorSaver;
    */


    // get Mobility Spectrum

    //cout << "calculateMobilitySpectrum\n";

    /*
    # v1
    #   _______________________________________________________       ______________________
    #   |T I CBRatio Thickness B1 .. Bn Us1 .. Usn Uy1 .. Uyn|  ---> |n1 mu1 n2 mu2 n3 mu3|
    #   _______________________________________________________       ______________________
    */
        cout << "NewSectionBeginHere\n";
        cout << Temperature <<"\t"  << current <<"\t"  << CBRatio <<"\t"  << sampleThickness <<"\t";
        TSignal const * B = p->getB();
        TSignal const * Hall = p->getHallEffect();
        TSignal const * MagnetoRes = p->getMagnetoResistance();
        for(auto i = 0u;i< B->size();i++)
        {
            cout << (*B)[i] << "\t";
        }
        for(auto i = 0u;i< Hall->size();i++)
        {
            cout << (*Hall)[i] << "\t";
        }
        for(auto i = 0u;i< MagnetoRes->size();i++)
        {
            cout << (*MagnetoRes)[i] << "\t";
        }
        cout << "\n";
        for (auto i = 0; i < 3; ++i)
        {
            cout << p->getTheorMobility(i) << "\t" << p->getTheorConcentration(i) << "\t";
        }




    /*
    # v2
    #  ___________________________________         ___________________________
    #  |Критерий1 Критерий2 ... КритерийN|  --->   |Относительная погрешность|
    #
    */
        cout << "\nNewSectionBeginHere\n";
        cout << (*Hall)[Hall->size()-1] << "\t";

        // Расчет спектр подвижности здесь:)
        // Сюда же хочу передавать параметры округления, т.к. оно должно происходить на последнем этапе
        p->calculateMobilitySpectrum(NumberOfDecPlaces);


        cout << "\n";



        for (auto i=0u;i<p->getHoleConcentration()->size();i++)
        {
            cout << p->getHoleMobility()->operator[](i) << "\t";
            cout << p->getHoleConcentration()->operator[](i) << "\t";

        }

        if(p->getElectronConcentration()->size()>=1)
        {
            cout << p->getElectronMobility()->operator[](0) << "\t";
            cout << p->getElectronConcentration()->operator[](0) << "\t";
        }

        cout << "\n";
        for (auto i = 0; i < 3; ++i)
        {
            cout << p->getTheorMobility(i) << "\t" << p->getTheorConcentration(i) << "\t";
        }

    TSignal ex(p->getMobility()->begin(),p->getMobility()->end());
    TSignal eY(p->getElectronConductivity()->begin(),p->getElectronConductivity()->end());
    TSignal hx(p->getMobility()->begin(),p->getMobility()->end());
    TSignal hY(p->getHoleConductivity()->begin(),p->getHoleConductivity()->end());

    // run MultiCarrierFit

    //long double VesGxx=1;
    //long double VesGxy=1;

    //cout << "runSmartMultiCarrierFit\n";

    //p->runMultiCarrierFit(VesGxx,VesGxy);
    //getAndDisplayMultiCarrierFitResults(p);

    //cout << "runMultiCarrierFit is ended\n";


    //fout << current << "\t" << CBRatio << "\t" << sampleThickness << "\n";

    


    //for (auto i =0; i< )

    /*if(p->getElectronConcentration()->size()>=1)
    {
    Form1->MobSpecResults->Cells[1][3]=FloatToStrF(p->getElectronConcentration()->operator[](0),ffExponent,5,2);
    Form1->MobSpecResults->Cells[2][3]=FloatToStrF(p->getElectronMobility()->operator[](0),ffFixed,5,5);
    }
    if(p->getHoleConcentration()->size()>=1)
    {
    Form1->MobSpecResults->Cells[1][1]=FloatToStrF(p->getHoleConcentration()->operator[](0),ffExponent,5,2);
    Form1->MobSpecResults->Cells[2][1]=FloatToStrF(p->getHoleMobility()->operator[](0),ffFixed,5,5);
    }
    if(p->getHoleConcentration()->size()>=2)
    {
    Form1->MobSpecResults->Cells[1][2]=FloatToStrF(p->getHoleConcentration()->operator[](1),ffExponent,5,2);
    Form1->MobSpecResults->Cells[2][2]=FloatToStrF(p->getHoleMobility()->operator[](1),ffFixed,5,5);
    }

    for(int i =0;i<ex.size();++i)
    {
    Form1->ChspElecComponent->AddXY(ex[i],eY[i],"",clBlue);
    Form1->ChSpHoleComponent->AddXY(hx[i],hY[i],"",clRed);
    }*/

        /*

        Допустим первый аргумент - имя файла с параметрами от которых мы начинаем работу (это чтобы можно было продолжить:))


        Второй файл - имя файла для результатов (или имя-шаблон такого файла)

        */
    //}


/*
# v3
#                               ______________________
#  |Спектр подвижности|   --->  |n1 mu1 n2 mu2 n3 mu3|
#                               ______________________
*/
    cout << "\nNewSectionBeginHere\n";
    for (auto& i: ex) {
        cout << i << "\t";
    }
    for (auto& i: eY) {
        cout << i << "\t";
    }
    for (auto& i: hx) {
        cout << i << "\t";
    }
    for (auto& i: hY) {
        cout << i << "\t";
    }
/*
    for (auto i = 0; i < ex.size(); ++i)
    {
        fout << ex[i] << "\t";
    }
    fout <<"\n";
    for (auto i = 0; i < eY.size(); ++i)
    {
        fout << eY[i] << "\t";
    }
    fout <<"\n";
    for (auto i = 0; i < hx.size(); ++i)
    {
        fout << hx[i] << "\t";
    }
    fout <<"\n";
    for (auto i = 0; i < hY.size(); ++i)
    {
        fout << hY[i] << "\t";
    }
    fout <<"\n";*/

    cout << "\n";

    for (auto i=0u;i<p->getHoleConcentration()->size();i++)
    {
        cout << p->getHoleMobility()->operator[](i) << "\t";
        cout << p->getHoleConcentration()->operator[](i) << "\t";
    }

    if(p->getElectronConcentration()->size()>=1)
    {
        cout << p->getElectronMobility()->operator[](0) << "\t";
        cout << p->getElectronConcentration()->operator[](0) << "\t";
    }

    cout << "\n";
    for (auto i = 0; i < 3; ++i)
    {
        cout << p->getTheorMobility(i) << "\t" << p->getTheorConcentration(i) << "\t";
    }

    fout.close();
    return 0;
}

void getAndDisplayMultiCarrierFitResults(MagneticFieldDependence * p)
{
    auto numberOfCarrierTypes=3u;
    InDataSpectr nMagSpectr;
    InDataSpectr nGxxIn;
    InDataSpectr nGxyIn;
    // Сюда сохраняются выходные значения.
    MyData_spektr outGxx;
    MyData_spektr outGxy;
    TStatistic outValues;

    p->getMultiCarrierFitResults(nMagSpectr, nGxxIn, nGxyIn,
    outGxx,  outGxy, outValues);

    // Вывод результатов на графики и форму.
    /*
    ChSxx_theor->Clear();
    ChSxy_theor->Clear();
    ChSxxExper->Clear();
    ChSxyExper->Clear();


    long double koefSxx=1, koefSxy=1;

    if (CheckBox3->Checked) // относительный график?
    {
        koefSxx=max_elem(outGxx);
        koefSxy=max_elem(outGxy);

        if (max_elem(nGxxIn)>koefSxx)
        {
            koefSxx=max_elem(nGxxIn);
        }
        if (max_elem(nGxyIn)>koefSxy)
        {
            koefSxy=max_elem(nGxyIn);
        }
    }

    for (int i = 0; i < outGxx.size(); ++i)
    {
        ChSxx_theor->AddXY(nMagSpectr[i],outGxx[i]/koefSxx,"",clRed);
        ChSxy_theor->AddXY(nMagSpectr[i],outGxy[i]/koefSxy,"",clRed);
    }

    for (int i = 0; i < nMagSpectr.size(); ++i)
    {
        ChSxxExper->AddXY(nMagSpectr[i],nGxxIn[i]/koefSxx,"",clBlue);
        ChSxyExper->AddXY(nMagSpectr[i],nGxyIn[i]/koefSxy,"",clBlue);
    }*/

    // Результаты идут в формате:
    // По первому индексу:
    // Подвижность электронов, подвижность легких дырок, подвижность тяжелых дырок
    // Концентрация электронов, концентрация легких дырок, концентрация тяжелых дырок
    // По второму индексу: минимальные, средние, СКО, СКО %

    TStringList * tsl = new TStringList();

    for (auto j = 0u; j < 4u; ++j)
    {

        for(auto i=0u;i<2*numberOfCarrierTypes;++i)
        {
            tsl->push_back(to_string(outValues[i][j]));
        }
    }
    SaveToFile(tsl,"MZFitRes.txt");
/*
    // min
    for(int i=0;i<2*numberOfCarrierTypes;++i)
    {

        if(i<numberOfCarrierTypes)
            FitResults->Cells[i+1][1]=FloatToStrF(outValues[i][0],ffFixed,6,6);
        else
      FitResults->Cells[i+1][1]=FloatToStrF(outValues[i][0],ffExponent,5,2);
    }
    //ErrorLog->Lines->Add("middle");
    for(int i=0;i<2*numberOfCarrierTypes+1;++i)
    {
        if(i<numberOfCarrierTypes)
            FitResults->Cells[i+1][2]=FloatToStrF(outValues[i][1],ffFixed,6,6);
        else
       FitResults->Cells[i+1][2]=FloatToStrF(outValues[i][1],ffExponent,5,2);
    }
    //ErrorLog->Lines->Add("SKO");
    for(int i=0;i<2*numberOfCarrierTypes+1;++i)
    {
        if(i<numberOfCarrierTypes)
            FitResults->Cells[i+1][3]=FloatToStrF(outValues[i][2],ffFixed,6,6);
        else
       FitResults->Cells[i+1][3]=FloatToStrF(outValues[i][2],ffExponent,5,2);
    }
    //ErrorLog->Lines->Add("SKOPers");
    for(int i=0;i<2*numberOfCarrierTypes+1;++i)
    {
        if(i<numberOfCarrierTypes)
            FitResults->Cells[i+1][4]=FloatToStrF(outValues[i][3],ffFixed,6,6);
        else
       FitResults->Cells[i+1][4]=FloatToStrF(outValues[i][3],ffExponent,5,2);
    }
    */
    delete tsl;
}
