#include <iostream>
#include "main.h"
using namespace std;
void getAndDisplayMultiCarrierFitResults(MagneticFieldDependence * p);

int main(int argc, char *argv[])
{

    /*
     * если kNoise=0 - отдавать чистые данные и шум даже не считать.
     * ./main T kNoise current sampleLength sampleWidth sampleThickness eHeavyHoleConcentration molarCadmiunValue
     * ./main 77 1 10e-3 30e-3 10e-3 1e-5 1e22 0.21
     *
     *



    А нужен ли нам автомейшн здесь? Мб просто выдавать значения для заданных параметров?
    А остальное из питона прикрутить....
    */

    int Temperature; // начальная температура
    int coef;// коэффициент шума

    const MyDataType shagB=0.001;

    MyDataType eMolarCompositionCadmium;
    size_t NumberOfCarrierTypes = 3;
    std::string inventoryNumber="007";

    MyDataType current;

    MyDataType sampleLength,
            sampleWidth,
            sampleThickness;

    MyDataType eHeavyHoleConcentration;
    MyDataType electronMobility;

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


        }
    else
    {
        Temperature = stoi(argv[1]); // Температура
        coef = stoi(argv[2]); // Коэффициент шума
        current = stold(argv[3]); // Величина силы тока
        sampleLength = stold(argv[4]); // Длина образца
        sampleWidth = stold(argv[5]); // Ширина образца
        sampleThickness = stold(argv[6]); // Толщина образца
        eHeavyHoleConcentration = stold(argv[7]); // Концентрация тяжелых дырок
        eMolarCompositionCadmium = stold(argv[8]); // Молярный состав кадмия
        electronMobility = stold(argv[9]); // Подвижность электронов
    }


    MagneticFieldDependence * p = new MagneticFieldDependence(current,Temperature,inventoryNumber,
                                                              sampleLength,sampleWidth,sampleThickness);

    size_t NumberOfPoints = 2 / shagB + 1;


    MyDataType eAFactor = 5; //(5-8)
    MyDataType eKFactor = 1.3; //(1.3-1.5)

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

    std::string eSamplingFRes,eBandWidthFRes,eAttenuationFRes,eLengthFilterRes;
    std::string eSamplingFHall,eBandWidthFHall,eAttenuationFHall,eLengthFilterHall;

    eSamplingFRes = "100";
    eBandWidthFRes= "10";
    eAttenuationFRes= "20";
    eLengthFilterRes= "50";

    eSamplingFHall= "100";
    eBandWidthFHall= "10";
    eAttenuationFHall= "20";
    eLengthFilterHall= "50";



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

    DataKind dataflag=(coef==0?CURRENT_DATA:FILTERED_DATA);
    //cout << "CalculateTenzor\n";
    p->calculateTenzor(dataflag);

    // save Tenzor
    DataSaver * tenzorSaver=new DataSaver(to_string(Temperature),
    to_string(current), inventoryNumber,to_string(sampleLength),to_string(sampleWidth),to_string(sampleThickness));

    std::string fileName="tempSaveFileName";
    tenzorSaver->SaveData(CURRENT_DATA,p->getAveragedB(),
    p->getSxy(), p->getSxx(), ALL_POINTS,fileName);

    // 11 точек пока не нужны.
    //tenzorSaver->SaveData(CURRENT_DATA,p->getAveragedB(),
    //p->getSxy(), p->getSxx(), POINTS_11,fileName);

    delete tenzorSaver;



    // get Mobility Spectrum

    //cout << "calculateMobilitySpectrum\n";
    p->calculateMobilitySpectrum();

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


    fout << current << "\t" << CBRatio << "\t" << sampleThickness << "\n";

    


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
# v1
#   _______________________________________________________       ______________________
#   |T I CBRatio Thickness B1 .. Bn Us1 .. Usn Uy1 .. Uyn|  ---> |n1 mu1 n2 mu2 n3 mu3|
#   _______________________________________________________       ______________________ 
*/
    cout <<"\t" << Temperature <<"\t"  << current <<"\t"  << CBRatio <<"\t"  << sampleThickness <<"\t";
    TSignal const * B = p->getB();    
    TSignal const * Hall = p->getHallEffect();
    TSignal const * MagnetoRes = p->getMagnetoResistance();
    for(auto i = 0;i< B.size();i++)
    {
        cout << B[i] << "\t"
    }
    for(auto i = 0;i< Hall.size();i++)
    {
        cout << Hall[i] << "\t"
    }
    for(auto i = 0;i< MagnetoRes.size();i++)
    {
        cout << MagnetoRes[i] << "\t"
    }
    cout << "\n";
    for (auto i = 0; i < 3; ++i)
    {
        cout << p->getTheorMobility(0) << "\t" << p->getTheorConcentration(0) << "\t";
    }

/*
# v2
#  ___________________________________         ___________________________
#  |Критерий1 Критерий2 ... КритерийN|  --->   |Относительная погрешность|
#
*/

/*
# v3
#                               ______________________
#  |Спектр подвижности|   --->  |n1 mu1 n2 mu2 n3 mu3|
#                               ______________________
*/

    for (auto i = 0; i < ex.size(); ++i)
    {
        cout << ex[i] << "\t";
    }
    for (auto i = 0; i < eY.size(); ++i)
    {
        cout << eY[i] << "\t";
    }
    for (auto i = 0; i < hx.size(); ++i)
    {
        cout << hx[i] << "\t";
    }
    for (auto i = 0; i < hY.size(); ++i)
    {
        cout << hY[i] << "\t";
    }

    cout << "\n";
    for (auto i = 0; i < 3; ++i)
    {
        cout << p->getTheorMobility(0) << "\t" << p->getTheorConcentration(0) << "\t";
    }

    return 0;
}

void getAndDisplayMultiCarrierFitResults(MagneticFieldDependence * p)
{
    int numberOfCarrierTypes=3;
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

    for (int j = 0; j < 4; ++j)
    {

        for(int i=0;i<2*numberOfCarrierTypes;++i)
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
