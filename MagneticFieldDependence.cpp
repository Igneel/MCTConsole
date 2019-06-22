#include "MagneticFieldDependence.h"
#include "multizoneFit.h"
#include "mobilityspectrum.h"
#include "smartCalculation.h"
#include <algorithm>

void MagneticFieldDependence::addNoiseToSignals(long double coef)
{
    noiseCoef = coef;
    /*
     * Разумный порядок добавления такой:
     * Обычный шум (нормальный)
     * Квантовый шум
     * Квантование по уровню
     *
     */
    if (coef == 0.0L)
    {
        calculateTenzor(CURRENT_DATA);
    }
    else
    {
        coef = coef*HallEffect[HallEffect.size()-1]/100.0L;
        GaussianNoiseGenerator(coef);
        QuativativeNoiseGenerator();
        calculateTenzor(CURRENT_DATA); // После добавления шума надо пересчитать тензора
    }

}

void MagneticFieldDependence::QuativativeNoiseGenerator()
{
    //QuantumShumAdding(std::vector<long double> const &x,std::vector<long double> &out,
    //    std::vector<long double> &ret, long double koeff,const int l,long double mzr, int amplitude)


    TSignal NoisyB;
    TSignal NoisyHallEffect;
    TSignal NoisyMagnetoResistance;

    TSignal returnDataB;
    TSignal returnDataMR;
    TSignal returnDataHE;
    int amplitude = 2; // Максимальная величина шума в количестве МЗР.

    QuantumShumAdding(B,NoisyB,returnDataB,B.size(),mzr[2]*10,amplitude);
    QuantumShumAdding(MagnetoResistance,NoisyMagnetoResistance,returnDataMR,MagnetoResistance.size(),mzr[1],amplitude);
    QuantumShumAdding(HallEffect,NoisyHallEffect,returnDataHE,HallEffect.size(),mzr[1],amplitude);
}
void MagneticFieldDependence::GaussianNoiseGenerator(long double coef)
{
    //ShumAdding(std::vector<long double> const &x,std::vector<long double> &out,
    //    std::vector<long double> &ret, long double koeff,const int l)

    // Ахтунг, здесь не добавлется шум к магнитному полю.

    //TSignal NoisyB;
    TSignal NoisyHallEffect;
    TSignal NoisyMagnetoResistance;

    TSignal returnDataMR;
    TSignal returnDataHE;

    NoisyHallEffect.resize(HallEffect.size());
    NoisyMagnetoResistance.resize(MagnetoResistance.size());

    ShumAdding(MagnetoResistance,NoisyMagnetoResistance,returnDataMR,coef,MagnetoResistance.size());
    ShumAdding(HallEffect,NoisyHallEffect,returnDataHE,coef,HallEffect.size());

    copy(NoisyMagnetoResistance.begin(),NoisyMagnetoResistance.end(),MagnetoResistance.begin());
    copy(NoisyHallEffect.begin(),NoisyHallEffect.end(),HallEffect.begin());
}


long double MagneticFieldDependence::getTheorMobility(int n)
{
    return theoreticalDependences->getMobility(n);
}

long double MagneticFieldDependence::getTheorConcentration(int n)
{
    return theoreticalDependences->getConcentration(n);
}


void MagneticFieldDependence::calculateDependencesFromFilm(size_t NumberOfPoints, MyDataType h,
                                                           MyDataType eMolarCompositionCadmium,
                                                           int eTemperature,MyDataType eHeavyHoleConcentration,
                                                           MyDataType eAFactor,MyDataType eKFactor,
                                                           MyDataType eSampleThickness,MyDataType eCBRatio,
                                                           MyDataType eCurrentIntensity,size_t NumberOfCarrierTypes)
{
    /*.
    Заполняем магнитное поле.
    Рассчитываем компоненты тензора по значениям параметров плёнки.
    Рассчитываем эффективные значения.
    Получаем измеряемые сигналы.
    */

    //h=0.001;
    //NumberOfPoints=2/h+1;

    //--------------------------Классы------------------------------------------
    theoreticalDependences=new clMagneticFieldDependences(NumberOfPoints,h,
    eMolarCompositionCadmium,
    eTemperature,eHeavyHoleConcentration,
    eAFactor,eKFactor,
    eSampleThickness,eCBRatio,
    eCurrentIntensity,NumberOfCarrierTypes);


/*
    theoreticalDependences->getConcentration(2);// электроны
    theoreticalDependences->getMobility(2);
    theoreticalDependences->getConcentration(1);
    theoreticalDependences->getMobility(1);
    theoreticalDependences->getConcentration(0);
    theoreticalDependences->getMobility(0); // тяжелые дырки
*/
    //--------------------------Классы------------------------------------------
    theoreticalDependences->calculateMagneticFieldDependences();

    const TSignal & tB   = theoreticalDependences->getSignalB();
    const TSignal & tUs  = theoreticalDependences->getSignalUs();
    const TSignal & tUy  = theoreticalDependences->getSignalUy();
    const TSignal & tSxx = theoreticalDependences->getSxx();
    const TSignal & tSxy = theoreticalDependences->getSxy();


    theorB.resize(tB.size());
    theorUs.resize(tUs.size());
    theorUy.resize(tUy.size());
    theorSxx.resize(tSxx.size());
    theorSxy.resize(tSxy.size());
    std::copy(tB.begin(),tB.end(),theorB.begin());
    std::copy(tUs.begin(),tUs.end(),theorUs.begin());
    std::copy(tUy.begin(),tUy.end(),theorUy.begin());
    std::copy(tSxx.begin(),tSxx.end(),theorSxx.begin());
    std::copy(tSxy.begin(),tSxy.end(),theorSxy.begin());

}

void MagneticFieldDependence::CopyTheorToSignals()
{
    B.resize(theorB.size());
    BHall.resize(theorB.size());
    BMagnetoResistance.resize(theorB.size());
    MagnetoResistance.resize(theorUs.size());
    HallEffect.resize(theorUy.size());
    sxx.resize(theorSxx.size());
    sxy.resize(theorSxy.size());




    std::copy(theorB.begin(),theorB.end(),B.begin());
    std::copy(theorB.begin(),theorB.end(),BMagnetoResistance.begin());
    std::copy(theorB.begin(),theorB.end(),BHall.begin());
    std::copy(theorUy.begin(),theorUy.end(),HallEffect.begin());
    std::copy(theorUs.begin(),theorUs.end(),MagnetoResistance.begin());
    std::copy(theorSxx.begin(),theorSxx.end(),sxx.begin());
    std::copy(theorSxy.begin(),theorSxy.end(),sxy.begin());

}

void MagneticFieldDependence::CopyTheorToNoisySignals()
{

}

MagneticFieldDependence::MagneticFieldDependence(MyDataType current, MyDataType temperature,
                                                 std::string SampleInventoryNumber,
                                                 MyDataType length, MyDataType width, MyDataType Thickness)
{
    filterParamsHall=new FilterParams(); // по идее нужно бы и инциализировать их тут, дабы не было проблем в случае чего:).
    filterParamsResistance=new FilterParams();
    saver =new DataSaver(to_string(temperature),to_string(current),
                         SampleInventoryNumber, to_string(length), to_string(width), to_string(Thickness));
    MobilitySpectrumObj=nullptr;
    paramsType=DIRECT;
    leftBound.resize(3);
    leftBound[DIRECT]=0;
    leftBound[REVERSE]=-2;
    leftBound[COMBINE]=-2;
    rightBound.resize(3);
    rightBound[DIRECT]=2;
    rightBound[REVERSE]=0;
    rightBound[COMBINE]=2;

    PowPolinomRes=2;
    PowPolinomHall=1;

    mzr.push_back(1.472331E-03L);
    mzr.push_back(4.206659E-04L); // обычно использую для холла и магнитосопротивления, диапазон значений 2.5В
    mzr.push_back(1.942127E-04L); // обычно использую для магнитного поля. диапазон значений 0,625В
    mzr.push_back(1.601563E-04L);
}

MagneticFieldDependence::~MagneticFieldDependence()
{
    if (saver)
    {
        delete saver;
        saver=nullptr;
    }
    
    if (filterParamsHall)
    {
        delete filterParamsHall;
        filterParamsHall=nullptr;
    }
    
    if (filterParamsResistance)
    {
        delete filterParamsResistance;
        filterParamsResistance=nullptr;
    }    

    if(MobilitySpectrumObj)
    {
        delete MobilitySpectrumObj;
        MobilitySpectrumObj=nullptr;
    }

    if(theoreticalDependences)
    {
        delete theoreticalDependences;
        theoreticalDependences = nullptr;
    }

}

void MagneticFieldDependence::loadDataHelper(TSignal &temp, string AnsiS, const std::string delimiter)
{
    temp.clear();
    std::string s1 = AnsiS.c_str();
    std::string s; 
    std::string::size_type start = 0;
    std::string::size_type finish = 0;
    ReplaceDotsToComma(s1,s);
    do 
        {
            finish = s.find(delimiter, start);
            std::string word = s.substr(start, finish-start);
            temp.push_back(stold(word.c_str()));
            start=finish+1;
        } while (finish != std::string::npos);
}

void MagneticFieldDependence::loadSampleDescription(TStringList *Names,TStringList *Values,std::string FileName)
{
    // Откусывает от имени файла часть, чтобы получить имя файла с описанием образца
    std::string temp = FileName.c_str();

    unsigned long rBound = temp.find_last_of("/\\");
    unsigned long rBound2 = temp.find_first_of("DCRT",rBound);
    std::string NewFileName = FileName.substr(0,rBound2)+"Description"+".txt";

    TStringList * tts = new TStringList();
    if(LoadFromFile(tts,NewFileName)) // загружаем это содержимое сюда.
    {

    const std::string delimiter = "\t"; // и парсим

    std::string s;
    for (auto i=tts->begin();i!=tts->end();++i)
    {

        ReplaceDotsToComma(*i,s); // если что - приводим всё к единому формату
        *i=s;
    }

    Names->clear();
    Values->clear();

    for (auto i=0u;i<tts->size();++i)
    {
        s=(*tts)[i];
        std::string::size_type start = 0;
        std::string::size_type finish = 0;

        if(s.empty())
            continue;

        finish = s.find(delimiter, start);
        Names->push_back(s.substr(start, finish-start).c_str());
        Values->push_back(s.substr(finish-start+1, s.length()).c_str());
    }

    setSampleDescription((*Values)[1],(*Values)[2],(*Values)[0],
    (*Values)[3],(*Values)[4],(*Values)[5]);
    }
    delete tts;
}

void MagneticFieldDependence::loadData(TStringList * tts)
{
    /*
    Загрузка самих данных.
    */
    const std::string delimiterTab = "\t";
    const std::string delimiterSpace = " ";
    const int numberOfColls=3;

    TSignal B;
    TSignal Hall;
    TSignal Resistance;
    TSignal temp;

    for(auto i=0u;i<tts->size();++i) // по количеству строк
    {
        if((*tts)[i].empty()) // пустые строки пропускаем
            continue;
            
        loadDataHelper(temp,(*tts)[i],delimiterTab);

        if(temp.size()!=numberOfColls)
        {          
        loadDataHelper(temp,(*tts)[i],delimiterSpace);
        if(temp.size()!=numberOfColls)
            return;
        }

        B.push_back(temp[0]);
        Resistance.push_back(temp[1]);
        Hall.push_back(temp[2]);
    }

    setDependence(B.begin(),B.end(),Hall.begin(),Resistance.begin());
    //testB= new MagneticField(B.begin(),B.end());
    //testHall= new THallEffect(Hall.begin(),Hall.end());
}

void MagneticFieldDependence::SaveAllData(std::string FileName)
{
    // ахтунг, это страшный костыль -
    // там где эта функция вызывается задаются имена файлов
    // Com Dir и Rev
    // по ним я и определяю тип сигнала.!!!!
    ParamsType pt;
    bool isCombinedParams=false;

    if(FileName[FileName.length()]=='m')
        {
            pt=COMBINE;
            isCombinedParams=true;
        }
    if(FileName[FileName.length()]=='v')  pt=REVERSE;
    if(FileName[FileName.length()]=='r')  pt=DIRECT;

        saver->setParamsType(pt);
        saver->SaveData(CURRENT_DATA,&B,&HallEffect,&MagnetoResistance,(isCombinedParams?POINTS_21:POINTS_11),FileName);
        saver->SaveData(CURRENT_DATA,&B,&HallEffect,&MagnetoResistance,ALL_POINTS,FileName);

        if(FilteredB.size()!=0)
        {
            saver->SaveData(FILTERED_DATA,&FilteredB,&FilteredHallEffect,&FilteredMagnetoResistance,(isCombinedParams?POINTS_21:POINTS_11),FileName);
            saver->SaveData(FILTERED_DATA,&FilteredB,&FilteredHallEffect,&FilteredMagnetoResistance,ALL_POINTS,FileName);
        }
        //saver->SaveData(EXTRAPOLATED_DATA,&ExtrapolatedB,&ExtrapolatedHallEffect,&ExtrapolatedMagnetoResistance,(isCombinedParams?POINTS_21:POINTS_11),FileName);
        if(ExtrapolatedB.size()!=0)
        saver->SaveData(EXTRAPOLATED_DATA,&ExtrapolatedB,&ExtrapolatedHallEffect,&ExtrapolatedMagnetoResistance,ALL_POINTS,FileName);

        if(sxx.size()!=0)
        {
            saver->SaveData(TENZOR_DATA,&AveragedB,&sxx,&sxy,POINTS_11,
                FileName+"R"+ filterParamsResistance->getFilterParams()+ "H"+ filterParamsHall->getFilterParams());
            saver->SaveData(TENZOR_DATA,&AveragedB,&sxx,&sxy,ALL_POINTS,
                FileName+"R"+ filterParamsResistance->getFilterParams()+ "H"+ filterParamsHall->getFilterParams());

            saver->SaveData(MOBILITY_DATA,&mobility,&holeConductivity,&electronConductivity,ALL_POINTS,FileName);
        }

        if (pt==COMBINE)
        {
            saver->SaveData(AVERAGED_DATA,&AveragedB,&AveragedHallEffect,&AveragedMagnetoResistance,ALL_POINTS,FileName);
        }

        if(MobilitySpectrumObj)
        {
            std::string t1=FileName.c_str();
            MobilitySpectrumObj->saveEigenValues(t1+"Eigen.txt");
            MobilitySpectrumObj->saveResults(t1+"MobSpecRes.txt");
            MobilitySpectrumObj->savePeakWeigth(t1+"PeakWeigth.txt");
        }
        saver->SaveSampleDescription(FileName);
}

//------------Подгонка данных-------------------------------------------------

void MagneticFieldDependence::averageData(TSignal & inY, TSignal &outY, FeatType featType, unsigned long centerOfSimmetry)
{
    /*
    Усредненный сигнал на данный момент существует лишь для положительного магнитного поля.
   
    */
    auto size=inY.size();

    // Размер рассчитывается в зависимости от точки симметрии.
    outY.resize(std::min(centerOfSimmetry + 1,size - centerOfSimmetry));
    
    
        switch(featType)
        {
        case ODD_FEAT: // нечетная подгонка
        /*
        Нечетное усреднение
        имеем
        f(-x)=-f(x)
        Тогда среднее есть сумма двух одинаковых значений поделенная пополам.
        В итоге получаем:

        (f(-x)+(-f(x)))/2.0
        Т.е.
        (f(-x)-f(x))/2.0

        Но возникает вопрос со знаком.
        По идее сигнал сюда приходит комбинированный, а значит
        нулевой элемент соответствует индукции -2 Тл.
        Получается нам надо взять знак этого элемента и изменить на противоположный -
        это и будет тот знак, который нам так нужен.

        Если пытаться использовать элементы из середины - велика вероятность ошибки,
        т.к. из-за неоднородностей образца для нулевого магнитного поля идеально
        нулевых значений напряжения Холла получить не удаётся.
        Оно будет близким к нулю, при этом справа или слева - предсказать не получится.
        
        С магнитным полем это действует достаточно хорошо, но вот если речь идет
        о сигнале - возникают некоторые проблемы.
        Поэтому лучше брать знак той точки, которую усредняем. Выглядит логично.


        И вот я нашел ещё один косят - а что если функция не меняет знака?
        На всей протяженности знак одинаковый и тут понеслась....

        Так, мы симметричны относительно точки B=0 и некоторого значения сигнала,
        которое соответствует полю равному 0.
        Вопрос - как нам узнать это значение?
        В общем случае один из краев или центр входного сигнала должен быть
        близким к нулю, но он может быть близким к нему, но так его и не пересечь.

        Возникли проблемы с усреднением фильтрованного графика холловского напряжения.


        */
        for(auto i=0u;i<size/2+1;++i)
            {
                /*
                Есть такая точка - центр симметрии.
                Мы усредняем график относительно этого центра.
                Таким образом получается
                */

                if (centerOfSimmetry+i<size && centerOfSimmetry>=i) // дабы не вывалиться за границы
                {
                    outY[i]=(inY[centerOfSimmetry+i]-inY[centerOfSimmetry-i])/2.0L;

                }
                
                else // Пока значений хватает всё вроде как хорошо
                    // Однако когда они заканчиваются тут начинается простое переписывание оставшихся
                    // Это не очень правильно. Нужно экстраполировать.
                { 
                    /*
                        Нужно определить шаг.
                        Определить сколько точек нужно экстраполировать.
                        Видимо до 2+h Тл
                        Вызывать экстраполяцию.
                        Дописать результаты сюда.
                        Этот подход будет правильней. И без глюков с прыжками.
                    */               

                    /*
                    if (centerOfSimmetry-i<0 && centerOfSimmetry+i<size)
                    {
                        outY[i]=inY[centerOfSimmetry+i];
                    }
                    else
                    {
                        if (centerOfSimmetry-i>=0)
                        {
                            outY[i]=-inY[centerOfSimmetry-i];
                        }                   
                    }*/
                }
/*

            if( sign(inY[0]*sign(inY.back()))<0) // тогда если по краям знаки разные
                outY[i]=-sign(inY[i])*fabs((inY[i]-inY[size-1-i]))/2.0;
            else
                outY[i]=sign(inY[i])*fabs((inY[i]-inY[size-1-i]))/2.0;
*/
                //outY[i]=sign(inY[size-1-i])*fabs((inY[i]-inY[size-1-i]))/2.0;
                //outY[i]=-sign(inY[i])*fabs((inY[i]-inY[size-1-i]))/2.0;
                //outY[size-1-i]=-inY[i];
            }
            // Поскольку начинали с дальнего края - реверсим вектор.
      //  std::reverse(outY.begin(),outY.end());
            break;
        case EVEN_FEAT: // четная подгонка
        /*
        Так, теперь выясняется что и тут не всё так просто.
        Есть у нас четная симметрия.
        Относительно оси Y
        Сие означает, что для противоположных значений индукции магнитного поля значения данной функции равны
        f(-B)=f(B)
        Для обычных входных данных нет никаких проблем и эта простая формула работает.
        Но для фильтрованных сигналов, с увеличением длины фильтра возрастает погрешность - мы усредняем не те точки.
        Нужно усреднять относительно точки В=0
        Т.е. надо:
        1. найти точку В=0
        2. получить её индекс
        3. усреднить относительно её те точки что получится
        4. недостающие (для высоких значений магнитных полей, при длинном фильтре они экстраполируются) взять из исходного сигнала, без усреднения
        P.S. Не уверен насчет поведения недостающих точек, но попробовать стоит.
        5. Поскольку мы живем в мире реальных сигналов, то минимум сигнала (сопротивления) для четной подгонки запросто не совпадает 
        с нулевым значением магнитного поля.
        Для этого вводим новый параметр - индекс элемента, относительно которого строим симметрию.
    

        */

                
        //ShowMessage(FloatToStr(inY[min])+" "+IntToStr((int)min));
            
        for(auto i=0u;i<size/2+1;++i) // Элементов станет в два раза меньше.
        {
            if (centerOfSimmetry+i<size && centerOfSimmetry>=i) // дабы не вывалиться за границы
            {
                outY[i]=(inY[centerOfSimmetry-i]+inY[centerOfSimmetry+i])/2.0L;
                if(fabs(outY[i])<0.000001L)
                    cout << "O No!!!! MagneticFieldDependence.cpp\n";
            }
            /*
            else
            {                
                if (centerOfSimmetry-i<0 && centerOfSimmetry+i<size)
                {
                    outY[i]=inY[centerOfSimmetry+i];
                }
                else
                {
                    if (centerOfSimmetry-i>=0)
                    {
                        outY[i]=inY[centerOfSimmetry-i];
                    }                   
                }
            }*/
            //outY[i]=(inY[i]+inY[size-1-i])/2.0;
            //tempY[size-1-i]=tempY[i];
        }
        break;
    }
}

void MagneticFieldDependence::featData(DataKind dataKind)
{
    TSignal tempInB;
    TSignal tempInBHall;
    TSignal tempInBResistance;
    TSignal tempInHall;
    TSignal tempInResistance;

    // Получаем указатели на нужные нам данные.

    tempInBHall=(*getPointerBHall(dataKind));
    tempInBResistance=(*getPointerBResistance(dataKind));
    tempInB=(*getPointerB(dataKind));
    tempInHall=(*getPointerHall(dataKind));
    tempInResistance=(*getPointerMagnetoResistance(dataKind));

    if(tempInB.size()==0 || tempInBHall.size()==0 || tempInBResistance.size()==0)
    {
        //"Получен пустой массив данных!!!");
        return;
    }
    
    /*
    Ещё один этап весьма важной проверки - если сигнал не является Совмещенным, т.е. существует только для положительного или только для отрицательного поля
    То усреднять вообще говоря нечего.
    */
    if(indexOfElemClosestToValue(tempInBHall,0)<tempInBHall.size()/10  || indexOfElemClosestToValue(tempInBHall,0)==tempInBHall.size()-1)
    {
        AveragedBHall=tempInBHall;
        AveragedBMagnetoResistance=tempInBResistance;
        AveragedB=tempInB;
        AveragedHallEffect=tempInHall;
        AveragedMagnetoResistance=tempInResistance;

        return;
    }

    averageData(tempInBHall,AveragedBHall,ODD_FEAT,indexOfElemClosestToValue(tempInBHall,0));
    averageData(tempInBResistance,AveragedBMagnetoResistance,ODD_FEAT,indexOfElemClosestToValue(tempInBResistance,0));

    averageData(tempInB,AveragedB,ODD_FEAT,indexOfElemClosestToValue(tempInB,0));
    averageData(tempInHall,AveragedHallEffect,ODD_FEAT,indexOfElemClosestToValue(tempInBHall,0));
    averageData(tempInResistance,AveragedMagnetoResistance,EVEN_FEAT,indexOfElemClosestToValue(tempInBResistance,0));
}

void MagneticFieldDependence::GetEqualNumberOfPoints(TSignal & B,
TSignal & BHall,TSignal & BRes, TSignal & Hall,
TSignal & Res)
{
    MyDataType left, right;
    size_t length;

    TSignal tempB;
    TSignal tempBHall;
    TSignal tempBRes;

    TSignal tempHall;
    TSignal tempResistance;
    
    std::vector<size_t> s;
    s.push_back(B.size());
    s.push_back(BHall.size());
    s.push_back(BRes.size());
    s.push_back(Hall.size());
    s.push_back(Res.size());

    std::sort(s.begin(), s.end());

    if (s[0]==0)
    {
        return;
    }

    length=s[0];
    left=B[0];
    right=B.back();

    // это должно уравнивать количество точек, после фильтрации с разной длиной.
    thiningSignal(BHall, Hall, tempBHall, tempHall, left, right, length);
    thiningSignal(BRes, Res, tempBRes, tempResistance, left, right, length);
    
    BHall.clear();
    Hall.clear();
    BRes.clear();
    Res.clear();
    B.clear();
    
    BHall=tempBHall;
    Hall=tempHall;
    BRes=tempBRes;
    Res=tempResistance;
    B=BHall;
}
//-------------------------------------------------------------------------------
void MagneticFieldDependence::filterData()
{

    // очищаем результаты предыдущей фильтрации
    clearFilteredParams();
    // уравниваем количество точек

    GetEqualNumberOfPoints(B,BHall,BMagnetoResistance,
    HallEffect,MagnetoResistance);

    DataKind dataKind;
    dataKind=CURRENT_DATA;

    //cout << B.size() << endl;
    // << BHall.size() << endl;
    //cout << MagnetoResistance.size() << endl;
    if(B.size()>0)
    {
    if( B[0]<-1.0L && B.back() >1.0L) // Если это комбинированный сигнал
    {
        cout <<"inside\n";
     /*   Dependence h(BHall, HallEffect,true);
        Dependence mr(BMagnetoResistance, MagnetoResistance,false);
        h.Average();
        mr.Average();
        //HallEffect.clear();
        //BMagnetoResistance.clear();
        //BHall.clear();
        TSignal temp(h.getY()->begin(),h.getY()->end());
        TSignal tempB(h.getX()->begin(),h.getX()->end());
        //HallEffect=temp;
        //BMagnetoResistance=tempB;
        //MagnetoResistance.clear();
        TSignal temp2(mr.getY()->begin(),mr.getY()->end());
        TSignal tempB2(mr.getX()->begin(),mr.getX()->end());
        //MagnetoResistance=temp2;
        //BHall=tempB2;*/
    featData(CURRENT_DATA);// его надо усреднить
      /*  AveragedBHall=tempB;
        AveragedBMagnetoResistance=tempB2;
        AveragedB=tempB; // надо потом подумать над компромисным вариантом, либо отказаться от него вообще
        AveragedHallEffect=temp;
        AveragedMagnetoResistance=temp2;
*/
    // и фильтровать будем усредненный сигнал.
    dataKind = AVERAGED_DATA;
    // уравниваем количество точек
    GetEqualNumberOfPoints(AveragedB,AveragedBHall,AveragedBMagnetoResistance,
    AveragedHallEffect,AveragedMagnetoResistance);
	}
    //cout << "Helper!\n";
    filterDataHelper((*filterParamsHall),HALL_EFFECT,dataKind);

    filterDataHelper((*filterParamsResistance),MAGNETORESISTANCE,dataKind);
    FilteredB=FilteredBHall;

    GetEqualNumberOfPoints(FilteredB,FilteredBHall,FilteredBMagnetoResistance,
    FilteredHallEffect,FilteredMagnetoResistance);
    }
    else
    {
        cout << "B.size() is " << B.size() << endl;
    }
}
//-------------------------------------------------------------------------------
void MagneticFieldDependence::filterData(FilterParams &fPHall, FilterParams &fPResistance)
{
    /*
    Вызывает предыдущую, после установки новых параметров фильтрации.
    */
    setFilterParamsHall(fPHall.SamplingFrequecy, fPHall.BandwidthFrequency, fPHall.AttenuationFrequency, fPHall.filterLength);
    setFilterParamsResistance(fPResistance.SamplingFrequecy, fPResistance.BandwidthFrequency, fPResistance.AttenuationFrequency, fPResistance.filterLength);
    filterData();

}
//-------------------------------------------------------------------------------
void MagneticFieldDependence::filterDataHelper(FilterParams &fP,
    SignalType dependenceType, DataKind dataKind)
{
    /*
    Так, функция у нас по заданному типу зависимостей выполняет фильтрацию.
    Стало быть есть указатели на 
    Магнитное поле и сигнал
    */

    // Указатели на поле и сигнал.
    TSignal * inB = nullptr;
    TSignal * inSignal = nullptr;

    // Заботу о том комбинированный сигнал или нет, перекладываем на вышестоящую функцию.
    // Получаем нужные зависимости, с учетом типа и вида сигнала.
    switch(dependenceType)
    {
    case HALL_EFFECT:
        inB= getPointerBHall(dataKind);
        inSignal = getPointerHall(dataKind);
        //cout << "signal Hall size " << inSignal->size() << endl;
        break;

    case MAGNETORESISTANCE:
        inB= getPointerBResistance(dataKind);
        inSignal = getPointerMagnetoResistance(dataKind);
        //cout << "signal MR size " << inSignal->size() << endl;
        break;

    default:
        break;
    }

    unsigned int NumberOfPoints; // Количество точек, пригодится для скорости.

    if(inSignal->size()!=inB->size())
    {
        cout << "Количество точек магнитного поля и сигнала не совпадает! filterdataHelper.\n";
        cout << inSignal->size() << " " << inB->size() << endl;
        return;
    }

    NumberOfPoints=inSignal->size();

    if(NumberOfPoints==0)
    {
        cout << "Получен пустой массив данных! filterdataHelper.\n";
        return;
    }

    if(2*NumberOfPoints<fP.filterLength)
    {
        cout << "Слишком длинный фильтр! filterdataHelper.\n";
        return;
    }

    // В эти массивы будут достраиваться данные для отрицательных магнитных полей.
    // Это очень мило, а если это сигнал для отрицательного магнитного поля?
    // То теоретически достраивается положительная часть.
    // Надо пофиксить тут комменты на адекватные действительности.
    TSignal tempInB(2*NumberOfPoints);
    TSignal tempInSignal(2*NumberOfPoints);

    TSignal tempOutB(2*NumberOfPoints);
    TSignal tempOutSignal(2*NumberOfPoints);

    switch(dependenceType)
    {
    case HALL_EFFECT:


    // формируем сигнал для фильтра.
    // достраивая его в отрицательные магнитные поля.
    for (auto i = 0u; i < NumberOfPoints; ++i)
    {
    /*
    Давайте внимательно сюда посмотрим.
    У эффекта Холла отрицательная симметрия, относительно точки
    B==0;
    С чего вообще я взял, что это нулевой элемент? /(О_о)\

    Получается для сигнала с положительным магнитным полем - это выполняется.
    Для сигнала комбинированного, т.е. уже объединенного - это выполняется,
    потому что фильтруется усредненный сигнал (по сути имеющий только значения
    положительного магнитного поля.

    Для отрицательного магнитного поля сие равенство, насколько мне ясно - не выполняется.
    

    */
        tempInSignal[i]=-(*inSignal)[NumberOfPoints-i-1]+2*(*inSignal)[0];
        tempInB[i]=-(*inB)[NumberOfPoints-i-1];
        tempInSignal[i+NumberOfPoints]=(*inSignal)[i];
        tempInB[i+NumberOfPoints]=(*inB)[i];        
    }
    break;
    case MAGNETORESISTANCE:


    for (auto i = 0u; i < NumberOfPoints; ++i)
    {
        tempInSignal[i]=(*inSignal)[NumberOfPoints-i-1];   // чет
        tempInB[i]=-(*inB)[NumberOfPoints-i-1];        
        tempInSignal[i+NumberOfPoints]=(*inSignal)[i];
        tempInB[i+NumberOfPoints]=(*inB)[i];        
    }
    break;

    default:
        break;
    }

    /*
    В случае отрицательного магнитного поля надо инвертировать порядок элементов
    Потому что впереди выстраиваются значения для положительного магнитного поля.
    */
    if(tempInB[0]>1.0L)
    {
        std::reverse(tempInB.begin(),tempInB.end());
        std::reverse(tempInSignal.begin(),tempInSignal.end());
    }
    // фильтруем 
    TrForMassiveFilter(tempInB,tempInSignal,tempOutB,tempOutSignal,
                fP.filterLength,fP.SamplingFrequecy,fP.BandwidthFrequency,fP.AttenuationFrequency);


    // Размер внутри фильтра меняется, в зависимости от длины фильтра.
    NumberOfPoints=tempOutB.size();
    for(unsigned int i=fP.filterLength;i<NumberOfPoints;++i)
    {
        
    switch(dependenceType)
    {
    case HALL_EFFECT:
        FilteredHallEffect.push_back(tempOutSignal[i]);
        FilteredBHall.push_back(tempOutB[i]);
        break;
    case MAGNETORESISTANCE:
        FilteredMagnetoResistance.push_back(tempOutSignal[i]);
        FilteredBMagnetoResistance.push_back(tempOutB[i]);
        break;
    default:
        break;
    }
    }
}
//-------------------------------------------------------------------------------
static bool abs_compare(long double a, long double b)
{
    return (std::fabs(a) < std::fabs(b));
}
//-------------------------------------------------------------------------------

bool MagneticFieldDependence::extrapolateData( DataKind dataKind, const int polinomPowForMagnetoResistance, const int polinomPowForHallEffect)
{
	/*
	Раз уж пошли такие баги - буду разбираться комментированием кода.
	Имеются три указателя на экстраполируемые данные.
	В зависимости от того какие нужно, это будут либо фильтрованные, либо усредненные данные.
	Если вдруг что-то пошло не так - мы убегаем.
	Экстраполяция проходит только по части сигнала - на данный момент по последним двум третям.

	Потом проходит увеличение веса точки 0, 0 для холловского напряжения.

	Выделяем память для коэффициентов и вызываем функцию поиска решения ( curveFittingUniversal ).
	Эта функция кстати тоже не безгрешна, но там баги искать это совсем страшное зло.

	Магнитное поле заполняется заданным шагом (кстати, возможно шаг нужно рассчитывать как-то по другому), начиная с последнего текущего значения.

	*/
    bool returnValue=true;

    TSignal * BHallToExtrapolate=nullptr;
    TSignal * BResToExtrapolate=nullptr;
    TSignal * BToExtrapolate= nullptr;
    TSignal * HallToExtrapolate= nullptr;
    TSignal * ResistanceToExtrapolate= nullptr;

    if (dataKind==FILTERED_DATA)
    {
        BToExtrapolate= & FilteredB;
        BHallToExtrapolate=& FilteredBHall;
        BResToExtrapolate=& FilteredBMagnetoResistance;
        HallToExtrapolate= & FilteredHallEffect;
        ResistanceToExtrapolate= & FilteredMagnetoResistance;
    }
    if(dataKind ==AVERAGED_DATA)
    {
        BToExtrapolate= & AveragedB;
        BHallToExtrapolate=& AveragedBHall;
        BResToExtrapolate=& AveragedBMagnetoResistance;
        HallToExtrapolate= & AveragedHallEffect;
        ResistanceToExtrapolate= & AveragedMagnetoResistance;
    }   
    // Если вдруг что-то пошло не так - мы убегаем.
    if(BToExtrapolate==nullptr)
    {
    	return false;
    } 

    TSignal inBHall;
    TSignal inBMagnetoResistance;

    TSignal inHallEffect;
    TSignal inMagnetoResistance;

    /*
    Экстраполяция данных ведется по последней четверти фильтрованных значений.

    Вот тут собственно и косяк.
    Экстраполировать нужно либо по фильтрованным, либо по усредненным значениям.
    */
    for (size_t i = (16*BToExtrapolate->size())/20; i < BToExtrapolate->size(); ++i)
    {
        inBHall.push_back( (*BToExtrapolate)[i]);
        inBMagnetoResistance.push_back( (*BToExtrapolate)[i]);
        inHallEffect.push_back( (*HallToExtrapolate)[i]);
        inMagnetoResistance.push_back( (*ResistanceToExtrapolate)[i]);
    }

    auto NumberOfPoints=inBHall.size();
    if(NumberOfPoints==0)
    {
    cout <<"Количество точек равно нулю! Я не хочу делить на ноль:) extrapolateData\n";
    return false;
    }

    MyDataType h=(*max_element(inBHall.begin(),inBHall.end(),abs_compare)
                  - *min_element(inBHall.begin(),inBHall.end()))
            /NumberOfPoints;

    TSignal koefMagnetoResistance(polinomPowForMagnetoResistance+1);
    TSignal koefHallEffect(polinomPowForHallEffect+1);

    if(!curveFittingUniversal(&inBMagnetoResistance,&inMagnetoResistance, &koefMagnetoResistance,polinomPowForMagnetoResistance))
    return false;
    if(!curveFittingUniversal(&inBHall,&inHallEffect, &koefHallEffect,polinomPowForHallEffect))
    return false;


    TSignal newB;
    TSignal newBHall;
    TSignal newBRes;
    TSignal newHallEffect;
    TSignal newMagnetoResistance;

    // Допустим нам всегда надо доводить поле до rightFieldBoundary
    long double rightFieldBoundary = 2.5L; // Тл
    // Тогда нам надо расчитать кол-во точек, исходя из среднего шага.
    if (rightFieldBoundary>(BToExtrapolate->back()-h))
    {
        NumberOfPoints=(rightFieldBoundary-(BToExtrapolate->back()-h))/h+1;
    }
    

    newB.push_back(BToExtrapolate->back()-h); // заполняем магнитное поле.
    newBHall.push_back(BHallToExtrapolate->back()-h);
    newBRes.push_back(BResToExtrapolate->back()-h);
    for (unsigned int i = 1; i<NumberOfPoints; ++i) {
        newB.push_back(newB[i-1]+h);
    }
    for (unsigned int i = 1; i<inBHall.size(); ++i) {
        newBHall.push_back(newBHall[i-1]+h);
    }
    for (unsigned int i = 1; i<inBMagnetoResistance.size(); ++i) {
        newBRes.push_back(newBRes[i-1]+h);
    }

    // вычисляем экстраполированные зависимости.

    calculatePolinomByKoef(newBRes,koefMagnetoResistance,newMagnetoResistance);
    calculatePolinomByKoef(newBHall,koefHallEffect,newHallEffect);

    ExtrapolatedB.clear();
    ExtrapolatedMagnetoResistance.clear();
    ExtrapolatedHallEffect.clear();

    ExtrapolatedB=newB;
    ExtrapolatedMagnetoResistance=newMagnetoResistance;
    ExtrapolatedHallEffect=newHallEffect;
    
    //----------А вот тут прикручиваем недостающий кусочек в сигналы----
    
    unsigned int i=0;
    while(i<BToExtrapolate->size() && newB[i]<=BToExtrapolate->back())
    ++i; // находим элемент на котором заканчивается фильтрованное магнитное поле.

    for(unsigned int j=i;j<NumberOfPoints;j++)
    {     // в конце дописываем экстраполированные значения.
        BToExtrapolate->push_back(newB[j]);
        BHallToExtrapolate->push_back(newBHall[j]);
        BResToExtrapolate->push_back(newBRes[j]);
        ResistanceToExtrapolate->push_back(newMagnetoResistance[j]);
        HallToExtrapolate->push_back(newHallEffect[j]);
    }
    //------------------------------------------------------------------
    inBHall.clear();
    inBMagnetoResistance.clear();
    inHallEffect.clear();
    inMagnetoResistance.clear();


    // Эта часть повторяется ради экстраполяции слева
    /*
    Экстраполяция данных ведется по последней четверти фильтрованных значений.

    Вот тут собственно и косяк.
    Экстраполировать нужно либо по фильтрованным, либо по усредненным значениям.
    */
    
    for (size_t i = 0; i < 4*BToExtrapolate->size()/20; ++i)
    {
        inBHall.push_back( (*BToExtrapolate)[i]);
        inBMagnetoResistance.push_back( (*BToExtrapolate)[i]);
        inHallEffect.push_back( (*HallToExtrapolate)[i]);
        inMagnetoResistance.push_back( (*ResistanceToExtrapolate)[i]);
    }

    NumberOfPoints=inBHall.size();
    if(NumberOfPoints==0)
    {
    cout << "Количество точек равно нулю! Я не хочу делить на ноль:) extrapolateData\n";
    return false;
    }
    h = (*max_element(inBHall.begin(),inBHall.end())
         - *min_element(inBHall.begin(),inBHall.end()))
            /NumberOfPoints;

    koefMagnetoResistance.clear();
    koefHallEffect.clear();

    koefMagnetoResistance.resize(polinomPowForMagnetoResistance+1);
    koefHallEffect.resize(polinomPowForHallEffect+1);

    if(!curveFittingUniversal(&inBMagnetoResistance,&inMagnetoResistance, &koefMagnetoResistance,polinomPowForMagnetoResistance))
    return false;
    if(!curveFittingUniversal(&inBHall,&inHallEffect, &koefHallEffect,polinomPowForHallEffect))
    return false;

    newB.clear();
    newBHall.clear();
    newBRes.clear();
    newHallEffect.clear();
    newMagnetoResistance.clear();


    // Допустим нам всегда надо доводить поле до rightFieldBoundary Тл
    // Тогда нам надо расчитать кол-во точек, исходя из среднего шага.
    if (rightFieldBoundary>fabs(BToExtrapolate->front()+h))
    {
        NumberOfPoints=(rightFieldBoundary+(BToExtrapolate->front()+h))/h+1;
    }

    newB.push_back(BToExtrapolate->front()+h); // заполняем магнитное поле.
    newBHall.push_back(BHallToExtrapolate->front()+h);
    newBRes.push_back(BResToExtrapolate->front()+h);
    for (auto i = 1u; i<NumberOfPoints; ++i) {
        newB.push_back(newB[i-1]-h);
    }
    for (auto i = 1u; i<NumberOfPoints; ++i) {
        newBHall.push_back(newBHall[i-1]-h);
    }
    for (auto i = 1u; i<NumberOfPoints; ++i) {
        newBRes.push_back(newBRes[i-1]-h);
    }

    // вычисляем экстраполированные зависимости.
    calculatePolinomByKoef(newBRes,koefMagnetoResistance,newMagnetoResistance);
    calculatePolinomByKoef(newBHall,koefHallEffect,newHallEffect);


   // ExtrapolatedB=newB;
   // ExtrapolatedMagnetoResistance=newMagnetoResistance;
   // ExtrapolatedHallEffect=newHallEffect;

    //----------А вот тут прикручиваем недостающий кусочек в сигналы----
    std::reverse(newB.begin(),newB.end());
    std::reverse(newBHall.begin(),newBHall.end());
    std::reverse(newBRes.begin(),newBRes.end());
    std::reverse(newHallEffect.begin(),newHallEffect.end());
    std::reverse(newMagnetoResistance.begin(),newMagnetoResistance.end());

    //i=0;
    //while(i<NumberOfPoints && newB[i]<=BToExtrapolate->front())
    //++i; // находим элемент на котором начинается фильтрованное магнитное поле.

    for(unsigned int k=0;k<BToExtrapolate->size();++k)
    {     // в конце дописываем экстраполированные значения.
        newB.push_back(BToExtrapolate->operator[](k));
        newBHall.push_back(BHallToExtrapolate->operator[](k));
        newBRes.push_back(BResToExtrapolate->operator[](k));
        newHallEffect.push_back(HallToExtrapolate->operator[](k));
        newMagnetoResistance.push_back(ResistanceToExtrapolate->operator[](k));
    }

    BToExtrapolate->clear();
    BHallToExtrapolate->clear();
    BResToExtrapolate->clear();
    ResistanceToExtrapolate->clear();
    HallToExtrapolate->clear();

    for (auto i = 0u; i < newB.size(); ++i)
    {
        BToExtrapolate->push_back(newB[i]);
        BHallToExtrapolate->push_back(newBHall[i]);
        BResToExtrapolate->push_back(newBRes[i]);
        ResistanceToExtrapolate->push_back(newMagnetoResistance[i]);
        HallToExtrapolate->push_back(newHallEffect[i]);
    }
        
    /* TStringList * tsl= new TStringList;
    for (int i=0; i<FilteredMagnetoResistance.size(); ++i)
    {
        tsl->Add(FloatToStr(FilteredB[i]) + "\t"+ FloatToStr(FilteredMagnetoResistance[i]));
    }
    tsl->SaveToFile("E:\\Дела\\Институт физики полупроводников\\Lcard\\MCTParamsExplorer\\filterOutPutAftExtr.txt");
*/

return returnValue;   
}

//-------------------------------------------------------------------------------
void MagneticFieldDependence::multiplyB(DataKind dataKind)
{
    TSignal * temp;
    temp=getPointerB(dataKind);

    TSignal::iterator pos;
    for(pos=temp->begin();pos!=temp->end();++pos)
    *pos*=10;
}  
//-------------------------------------------------------------------------------
void MagneticFieldDependence::multiplySignal(SignalType s, MyDataType x)
{
TSignal * p = nullptr;
if(s==HALL_EFFECT) p=&HallEffect;
if(s==MAGNETORESISTANCE) p=&MagnetoResistance;
if(p)
{
    TSignal::iterator begin=p->begin();
    TSignal::iterator end=p->end();

    for(;begin!=end;++begin)
        *begin*=x;
}
}

//-------------------------------------------------------------------------------
inline void MagneticFieldDependence::ReplaceDotsToComma(std::string &in, std::string & out)
{
    unsigned int findIndex=0;
	std::string s=in;
	std::string strToReplaceWhich=","; // на что меняем
	std::string strToSearch=".";   // что ищем
	while ((findIndex=s.find(strToSearch,findIndex))!=std::string::npos)
    {
        s.replace(s.begin()+findIndex,s.begin()+findIndex+strToSearch.length(),
        strToReplaceWhich.begin(),strToReplaceWhich.end());
    } 
    out=s;
}

//-------------------------------------------------------------------------------
bool MagneticFieldDependence::setFilterParamsHall(string samplingFrequecy, string bandwidthFrequency,
                                                  string attenuationFrequency, string length)
{
    filterParamsHall->setFilterParams(stold(samplingFrequecy), stold(bandwidthFrequency),
                                      stold(attenuationFrequency), stoll(length));
    return true;
}
//-------------------------------------------------------------------------------
bool MagneticFieldDependence::setFilterParamsHall(MyDataType samplingFrequecy,
                                                  MyDataType bandwidthFrequency,MyDataType attenuationFrequency, long length)
{
    filterParamsHall->setFilterParams(samplingFrequecy, bandwidthFrequency, attenuationFrequency, length);
    return true;
}
//-------------------------------------------------------------------------------
bool MagneticFieldDependence::setFilterParamsResistance(std::string samplingFrequecy,std::string bandwidthFrequency,std::string attenuationFrequency, std::string length)
{
    filterParamsResistance->setFilterParams(stold(samplingFrequecy), stold(bandwidthFrequency), stold(attenuationFrequency), stoll(length));
    return true;
}
//-------------------------------------------------------------------------------
bool MagneticFieldDependence::setFilterParamsResistance(MyDataType samplingFrequecy,
                                                        MyDataType bandwidthFrequency,
                                                        MyDataType attenuationFrequency, long length)
{
    filterParamsResistance->setFilterParams(samplingFrequecy, bandwidthFrequency, attenuationFrequency, length);
    return true;
}
//-------------------------------------------------------------------------------
void MagneticFieldDependence::setRoundNeeded(bool needRound)
{
    saver->setRoundNeeded(needRound);
}

//-------------------------------------------------------------------------------


void MagneticFieldDependence::setDependence(TSignal::iterator beginB, 
        TSignal::iterator endB, TSignal::iterator beginHall, 
        TSignal::iterator beginResistance)
{
    clearCurrentParams();

    for ( ; beginB != endB; ++beginB, ++beginResistance, ++beginHall)
    {
        B.push_back(*beginB);
        MagnetoResistance.push_back(*beginResistance);
        HallEffect.push_back(*beginHall); 
    }

    if(fabs(*max_element(B.begin(),B.end(),abs_compare))<0.5L)
    {
        multiplyB(CURRENT_DATA);
    }
    BHall=B;
    BMagnetoResistance=B;


    filterData();
    extrapolateData(FILTERED_DATA,PowPolinomRes,PowPolinomHall);
    // в перспективе степень будет зависеть от температуры и возможно чего-нибудь ещё.      
}

void MagneticFieldDependence::setSampleDescription(MyDataType Temperature, MyDataType Current, std::string SampleInventoryNumber, MyDataType length, MyDataType width, MyDataType Thickness)
{
    saver->setSampleDescription(to_string(Temperature), to_string(Current), SampleInventoryNumber, to_string(length), to_string(width), to_string(Thickness));
}

void MagneticFieldDependence::setSampleDescription(std::string Temperature, std::string Current, std::string SampleInventoryNumber,
    std::string length, std::string width, std::string Thickness)
{
    saver->setSampleDescription(Temperature, Current, SampleInventoryNumber, length, width, Thickness);
}
//-------------------------------------------------------------------------------

FilterParams const * MagneticFieldDependence::getFilterParamsHall()
{
    return filterParamsHall;
}
//------------------------------------------------------------------------------
FilterParams const * MagneticFieldDependence::getFilterParamsResistance()
{
    return filterParamsResistance;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getB()
{
    return &B;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getFilteredB()
{
    return &FilteredB;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getExtrapolatedB()
{
    return &ExtrapolatedB;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getHallEffect()
{
    return &HallEffect;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getMagnetoResistance()
{
    return &MagnetoResistance;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getFilteredHallEffect()
{
    return &FilteredHallEffect;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getFilteredMagnetoResistance()
{
    return &FilteredMagnetoResistance;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getExtrapolatedHallEffect()
{
    return &ExtrapolatedHallEffect;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getExtrapolatedMagnetoResistance()
{
    return &ExtrapolatedMagnetoResistance;
}
//-------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getSxx()
{
    return &sxx;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getSxy()
{
    return &sxy;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getAveragedB()
{
    return &AveragedB;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getRh_eff()
{
    return &Rh_eff;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getS_eff()
{
    return &s_eff;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getMobility()
{
    return &mobility;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getHoleConductivity()
{
    return &holeConductivity;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getElectronConductivity()
{
    return &electronConductivity;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getHoleConcentration()
{
    return &holeConcentration;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getHoleMobility()
{
    return &holeMobility;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getElectronConcentration()
{
    return &electronConcentration;
}
//----------------------------------------------------------------------------------
TSignal const * MagneticFieldDependence::getElectronMobility()
{
    return &electronMobility;
}
//----------------------------------------------------------------------------------

void MagneticFieldDependence::calculateEffectiveParamsFromSignals()
{
    /*
    Расчет производится по усредненным значениям.
    */
    if (AveragedB.size() && AveragedMagnetoResistance.size() && AveragedHallEffect.size())
    {
        std::string T;
        std::string Current;
        std::string SampleLength;
        std::string SampleWidth;
        std::string SampleThickness;
        std::string InventoryNumber;

        saver->getSampleDescription(T, Current, InventoryNumber, SampleLength, SampleWidth, SampleThickness);
        MyDataType I=stold(Current)*1E-6L;
        MyDataType SLength=stold(SampleLength)*1E-3L;
        MyDataType SWidth=stold(SampleWidth)*1E-3L;
        MyDataType SThickness=stold(SampleThickness)*1E-6L;

        auto NumberOfPoints = AveragedB.size();
        s_eff.resize(NumberOfPoints);
        Rh_eff.resize(NumberOfPoints);

        for (auto i = 0u; i < NumberOfPoints ; ++i)
        {
            if(fabs(AveragedMagnetoResistance[i])<THEALMOSTZERO)
                s_eff[i]=0;
            else
            {
                s_eff[i]=SLength/SWidth/SThickness*I/AveragedMagnetoResistance[i];
            }
            if(AveragedB[i] == 0.0L)
                Rh_eff[i] = 0;
            else
                Rh_eff[i]=SThickness*AveragedHallEffect[i]/I;
        }
    }
}

//--------------------------------------------------------------------------------

void MagneticFieldDependence::calculateTenzorFromEffectiveParams()
{
    if (s_eff.size() && Rh_eff.size() && s_eff.size() == Rh_eff.size())
    {
        auto NumberOfPoints = s_eff.size();
        sxx.resize(NumberOfPoints);
        sxy.resize(NumberOfPoints);

        for (auto i = 0u; i < NumberOfPoints ; ++i)
        {
            sxx[i]=s_eff[i]/
                (Rh_eff[i]*Rh_eff[i]*s_eff[i]*s_eff[i]+1.0L);
            sxy[i]=s_eff[i]*s_eff[i]*Rh_eff[i]/
               (Rh_eff[i]*Rh_eff[i]*s_eff[i]*s_eff[i]+1.0L);
        }
    }
}
//-------------------------------------------------------------------------
TSignal * MagneticFieldDependence::getPointerB(DataKind dataKind)
{
    switch(dataKind)
    {
    case CURRENT_DATA:
        return &B;
    case FILTERED_DATA:
        return &FilteredB;
    case EXTRAPOLATED_DATA:
        return &ExtrapolatedB;
    case AVERAGED_DATA:
        return &AveragedB;
    default:
        break;
    }
    return nullptr;
}

TSignal * MagneticFieldDependence::getPointerBHall(DataKind dataKind)
{
    switch(dataKind)
    {
    case CURRENT_DATA:
        return &BHall;
    case FILTERED_DATA:
        return &FilteredBHall;
    case EXTRAPOLATED_DATA:
        return &ExtrapolatedBHall;
    case AVERAGED_DATA:
        return &AveragedBHall;
    default:
        break;
    }
    return nullptr;
}

TSignal * MagneticFieldDependence::getPointerBResistance(DataKind dataKind)
{
    switch(dataKind)
    {
    case CURRENT_DATA:
        return &BMagnetoResistance;
    case FILTERED_DATA:
        return &FilteredBMagnetoResistance;
    case EXTRAPOLATED_DATA:
        return &ExtrapolatedBMagnetoResistance;
    case AVERAGED_DATA:
        return &AveragedBMagnetoResistance;
    default:
        break;
    }
    return nullptr;
}
//-------------------------------------------------------------------------
TSignal * MagneticFieldDependence::getPointerHall(DataKind dataKind)
{
    switch(dataKind)
    {
    case CURRENT_DATA:
        return &HallEffect;
    case FILTERED_DATA:
        return &FilteredHallEffect;
    case EXTRAPOLATED_DATA:
        return &ExtrapolatedHallEffect;
    case AVERAGED_DATA:
        return &AveragedHallEffect;
    default:
        break;
    }
    return nullptr;
}
//-------------------------------------------------------------------------
TSignal * MagneticFieldDependence::getPointerMagnetoResistance(DataKind dataKind)
{
    switch(dataKind)
    {
    case CURRENT_DATA:
        return &MagnetoResistance;
    case FILTERED_DATA:
        return &FilteredMagnetoResistance;
    case EXTRAPOLATED_DATA:
        return &ExtrapolatedMagnetoResistance;
    case AVERAGED_DATA:
        return &AveragedMagnetoResistance;
    default:
        break;
    }
    return nullptr;
}
//-------------------------------------------------------------------------
TSignal * MagneticFieldDependence::getPointerSxx(/*DataKind dataKind*/)
{
    return &sxx;
}
//-------------------------------------------------------------------------
TSignal * MagneticFieldDependence::getPointerSxy(/*DataKind dataKind*/)
{
    return &sxy;
}
//-------------------------------------------------------------------------

void MagneticFieldDependence::calculateTenzor(DataKind dataKind)
{
    featData(dataKind);
    calculateEffectiveParamsFromSignals();
    calculateTenzorFromEffectiveParams();
}

//-------------------------------------------------------------------------
void  MagneticFieldDependence::cutData(/*DataKind dataKind*/)
{
    // А если тут ничего не происходит... зачем я её делал?))))
    // оставляет только положительные значения магнитного поля
    // Хм, надо будет разобраться надо ли оно нам. По идее - её нужно бы использовать для обработки
    // сигналов положительных и отрицательных.
   ;
}

void MagneticFieldDependence::clearCurrentParams()
{
    B.clear();
    BHall.clear();
    BMagnetoResistance.clear();
    HallEffect.clear();
    MagnetoResistance.clear();    
}

void MagneticFieldDependence::clearFilteredParams()
{
    FilteredB.clear();
    FilteredBHall.clear();
    FilteredBMagnetoResistance.clear();
    FilteredHallEffect.clear();
    FilteredMagnetoResistance.clear();
}

void MagneticFieldDependence::setParamsType(ParamsType pt)
{
    paramsType=pt;
}

void MagneticFieldDependence::shiftCurve(DataKind dataKind,SignalType dependenceType,MyDataType shiftValue,MyDataType leftBound, MyDataType rightBound)
{
    TSignal * pointToY=nullptr;
    TSignal * pointToX=nullptr;
    switch (dependenceType)
    {
        case HALL_EFFECT:
            pointToY=getPointerHall(dataKind);
            pointToX=getPointerB(dataKind);
            break;
        case MAGNETORESISTANCE:
            pointToY=getPointerMagnetoResistance(dataKind);
            pointToX=getPointerB(dataKind);
            break;
        default:
        break;
    }
    if (pointToY==nullptr)
    {
        return;
    }

    for (auto i = pointToY->begin(), j=pointToX->begin(); i != pointToY->end() && j!=pointToX->end(); ++i,++j)
    {
        if (*j>=leftBound && *j<=rightBound)
        {
            *i+=shiftValue;    
        }        
    }
}

void MagneticFieldDependence::setExtrapolateParams(int powPolinowHall,int powPolinomRes)
{
    PowPolinomRes=powPolinowHall;
    PowPolinomHall=powPolinomRes;
}

// Спектр подвижности
bool MagneticFieldDependence::calculateMobilitySpectrum(unsigned int NumberOfDecPlaces)
{
    unsigned int NumberOfDecimalPlaces = 5;
    TSignal nB;
    TSignal nSxx;
    TSignal nSxy;

    thiningSignal(AveragedB, sxx, nB, nSxx,0, 2, 11);
    thiningSignal(AveragedB, sxy, nB, nSxy,0, 2, 11);

    roundM(nB, NumberOfDecimalPlaces);
    roundM(nSxx, NumberOfDecimalPlaces);
    roundM(nSxy, NumberOfDecimalPlaces);

    if (MobilitySpectrumObj!=nullptr)
    {
        delete MobilitySpectrumObj;
    }

    MobilitySpectrumObj= new mobilitySpectrum(nB,nSxx,nSxy,nB.size());
    MobilitySpectrumObj->getResults(mobility, holeConductivity, electronConductivity);
    MobilitySpectrumObj->getExtremums(holeConcentration, holeMobility, electronConcentration, electronMobility);
    MobilitySpectrumObj->calculatePeaksWeigth("mobSpecCriteria.txt");
    return true;
}
bool MagneticFieldDependence::runSmartCalcutation()
{
    smartCalculation sC(this);
    sC.processData();
    sC.saveResults("smartCalculationResults.txt");
    upBound=sC.getUpBound();
    lowBound=sC.getLowBound();
    return true;
}

// Многозонная подгонка
bool MagneticFieldDependence::runSmartMultiCarrierFit(long double VesGxx, long double VesGxy)
{
    srand(time(nullptr));
    // тут границы для поиска кривые... надо исправлять.
    std::vector<long double> LowBound=lowBound; // сюда - границы для поиска
    std::vector<long double> UpBound=upBound;

    // Получаем сам спектр и компоненты тензора
    InDataSpectr MagSpectr(AveragedB.begin(),AveragedB.end());
    InDataSpectr GxxIn(sxx.begin(),sxx.end());
    InDataSpectr GxyIn(sxy.begin(),sxy.end());

    thiningSignal(MagSpectr, GxxIn, nMagSpectr, nGxxIn,0, 2, 11);
    thiningSignal(MagSpectr, GxyIn, nMagSpectr, nGxyIn,0, 2, 11);

    // auto GxxSize=nMagSpectr.size();
    int numberOfCarrierTypes=3;

    MultiZoneFit * mzf=new MultiZoneFit(100,VesGxx, VesGxy);

    mzf->RunMultizoneFeat(LowBound,  UpBound,
                          nMagSpectr,  nGxxIn, nGxyIn,
                           outGxx,  outGxy,
                           outValues,
                           numberOfCarrierTypes);
    
    delete mzf;
    return true;
}

// Многозонная подгонка
bool MagneticFieldDependence::runMultiCarrierFit(long double VesGxx, long double VesGxy)
{
    srand(time(nullptr));
    // важный вопрос - в каком порядке нужно помещать сюда данные.
    // Порядок такой: подвижность электронов, легких дырок и тяжелых дырок.
    // Далее - концентрации в том же порядке.
    std::vector<long double> ExactBound; // сюда нужны пики из спектра

    ExactBound.push_back(-electronMobility[0]); // Подвижность электронов
    ExactBound.push_back(holeMobility[1]); // Подвижность легких дырок
    ExactBound.push_back(holeMobility[0]); // Подвижность тяжелых дырок

    ExactBound.push_back(-electronConcentration[0]); // Концентрация электронов
    ExactBound.push_back(holeConcentration[1]); // Концентрация легких дырок
    ExactBound.push_back(holeConcentration[0]); // Концентрация тяжелых дырок


    /* Необходимо проверить все ли границы установлены корректно.
    Если подвижность пика равна минимальной исследуемой подвижности - пик не был найден.
    В этом случае необходимо задать границы поиска исходя из другой априорной информации.*/

    /*
    В теории на данный момент спектры не имеющие пика не допускаются к расчету.
    */
    long double leftBoundMobility=0.001L;

    if (ExactBound[0]==leftBoundMobility ||
        ExactBound[1]==leftBoundMobility ||
        ExactBound[2]==leftBoundMobility)
    {
        /* code */
    }

    // Нужно выяснить отчего зависит разброс.
    // Разброс нужно сделать зависимым от статистики спектра подвижности.
    // Полагаю это какие-то коэффициенты
    std::vector<long double> coef = {
        0.7L, // Подвижность электронов
        0.3L, // Подвижность легких дырок
        0.2L, // Подвижность тяжелых дырок
        0.9L, // Концентрация электронов
        0.3L, // Концентрация легких дырок
        0.2L, // Концентрация тяжелых дырок
    };

    std::vector<long double> LowBound; // сюда - границы для поиска
    std::vector<long double> UpBound;

    for (auto i = 0u; i < ExactBound.size(); ++i)
    {
        LowBound.push_back(ExactBound[i]-ExactBound[i]*coef[i]);
        UpBound.push_back(ExactBound[i]+ExactBound[i]*coef[i]);
    }

    // Получаем сам спектр и компоненты тензора
    InDataSpectr MagSpectr(AveragedB.begin(),AveragedB.end());
    InDataSpectr GxxIn(sxx.begin(),sxx.end());
    InDataSpectr GxyIn(sxy.begin(),sxy.end());

    


    thiningSignal(MagSpectr, GxxIn, nMagSpectr, nGxxIn,0, 2, 11);
    thiningSignal(MagSpectr, GxyIn, nMagSpectr, nGxyIn,0, 2, 11);

    

    // int GxxSize=nMagSpectr.size();
    int numberOfCarrierTypes=3;

    MultiZoneFit * mzf=new MultiZoneFit(100,VesGxx, VesGxy);


    int result=mzf->RunMultizoneFeat(LowBound,  UpBound,
                          nMagSpectr,  nGxxIn, nGxyIn,
                           outGxx,  outGxy,
                           outValues,
                           numberOfCarrierTypes);
    
    delete mzf;
    return true;
}

void MagneticFieldDependence::getMultiCarrierFitResults(InDataSpectr & nMagSpectr, InDataSpectr & nGxxIn, InDataSpectr & nGxyIn, 
    MyData_spektr & outGxx, MyData_spektr &  outGxy, TStatistic & outValues)
{
    nMagSpectr=this->nMagSpectr;
    nGxxIn=this->nGxxIn;
    nGxyIn=this->nGxyIn;
    outGxx=this->outGxx;
    outGxy=this->outGxy;
    outValues=this->outValues;
}

void MagneticFieldDependence::getExtrapolateParams(int & powPolinomH,int & powPolinomR)
{
    powPolinomH=PowPolinomHall;
    powPolinomR=PowPolinomRes;
}
