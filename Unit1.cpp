
/*

//-------------------Открытие файла------------------------------------------
// http://www.cplusplus.com/doc/tutorial/files/
void  openFileWithSignal(AnsiString filename)
{
        TStringList *tts=new TStringList();  // сюда будем загружать из файла
        tts->LoadFromFile(filename);// загрузили

        MagneticFieldDependence * p=InitParams();
        if (p)
        {
            StatusBar->Panels->Items[2]->Text=filename;
            TStringList *Names=new TStringList();
            TStringList *Values=new TStringList();

            p->loadSampleDescription(Names,Values,OpenDialog1->Files->Strings[0]);
            UpdateSampleDescription(Names,Values);
            p->loadData(tts);
            delete Names;
            delete Values;
        }
        delete tts;
}


// Сохранить всё.
void  N11Click()
{
    if((*ActiveParams()))
    {
    (*ActiveParams())->setSampleDescription(uiTemperature->Text,uiCurrent->Text,
        uiInventoryNumber->Text,uiSampleLength->Text,uiSampleWidth->Text,uiSampleThickness->Text);
    (*ActiveParams())->setParamsType(ResCurveIndex->ItemIndex);
    } 
    params->SaveAllData(SaveDialog1->FileName+"Com"); 
    paramsDirect->SaveAllData(SaveDialog1->FileName+"Dir");
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void concatDependence()
{
    if(!(paramsDirect && paramsReverse))
    {
    ShowMessage("Для объединения нужно измерить две зависимости");
    return;
    }

    ResCurveIndex->ItemIndex=2;    
    MagneticFieldDependence * p=InitParams();
    if(p)

    {
        params->setFilterParamsResistance(eSamplingFRes->Text, eBandWidthFRes->Text,
         eAttenuationFRes->Text, eLengthFilterRes->Text);
        params->setFilterParamsHall(eSamplingFHall->Text, eBandWidthFHall->Text,
         eAttenuationFHall->Text, eLengthFilterHall->Text);
    }

    TSignal B; // создаем буфер для новых зависимостей.
    TSignal Hall;
    TSignal Resistance;

    TSignal B2; // создаем буфер для новых зависимостей.
    TSignal Hall2;
    TSignal Resistance2;

    StatusBar->Panels->Items[1]->Text="Объединение зависимостей.";
    Form1->Update();
    // вбрасываем в обратном порядке зависимости для отрицательного магнитного поля.
    for (TSignal::const_reverse_iterator i = paramsReverse->getB()->rbegin(); i != paramsReverse->getB()->rend(); ++i)
    {
        B.push_back(*i);    
    }

    for (TSignal::const_reverse_iterator i = paramsReverse->getHallEffect()->rbegin(); i != paramsReverse->getHallEffect()->rend(); ++i)
    {
        Hall.push_back(*i);    
    }

    for (TSignal::const_reverse_iterator i = paramsReverse->getMagnetoResistance()->rbegin(); i != paramsReverse->getMagnetoResistance()->rend(); ++i)
    {
        Resistance.push_back(*i);    
    }
    // вбрасываем в прямом порядке зависимости для положительного поля.
    for (TSignal::const_iterator i = paramsDirect->getB()->begin(); i != paramsDirect->getB()->end(); ++i)
    {
        B.push_back(*i);   
        B2.push_back(*i);  
    }

    for (TSignal::const_iterator i = paramsDirect->getHallEffect()->begin(); i != paramsDirect->getHallEffect()->end(); ++i)
    {
        Hall.push_back(*i);  
        Hall2.push_back(*i);   
    }

    for (TSignal::const_iterator i = paramsDirect->getMagnetoResistance()->begin(); i != paramsDirect->getMagnetoResistance()->end(); ++i)
    {
        Resistance.push_back(*i);  
        Resistance2.push_back(*i);   
    }

    TSignal outB;
    TSignal outHall;
    TSignal outResistance;

    size_t minimalLength=paramsReverse->getB()->size()>paramsDirect->getB()->size()?
    paramsDirect->getB()->size() : paramsReverse->getB()->size();

    // что-то похоже что эта функция работает немного не так как надо...
    //thiningSignal(B, Hall, outB, outHall, -2, 0, 2*minimalLength);
    //thiningSignal(B, Resistance, outB, outResistance, -2, 2, 2*minimalLength);
    // и у меня есть некоторое подозрение почему оно так.
    // возможно более правильным будет сначала определить сигнал с наименьшим количество точек (прямой или обратный)
    // потом проредить оба сигнала отдельно и после этого объединить их.
    thiningSignal(B, Hall, outB, outHall, -2, 2, 2*minimalLength);
    thiningSignal(B, Resistance, outB, outResistance, -2, 2, 2*minimalLength);

    params->setDependence(outB.begin(),outB.end(),outHall.begin(),outResistance.begin());

}


long double calcConcentrationFromGp(long double G_p, long double Mu)
{
    long double electronCharge=1.602e-19;
    return G_p/(Mu*electronCharge);
}


void addPeak(TChartSeries *Sender,int ValueIndex)
{
    long double electronCharge=1.602e-19;
    long double G_p=Sender->YValues->Value[ValueIndex];
    long double Mu= Sender->XValues->Value[ValueIndex];
    long double Concentration=G_p/(Mu*electronCharge);

    //Form1->MobSpecResults->Cells[2][Form1->MobSpecResults->Selection.Top]= FloatToStrF( Mu, ffFixed, 5, 5);
    //Form1->MobSpecResults->Cells[1][Form1->MobSpecResults->Selection.Top]= FloatToStrF(Concentration, ffExponent, 5, 2);

    TGridRect tgr=Form1->MobSpecResults->Selection;
    tgr.Top++;
    tgr.Bottom++;
    if(tgr.Top>3)
    {
        tgr.Top=1;
        tgr.Bottom=1;
    }
    Form1->MobSpecResults->Selection =tgr;

    // Сюда нужно добавить обновление значений пиков в объект магнитополневых зависимостей!!!
}


void TForm1::getAndDisplayMultiCarrierFitResults(MagneticFieldDependence * p)
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
    }

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
            tsl->Add(FloatToStr(outValues[i][j]));
        }
    }

    tsl->SaveToFile(SaveDialog1->FileName+eBandWidthFRes->Text+eBandWidthFHall->Text+"MZFitRes.txt");

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
    delete tsl;
}

void  TForm1::btnMultiCarrierFitClick()
{
    MagneticFieldDependence ** par=ActiveParams();
    MagneticFieldDependence * p;
    p=*par;
    long double VesGxx=StrToFloat(LabeledEdit1->Text);
    long double VesGxy=StrToFloat(LabeledEdit2->Text);

    if(p->runMultiCarrierFit(VesGxx,VesGxy))
    {
        getAndDisplayMultiCarrierFitResults(p);
    }
}

std::vector<std::string> getAllFileNamesWithinFolder(std::string folder)
{
    std::vector<std::string> names;
    char search_path[200];
    if (folder.length() <200)
    {
        sprintf(search_path, "%s/*.*", folder.c_str());
        WIN32_FIND_DATA fd; 
        HANDLE hFind = ::FindFirstFile(search_path, &fd); 
        if(hFind != INVALID_HANDLE_VALUE) { 
            do { 
                // read all (real) files in current folder
                // , delete '!' read other 2 default folder . and ..
                if(! (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ) {
                    names.push_back(fd.cFileName);
                }
            }while(::FindNextFile(hFind, &fd)); 
            ::FindClose(hFind); 
        } 
        return names;
    }
    return names;
}

void filterFileNames(std::vector<std::string> & files, std::string rem)
{
    for (int i = 0; i < files.size(); ++i)
    {
        size_t n = files[i].find(rem);
        if(n!= std::string::npos)
            files[i]="";
    }
}
*/

//void  TForm1::Button3Click()
//{
/*
Функция должна обеспечивать:
1. Открытие файла.
2. Фильтрацию, расчет тензоров, спектра и подгонки.
3. Сохранение всех результатов.
*/
/*
// открыть папку и обработать все файлы в ней
if (OpenDialog1->Execute())
{

    std::string folder = ExtractFilePath(OpenDialog1->Files->Strings[0]).c_str();

    std::vector<std::string> files=getAllFileNamesWithinFolder(folder);

    // Нужно удалить из списка файлы, содержащие в имени "Description"

    filterFileNames(files, "Description");

    // предусмотреть объединение сигналов
    // пока предполагаем что сигнал уже объединен - ведь на модели испытывать будем.

    int temp=files.size();

    for (int j = 0; j < files.size(); ++j) // По всем именам файлов
    {
        if (files[j]!="")
        {
        openFileWithSignal(files[j].c_str());


        char folderPath[500];
        if (files[j].length() <500)
        {
            sprintf(folderPath, "%s%s1", folder.c_str(), files[j].c_str());

            if(CreateDirectory(folderPath,NULL)) // Создаём отдельную папку для сохранения данных по каждому файлу.
            {
                std::string name=folder+files[j]+"1\\1"; // Дабы сохранять в отдельную папку.
                long double step=0.1;
                for (long double BandWidthFHall=0.01, AttenuationFHall=0.1; BandWidthFHall < 10; BandWidthFHall+=step, AttenuationFHall+=step)
                {
                    if(step==0.1 && BandWidthFHall>1.0)
                    {
                        step=1.0;
                    }
                    eBandWidthFHall->Text=FloatToStr(BandWidthFHall);
                    eAttenuationFHall->Text=FloatToStr(AttenuationFHall);

                    AnsiString t=name.c_str()+FloatToStr(BandWidthFHall);
                    //eLengthFilterRes->Text="300";
                    long double step=0.1;
                    for (long double BandWidthFRes=0.01, AttenuationFRes=0.1; BandWidthFRes < 10; BandWidthFRes+=step, AttenuationFRes+=step) // Это будет страшный цикл по разным параметрам фильтрации
                    { // Но для начала буду менять только частоту среза
                        if(step==0.1 && BandWidthFRes>1.0)
                        {
                            step=1.0;
                        }
                        AnsiString t2=t+FloatToStr(BandWidthFRes);
                        SaveDialog1->FileName=t2;

                        eBandWidthFRes->Text=FloatToStr(BandWidthFRes);
                        eAttenuationFRes->Text=FloatToStr(AttenuationFRes);

                        Application->ProcessMessages();

                        if(CheckBox2->Checked==false)
                        {
                        N26->Click(); // Выполняет Фильтрацию, расчет тензоров, спектра и подгонки.
                        }
                        else
                        {
                        N25->Click();
                        }

                        if(params)
                        {
                            params->SaveAllData(SaveDialog1->FileName+"Com");
                        }
                        //N11->Click(); // Сохранить всё
                    }

                }

            }
            else
                ;//ShowMessage(IntToStr(GetLastError()));
        }
        }
    }
}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void  TForm1::Button5Click()
{
    const TSignal * bv=(*ActiveParams())->getB();
    const TSignal * hallv=(*ActiveParams())->getHallEffect();
    const TSignal *mv=(*ActiveParams())->getMagnetoResistance();
    std::vector<long double> b=*bv;
    std::vector<long double> hall=*hallv;
    std::vector<long double> m=*mv;
    long double h = fabs(b[1]-b[0]);
    TSignal d1hall, d1m;
    TSignal d2hall, d2m;

    d1hall=calculateFirstDerivative(hall,h);
    d2m=calculateSecondDerivative(m,h);
    d1m=calculateFirstDerivative(m,h);
    d2hall=calculateSecondDerivative(hall,h);

    TStringList * tsl=new TStringList();

    for (int i = 0; i < d2hall.size(); ++i)
    {
        tsl->Add(FloatToStr(d1hall[i])+"\t"+FloatToStr(d2hall[i])+"\t"+FloatToStr(d1m[i])+"\t"+FloatToStr(d2m[i]));
    }

    tsl->SaveToFile("derivativeOfSignals.txt");   

    delete tsl;
}


void  TForm1::Button6Click()
{
   MagneticFieldDependence ** par=ActiveParams();
    MagneticFieldDependence * p;
    if (*par==NULL)
    {
        ShowMessage("Вероятно выбран не тот график.");
        return;
    }
    else
    {
        long double VesGxx=StrToFloat(LabeledEdit1->Text);
        long double VesGxy=StrToFloat(LabeledEdit2->Text);
        p=*par;
        p->runSmartCalcutation();
        p->runSmartMultiCarrierFit(VesGxx,VesGxy);    
        getAndDisplayMultiCarrierFitResults(p);
    }
}*/
//---------------------------------------------------------------------------
