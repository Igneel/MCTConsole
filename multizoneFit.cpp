#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <cmath>  /* round*/

#include "multizoneFit.h"




void MultiZoneFit::memoryAlloc()
{
    MagField_spektr.resize(MaxPoints+1);
    Gxx.resize(MaxPoints+1);
    Gxy.resize(MaxPoints+1);
    GxxExp.resize(MaxPoints+1);
    GxyExp.resize(MaxPoints+1);

    Data.resize(SIZE);
    Data0.resize(SIZE);
    Min_Value.resize(SIZE);
    Max_Value.resize(SIZE);
    D_Step.resize(SIZE);
    D_StepOld.resize(SIZE);
    DataOld.resize(SIZE);
    DataBef.resize(SIZE);

    d1.resize(MaxParameters+1);
    for (auto i = 0u; i < MaxParameters+1; ++i) {
        d1[i].resize(MaxRepeat+1);
    }
}

long double MultiZoneFit::func_hall8(const Image & Data, const Data_spektr & Magfield_spektr, Data_spektr &  Gxx,
  Data_spektr &  Gxy, const Data_spektr &  GxxExp, const Data_spektr &  GxyExp,
 const long double NPoint)
{
  unsigned int myi;
  long double g,g1,cond1,cond2,cond3;
for (myi=0; myi<NPoint; ++myi)
  {
  cond1=Data[4]*Data[1]/(1+powl(Data[1]*Magfield_spektr[myi],2));
  cond2=Data[5]*Data[2]/(1+powl(Data[2]*MagField_spektr[myi],2));
  cond3=Data[6]*Data[3]/(1+powl(Data[3]*MagField_spektr[myi],2));
  Gxx[myi]=ElectronCharge*(cond1+cond2+cond3);
  if (myi==0 )
   Gxy[0]=0;
  else
   Gxy[myi]=(cond1*Data[1]+cond2*Data[2]+cond3*Data[3])
        *MagField_spektr[myi]*ElectronCharge;
  }
   g=0;
   g1=0;
  for (myi=0; myi<NPoint; ++myi)
  {
     //g+=powl(Gxx[myi]-GxxExp[myi],2)/(GxxExp[myi]+Gxx[myi]);
    g+=powl(Gxx[myi]-GxxExp[myi],4)/(GxxExp[myi]+Gxx[myi]);
   if (myi==0)
    {
      g1=0;
   }
       else
    {
     g1+=powl(Gxy[myi]-GxyExp[myi],2)/(fabsl(GxyExp[myi])+fabsl(Gxy[myi]));
    }
  }

  return 100.0L *(sqrtl(g*VesGxx)+sqrtl(g1*VesGxy))/(NPoint);
} /* of func */



  void MultiZoneFit::BegRand(const int N_Make)
{
    Image Params_Best;
    Params_Best.resize(SIZE);
    Data.resize(SIZE);
    Fbefore=1.e8L;
    for ( int ii=1 ; ii<= N_Make ;++ii)  // {до N_Make} { Number of steps }
    {
        for ( int j=1 ; j<= 100 ;++j)
        {   // присваеваем каждому параметру случайное значение
            // в заданных пределах
            for (auto i=1u; i<= N_Data ;++i)
              Data[i]=Min_Value[i]+(static_cast<long double >(rand())/RAND_MAX)*(Max_Value[i]-Min_Value[i]);
            // вычисляем значение целевой функции в данной точке
            Fnew=func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint);
             // случайно сгенерированных параметров
            if (Fnew<Fbefore)
            {      // запоминаем лучшие значения параметров
                Params_Best=Data;
              Fbefore=Fnew; // запоминаем значение соотв. им функции
            }

        }
        // после 100 случайных забросов записываем лучшие значения параметров
        // в Data и рисуем по ним график
        Fnew=Fbefore;
        Data=Params_Best;
        // GraphikH;
    }


}   // в итоге мы получаем лучшую из 100*N_Make точку,
       // используя её в качестве базовой.

 void MultiZoneFit::InitVar()
 {
  const int magicSix = 6;
  //int fi;
   PP=1;
   DataOld=Data;
   DataBef=Data;
   Fbefore=Fnew;
   for (auto fi=1u; fi<= magicSix ;++fi)          // выбор шага приращения
    if ( Data0[fi]!=0.0L ) D_Step[fi]=D_StepOld[fi]*Data[fi]/Data0[fi];
    else D_Step[fi]=0;
 }

void MultiZoneFit::CheckLimits()   // проверяет, не принял ли параметр значение,
// лежащее за заданными пределами
  {

   for (auto fi=1u; fi<= N_Data ;++fi)
   {       // если true -  параметрам присваиваются граничные значения
     if ( Data[fi]<Min_Value[fi] ) Data[fi]=Min_Value[fi];
     if ( Data[fi]>Max_Value[fi] ) Data[fi]=Max_Value[fi];
    }
  }

void MultiZoneFit::Research()
{
    unsigned int k=1;
    Image   Params_Old;
    int as0[7]={0}; // 1..6

 Fnew=func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint);
 Params_Old=Data;
 do
 {
     ++as0[k];
     Fold=Fnew;
     Data[k]=Params_Old[k]+as0[k]*D_Step[k];
     CheckLimits();
     Fnew=func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint);
     if ( Fnew>Fold )   // удачен ли сдвиг?
       {    // неудачен - шаг назад
         Data[k]=Params_Old[k]-as0[k]*D_Step[k];
         CheckLimits();
         Fnew=func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint);
       }
     if ( Fnew>Fold )   // удачен ли сдвиг?
       {    // если снова неудачен - оставляем старое значение
         Fnew=Fold;
         Data[k]=Params_Old[k];
        }
     else  // удачен - запоминаем и переходим к следующему параметру
       Params_Old[k]=Data[k];
     ++k;
}
 while (!(k>N_Data));        // делает шаг по всем параметрам, запоминает
                        // лучшие за шаг
 Data=Params_Old; // переписываем лучшие в Data
} //{of void}

void MultiZoneFit::Hook()
{
    for(auto i=1u; i<= N_Data ;++i)
    {
        D_StepOld[i]=Data[i]*0.1L;
    }
    Data0=Data;
    Fnew=func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint);
    FlagEnd=false;
    FlagDipl=false;
   do
   {
     if (!FlagDipl) InitVar();
     Research(); // исследование
      if ( fabsl(Fnew)<eps1)
        {
          FlagEnd=true;
        }
      // улучшено ли значение функции && достаточно ли мало изменение функции
      while ( (Fnew<Fbefore) &&(fabsl(Fbefore-Fnew)>epsilon))
       {
          // запоминаем параметры старой базовой точки
          DataOld=Data;
         Fbefore=Fnew;             //  {function in the new base point }
         for (auto i=1u; i<= N_Data ;++i)
         {    // новая базовая точка
            Data[i]=Data[i]+a[i]*(Data[i]-DataBef[i]); // экстраполяция
            if ( Data[i]<Min_Value[i] ) Data[i]=Min_Value[i];
            if ( Data[i]>Max_Value[i] ) Data[i]=Max_Value[i];
         }
         DataBef=Data;
         for (auto i=1u; i<= N_Data ;++i)  // увеличиваем шаг приращения
          if ( Data0[i]!=0.0L )
          {
              D_Step[i]=D_StepOld[i]*Data[i]/Data0[i];
          }
          else
          {
              D_Step[i]=0;
          }
         FlagDipl=true;
         Research(); // исследование
         if ( Fnew<Fbefore ) // если значение функции уточнено
         {   // рисуем график
          // GraphikH;
         }
         if ( (fabsl((Fnew-Fbefore))>epsilon) )
            {
              NP=0;
            }
         if ( (fabsl((Fnew-Fbefore))<=epsilon) )
         {
          ++NP;
         break;
         }
         if ( Fnew<eps1 )  // если значение функции меньше заданной величины
         {   // конец поиска.
           FlagEnd=true;
         break;
         }
            
       }
      // если изменения функции незначительные
      if ( FlagDipl && (fabsl((Fnew-Fbefore))<=epsilon) )
         //  FlagEnd=True
           FlagDipl=false;  // не заверншаем, а увеличиваем шаг приращения
      else  // если новое значение функции стало больше предыдущего
        if ( FlagDipl && (Fnew>Fbefore) )
          {  // возвращаемся к предыдущему значению функции
            Data=DataOld;
           Fnew=Fbefore;
           FlagDipl=false;
          }
        else
          {  // данный вариант возможен лишь при  FlagDipl=false.

          for (auto i=1u; i<= N_Data ;++i)  // уменьшаем шаг
          {
            D_StepOld[i]=D_StepOld[i]/2.0L;
          }
          Data=DataOld;// к предыдущей точке
          Fnew=Fbefore;
          PP=0;
          ++NP;
          }
          if ( NP>=Nstep )
          {
            FlagEnd=true;
          }

} while (!FlagEnd);
}

////////////////////////////////////////////////////////////////////////////
///////////////// НАЧАЛО "ХОЛЛ. МНОГОЗОННАЯ ПОДГОНКА" //////////////////////
////////////////////////////////////////////////////////////////////////////

void MultiZoneFit::Optimiz_hall8()
{
//  srand(time(NULL));
  N_Data=6;
  BegRand(1000);   // в данном случае c 1000 гораздо лучше ищет
  Hook();
}

long long justRound(long double x)
{
  if ((static_cast<int>(x*10))%10<5)
  {
    return floorl(x);
  }
  else
  {
    return ceill(x);
  }
}

void MultiZoneFit::btnFeatMultiClick()
{
  long double funcMin=10e8L;
  //Form1->ProgressBar1->Max=repeatQuantity;

  for (auto l=1u; l<= repeatQuantity ; ++l)
  {
    //Form1->ProgressBar1->Position=l;


    Optimiz_hall8();
    for (auto i=1u; i<= 6 ; ++i)
    d1[i][l]=Data[i];
    d1[7][l]=func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint);

    if(  d1[7][l]<funcMin )
    {
      funcMin=d1[7][l];
      //VesGxx= justRound(funcMin); // Вот тут меняется параметр веса!
    }
  }
 //Gistogram(d1,SerArr,7,100);
}
/*
bool Interval(TStringGrid & StringGrid3)   // рассчитывает подвижность и концентрацию

{

int i;             // 3-х типов носителей в интервале (1%; 150%)
bool Stg3_pusto;
long double max,min;
  min=-0.5;
  max=0.5;
  Stg3_pusto=false;
  for ( i=1 ; <= 3 ; ++i)
      if( (StringGrid3.Cells[1,i]="") or (StringGrid3.Cells[2,i]="") )
         Stg3_pusto=true;
      if( Stg3_pusto )
      {
        Interval=true;
        Exit();
      }
      else
      {
        Max_Value[1]=RoundTo(StrToFloat(StringGrid3.Cells[2,3])+
        StrToFloat(StringGrid3.Cells[2,3])*min,-6);
        Max_Value[2]=RoundTo(StrToFloat(StringGrid3.Cells[2,2])+
        StrToFloat(StringGrid3.Cells[2,2])*max,-3);
        Max_Value[3]=RoundTo(StrToFloat(StringGrid3.Cells[2,1])+
        StrToFloat(StringGrid3.Cells[2,1])*max,-3);
        Max_Value[4]=RoundTo(StrToFloat(StringGrid3.Cells[1,3])+
        StrToFloat(StringGrid3.Cells[1,3])*min,14);
        Max_Value[5]=RoundTo(StrToFloat(StringGrid3.Cells[1,2])+
        StrToFloat(StringGrid3.Cells[1,2])*max,15);
        Max_Value[6]=RoundTo(StrToFloat(StringGrid3.Cells[1,1])+
        StrToFloat(StringGrid3.Cells[1,1])*max,19);
        Min_Value[1]=RoundTo(StrToFloat(StringGrid3.Cells[2,3])+
        StrToFloat(StringGrid3.Cells[2,3])*max,-3);
        Min_Value[2]=RoundTo(StrToFloat(StringGrid3.Cells[2,2])+
        StrToFloat(StringGrid3.Cells[2,2])*min,-6);
        Min_Value[3]=RoundTo(StrToFloat(StringGrid3.Cells[2,1])+
        StrToFloat(StringGrid3.Cells[2,1])*min,-6);
        Min_Value[4]=RoundTo(StrToFloat(StringGrid3.Cells[1,3])+
        StrToFloat(StringGrid3.Cells[1,3])*max,15);
        Min_Value[5]=RoundTo(StrToFloat(StringGrid3.Cells[1,2])+
        StrToFloat(StringGrid3.Cells[1,2])*min,14);
        Min_Value[6]=RoundTo(StrToFloat(StringGrid3.Cells[1,1])+
        StrToFloat(StringGrid3.Cells[1,1])*min,18);
        Interval=false;
      }
      
}*/

// расчет среднего и СКО
// вот сюда указатель на массив надо давать.
void MultiZoneFit::Statistic(const DataValue & mass, std::vector<long double> & ResulitsForStatistic,
                             unsigned int m, unsigned int n)
{
	long double Xsr, S;
    for (auto ll=1u; ll<= m ; ++ll)
    {
        Xsr=0;
        S=0;
        for (auto  l=1u; l<= n ; ++l)
        Xsr=Xsr+mass[ll][l];
        Xsr=Xsr/static_cast<long double>(n);
        for (auto  l=1u; l<= n ; ++l)
        S+=powl(mass[ll][l]-Xsr,2);
        S=sqrtl(S/(n-1));//



        ResulitsForStatistic[ll]= Xsr;
        ResulitsForStatistic[ll+m]= S;
        ResulitsForStatistic[ll+2*m]= (S/fabsl(Xsr))*100.0L;
    }

  
  for (auto ll=1u; ll<= m ; ++ll)
    {
      long double Xmiddle=0, s=0;
      long double count=0;
      for (auto  l=1u; l<= n ; ++l)
        if(fabsl(mass[ll][l])<=fabsl(ResulitsForStatistic[ll]+ResulitsForStatistic[ll+m])
           && fabsl(mass[ll][l])>=fabsl(ResulitsForStatistic[ll]-ResulitsForStatistic[ll+m]))
        {
          Xmiddle+=mass[ll][l];
          count++;
        }
      
      Xmiddle=Xmiddle/static_cast<long double>(count);
        
      for (auto  l=1u; l<= n ; ++l)
        if(fabsl(mass[ll][l])<=fabsl(ResulitsForStatistic[ll]+ResulitsForStatistic[ll+m])
           && fabsl(mass[ll][l])>=fabsl(ResulitsForStatistic[ll]-ResulitsForStatistic[ll+m]))
        {
          s+=powl(mass[ll][l]-Xmiddle,2);
        }
      
      s=sqrtl(s/(count-1));

      ResulitsForStatistic[ll]= Xmiddle;
      ResulitsForStatistic[ll+m]= s;
      ResulitsForStatistic[ll+2*m]= (s/fabsl(Xmiddle))*100.0L;
    }
}

////////////////////////////////////////////////////////////////////////////
///////////////// КОНЕЦ "ХОЛЛ. МНОГОЗОННАЯ ПОДГОНКА" //////////////////////
////////////////////////////////////////////////////////////////////////////

/*
Итак, картина следующая
На вход подаем зависимости:

Минимальные и максимальные значения. Т.е. коридор.
Дальше данные тензоров.

read(Data_File,Max_Value[1],Max_Value[2],Max_Value[3],Max_Value[4],
         Max_Value[5],Max_Value[6]);
         read(Data_File,Min_Value[1],Min_Value[2],Min_Value[3],Min_Value[4],
         Min_Value[5],Min_Value[6]);
         for ( i=0 ; <= NPoint-1 ; ++i)
         readln(Data_File,Magfield_spektr[i],GxxExp[i],GxyExp[i]);
         Memo6.Lines.LoadFromFile(OpenDialog3.FileName);

На выходе читаем данные из:

Минимальные значения функции и значения параметров соотв. минимальному значению функции.

Графики тензоров (экспериментальные и теоретические).

mmoFeatFunctionValue.Lines.Add(FloatToStr(func_hall8));
   mmoFeatMuElectronValue.Lines.Add(FloatToStr(DATA[1]));
   MemoMnog3.Lines.Add(FloatToStr(DATA[2]));
   MemoMnog4.Lines.Add(FloatToStr(DATA[3]));
   mmoConcentrationElectrons.Lines.Add(FloatToStr(DATA[4]));
   MemoMnog6.Lines.Add(FloatToStr(DATA[5]));
   MemoMnog7.Lines.Add(FloatToStr(DATA[6]));
   for ( i=0 ; <= NPoint-1  ; ++i)
   {
     Series7.AddXY(MagField_spektr[i], Gxx[i]);
     Series8.AddXY(MagField_spektr[i], GxxExp[i]);
     Series9.AddXY(MagField_spektr[i], Gxy[i]);
     Series10.AddXY(MagField_spektr[i], GxyExp[i]);

     Надо найти способ выудить отсюда ещё средние значения.
     СКО и СКО в %

}*/
/*
long double MultiZoneFit::someTestFunc(long double VesGxx)
{
    VesGxx=VesGxx+1;
    return VesGxx;
}*/

// Возможно стоит сделать инициализатор, дабы не отжирать тонну памяти?

// Возвращаем по традиции количество точек?
// Веса, нижняя и верхняя границы, Спектр, две компоненты тензора
int MultiZoneFit::RunMultizoneFeat (const std::vector<long double> LowBound,const  std::vector<long double> UpBound,
  const InDataSpectr MagSpectr,const  InDataSpectr GxxIn,const  InDataSpectr GxyIn,
  MyData_spektr & outGxx, MyData_spektr & outGxy,
  TStatistic & outValues,
  unsigned int numberOfCarrierTypes)
{
memoryAlloc();
    std::vector<long double> ResulitsForStatistic;
    unsigned long SizeForStatistic;
  ResulitsForStatistic.resize(3*(2*numberOfCarrierTypes+1)+1);
  // Это надо получить из управляющей программы:
  NPoint=MagSpectr.size();

  // Эти переменные объявлены в optim1 и их тоже надо получить.
    Max_Value.resize(numberOfCarrierTypes*2+1);
    Min_Value.resize(numberOfCarrierTypes*2+1);
  for (auto i=1u; i<= numberOfCarrierTypes*2 ; ++i)
  {
    Max_Value[i]=UpBound[i-1];
    Min_Value[i]=LowBound[i-1];
  }
    MagField_spektr.resize(NPoint);
    GxxExp.resize(NPoint);
    GxyExp.resize(NPoint);
  // ну и разумеется получить сами точки.
  for (auto i=0u; i<= NPoint-1 ; ++i)
  {
    MagField_spektr[i]=MagSpectr[i];
    GxxExp[i]=GxxIn[i];
    GxyExp[i]=GxyIn[i];
  }  

  // После того как все данные получены начинаем выполнение.
  btnFeatMultiClick();

  // программа отработала
  // Так, где лежат результаты?

  for (auto i=0u; i<= NPoint-1 ; ++i)
  {
    outGxx.push_back(Gxx[i]);
    outGxy.push_back(Gxy[i]);
  }
  // эта функция должна посчитать нам статистику.
  Statistic(d1, ResulitsForStatistic, MaxParameters-1,MaxRepeat);
  SizeForStatistic=MaxParameters-1;
  outValues.resize(MaxParameters);
  for (auto i=0u; i<= MaxParameters-1 ; ++i)
  {

    outValues[i].push_back(Data[i+1]);
    outValues[i].push_back(ResulitsForStatistic[i+1]);
    outValues[i].push_back(ResulitsForStatistic[i+SizeForStatistic+1]);
    outValues[i].push_back(ResulitsForStatistic[i+2*SizeForStatistic+1]);
  }

  outValues[7].push_back(func_hall8(Data,MagField_spektr,Gxx,Gxy,GxxExp,GxyExp,NPoint)); // минимальное значение функции

  return 0;

}  
