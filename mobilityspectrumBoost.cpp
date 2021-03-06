#include "mobilityspectrumBoost.h"
#include <float.h>
#include "commonFunctions.h"

using namespace std;


//-------------------------------------------------------------------
std::vector<myFloat> calculateSecondDerivative(std::vector<myFloat> & y, myFloat h)
{
  TSignal d2Y(y.size()-2);

      for(size_t i =0;i<y.size()-2;i++)
      {
        //d2Y[i]=1.0/(fabs(dY[i+1]-dY[i]))*(2.0*dY[i+1]-dY[i+2]/2.0-3.0/2.0*dY[i]);
        // формула f(x-h)-2f(x)+f(x+h)/h^2
        d2Y[i]=(y[i]-2*y[i+1]+y[i+2])/h/h;
        //d2Y[i]=(dY[i+2]-dY[i])/2.0/h;
      }
  return d2Y;
}

//-------------------------------------------------------------------
std::vector<myFloat> calculateFirstDerivative(std::vector<myFloat> & y, myFloat h)
{
  TSignal dY(y.size()-2);
  // Посчитаем производную методом конечных разностей
  // формула df/dx=(f(x+h)-f(x-h))/2/h;
  for(size_t i =0;i<y.size()-2;i++)
      {
        //dY[i]=1.0/(fabs(y[i+1]-y[i]))*(2.0*y[i+1]-y[i+2]/2.0-3.0/2.0*y[i]);
        dY[i]=(y[i+2]-y[i])/2.0/h;
      }
  return dY;
}



void mobilitySpectrum::SetLength (vector < vector <myFloat> > & v, const size_t size1, const size_t size2)
{
   v.resize(size1);
   for (size_t i = 0; i < size1; ++i) {
       v[i].resize(size2);
   }
}


////////////////////////////////////////////////////////////////////////////
/////////////////////// НАЧАЛО "ХОЛЛ. ПОДВИЖНОСТЬ" /////////////////////////
////////////////////////////////////////////////////////////////////////////

void  mobilitySpectrum::InitArray()
{
   IntMagField.resize(SizeData);
   IntGxx.resize(SizeData);
   IntGxy.resize(SizeData);
}

void  mobilitySpectrum::InitArray2()
{
   Spectr_e.resize(SizeData);
   Spectr_p.resize(SizeData);
   Mobility.resize(SizeData);
}

void  mobilitySpectrum::GetCoef(const Data_spektr &A, const Data_spektr &X, const myFloat  b,
    myFloat & p0, myFloat & p1, myFloat & p2)

{
    int i,j;
    myFloat r,r1,r2,q,s;
    p0 = 0;
    p1 = 0;
    p2 = 0;
    for (i = 0;i<=NumberOfPoints;++i)
    {
        r = 1;r1 = 0;r2 = 0;
        for (j = 0;j<=NumberOfPoints;++j)
        {
            if (i!=j)
            {
                q = X[i]-X[j];
                s = b-X[j];
                r2=(r2*s+2*r1)/q;
                r1=(r1*s+r)/q;
                r = r*s/q;
            }
        }
        p0 = p0+r*A[i];
        p1 = p1+r1*A[i];
        p2 = p2+r2*A[i];
    }
}

// Эта функция также нужна разве что мне для справки, т.к. вывод на график будет
// в основной программе.
// Или нет, не совсем понятно - возможно программа использует построенные графики
// для дальнейших расчетов, что тоже не совсем правильно.
// строит экспериментальные точки на графиках компонент тензора

void  mobilitySpectrum::AddExpPoints(TLineSeries& ExpXX, TLineSeries& ExpXY)
{
    int i;

    ExpXX.clear();
    //ExpXX.Pointer.HorizSize = 2;
    //ExpXX.Pointer.VertSize = 2;
    for (i = 0 ; i<= NumberOfPoints; ++i )
        if ( MagField_spektr[i]==0 ) ExpXX.push_back(make_pair(0.001,GxxExp[i])); // ExpXX.AddXY(0.001,GxxExp[i],"",clTeeColor);
        else ExpXX.push_back(make_pair(MagField_spektr[i],GxxExp[i])); // ExpXX.AddXY(MagField_spektr[i],GxxExp[i],"",clTeeColor);


    ExpXY.clear();
    //ExpXY.Pointer.HorizSize = 2;
    //ExpXY.Pointer.VertSize = 2;
    for (i = 0 ; i<= NumberOfPoints; ++i )
     if ( MagField_spektr[i]==0 ) ExpXY.push_back(make_pair(0.001,GxyExp[i])); // ExpXY.AddXY(0.001,GxyExp[i],"",clTeeColor);
       else ExpXY.push_back(make_pair(MagField_spektr[i],GxyExp[i])); // ExpXY.AddXY(MagField_spektr[i],GxyExp[i],"",clTeeColor);

}

void  mobilitySpectrum::GetLnLimits(myInt& Xmin, myInt& Xmax )
{
    myFloat f;
   if ( MagField_spektr[0]>0 )
    {
     f = log10(MagField_spektr[0]);
     Xmin = static_cast<int>(f);
    }
     else Xmin=-3;

    f = log10(MagField_spektr[NumberOfPoints]);
    Xmax = static_cast<int>(f);
    myFloat temp2;
    myFloat temp=modf(f,&temp2);
    if ( temp>0.001 )
        Xmax = Xmax+1;

}

void  mobilitySpectrum::CS(Data_spektr& X, Data_spektr& F, Data_spektr& C, const myFloat p1, const myFloat pn)
{
    int i,j,m;
    Data_spektr K;
    myFloat A,B,R;

    K.resize(MaxPoints);

   K[1]=0;
   C[1]=p1;
   A = X[1]-X[0];
   B = X[2]-X[1];
   K[2]=(3*((F[2]-F[1])/(X[2]-X[1])-(F[1]-F[0])/(X[1]-X[0]))-
       (X[1]-X[0])*C[1])/2/(A+B);
   C[2]=B/2/(A+B);
   for (i= 3; i<= NumberOfPoints; ++i )
    {
     j = i-1;m = j-1;
     A = X[i]-X[j];
     B = X[j]-X[m];
     R = 2*(A+B)-B*C[j];
     C[i]=A/R;
     K[i]=(3*((F[i]-F[j])/A-(F[j]-F[m])/B)-B*K[j])/R;
    }
   C[NumberOfPoints]=K[NumberOfPoints]-C[NumberOfPoints]*pn;

   for (i= NumberOfPoints-1 ; i>= 2; --i ) C[i]=K[i]-C[i]*C[i+1];
}

myFloat mobilitySpectrum::Sp(Data_spektr & X, Data_spektr &F, Data_spektr & C, const myFloat x1)

{
    int i,j;
    myFloat A,B,D,Q,R,P;

i = 1;
while (x1>X[i])
    ++i;

 j = i-1;
 A = F[j];
 B = X[j];
 Q = X[i]-B;
 R = x1-B;
 P = C[i];
 D = C[i+1];
 B=(F[i]-A)/Q-(D+2*P)*Q/3;
 D=(D-P)/3/Q;
 return A+R*(B+R*(P+D*R));
}

// Эта функция нигде никем не вызывается О_о, что крайне странно - надо искать в чем дело.
void  mobilitySpectrum::MakeInterpolate(TLineSeries &Gxx,TLineSeries& Gxy,
    TLineSeries& ExpXX, TLineSeries& ExpXY)


{
    Data_spektr temp_l,temp_t,AGxx,AGxy,AField;
    myFloat sf,lm,p,p1;
    int i,j;
    int k;
    myInt Lmin,Lmax;
temp_l.resize(MaxPoints);
temp_t.resize(MaxPoints);
AGxx.resize(MaxPoints);
AGxy.resize(MaxPoints);
AField.resize(MaxPoints);

  //Формируем новые матрицы для расчета производных в точке В=0
  AField[0]=-MagField_spektr[1];
  AGxx[0]=GxxExp[1];
  AGxy[0]=-GxyExp[1];
  for (i= 1; i<= NumberOfPoints;++i)
   {
    AField[i]=MagField_spektr[i-1];
    AGxx[i]=GxxExp[i-1];
    AGxy[i]=GxyExp[i-1];
   }

   //Вычисляем производные в точке В=0
   --NumberOfPoints;
   GetCoef(AGxx,AField,AField[1],p,p1,A1);
   p1 = 0;
   GetCoef(AGxy,AField,AField[1],p,p1,B1);
   B1 = 0;

   //Вычисляем производные в точке В=Вмах
   GetCoef(GxxExp,MagField_spektr,MagField_spektr[NumberOfPoints],p,p1,An);
   GetCoef(GxyExp,MagField_spektr,MagField_spektr[NumberOfPoints],p,p1,Bn);
   An = 0;
   CS(MagField_spektr,GxxExp,temp_l,A1,An);
   CS(MagField_spektr,GxyExp,temp_t,B1,Bn);

   AddExpPoints(ExpXX,ExpXY);

   Gxx.clear();
   Gxy.clear();

   GetLnLimits(Lmin,Lmax);
   SizeData=static_cast<long int>((Lmax-Lmin+1)*PointPerInt+1);
   // SizeData:=(Lmax-Lmin+1)*sizeof(ImageDat);
   // Надо будет проверить и посмотреть кто прав.
   InitArray();
   k = 0;
   for (i = 0;i<= (Lmax-Lmin); ++i )
    {
     lm = exp((static_cast<int>(Lmin)+i)*log(10));
     sf = lm;
     for (j = 1 ; j<= PointPerInt-1 ; ++j )
      {
       IntMagField[k]=sf;
       IntGxx[k]=Sp(MagField_spektr,GxxExp,temp_l,sf);
       IntGxy[k]=Sp(MagField_spektr,GxyExp,temp_t,sf);
       Gxx.push_back(make_pair(IntMagField[k],IntGxx[k]));
       Gxy.push_back(make_pair(IntMagField[k],IntGxy[k]));
       //Gxx.AddXY(IntMagField[k],IntGxx[k],"",clTeeColor);
       //Gxy.AddXY(IntMagField[k],IntGxy[k],"",clTeeColor);
       sf = lm*exp(static_cast<myFloat>(j)/static_cast<myFloat>(PointPerInt)*log(10));
       if ( sf>MagField_spektr[NumberOfPoints] ) break;
       ++k;
      }
     if ( sf>MagField_spektr[NumberOfPoints] ) break;
    }

}


void  mobilitySpectrum::MakeMNK( bool a,TLineSeries& Gxx, TLineSeries& Gxy, TLineSeries& ExpXX,TLineSeries& ExpXY)
{
    mat tmp_m;
        Data_spektr coef_t,coef_l;
        int Kind,k;
        myInt Lmin,Lmax;
        myFloat lm,sf;
  SetLength(tmp_m,MaxPoints,MaxPoints);

  coef_t.resize(MaxPoints);
  coef_l.resize(MaxPoints);


   Power_spektr = 3;  // какая-то степень, или мощность, сильно похоже на кол-во типов носителей
   Kind = 2;  // разновидности, тоже пока не ясно зачем и что
   // эти вызовы заполняют матрицу tmp_m
   // которая кстати уже не пустая О_о

   gram(NumberOfPoints,static_cast<int>(Power_spektr),Kind,MagField_spektr,GxxExp,tmp_m);
   // Гаусс - судя по всему решает матрицу методом Гаусса, сохраняет всё в coef-l

   gauss(static_cast<int>(Power_spektr),tmp_m,coef_l);

   gram(NumberOfPoints,static_cast<int>(Power_spektr),Kind,MagField_spektr,GxyExp,tmp_m);

   gauss(static_cast<int>(Power_spektr),tmp_m,coef_t);

   Gxx.clear(); // чистим графики компонент
   Gxy.clear();


   if ( a )
     {
      AddExpPoints(ExpXX,ExpXY); // добавляет точки на график, экспериментальные
      GetLnLimits(Lmin,Lmax); // получает пределы, логарифмические
      SizeData=static_cast<long>((Lmax-Lmin+1)*PointPerInt+1); // считаем размер данных
      InitArray();  // выделяем его---------------------------------------------------------------

      k = 0;
      for (int i= 0 ; i<= (Lmax-Lmin); ++i )
       {
       lm = exp((static_cast<int>(Lmin)+i)*logl(10));
       sf = lm;
       for (int j = 1 ; j<= PointPerInt-1; ++j )
        {
        IntMagField[k]=sf;

        // а тут происходит самое страшное - считаются два основных графика
        fi(NumberOfPoints,static_cast<int>(Power_spektr),Kind,coef_l,MagField_spektr,sf,IntGxx[k]);
        fi(NumberOfPoints,static_cast<int>(Power_spektr),Kind,coef_t,MagField_spektr,sf,IntGxy[k]);
        // и судя по всему ось х - действительно подвижность
        // а ось у - это величины компонент тензора проводимости
        // но они как-то модифицированы
        // видимо в предыдущей функции, надо уточнить
        Gxx.push_back(make_pair(IntMagField[k],IntGxx[k]));
        Gxx.push_back(make_pair(IntMagField[k],IntGxy[k]));
        //Gxx->AddXY(IntMagField[k],IntGxx[k],"",clTeeColor);
        //Gxy->AddXY(IntMagField[k],IntGxy[k],"",clTeeColor);
        sf = lm*exp(static_cast<myFloat>(j)/static_cast<myFloat>(PointPerInt)*logl(10));
        if ( sf>MagField_spektr[NumberOfPoints] )
            break;
        ++k;
       }
      if ( sf>MagField_spektr[NumberOfPoints] )
          break;
     }

    } else // тут это, такая картина - ежели идет эта ветка то sf никак не определен...
    { // но по факту эта ветка никогда не вызывается.
       sf=0;
     for (int i= 0 ;i <= NumberOfPoints; ++i )
      {
        fi(NumberOfPoints,static_cast<int>(Power_spektr),Kind,coef_l,MagField_spektr,sf,GxxExp[i]);
        fi(NumberOfPoints,static_cast<int>(Power_spektr),Kind,coef_t,MagField_spektr,sf,GxyExp[i]);
      }
    }

}



void  mobilitySpectrum::BAS(const int n, const int m, const int L, const myFloat x1, Data_spektr & x, Data_spektr & t)
{
    int k;
    myFloat z,r,denominator;

 denominator=(x[n]-x[0]);
 if ( denominator==0 )
 z = 2.0*(x1-x[0])-1.0;         //--------------------АХТУНГ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else
 z = 2.0*(x1-x[0])/(x[n]-x[0])-1.0;
 t[0]=1.0;
 t[1]=z;

 for (k = 1 ; k<= m-1; ++k )
 {
 r = z*t[k];
  switch (L) {
    case 1:
    r = r-t[k-1]/4.0;
    break;
    case 2:
    r = 2*r-t[k-1];
    break;
    case 3:
    r=((k+k+1)*r-k*t[k-1])/static_cast<myFloat>(k+1);
    break;
    default:
    ;//ShowMessage("Case Error");
    break;
  }  
  t[k+1]=r;
  }
 }

void  mobilitySpectrum::gram(const int N,const int m,const int l,
                             Data_spektr & x, Data_spektr & f, mat & a)
{
    // Работает не правильно ибо возвращает нам тонну нулей.
    int i,j,k;
    myFloat q,r,s;
    Data_spektr t;
    vector< vector <myFloat> > p;

  SetLength(p,MaxPoints,5*MaxPoints);
  t.resize(MaxPoints);
 for (i= 0 ; i<= N; ++i )
 {
   BAS(N,m,l,x[i],x,t);
   for (j = 0 ; j<= m; ++j )
     p[j][i]=t[j];
 }
 for (k = 0 ; k<= m; ++k )
 {
   r = 0;
  for (j = k ; j<= m; ++j )
  {
  s = 0;
  r = 0;
   for (i= 0 ; i<= N;++i )
   {
    q = p[k][i];
    s +=q*p[j][i];
     if ( j==m )
      r += q*f[i];
   }
   a[k][j]=s;
   a[j][k]=s;
  }
  a[k][m+1]=r;  // Warning you are not right.
 }
}

void  mobilitySpectrum::gauss(const int N, mat & a, Data_spektr & x)
 {
    int i,j,k,k1,n1;
    myFloat r,s;
   n1 = N+1;
   for (k = 0 ;k <= N;++k )
   {
        k1 = k+1;
        s = a[k][k];

        if ( s == 0 )
            s = 1;

        for (j = k1 ; j<= n1; ++j )
          a[k][j]=a[k][j]/s;
        for (i= k1 ; i<= N; ++i )
        {
          r = a[i][k];
          for (j = k1 ; j<= n1; ++j )
          a[i][j]=a[i][j]-a[k][j]*r;
        }
    }
     for (i= N ; i>= 0; --i )
     {
      s = a[i][n1];
      for (j = i+1 ; j<= N; ++j )
      s = s-a[i][j]*x[j];
      x[i]=s;
    }
  }

void  mobilitySpectrum::fi(const int n, const int m, const int l, Data_spektr & c, Data_spektr & x,
    const myFloat x1, myFloat& s)
{
    int i;
    Data_spektr t;
 t.resize(MaxPoints);
 s = c[0];
 BAS(n,m,l,x1,x,t);
 for (i= 1 ; i<= m; ++i )
     s = s+c[i]*t[i];
}


void  mobilitySpectrum::MakeLagranj()
{
  myFloat X,Y,Y1,Y2;
  X = MagField_spektr[NumberOfPoints]-(MagField_spektr[NumberOfPoints]-
        MagField_spektr[NumberOfPoints-1])/4;
  GetCoef(GxxExp,MagField_spektr,X,Y,Y1,Y2);
  GxxExp[NumberOfPoints+1]=GxxExp[NumberOfPoints];
  GxxExp[NumberOfPoints]=Y;
  GetCoef(GxyExp,MagField_spektr,X,Y,Y1,Y2);
  GxyExp[NumberOfPoints+1]=GxyExp[NumberOfPoints];
  GxyExp[NumberOfPoints]=Y;
  MagField_spektr[NumberOfPoints+1]=MagField_spektr[NumberOfPoints];
  MagField_spektr[NumberOfPoints]=X;
  ++NumberOfPoints;
}



void  mobilitySpectrum::Tred2(const int n, Dat1 & d,Dat1 & e,
                 Dat2 & a,Dat2 & z, bool & fail)
{
    int i,j,k,l;
  myFloat f,g,h,hh;
fail = false;
for (i = 1 ; i<= n; ++i )
  for (j = 1 ; j<= i; ++j )
    {
        z[i][j] = a[i][j];
    }

for (i = n ; i>= 2; --i )
{
  l = i-2;
  f = z[i][i-1];
  g = 0;
  for (k = 1 ; k<= l; ++k ) g = g+z[i][k]*z[i][k];

  h = g+f*f;

  ++l;
  if ( f>=0.0 )
    { e[i]=-sqrt(h);
        g=-sqrt(h);
    }
    else
        {
            e[i]=sqrt(h);
            g = sqrt(h);
        }
  h = h-f*g;
  z[i][i-1]=f-g;
  f = 0.0;

  for (j = 1 ; j<= l; ++j )
  {
    z[j][i]=z[i][j]/h;
    g = 0.0;
    for (k = 1 ; k<= j; ++k ) g = g+z[j][k]*z[i][k];
    for (k = j+1 ; k<= l; ++k ) g = g+z[k][j]*z[i][k];
    e[j]=g/h;
    f = f+g*z[j][i];
  }

  hh = f/(h+h);
//  matrix reduce
  for (j = 1 ; j<= l; ++j )
  {
    f = z[i][j];
g = e[j]-hh*f;
e[j]=e[j]-hh*f;
    for (k = 1 ; k<= j; ++k )
      z[j][k]=z[j][k]-f*e[k]-g*z[i][k];
  }

//{Skip:}
  d[i]=h;
} //{ i }

  d[1] = 0;
  e[1] = 0;
// { store transformation matrix }
  for (i = 1 ; i<= n; ++i )
  {
    l = i-1;
    if ( d[i]!=0. )
    {
      for (j = 1 ; j<= l; ++j )
      {
        g = 0.;
        for (k = 1 ; k<= l; ++k ) g = g+z[i][k]*z[k][j];
        for (k = 1 ; k<= l; ++k ) z[k][j]=z[k][j]-g*z[k][i];
      }
    }
      d[i]=z[i][i];
      z[i][i]=1;
      for (j = 1 ; j<= l; ++j ) {
        z[i][j]=0;
        z[j][i]=0;
      }
  }
} // { Tred2 }

/*
{ **************************************************************** }

{ J.H.Wilkinson, C.Reinch. Handbook for Automatic Computation.
  Linear Algebra. Springer Verlag: Heidelberg, NewYork, Berlin.
  algorithm II-4.
  void  calculate all eigenvalues and eigenvectors
  of real tridiagonal symmetric matrix by QL-algorithm
  with implicit shift.

  Data:
  n - matrix order;
  macheps - smallest possible value such thet 1+macheps>1;
  d - array [1..n] - containts diagonal elements of matrix;
  e - array [1..n] - containts underdiagonal elements of matrix;
   e[1] - arbitraray;
  z - array [1..n,1..n] - unit matrix or matrix of reduction
   ; <= tridiagonal form by Tred2.

  Results:
  d - eigenvalues in increasing order;
  z - eigenvectors (z[i,k],k=1..n);
  e - information is losed.

  Test:
  d=(1,1e2,1e4,1e6,1e8,1e10,1e12);
  e=(0,10,1e3,1e5,1e7,1e9,1e11);
   eigenvalues:         eigenvector for min eigenvalue

  -9.46347415648e8      9.99949513816e-1
  -9.46346919727e2     -1.00974711993e-5
   9.99899020189e-1    -9.99849548745e-3
   1.04633721478e3      9.99850548602e-4
   1.00989903020e6     -1.00974710982e-14
   1.04633771269e9     -9.99850548500e-6
   1.01000000980e12     9.99850548498e-7

  d=(1e12,1e10,1e8,1e6,1e4,1e2,1);
  e=(0,1e11,1e9,1e7,1e5,1e3,10);
   eigenvalues:         eigenvector for min eigenvalue

  -9.46347415707e8      9.99850548432e-7
  -9.46346919876e2     -9.99850548417e-6
   9.99899020157e-1    -1.12526974437e-14
   1.04633721466e3      9.99850548546e-4
   1.00989903020e6     -9.99849548703e-3
   1.04633771265e9     -1.00974740462e-5
   1.01000000980e12     9.99949513813e-1                   }
*/
void  mobilitySpectrum::Imtql2(const int n, const myFloat macheps, Dat1 & d, Dat1 & e,
                 Dat2 & z, bool & fail)
// label Test,Cont,Fail_exit;
{
int i,ia,j,k,m,its;
myFloat h,c,p,q,s,t,u;


fail = false;
for (i= 2 ; i<= n; ++i )
    e[i-1]=e[i];
e[n]=0.0;
k = n-1;
for (j = 1 ; j<= n; ++j )
{
  its = 0;

// { searching of negligibly small underdiagonal element}
Test:
  for (m = j ; m<= k; ++m )
    if ( fabs(e[m])<=macheps*(fabs(d[m])+fabs(d[m+1])) )
        goto Cont;
  m = n;

Cont:
  u = d[j];
  if ( m!=j )
  {
    if ( its==30 )
        {
            fail = true;
            return;
        }
    ++its;
// { formation of shift }
    q = (d[j+1]-u)/(2.0*e[j]);
    t = sqrt(1.0+q*q);
    if ( q<0.0 )
        q = d[m]-u+e[j]/(q-t);
    else q = d[m]-u+e[j]/(q+t);
    u = 0.;
    s = 1.0;
    c = 1.0;

    for (i= m-1 ; i>= j; --i )
    {
      p = s*e[i];
      h = c*e[i];
      if ( fabs(p)>=fabs(q))
      {
        c = q/p;
        t = sqrt(c*c+1.0);
        e[i+1]=p*t;
        s = 1./t;
        c = c*s;
      }
      else
      {
        s = p/q;
        t = sqrt(s*s+1.0);
        e[i+1]=q*t;
        c = 1./t;
        s = s*c;
      }

      q = d[i+1]-u;
      t=(d[i]-q)*s+2.0*c*h;
      u = s*t;
      d[i+1]=q+u;
      q = c*t-h;

// { calculation of eigenvector }
      for (ia = 1 ; ia<= n; ++ia )
      {
        p = z[ia][i+1];
        z[ia][i+1]=s*z[ia][i]+c*p;
        z[ia][i]=c*z[ia][i]-s*p;
      } //{ ia }
    } //{ i }
    d[j]=d[j]-u;
    e[j]=q;
    e[m]=0.;
    goto Test;
  } //{ m!=j }
} //{ j }

// { sorting of eigenvalues and eigenvectors }
for (i= 1 ; i<= n; ++i )
{
  k = i;
  p = d[i];
  for (j = i+1 ; j<= n; ++j )
    if ( d[j]<p )
        { k = j;
            p = d[j];
        }
  if ( k!=i )
  {
    d[k]=d[i];
    d[i]=p;
    for (j = 1 ; j<= n; ++j )
    { p = z[j][i];
        z[j][i]=z[j][k];
        z[j][k]=p;
         }
  }
} //{ i }
} //{ Imtql2 }

myFloat mobilitySpectrum::GetElem(const int j1,const int k1,const int i1)
 {
      myFloat s;
  int ii;
    s = 0;
    for  (ii = i1 ; ii<= (NumberOfPoints-(j1-2)); ++ii )
    {
     if ( j1==2 )
     {

       if ( ii!=k1 )
          s = s+B_spektr[ii]*B_spektr[ii];
     }
     else
        if ( ii!=k1 )
        s = s+B_spektr[ii]*B_spektr[ii]*GetElem(j1-1,k1,ii+1);
     }
    return s;
 }


void  mobilitySpectrum::MakeMatrC()
{

 for (int j = 1 ; j<= NumberOfPoints; ++j )
  for (int k = 1 ; k<= NumberOfPoints; ++k )
  {
    Cr[j][k] = 0;
    Cl[j][k] = 0;
    Cr_t[j][k] = 0;
    Cl_t[j][k] = 0;
   }

 for (int j = 1 ; j<= NumberOfPoints; ++j )
 for (int k = 1 ; k<= NumberOfPoints; ++k )
  {
   if ( j==1 )
   {
    Cr[j][k]=1;
   }

   else
   {
     Cr[j][k]=GetElem(j,k,1);
   }
     Cl[j][k]=-Cr[j][k]*B_spektr[k];
   }
}

void  mobilitySpectrum::MakeMatrA()
 {
  for (int j = 1 ; j<= NumberOfPoints; ++j)
  for (int k = 1 ; k<= NumberOfPoints; ++k)
  {
    Am[j][k]=0;
    }

  for (int i= 1 ; i<= NumberOfPoints; ++i )
   for (int j = 1 ; j<= NumberOfPoints; ++j )
    if ((i+j)%2==1)
      for (int k = 1 ; k<= NumberOfPoints; ++k )
           Am[i][j]=Am[i][j]+Gxy_sp[k]*Cl_t[k][(i+j-1) >> 1];
               else
       for (int k = 1 ; k<= NumberOfPoints; ++k )
           Am[i][j]=Am[i][j]+Gxx_sp[k]*Cr_t[k][(i+j) >> 1];
 }

//if su=0 then det(C)=0
void  mobilitySpectrum::InverseMatrC(Dat2 & Ci,Dat2 & C,myFloat & Su,const int NP)
{
    Dat3 at;
    myFloat sr;

    SetLength(at,2*MaxPoints+1,2*MaxPoints+1); // SetLength(at,MaxPoints+1,2*MaxPoints+1);
    for (int i= 1 ; i<= NP; ++i )
        for (int j = 1 ; j<= NP; ++j )
            at[i][j]=Ci[i][j];

    for (int i = 1 ; i<= NP; ++i )
    {
        for (int j = NP+1 ; j<= 2*NP; ++j )
            at[i][j]=0;
        at[i][i+NP]=1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    int j;
    for (int k = 1 ;k <= NP; ++k )
    {
        Su = at[k][k];
        j = k;
        for (int i= k+1 ; i<= NP; ++i)
        {
            sr = at[i][k];
            if ( fabs(sr)>fabs(Su) )
            {
                Su = sr;
                j = i;
            }
        }
        if ( Su==0 )
            return;//--------------------
        if ( j!=k )
            for (int i= k ; i<= 2*NP; ++i )
            {
                sr = at[k][i];
                at[k][i]=at[j][i];
                at[j][i]=sr;
            }
        for (int j = k+1 ; j<= 2*NP; ++j )
            at[k][j]=at[k][j]/Su;
        for (int i = k+1 ; i<= NP; ++i )
        {
            sr = at[i][k];
            for (int j = k+1 ; j<= 2*NP; ++j )
                at[i][j]=at[i][j]-at[k][j]*sr;
        }
    }
    if ( Su!=0 )
    for (int j = NP+1 ; j<= 2*NP; ++j )
        for (int i = NP-1 ; i>= 1; --i )
        {
            sr = at[i][j];
            for (int k = i+1 ; k<= NP; ++k )
                sr = sr-at[k][j]*at[i][k];
            at[i][j]=sr;
        }
    if ( Su!=0 )
        for (int i= 1 ; i<= NP; ++i )
            for(int j = 1 ; j<= NP; ++j )
                C[i][j]=at[i][j+NP];
}

myFloat mobilitySpectrum::S_s(const myFloat Mi)

   {
int im,js;
myFloat a,bb;
    // { make Vm }
     Mv[1]=1;
     if ( fabs(Mi)<1e-26 )
      for (im = 2 ; im<= NumberOfPoints; ++im ) Mv[im]=0;
     else
     for (im = 2 ; im<= NumberOfPoints; ++im )
      if ( (im)%2==1 )
       Mv[im]=exp((im-1)*log(fabs(Mi)));
      else
       if ( Mi<0 )
           Mv[im]=exp((im-1)*log(fabs(Mi)));
       else
           Mv[im]=-exp((im-1)*log(fabs(Mi)));
/*(*     { compute ¦Vm¦¤}
     vm = 0;
     for is = 1 ; <= NPoint )
       vm = vm+Mv[is]*Mv[is];
     { compute am/Go for Hp=0}
     a = 0;
     for is = 1 ; <= NPoint ) a = a+(Mv[is]*Mv[is]);*)*/
     a = 1;
     //{ compute am/Go }
     for (im = 1 ; im<= NumberOfPoints; ++im)
         a = a*(1+B_spektr[im]*Mi*B_spektr[im]*Mi);
     for (im = 1 ; im<= NumberOfPoints; ++im)
     {
         Vpr[im]=0;
     for (js = 1 ; js<= NumberOfPoints; ++js)
          Vpr[im]+=Qm[js][im]*Mv[js];
      Vpr[im]=Vpr[im]*Vpr[im];
     }
     bb = 0;
     for (im = 1 ; im<= NumberOfPoints; ++im ) {
//{      if ( (abs(Lv[im])<1e-30)and(Vpr[im]<1e-15) ) continue}
       if ( (fabs(Lv[im])==0) && (Vpr[im]>1e-15) ) { a = 0;break;}
      bb = bb+Vpr[im]/fabs(Lv[im]); }
     return (a/bb);
   }

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void  mobilitySpectrum::MobilitySpectrumFunc(TLineSeries& LineSeries1, TLineSeries& Series5)
{
   myFloat Sf,lm;
   int j,k,i;
   int Lmin,Lmax;

   for (i = 0 ; i<= NumberOfPoints; ++i )
    {
     B_spektr[i+1]=MagField_spektr[i];
     Gxx_sp[i+1]=GxxExp[i];
     Gxy_sp[i+1]=GxyExp[i];
    }
   ++NumberOfPoints;
   MakeMatrC();
   InverseMatrC(Cr,Cr_t,Sf,NumberOfPoints);
   if ( Sf==0 )
   {
     ;//ShowMessage("Определитель равен нулю!");
     --NumberOfPoints;
     return;
   }
   if ( MagField_spektr[0]==0 )
   {
        for (j = 1 ;j <= NumberOfPoints-1; ++j )
        for (k = 1 ; k<= NumberOfPoints-1; ++k )
            Cm[j][k]=Cl[j][k+1];

        InverseMatrC(Cm,Cm_t,Sf,NumberOfPoints-1);
        for (j = 1 ; j<= NumberOfPoints; ++j )
            for (k = 1 ; k<= NumberOfPoints;++k )
                if ( j==1 )
                    Cl_t[1][k]=0;
                else if ( k==NumberOfPoints )
                    Cl_t[j][NumberOfPoints]=0;
                else
                    Cl_t[j][k]=Cm_t[j-1][k];
   }
   else
   {
    InverseMatrC(Cl,Cl_t,Sf,NumberOfPoints);
   }

    MakeMatrA();
    // интересно что это за числа ниже?
    // параметр 6e-4913 в tred2 не используется О_о
    Tred2(NumberOfPoints,Lv,Xv,Am,Qm,bulua);

    Imtql2(NumberOfPoints,LDBL_EPSILON,Lv,Xv,Qm,bulua); // 5.42e-20 // 5.42e-80 // LDBL_EPSILON == 1e-19

    eigenValues=Lv;
    eigenVectors=Qm;
    // В Lv сохраняется собственные значения матрицы
    // В Qm - собственные вектора матрицы.
    LineSeries1.clear();
    Series5.clear();
   Lmin = MSLeft;
   Lmax = MSRight;
   SizeData=static_cast<int>((Lmax-Lmin+1)*PointPerInt+1); // возможно тут придется домножать на размер элемента
   // впрочем нет - там же сейчас другая функция.
   InitArray2();
   k = 0;
   for (i = 0 ;i <= (Lmax-Lmin); ++i )
    {
     lm = exp((Lmin+i)*logl(10));
     Sf = lm;
     for (j = 1 ; j<= PointPerInt-1; ++j )
      {
       Mobility[k]=Sf;
       Spectr_e[k]=S_s(-Sf);
       Spectr_p[k]=S_s(Sf);
       LineSeries1.push_back(make_pair(Mobility[k],Spectr_e[k]));
       Series5.push_back(make_pair(Mobility[k],Spectr_p[k]));
       //LineSeries1.AddXY(Mobility[k],Spectr_e[k],"",clTeeColor);
       //Series5.AddXY(Mobility[k],Spectr_p[k],"",clTeeColor);
       Sf = lm*exp(j/static_cast<myFloat>(PointPerInt)*logl(10));
       ++k;
      // if ( sf>10 ) break;
      }
    // if ( sf>10 ) break;
    }
    GridPoints = k-1;
   --NumberOfPoints;
}


 /*   {
      Итак, по сути, здесь происходит тыканье значений максимумов и
      отображение их в таблице под графиком спектра подвижности.
      Т.е. по факту этот код будет уже на сях и здесь он не нужен.
      В комментарии его.



void  chtSpectrClickSeries(chtSpectr: TCustomChart;
  Series: TChartSeries; ValueIndex: int ;
  Shift: TShiftState; X, Y: int ; Series5:TCustomSeries;
  LineSeries1:TLineSeries; StringGrid3:TStringGrid;RowInFocus:int );
//var Mu,G_e,G_p,con_p,con_e:myFloat;
 //     Ind, Ind2:longint;


{  }
    {
with chtSpectr )
   {

    ind = Series5.GetCursorValueIndex;
    Ind2:= LineSeries1.GetCursorValueIndex;
    if ( ind<>-1 )
      {                 // это по дыркам.
       G_p = Series5.YValue[ind]; // наверное проводимость
       Mu = Series5.XValue[ind];  // подвижность
       con_p = g_p/(mu*1.602e-19); // концентрация
       StringGrid3.Cells[1,RowInFocus]:=to_string(con_p); // показываем выбранные результаты
       StringGrid3.Cells[2,RowInFocus]:=to_string(Mu);
       if ( RowInFocus<3 )
         Inc(RowInFocus) // и да, переключаем фокус таблицы
       else
         RowInFocus = 1;
      }
      if ( ind2<>-1 )
      {
       G_e = LineSeries1.YValue[ind2]; // это по электронам
       Mu = Series5.XValue[ind2];
       con_e = g_e/(mu*1.602e-19);
       StringGrid3.Cells[1,RowInFocus]:=to_string(-con_e);
       StringGrid3.Cells[2,RowInFocus]:=to_string(-Mu);
       if ( RowInFocus<3 )
         Inc(RowInFocus)
       else
         RowInFocus = 1;
      }
   }

} }

// Запись результатов Спектра подвижности
// Количество точек
// Минимумы и максимумы по носителям зарядов
// И дальше неизвестно что:)
// Но нам надо то, что оно выводило на графики.
*/
////////////////////////////////////////////////////////////////////////////
/////////////////////// КОНЕЦ "ХОЛЛ. ПОДВИЖНОСТЬ"///////////////////////////
////////////////////////////////////////////////////////////////////////////

long double mobilitySpectrum::getResultEY(const int i)
{
  return static_cast<long double>(electronMobilitySpectrum.operator [](i).second); // YValue[i];
}

long double  mobilitySpectrum::getResultEX(const int i)
{
  return static_cast<long double>(electronMobilitySpectrum.operator [](i).first);
}

long double mobilitySpectrum::getResultHY(const int i)
{
  return static_cast<long double>(holeMobilitySpectrum.operator [](i).second);
}

long double mobilitySpectrum::getResultHX(const int i)
{
  return static_cast<long double>(holeMobilitySpectrum.operator [](i).first);
}


// входные параметры есть, надо подумать что мы возвращаем, где оно хранится
// и как мы это будем возвращать.

size_t mobilitySpectrum::getResultSize()
{
    return electronMobilitySpectrum.size();
}

mobilitySpectrum::mobilitySpectrum(Data_spektr &MagneticFieldP, Data_spektr &Exx, Data_spektr &Exy, int size)
{
    PointPerInt=50*50; // 50*50
    MaxPoints = size;
    if ( MaxPoints>1000 )
    MaxPoints = 11;


    /*imtql2 testing*/
/*
    Data_spektr d={0,1,1e2,1e4,1e6,1e8,1e10,1e12};
    Data_spektr e={0,0,10,1e3,1e5,1e7,1e9,1e11};
    mat z;
    z.resize(8);
    for (int i = 0; i < 8; ++i) {
        z[i].resize(8);
    }
    bool fail;
    Imtql2(7,5.42e-20,d,e,z,fail);
*/


    TLineSeries Gxx, Gxy;
    TLineSeries ExpXX, ExpXY;

    NumberOfPoints = MaxPoints-1;

    MagField_spektr.resize(MaxPoints);
    GxxExp.resize(MaxPoints);
    GxyExp.resize(MaxPoints);

    // Оказывается эти две переменные... внезапно пределы расчета подвижности
    // Ну т.е. их степени, типа 10^-2 и 10^1
    // И кстати, надо сделать функцию для их изменения.
    // Либо передавать их при вызове функции по расчету спектра.
    MSLeft=-3;
    MSRight = 1;

    B_spektr.resize(MaxPoints+1);
    Gxx_sp.resize(MaxPoints+1);
    Gxx_MC.resize(MaxPoints+1);
    Gxy_MC.resize(MaxPoints+1);
    Gxy_sp.resize(MaxPoints+1);
    Xr.resize(MaxPoints+1);
    Lv.resize(MaxPoints+1);
    Xv.resize(MaxPoints+1);
    Mv.resize(MaxPoints+1);
    Vpr.resize(MaxPoints+1);

    SetLength(Am,MaxPoints+1,MaxPoints+1);
    SetLength(Qm,MaxPoints+1,MaxPoints+1);
    SetLength(Cl,MaxPoints+1,MaxPoints+1);
    SetLength(Cr,MaxPoints+1,MaxPoints+1);
    SetLength(Cl_t,MaxPoints+1,MaxPoints+1);
    SetLength(Cr_t,MaxPoints+1,MaxPoints+1);
    SetLength(Cm,MaxPoints+1,MaxPoints+1);
    SetLength(Cm_t,MaxPoints+1,MaxPoints+1);

    for (int i = 0 ; i< MaxPoints; ++i )
    {
        MagField_spektr[i]=MagneticFieldP[i];
        GxxExp[i]=Exx[i];
        GxyExp[i]=Exy[i];
    }

    //MakeInterpolate(Gxx,Gxy,ExpXX,ExpXY);
    MakeMNK(true,Gxx,Gxy,ExpXX,ExpXY);
    MobilitySpectrumFunc(electronMobilitySpectrum,holeMobilitySpectrum);
    for (size_t i = 0; i < electronMobilitySpectrum.size(); ++i)
    {
      resultElectronConductivity.push_back(static_cast<long double>(electronMobilitySpectrum.operator [](i).second));
      resultHoleConductivity.push_back(static_cast<long double>(holeMobilitySpectrum.operator [](i).second));
      resultMobility.push_back(static_cast<long double>(electronMobilitySpectrum.operator [](i).first));
    }
    if(resultMobility.size()!=0)
    {
    //cout << "If inside it \n";
      

    size_t index = searchSignificantPeak(resultElectronConductivity, 0, resultMobility[1]-resultMobility[0]); // электроны
    //cout << "anyone can here me?\n";
      if(index!=resultMobility.size()) // если это не правый край графика
      {
        size_t index2 = searchSignificantPeak(resultElectronConductivity, index+40, resultMobility[1]-resultMobility[0]); // электроны
        if (index2!=resultMobility.size())
        {
          extremumElectronIndex.push_back(index2);
        }
        else
        {
          extremumElectronIndex.push_back(index);
        }
      }
      else
      {
        index = searchSignalSlowdown(resultElectronConductivity, 0, resultMobility[1]-resultMobility[0]); // электроны
        if(index!=resultMobility.size()) // если это не правый край графика
        {
          // Этот пик может быть зеркальным. Пока просто смотрим есть ли пик дальше - если есть, то ставим его.
          size_t index2 = searchSignalSlowdown(resultElectronConductivity, index+40, resultMobility[1]-resultMobility[0]); // электроны
          if(index2!=resultMobility.size())
          {
            extremumElectronIndex.push_back(index2);
          }
          else
          {
          extremumElectronIndex.push_back(index);
          }
        }
      }
      //cout << "middle\n";

      index = searchSignificantPeak(resultHoleConductivity, 0, resultMobility[1]-resultMobility[0]); // тяжелые дырки
      if(index!=resultMobility.size())
      {
          extremumHoleIndex.push_back(index);
      }

      index = searchSignificantPeak(resultHoleConductivity, index+4, resultMobility[1]-resultMobility[0]); // легкие дырки
      if(index!=resultMobility.size())
      {
          extremumHoleIndex.push_back(index);
      }
      else
      {
        if (extremumHoleIndex.size()!=0)
        {
        index = searchSignalSlowdown(resultHoleConductivity, extremumHoleIndex.back()+40, resultMobility[1]-resultMobility[0]); // легкие дырки
        if(index!=resultMobility.size())
        {
          extremumHoleIndex.push_back(index);
        }
        }
      }

      //cout << "okay\n";

      if (extremumElectronIndex.size()==0)
      {
        extremumElectronIndex.push_back(0);
      }
      if(extremumHoleIndex.size()==0)
      {
        extremumHoleIndex.push_back(0);
        extremumHoleIndex.push_back(0);
      }
      if(extremumHoleIndex.size()==1)
      {
        extremumHoleIndex.push_back(0);
      }

      
      }
}
//-----------------------------------------------------------------------

long double mobilitySpectrum::calcConcentrationFromGp(long double G_p, long double Mu)
{
    const long double electronCharge=1.60217656535E-19;// Кл
    return G_p/(Mu*electronCharge);
}

void mobilitySpectrum::getExtremums(TSignal & holeConcentration, TSignal & holeMobility, TSignal & electronConcentration, TSignal & electronMobility)
{
  holeConcentration.clear();
  holeMobility.clear();
  electronConcentration.clear();
  electronMobility.clear();

  for (size_t i = 0; i < extremumHoleIndex.size(); ++i)
  {
    holeConcentration.push_back(calcConcentrationFromGp(resultHoleConductivity[extremumHoleIndex[i]],resultMobility[extremumHoleIndex[i]]));
    holeMobility.push_back(resultMobility[extremumHoleIndex[i]]);
  }
  for (size_t i = 0; i < extremumElectronIndex.size(); ++i)
  {
    electronConcentration.push_back(calcConcentrationFromGp(resultElectronConductivity[extremumElectronIndex[i]],resultMobility[extremumElectronIndex[i]]));
    electronMobility.push_back(resultMobility[extremumElectronIndex[i]]);
  }
}

size_t mobilitySpectrum::searchSignalSlowdown(TSignal &y, size_t startPosition, long double inh)
{
    myFloat h = static_cast<myFloat>(inh);
  size_t size=y.size();
    /*
    Функция должна реагировать на замедление изменения сигнала.
    Т.е. на малозаметный пик.
    */
    if(startPosition>=size)
        return size;

      int dsize=size-2;
      size_t d2size=dsize-2;
      // Посчитаем производную методом конечных разностей
      // формула df/dx=1/h*(2*f(x+h)-f(x+2h)/2-3/2*f(x))

      // формула df/dx=(f(x+h)-f(x-h))/2/h;
      std::vector<myFloat> dY(size);

      for(int i =0;i<dsize;i++)
      {
        dY[i]=(y[i+2]-y[i])/2.0/h;
      }

      std::vector<myFloat> d2Y(size);

      for(size_t i =0;i<d2size;i++)
      {
        //d2Y[i]=(dY[i+2]-dY[i])/2.0/h;
        d2Y[i]=(y[i]-2*y[i+1]+y[i+2])/h/h;
      }

      for (int i = startPosition; i < dsize-1; ++i)
        {
            /*
            Поиск такой:
            1. ищем участок на котором первая производная отрицательная.
            */
            while (i<dsize-1 && (dY[i]>0)) ++i;

            while (i<dsize-1 && (dY[i]<0 && dY[i+1]-dY[i]<0)) ++i; // первая производная отрицательная и уменьшается

            while (i<dsize-1 && (dY[i]<0 && dY[i+1]-dY[i]>0)) ++i; // первая производная отрицательная и увеличивается

            --i;
            // Теперь проверим остальные условия - вторая производная должна изменить знак.
            // Вероятно стоит расширить диапазон поиска до +-10 значений
            if (i+10<dsize-1 && i-10>=0 && (d2Y[i-10]>0 && d2Y[i+10]<0) )
            {
                return i;
            }
            if( i<(int)startPosition)
            {
            return size;
            }
        }
        return size;

}

size_t mobilitySpectrum::searchSignificantPeak(TSignal &y, size_t startPosition, long double h)
{
    size_t size=y.size();
    if(startPosition>=size)
        return size;

    size_t dsize=size-2;
    //size_t d2size=dsize-2;
    /*for(auto i = y.begin(); i!=y.end();++i)
    {
        cout << *i << endl;
    }*/

    std::vector<long double> dY=calculateFirstDerivative(y,h);

    std::vector<long double> d2Y=calculateSecondDerivative(y,h);

    /*
      Поиск пиков. Считаем производные первого и второго порядков.
      Самые ярко выраженные пики должны иметь такие признаки:
      1. Первая производная сначала растет и положительна.
      2. После пика - убывает и отрицательна 

      3. Вторая производная до пика - положительна.
      4. После пика - отрицательна.
      5. Три-четыре точки, там где находится пик - скачки производных, что логично, т.к. там должны быть точки разрыва.
    */
    cout << "start search\n"<< startPosition<< " " << dsize << endl;
    for (size_t i = startPosition; i < dsize-1; ++i)
    {
        /*
        Поиск такой:
        1. ищем участок на котором первая производная положительна и растет.
        */
        while (i<dsize-1 && (dY[i]<0)) ++i;

        // ищем точку, с которой рост первой производной прекращается
        while (i<dsize-1 && (dY[i]>0)) ++i;

        --i;
        // Теперь проверим остальные условия - вторая производная отрицательна
        // Вероятно стоит расширить диапазон поиска до +-10 значений
        if (d2Y[i]<0 )
        {
            return i;
        }
    }
    return size;
}



size_t mobilitySpectrum::searchPeakLeftBorder(std::vector<myFloat> dh,std::vector<myFloat> d2h, size_t index)
{
  // Ищет границу пика слева, т.е. такую точку, после которой функция меняет направление изменения (начинает увеличиваться).
  size_t i=1;
  while(index-i<d2h.size() && d2h[index-i]>d2h[index-i+1])
  {
    ++i;
  }
  --i;
  i=index-i;
  return i;
}
size_t mobilitySpectrum::searchPeakRigthBorder(std::vector<myFloat> dh,std::vector<myFloat> d2h, size_t index)
{
  // Ищет границу пика справа. 
  size_t j=1;
  while(index+j<d2h.size() && d2h[index+j]<d2h[index+j+1])
  {
    ++j;
  }
  j=index+j;

  if(j==index+1)
  {
  j=1;
      while(index+j<d2h.size() && d2h[index+j]>d2h[index+j+1])
    {
      ++j;
    }
    while(index+j<d2h.size() && d2h[index+j]<d2h[index+j+1])
    {
      ++j;
    }
    j=index+j;
  }

  return j;
}

void mobilitySpectrum::constructPeakCriteria(PeaksCriteria & peaksCriteria, TStringList * tsl, const std::vector<myFloat> & resMob, const std::vector<myFloat> & resCond, int index, int i, int j)
{
  myFloat peakHeigh = resCond[index]-(resCond[i]+resCond[j])/2.0;
  myFloat peakWidth = resultMobility[j]-resultMobility[i];
  myFloat symmetri=0;
  myFloat peakVelocityR=0;
  myFloat peakVelocity2R=0;
  myFloat peakVelocityL=0;
  myFloat peakVelocity2L=0;
  if(index+(i+j)/4+1 < (int)resultMobility.size() && index-(i+j)/4 >=0)
  {
    symmetri = resCond[i]-resCond[j];
    peakVelocityR = fabs(resCond[index+(i+j)/4]-resCond[index+(i+j)/4+1]);
    peakVelocity2R = peakVelocityR/(resultMobility[index+(i+j)/4]-resultMobility[index+(i+j)/4+1]);

    peakVelocityL = fabs(resCond[index-(i+j)/4]-resCond[index-(i+j)/4+1]);
    peakVelocity2L = peakVelocityL/(resultMobility[index-(i+j)/4]-resultMobility[index-(i+j)/4+1]);

  }
  if(index+(i+j)/4+1 >= (int)resultMobility.size() && index > 10)
  {
    size_t t=resultMobility.size()-index-2;
    peakVelocityR = fabs(resCond[index+t]-resCond[index+t+1]);
    peakVelocity2R = peakVelocityR/(resultMobility[index+t]-resultMobility[index+t+1]);
    peakVelocityL = fabs(resCond[index-t]-resCond[index-t+1]);
    peakVelocity2L = peakVelocityL/(resultMobility[index-t]-resultMobility[index-t+1]);
  }

  peaksCriteria.peakHeigh=static_cast<long double>(peakHeigh);
  peaksCriteria.peakWidth=static_cast<long double>(peakWidth);
  peaksCriteria.peakWidthOrd=j-i;
  peaksCriteria.peakVelocityR=static_cast<long double>(peakVelocityR);
  peaksCriteria.peakVelocity2R=static_cast<long double>(peakVelocity2R);
  peaksCriteria.symmetri=static_cast<long double>(symmetri);// симметричность пиков
  peaksCriteria.peakVelocityL=static_cast<long double>(peakVelocityL); // модуль скорости изменения слева от пика
  peaksCriteria.peakVelocity2L=static_cast<long double>(peakVelocity2L); // скорость изменения слева от пика

  tsl->push_back(to_string(static_cast<long double>(peakHeigh))+"\t"+to_string(static_cast<long double>(peakWidth))+
                 "\t"+to_string(j-i)+"\t"+to_string(static_cast<long double>(peakVelocityR))+"\t"+
                 to_string(static_cast<long double>(peakVelocity2R)));
  tsl->push_back(to_string(static_cast<long double>(symmetri))+"\t"+
                 to_string(static_cast<long double>(peakVelocityL))+"\t"+to_string(static_cast<long double>(peakVelocity2L))); // симметричность пиков и скорость изменения слева от него

}



void mobilitySpectrum::calculatePeakWeigth(PeaksCriteria & peaksCriteria,TStringList * tsl, 
  TSignal & resultMobility, TSignal & resultConductivity,
  std::vector<myFloat> & d, std::vector<myFloat> & d2, size_t extremumIndex)
{
  size_t index=extremumIndex;
  int i=searchPeakLeftBorder(d,d2,index);
  size_t j=searchPeakRigthBorder(d,d2,index);
  constructPeakCriteria(peaksCriteria,tsl, resultMobility, resultConductivity,index,i,j);

  peaksCriteria.peakLeftBorderFirstDerivative=static_cast<long double>(d[i+1]);
  peaksCriteria.peakFirstDerivative=static_cast<long double>(d[index]);
  peaksCriteria.peakRightBorderFirstDerivative=static_cast<long double>(d[j-1]);
  peaksCriteria.peakLeftMiddleFirstDerivative=static_cast<long double>(d[(i+index)/2]);
  peaksCriteria.peakRightMiddleFirstDerivative=static_cast<long double>(d[(j+index)/2]);

  peaksCriteria.peakLeftBorderSecondDerivative=static_cast<long double>(d2[i+1]);
  peaksCriteria.peakSecondDerivative=static_cast<long double>(d2[index]);
  peaksCriteria.peakRightBorderSecondDerivative=static_cast<long double>(d2[j-1]);
  peaksCriteria.peakLeftMiddleSecondDerivative=static_cast<long double>(d2[(i+index)/2]);
  peaksCriteria.peakRightMiddleSecondDerivative=static_cast<long double>(d2[(j+index)/2]);

  tsl->push_back(to_string(static_cast<long double>(d[i+1]))+"\t"+to_string(static_cast<long double>(d[index]))+"\t"+to_string(static_cast<long double>(d[j-1])));
  tsl->push_back(to_string(static_cast<long double>(d[(i+index)/2]))+"\t"+to_string(static_cast<long double>(d[(j+index)/2])));

  tsl->push_back(to_string(static_cast<long double>(d2[i+1]))+"\t"+to_string(static_cast<long double>(d2[index]))+"\t"+
          to_string(static_cast<long double>(d2[j-1])));
  tsl->push_back(to_string(static_cast<long double>(d2[(i+index)/2]))+"\t"+to_string(static_cast<long double>(d2[(j+index)/2])));
}

long double mobilitySpectrum::calculatePeaksWeigth()
{
  PeaksCriteria peaksCriteria;
  vPC.clear();
  // есть пик.
  // в две стороны от него запускаем поиск точки перегиба
  // как только нашли - это конец пика
  // Вычисляем высоту пика относительно этой точки
  // Ширину пика, относительно 0,606 высоты
  TStringList *tsl = new TStringList;

  long double hMobility = static_cast<long double>(resultMobility[1]-resultMobility[0]); // шаг по подвижности

  std::vector<myFloat> de=calculateFirstDerivative(resultElectronConductivity, hMobility); // электроны
  std::vector<myFloat> dh=calculateFirstDerivative(resultHoleConductivity, hMobility); // тяжелые дырки

  std::vector<myFloat> d2e=calculateSecondDerivative(resultElectronConductivity, hMobility); // электроны
  std::vector<myFloat> d2h=calculateSecondDerivative(resultHoleConductivity, hMobility); // тяжелые дырки

for (size_t k = 0; k < extremumHoleIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,resultMobility, resultHoleConductivity,dh,d2h,extremumHoleIndex[k]);
  vPC.push_back(peaksCriteria);
}

for (size_t k = 0; k < extremumElectronIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,resultMobility, resultElectronConductivity,de,d2e,extremumElectronIndex[k]);
  vPC.push_back(peaksCriteria);
  }

  additionalData.leftPointElectronConductivity=resultElectronConductivity[0];
  additionalData.leftPointHoleConductivity=resultHoleConductivity[0];

  tsl->push_back(to_string(resultElectronConductivity[0])+"\t"+to_string(resultHoleConductivity[0])); // первые точки

  // добавить минимальные значения спектров

  // добавить максимальную скорость изменения (хотя с этим могут быть проблемы)

  // расстояние между пиками?

  additionalData.leftPointMobility=resultMobility[0];
  additionalData.rightPointMobility=resultMobility.back();

  additionalData.rightPointElectronConductivity=resultElectronConductivity.back();
  additionalData.rightPointHoleConductivity=resultHoleConductivity.back();

  tsl->push_back(to_string(resultElectronConductivity.back())+"\t"+to_string(resultHoleConductivity.back())); // последние точки


  /*
  После всего этого - пробую сделать всё тоже самое, только для логарифмов нашего спектра.
  */

  TSignal log10resultElectronConductivity;
  TSignal log10resultHoleConductivity;
  TSignal log10resultMobility;

  for (size_t i = 0; i < resultElectronConductivity.size(); ++i)
  {
    log10resultElectronConductivity.push_back(log10l(resultElectronConductivity[i]));
    log10resultHoleConductivity.push_back(log10l(resultHoleConductivity[i]));
    log10resultMobility.push_back(log10l(resultMobility[i]));
  }

  myFloat loghMobility = log10resultMobility[1]-log10resultMobility[0];

  std::vector<myFloat> del=calculateFirstDerivative(log10resultElectronConductivity, loghMobility); // электроны
  std::vector<myFloat> dhl=calculateFirstDerivative(log10resultHoleConductivity, loghMobility); // тяжелые дырки

  std::vector<myFloat> d2el=calculateSecondDerivative(log10resultElectronConductivity, loghMobility); // электроны
  std::vector<myFloat> d2hl=calculateSecondDerivative(log10resultHoleConductivity, loghMobility); // тяжелые дырки

for (size_t k = 0; k < extremumHoleIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,log10resultMobility, log10resultHoleConductivity,dhl,d2hl,extremumHoleIndex[k]);
  vPC.push_back(peaksCriteria);
}

for (size_t k = 0; k < extremumElectronIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,log10resultMobility, log10resultElectronConductivity,del,d2el,extremumHoleIndex[k]);
  vPC.push_back(peaksCriteria);
}

  additionalData.leftPointElectronConductivityLog=log10resultElectronConductivity[0];
  additionalData.leftPointHoleConductivityLog=log10resultHoleConductivity[0];

  tsl->push_back(to_string(log10resultElectronConductivity[0])+"\t"+to_string(log10resultHoleConductivity[0])); // первые точки

  // добавить минимальные значения спектров

  // добавить максимальную скорость изменения (хотя с этим могут быть проблемы)

  // расстояние между пиками?

  additionalData.rightPointElectronConductivityLog=log10resultElectronConductivity.back();
  additionalData.rightPointHoleConductivityLog=log10resultHoleConductivity.back();

  tsl->push_back(to_string(log10resultElectronConductivity.back())+"\t"+to_string(log10resultHoleConductivity.back())); // последние точки
  delete tsl;
  return 1;
}

long double mobilitySpectrum::calculatePeaksWeigth(std::string filename)
{

  PeaksCriteria peaksCriteria;

  vPC.clear();
  // есть пик.
  // в две стороны от него запускаем поиск точки перегиба
  // как только нашли - это конец пика
  // Вычисляем высоту пика относительно этой точки
  // Ширину пика, относительно 0,606 высоты
  TStringList *tsl = new TStringList;

  long double hMobility = static_cast<long double>(resultMobility[1]-resultMobility[0]); // шаг по подвижности

  std::vector<myFloat> de=calculateFirstDerivative(resultElectronConductivity, hMobility); // электроны
  std::vector<myFloat> dh=calculateFirstDerivative(resultHoleConductivity, hMobility); // тяжелые дырки

  std::vector<myFloat> d2e=calculateSecondDerivative(resultElectronConductivity, hMobility); // электроны
  std::vector<myFloat> d2h=calculateSecondDerivative(resultHoleConductivity, hMobility); // тяжелые дырки

for (size_t k = 0; k < extremumHoleIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,resultMobility, resultHoleConductivity,dh,d2h,extremumHoleIndex[k]);
  vPC.push_back(peaksCriteria);
  /*
  size_t index=extremumHoleIndex[k];
  int i=searchPeakLeftBorder(dh,d2h,index);
  size_t j=searchPeakRigthBorder(dh,d2h,index);
  constructPeakCriteria(peaksCriteria,tsl, resultMobility, resultHoleConductivity,index,i,j);  
  */
}

for (size_t k = 0; k < extremumElectronIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,resultMobility, resultElectronConductivity,de,d2e,extremumElectronIndex[k]);
  vPC.push_back(peaksCriteria);
  /*
  int index=extremumElectronIndex[k];

  int i=searchPeakLeftBorder(de,d2e,index);

  size_t j=searchPeakRigthBorder(de,d2e,index);
  
  constructPeakCriteria(peaksCriteria,tsl, resultMobility, resultElectronConductivity,index,i,j);

  */
  }

  additionalData.leftPointElectronConductivity=resultElectronConductivity[0];
  additionalData.leftPointHoleConductivity=resultHoleConductivity[0];

  tsl->push_back(to_string(resultElectronConductivity[0])+"\t"+to_string(resultHoleConductivity[0])); // первые точки

  // добавить минимальные значения спектров

  // добавить максимальную скорость изменения (хотя с этим могут быть проблемы)

  // расстояние между пиками?

  additionalData.rightPointElectronConductivity=resultElectronConductivity.back();
  additionalData.rightPointHoleConductivity=resultHoleConductivity.back();

  tsl->push_back(to_string(resultElectronConductivity.back())+"\t"+to_string(resultHoleConductivity.back())); // последние точки


  /*
  После всего этого - пробую сделать всё тоже самое, только для логарифмов нашего спектра.
  */

  TSignal log10resultElectronConductivity;
  TSignal log10resultHoleConductivity;
  TSignal log10resultMobility;

  for (size_t i = 0; i < resultElectronConductivity.size(); ++i)
  {
    log10resultElectronConductivity.push_back(log10l(resultElectronConductivity[i]));
    log10resultHoleConductivity.push_back(log10l(resultHoleConductivity[i]));
    log10resultMobility.push_back(log10l(resultMobility[i]));
  }

  myFloat loghMobility = log10resultMobility[1]-log10resultMobility[0];

  std::vector<myFloat> del=calculateFirstDerivative(log10resultElectronConductivity, loghMobility); // электроны
  std::vector<myFloat> dhl=calculateFirstDerivative(log10resultHoleConductivity, loghMobility); // тяжелые дырки

  std::vector<myFloat> d2el=calculateSecondDerivative(log10resultElectronConductivity, loghMobility); // электроны
  std::vector<myFloat> d2hl=calculateSecondDerivative(log10resultHoleConductivity, loghMobility); // тяжелые дырки

for (size_t k = 0; k < extremumHoleIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,log10resultMobility, log10resultHoleConductivity,dhl,d2hl,extremumHoleIndex[k]);
  vPC.push_back(peaksCriteria);
  /*
  size_t index=extremumHoleIndex[k];

  int i=searchPeakLeftBorder(dhl,d2hl,index);

  size_t j=searchPeakRigthBorder(dhl,d2hl,index);

  constructPeakCriteria(peaksCriteria,tsl, log10resultMobility, log10resultHoleConductivity,index,i,j);

  */
}

for (size_t k = 0; k < extremumElectronIndex.size(); ++k)
{
  calculatePeakWeigth(peaksCriteria,tsl,log10resultMobility, log10resultElectronConductivity,del,d2el,extremumHoleIndex[k]);
  vPC.push_back(peaksCriteria);
  /*
  int index=extremumElectronIndex[k];

  int i=searchPeakLeftBorder(del,d2el,index);

  size_t j=searchPeakRigthBorder(del,d2el,index);
  
  constructPeakCriteria(peaksCriteria,tsl, log10resultMobility, log10resultElectronConductivity,index,i,j);

  */
  }

  additionalData.leftPointElectronConductivityLog=log10resultElectronConductivity[0];
  additionalData.leftPointHoleConductivityLog=log10resultHoleConductivity[0];

  tsl->push_back(to_string(log10resultElectronConductivity[0])+"\t"+to_string(log10resultHoleConductivity[0])); // первые точки

  // добавить минимальные значения спектров

  // добавить максимальную скорость изменения (хотя с этим могут быть проблемы)

  // расстояние между пиками?

  additionalData.rightPointElectronConductivityLog=log10resultElectronConductivity.back();
  additionalData.rightPointHoleConductivityLog=log10resultHoleConductivity.back();

  tsl->push_back(to_string(log10resultElectronConductivity.back())+"\t"+to_string(log10resultHoleConductivity.back())); // последние точки


  /*for (int i = 0; i < d2e.size(); ++i)
  {
    tsl->Add(to_string(resultMobility[i])+"\t"+to_string(de[i])+"\t"+to_string(dh[i])+"\t"+to_string(d2e[i])+"\t"+to_string(d2h[i]));
  }*/

  SaveToFile(tsl,filename);


  delete tsl;
  return 1;
}


void mobilitySpectrum::getExtremumHoleIndex(std::vector<size_t> &v)
{
  v=extremumHoleIndex;
}
void mobilitySpectrum::getExtremumElectronIndex(std::vector<size_t> &v)
{
  v=extremumElectronIndex;
}

void mobilitySpectrum::getResults(TSignal & mobility, TSignal & HoleConductivity, TSignal & ElectronConductivity)
{
  mobility=resultMobility;
  HoleConductivity=resultHoleConductivity;
  ElectronConductivity=resultElectronConductivity;
}

std::vector<myFloat> mobilitySpectrum::getEigenValues()
{
  return eigenValues;
}
std::vector< std::vector<myFloat> > mobilitySpectrum::getEigenVectors()
{
  return eigenVectors;
}

void mobilitySpectrum::savePeakWeigth(std::string filename)
{
  calculatePeaksWeigth(filename);
}


void mobilitySpectrum::saveEigenValues(std::string filename)
{
  TStringList *tsl = new TStringList();

  for (size_t i = 0; i < eigenVectors.size(); ++i)
  {
    string t;
    for (size_t j = 0; j < eigenVectors[i].size(); ++j)
    {
      t=t+to_string(static_cast<long double>(eigenVectors[i][j]))+"\t";
    }
    tsl->push_back(t);
  }

  tsl->push_back("\n");

  for (size_t i = 0; i < eigenValues.size(); ++i)
  {
    tsl->push_back(to_string(static_cast<long double>(eigenValues[i])));
  }

  SaveToFile(tsl,filename);

  delete tsl;
}


void mobilitySpectrum::saveResults(std::string filename)
{
  TStringList *tsl = new TStringList();

  TSignal holeConcentration;
  TSignal holeMobility;
  TSignal electronConcentration;
  TSignal electronMobility;

  this->getExtremums(holeConcentration, holeMobility, electronConcentration, electronMobility);

  for (size_t i = 0; i < holeMobility.size(); ++i)
  {
    tsl->push_back(to_string(holeConcentration[i])+"\t"+to_string(holeMobility[i]));
  }

  for (size_t i = 0; i < electronMobility.size(); ++i)
  {
    tsl->push_back(to_string(electronConcentration[i])+"\t"+to_string(electronMobility[i]));
  }

  SaveToFile(tsl,filename);


  delete tsl;
}


std::vector <PeaksCriteria> mobilitySpectrum::getPeaksCriteria()
{
  return vPC;
}

AdditionalData mobilitySpectrum::getAdditionalData()
{
  return additionalData;
}
