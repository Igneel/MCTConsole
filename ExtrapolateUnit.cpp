#include "ExtrapolateUnit.h"
//---------------------------------------------------------------------------

#define DEBUG 0

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Модуль экстраполяции.
// Вычисление определителя.

bool testExtrapolateUnit()
{
	// Попытаюсь протестировать все функции по очереди
	// Расчет определителя
	int size=5;
	long double ** A;
	A= new long double * [size];
	for (int i = 0; i < size; ++i)
	{
		A[i]= new long double [size];
	}

	A[0][0]=4;
	A[0][1]=3;
	A[0][2]=-4;
	A[0][3]=2;
	A[0][4]=2;
	A[1][0]=2;
	A[1][1]=-1;
	A[1][2]=-3;
	A[1][3]=-4;
	A[1][4]=2;
	A[2][0]=3;
	A[2][1]=1;
	A[2][2]=1;
	A[2][3]=2;
	A[2][4]=-1;
	A[3][0]=1;
	A[3][1]=2;
	A[3][2]=3;
	A[3][3]=4;
	A[3][4]=-1;
	A[4][0]=-1;
	A[4][1]=1;
	A[4][2]=-1;
	A[4][3]=-2;
	A[4][4]=3;

	if (determinant(A, size)!=-48.0)
	{
		return false;
	}

	TMatrix B(5);

	B[0].push_back(4);
	B[0].push_back(3);
	B[0].push_back(-4);
	B[0].push_back(2);
	B[0].push_back(2);
	B[1].push_back(2);
	B[1].push_back(-1);
	B[1].push_back(-3);
	B[1].push_back(-4);
	B[1].push_back(2);
	B[2].push_back(3);
	B[2].push_back(1);
	B[2].push_back(1);
	B[2].push_back(2);
	B[2].push_back(-1);
	B[3].push_back(1);
	B[3].push_back(2);
	B[3].push_back(3);
	B[3].push_back(4);
	B[3].push_back(-1);
	B[4].push_back(-1);
	B[4].push_back(1);
	B[4].push_back(-1);
	B[4].push_back(-2);
	B[4].push_back(3);

	if (determinant(B, size)!=-48.0)
	{
		return false;
	}

//curveFittingUniversal

	TSignal inX;
	TSignal koef;
	koef.push_back(2);
	koef.push_back(3);
	koef.push_back(4);
	koef.push_back(5);
	koef.push_back(6);
	koef.push_back(7);

	TSignal outF;
	inX.push_back(0);
	for (int i = 0; i < 101; ++i) // надо много где проверить количество точек.
	{
		inX.push_back(inX.back()+0.02);
	}
	// должно посчитаться выражение вида 
	// 2x^5+3x^4+4x^3+5x^2+6x+7
	calculatePolinomByKoef(inX, koef, outF);

	TSignal ans;
ans.push_back(7.00000000000000);
ans.push_back(7.12203248640000);
ans.push_back(7.24826388480000);
ans.push_back(7.37890443520000);
ans.push_back(7.51417743360000);
ans.push_back(7.65432000000000);
ans.push_back(7.79958384640000);
ans.push_back(7.95023604480000);
ans.push_back(8.10655979520000);
ans.push_back(8.26885519360000);
ans.push_back(8.43744000000000);
ans.push_back(8.61265040640000);
ans.push_back(8.79484180480000);
ans.push_back(8.98438955520000);
ans.push_back(9.18168975360000);
ans.push_back(9.38716000000000);
ans.push_back(9.60124016640000);
ans.push_back(9.82439316480000);
ans.push_back(10.0571057152000);
ans.push_back(10.2998891136000);
ans.push_back(10.5532800000000);
ans.push_back(10.8178411264000);
ans.push_back(11.0941621248000);
ans.push_back(11.3828602752000);
ans.push_back(11.6845812736000);
ans.push_back(12.0000000000000);
ans.push_back(12.3298212864000);
ans.push_back(12.6747806848000);
ans.push_back(13.0356452352000);
ans.push_back(13.4132142336000);
ans.push_back(13.8083200000000);
ans.push_back(14.2218286464000);
ans.push_back(14.6546408448000);
ans.push_back(15.1076925952000);
ans.push_back(15.5819559936000);
ans.push_back(16.0784400000000);
ans.push_back(16.5981912064000);
ans.push_back(17.1422946048000);
ans.push_back(17.7118743552000);
ans.push_back(18.3080945536000);
ans.push_back(18.9321600000000);
ans.push_back(19.5853169664000);
ans.push_back(20.2688539648000);
ans.push_back(20.9841025152000);
ans.push_back(21.7324379136000);
ans.push_back(22.5152800000000);
ans.push_back(23.3340939264000);
ans.push_back(24.1903909248000);
ans.push_back(25.0857290752000);
ans.push_back(26.0217140736000);
ans.push_back(27.0000000000000);
ans.push_back(28.0222900864000);
ans.push_back(29.0903374848000);
ans.push_back(30.2059460352000);
ans.push_back(31.3709710336000);
ans.push_back(32.5873200000000);
ans.push_back(33.8569534464000);
ans.push_back(35.1818856448000);
ans.push_back(36.5641853952000);
ans.push_back(38.0059767936000);
ans.push_back(39.5094400000000);
ans.push_back(41.0768120064000);
ans.push_back(42.7103874048000);
ans.push_back(44.4125191552000);
ans.push_back(46.1856193536000);
ans.push_back(48.0321600000000);
ans.push_back(49.9546737664000);
ans.push_back(51.9557547648000);
ans.push_back(54.0380593152000);
ans.push_back(56.2043067136000);
ans.push_back(58.4572800000000);
ans.push_back(60.7998267264000);
ans.push_back(63.2348597248000);
ans.push_back(65.7653578752000);
ans.push_back(68.3943668736000);
ans.push_back(71.1250000000000);
ans.push_back(73.9604388864000);
ans.push_back(76.9039342848000);
ans.push_back(79.9588068352000);
ans.push_back(83.1284478336000);
ans.push_back(86.4163200000000);
ans.push_back(89.8259582464000);
ans.push_back(93.3609704448000);
ans.push_back(97.0250381952000);
ans.push_back(100.821917593600);
ans.push_back(104.755440000000);
ans.push_back(108.829512806400);
ans.push_back(113.048120204800);
ans.push_back(117.415323955200);
ans.push_back(121.935264153600);
ans.push_back(126.612160000000);
ans.push_back(131.450310566400);
ans.push_back(136.454095564800);
ans.push_back(141.627976115200);
ans.push_back(146.976495513600);
ans.push_back(152.504280000000);
ans.push_back(158.216039526400);
ans.push_back(164.116568524800);
ans.push_back(170.210746675200);
ans.push_back(176.503539673600);
ans.push_back(183.000000000000);

	for (int i = 0; i < 101; ++i)
	{
		if (fabs(outF[i]-ans[i])>0.00001 )
		{
			return false;
		}
	}

	return true;
}


long double determinant (long double ** Arr,int size)
{
	int i,j;
	long double det=0;
	long double** matr;
	switch (size) {
	case 1: return Arr[0][0];
	case 2: return Arr[0][0]*Arr[1][1]-Arr[0][1]*Arr[1][0];
	default:
	  matr=new long double*[size-1];
					for(i=0;i<size;++i)
					{
							for(j=0;j<size-1;++j)
							{
									if(j<i)
											matr[j]=Arr[j];
									else
											matr[j]=Arr[j+1];
							}
							det+=pow((long double)-1, (i+j))*determinant(matr, size-1)*Arr[i][size-1];
					}
					delete[] matr;
	return det;
	}
}

long double determinant (std::vector<std::vector<long double> > & Arr,int size)
{
	int i,j;
	long double det=0;
	std::vector<std::vector<long double> > matr;
	switch (size) {
	case 1: return Arr[0][0];
	case 2: return Arr[0][0]*Arr[1][1]-Arr[0][1]*Arr[1][0];
	default:
		matr.resize(size-1);
					for(i=0;i<size;++i)
					{
							for(j=0;j<size-1;++j)
							{
									if(j<i)
											matr[j]=Arr[j];
									else
											matr[j]=Arr[j+1];
							}
							det+=pow((long double)-1, (i+j))*determinant(matr, size-1)*Arr[i][size-1];
					}
	  return det;
	}
}


// Функция вычисляет подгоночные коэффициенты.
// Для поперечного сигнала.
// Напрямую не вызывать, пользоваться extrapolate.

// Возвращает кол-во точек
//Принимает указатель на входные данные, минимальное и максимальное значения по х, величину шага, указатель на выходные данные.

//---------------------------------------------------------------------------
//-----------curveFittingUniversal----------------------------------------------
int curveFittingUniversal(long double * inX, long double *inY, const int length,long double *outKoef,const int powPolinom)
{

	// как много кода.
	// по сути - формируем матрицу
	// перемножаем её транпонированную на её изначальную
	// дальше решаем систему методом Крамера (это тут), а вообще любый другим методом тоже.

	int a=powPolinom+1;
	// выражение f(x)=x^2*p1+x*p2+p3

	long double ** fullMatrix=new long double*[length];

	for (int i = 0; i < length; i++) {
		fullMatrix[i]=new long double[a+1];

	}

	long double **K=new long double *[a+1];
	for (int i = 0; i < a+1; i++) {
        K[i]=new long double [a+1];
	}

	long double **K5= new long double* [a];
	long double *Ks=new long double [a];

	long double ***delta=new long double **[a];  // матрицы для вычисления корней, один из столбцов - сободные члены.

		for (int i = 0; i < a; i++) {
			 delta[i]=new long double*[a];
			for (int j = 0; j < a; j++) {
				delta[i][j]=new long double[a];

		}}

	for (int i = 0; i < a; i++)
		K5[i]= new long double[a];
	long double *p=new long double [a];
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------

	/*
	Итак, много временных переменных:
	*/
    long double temp1,temp2,temp3;
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	/*
	формируем матрицу вида

	x^a x^a-1 ...  x^0 y
	*/
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < a; j++) {
		// копируем значения в матрицу
		// по столбцам, шестой - столбец единиц

			fullMatrix[i][j]=((a-j-1)==0?1:pow(inX[i],a-j-1));

		}
		// последний столбец формируется здесь из значений игрек.

		fullMatrix[i][a]=inY[i];
	}

	for (int i = 0; i <= a; i++) {
		for (int j = 0; j <= a; j++) {
			for (int k=0; k < length; k++) {
				K[i][j]+=fullMatrix[k][i]*fullMatrix[k][j];
	}}}
	// перемножение работает, проверено матлабом
	// K5 - содержит a строк и a столбцов, т.е. без свободных членов.


	// Ks - столблец из 2*a строк, свободные члены.


	for (int i = 0; i < a-1; i++) {
		for (int j = 0; j < a; j++) {

		K5[i][j]=K[i][j]; // копируем первые 2a-1 строк
		}
	}

	// а 6ая строка - среднее арифметическое 6 и 7 строк
	for (int i = 0; i < a; i++) {

		K5[a-1][i]=(K[a-1][i]+K[a][i])/2;
	}

	for (int i = 0; i < a-1; i++) {

		Ks[i]=K[i][a];
	}
	Ks[a-1]=(K[a-1][a]+K[a][a])/2;

	long double d0=determinant(K5,a); // детерминант правилен +

	for (int i = 0; i < a; i++) {
		for (int j = 0; j < a; j++) {
			for (int k = 0; k < a; k++) {
				if (k!=i)
				{
					// составляем матрицы по которым будем находить детерминанты для переменных.
					 delta[i][j][k]=K5[j][k];
				}
				else
				{
					// столбец нужной переменной заменяется свободными членами.
					delta[i][j][k]=Ks[j];
				}
	}}}

	for (int i = 0; i < a; i++) {
		outKoef[i]=p[i]=determinant(delta[i],a)/d0;
	}


	delete[] p;
	// прибираемся.
	for (int i = 0; i < a; i++) {
		delete[] K5[i];
	}

	delete[] K5;


	for (int i = 0; i < a+1; i++) {
        delete[] K[i];
	}
	delete[] K;

	delete[] Ks;

	for (int i = 0; i < length; i++) {
		delete[] fullMatrix[i];
	}
	delete[] fullMatrix;
	for (int i = 0; i < a; i++) {
			for (int j = 0; j < a; j++) {
				delete[] delta[i][j];
		}
				 delete[] delta[i];
		}

	delete[] delta;

	return a;
}



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//-----------curveFittingUniversal----------------------------------------------
int curveFittingUniversal(std::vector<long double> * inX, std::vector<long double> *inY, std::vector<long double> *outKoef,const int powPolinom)
{
	unsigned int lenght=(*inY).size();
    if(lenght!=inX->size())
    {
    //("Длина массивов в экстраполяции не совпадает!");
    return 0;
    }
	// как много кода.
	// по сути - формируем матрицу
	// перемножаем её транпонированную на её изначальную
	// дальше решаем систему методом Крамера (это тут), а вообще любый другим методом тоже.

	int a=powPolinom+1;
	// выражение f(x)=x^2*p1+x*p2+p3

	std::vector<std::vector<long double> > fullMatrix(lenght);
    std::vector<std::vector<long double>> K;
    K.resize(a+1);
    for (auto i=K.begin();i!=K.end();++i)
    {
        i->resize(a+1);
    }
    std::vector<std::vector<long double>> K5;
    K5.resize(a);
    for (auto i=K5.begin();i!=K5.end();++i)
    {
        i->resize(a);
    }
	std::vector<long double> Ks(a);
    std::vector<std::vector<std::vector<long double> > > delta;
    delta.resize(a);
    for (auto i=delta.begin();i!=delta.end();++i)
    {
        i->resize(a);
    }
	std::vector<long double> p(a);

	for(int i=0;i<a; ++i)
	{

		for (int j=0; j < a; ++j) {
            delta[i][j].resize(a);
		}
	}

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------

	/*
	Итак, много временных переменных:
	*/
    long double temp1,temp2,temp3;
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	/*
	формируем матрицу вида

	x^a x^a-1 ...  x^0 y
	*/
	for (int i = 0; i < lenght; ++i) {
		for (int j = 0; j < a; ++j) {
		// копируем значения в матрицу
		// по столбцам, предпоследний - столбец единиц

			fullMatrix[i].push_back(((a-j-1)==0?1:powl((*inX)[i],a-j-1)));
		}
		// последний столбец формируется здесь из значений игрек.
		fullMatrix[i].push_back((*inY)[i]);
	}

	for (int i = 0; i <= a; ++i) {
		for (int j = 0; j <= a; ++j) {
			for (int k=0; k < lenght; ++k) {
				K[i][j]+=fullMatrix[k][i]*fullMatrix[k][j];
	}}}
	// перемножение работает, проверено матлабом
	// K5 - содержит a строк и a столбцов, т.е. без свободных членов.


	// Ks - столблец из 2*a строк, свободные члены.


	for (int i = 0; i < a-1; ++i) {
		for (int j = 0; j < a; ++j) {
		K5[i][j]=K[i][j]; // копируем первые 2a-1 строк
		}
	}

	// а 6ая строка - среднее арифметическое 6 и 7 строк
	for (int i = 0; i < a; ++i) {
		K5[a-1][i]=(K[a-1][i]+K[a][i])/2;
	}

	for (int i = 0; i < a-1; ++i) {
		Ks[i]=K[i][a];
	}
	Ks[a-1]=(K[a-1][a]+K[a][a])/2;

	long double d0=determinant(K5,a); // детерминант правилен +

	for (int i = 0; i < a; ++i) {
		for (int j = 0; j < a; ++j) {
			for (int k = 0; k < a; ++k) {
				if (k!=i)
				{
					// составляем матрицы по которым будем находить детерминанты для переменных.
					 delta[i][j][k]=K5[j][k];
				}
				else
				{
					// столбец нужной переменной заменяется свободными членами.
					delta[i][j][k]=Ks[j];
				}
	}}}

	for (int i = 0; i < a; ++i) {
		(*outKoef)[i]=p[i]=determinant(delta[i],a)/d0;
	}

	return a;
}
//---------------------------------------------------------------------------




void LagrangeExtrapolation(TSignal & inX, TSignal & inY, TSignal & newX, TSignal & outY)
{
	outY.clear();
	
	size_t size=inX.size();
	for (int k = 0; k < newX.size(); ++k)
	{
	long double out=0;
	for (int i = 0; i < size; ++i)
	{
		long double l=1;
		for (int j = 0; j < size; ++j)
		{
			if (i!=j)
			{
				l=l*(newX[k]-inX[j])/(inX[i]-inX[j]);
			}
		}
		out+=inY[i]*l;
	}
	outY.push_back(out);
	}
}
