//---------------------------------------------------------------------------
#include "NoiseUnit.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// преобразование Бокса-Мюллера
void BoxMull(std::vector<long double> &z,unsigned long l)
{
	long double x,y;
    long double s;

    int ibound = 1000000;
    long double dbound = ibound;
	// вся эпичность этой ошибки потрясает%)
	// это ж надо было не подумать, что
	// если l - число нечетное, то мы вылезем за границу)
    for(auto i=0u;i<l;i+=2)
	{
		do
		{
            x=(rand()%(2*ibound)-ibound)/dbound;
            y=(rand()%(2*ibound)-ibound)/dbound;
			// шум в пределах от -1 до 1
			s=x*x+y*y;
        }while (s>1.0L || s<=0.0L);
        z[i]=(fmodl((x*sqrtl(-2.0L*logl(s)/s)),2*ibound)-ibound)/dbound;
		if((i+1)<l)
            z[i+1]=(fmodl(y*sqrtl(-2.0L*logl(s)/s),2*ibound)-ibound)/dbound;
	}
}


//--------------------------------------------------------------------------
// расчет среднего квадратичного отклонения
//--------------------------------------------------------------------------

long double Sko (std::vector<long double> const &x0,std::vector<long double> const &x,unsigned long l)
{
	long double z=0;
    for(auto i=0u;i<l;i++)
    z+= powl(fabsl((fabsl(x[i])-fabsl(x0[i]))),2);
    z/=static_cast<long double>(l);
    return sqrtl(z);
}

//---------------------------------------------------
// расчет математического ожидания
//-------------------------------------------------

 long double Mo (std::vector<long double> const &x,unsigned long l)
 {
	 long double M=0;
     for(auto i=0u;i<l;i++)
	 M+=x[i];
	 return M/l;
 }

//---------------------------------------------------------------------------
// возвращает величину математического ожидания шума, ско и ско в % в параметр ret

void ShumAdding(std::vector<long double> const &x,std::vector<long double> &out,
    std::vector<long double> &ret, long double koeff,const unsigned long l)
{
std::vector<long double> y;
y.resize(l);

BoxMull(y,l);
for(auto i=0u;i<l;i++)
{
    out[i]=x[i]+y[i]*koeff;
}
	// амплитуда шума определяется как 1/koef
	ret.push_back(Mo(y,l)); // математическое ожидание
	ret.push_back(Sko(x,out,l));// СКО
	ret.push_back(ret[0]/Mo(x,l)*100);// СКО в %  от математического ожидания исходного графика
}

// Генератор шума, с квантованием по уровню
// Также тут реализована функция, которая квантует сам сигнал по уровню.
void QuantumNoiseGenerator(std::vector<long double> &buf,long double mzr, int amplitude, unsigned long n)
{
/*
Стоит вопрос в частоте нулей, по сути при измерениях у нас нулей скажем 81%
Видимо тут надо сделать аналогично.
*/
    for (auto i = 0u; i < n; ++i) {
		int val=rand()%(50*amplitude+1);
		if(val>amplitude) val=0;
		int sign=rand()%2;
		long double resPoint=(sign==0?1:-1)*val*mzr;
		buf[i]=resPoint;
	}
}

void SignalQuantifier(std::vector<long double> & signal,long double mzr,unsigned long n)
{
     for (auto i = 0u; i < n; ++i) {
        long long temp=floorl(signal[i]/mzr+0.5L);
		signal[i]=temp*mzr;
	}
}
/* -----------------------------------------------------------------------------------
 *
 * Функция QuantumShumAdding
 * Вызывает генератор шума квантования по уровню
 * Добавляет шум к исходному сигналу
 * После чего проводить квантование сигнала по уровню
 * -----------------------------------------------------------------------------------

*/
void QuantumShumAdding(std::vector<long double> const &x,std::vector<long double> &out,
    std::vector<long double> &ret, const unsigned long l,long double mzr, int amplitude)
{
	std::vector<long double> y;
	y.resize(l);
	std::vector<long double> ytemp;
	ytemp.resize(l);
	QuantumNoiseGenerator(y,mzr,amplitude,l);

	std::vector<long double> qx;
	qx.resize(l);
    for(auto i=0u;i<l;i++)
	 {
		qx[i]=x[i];
	 }

	SignalQuantifier(qx,mzr,l);
    out.resize(l);
    for(auto i=0u;i<l;i++)
    {
		out[i]=qx[i]+y[i];
    }

    ret.push_back(Mo(y,l)); // математическое ожидание
    ret.push_back(Sko(x,out,l));// СКО
    ret.push_back(ret[0]/Mo(x,l)*100);// СКО в %  от математического ожидания исходного графика
}
