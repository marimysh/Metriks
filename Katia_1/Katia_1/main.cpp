#include <iostream>
#include <math.h>
#include <stdio.h>

const int U2F = 100;
const double PI = 3.1415926, T = 50, f = 0.02, A0 = 0.0005;

double Ua, Ub, Uc, Ud, Ue, Uf, Udd;
int t1,t2;                                      //какой номер функции
int kglob;                                      //для удобства подсчета Фурье

inline double fU (double t, int var1, int var2)         //гибкий вариант метрики. возможны все 6вариантов. var1знак, var2 период
{
    return var2*U2F*sqrt(2)*sin(2*PI*f*t + var1*2*PI/3);
}

inline double abs(double x)                             //тупо модуль числа
{
    if (x>0)
        return x;
    else
        return -x;
}

void Max6 (double u1, double u2, double u3, double u4, double u5, double u6, double &res, int &t1, int &t2);    //максимальная из 6 метрик. выводчисло и параметры максимальнойметрики
double MSimpson (double a, double b, double (*f)(double, int, int), bool fl2, int n);       //метод симпсона, n - кол-во шагов
double MSimpsonToch (double a, double b, double (*f)(double, int, int), double toch);       //точность для использования м. симпсона
void GraphU (double a, double b, int var1, int var2, bool Ub);                      //вывод значений метрики
void GraphUMax (double a, double b);                                                //вывод значений максимальной по всем метрикам
// a, b - границы, откуда - до куда считать

void Max6 (double u1, double u2, double u3, double u4, double u5, double u6, double &res, int &t1, int &t2)
{							// & - чтобы полученное значение можно было использовать
    res = u1;
    res = std::max(res, u2);
    res = std::max(res, u3);
    res = std::max(res, u4);
    res = std::max(res, u5);
    res = std::max(res, u6);
    if (res == u1)
    {
        t1 = 0;
        t2 = 1;
        return;
    }
    if (res == u2)
    {
        t1 = -1;
        t2 = 1;
        return;
    }
    if (res == u3)
    {
        t1 = 1;
        t2 = 1;
        return;
    }
    if (res == u4)
    {
        t1 = 0;
        t2 = -1;
        return;
    }
    if (res == u5)
    {
        t1 = -1;
        t2 = -1;
        return;
    }
    if (res == u6)
    {
        t1 = 1;
        t2 = -1;
        return;
    }
}

double MSimpson (double a, double b, double (*f)(double, int, int), bool fl2, int n)
{
    if (fl2)
        n /= 2;                             //для удобства подсчета с точностью
    double h = (b-a)/n;                     //размер шага
    double sum1 = 0, sum2 = 0;              //чет/нечет суммы
    bool odd = false;                       //отвечает за чет-нечет ход
    for (int i=1; i<n; ++i)
    {
        double x = a + i*h;
        if (odd)			//Сумма значений подынтегральной функции на концах четных внутренних отрезков
            sum2 += f(x, t1, t2);
        else                //... на концах нечетных внутренних отрезков
            sum1 += f(x, t1, t2);
        odd = not odd;
    }
    return h*(f(a, t1, t2)+f(b, t1, t2)+4*sum1+2*sum2)/3;
}

double MSimpsonToch (double a, double b, double (*f)(double, int, int), double toch)
{
    double n = 12,              //начальный шаг
           MS1=15, MS2=0;       //I^h, I^(h/2)
    while (toch < abs(MS1-MS2)/15)
    {
        MS1 = MSimpson(a,b,f,false,n);
        MS2 = MSimpson(a,b,f,true,n);       //шаг в 2 раза меньше
        n++;
//      std::cout << MS1 << " " << MS2 << " " << n << " " << (MS1-MS2)/15 << std::endl;
    }
//    std::cout << "n= " << n << std::endl;
    return MS1;
}

void GraphU (double a, double b, int var1, int var2, bool Ub)
{
    int n = 300;                    //кол-во шагов
    double h = (b-a)/n;             //длина шага
    if (!Ub)                        //просто вывести метрику
        for (int i=0; i<=n; ++i)
        {
            double x = a +i*h;
            std::cout << x << " " << fU(x,var1,var2) << std::endl;
        }
/*    else                            //которая зануляется на 4м диоде
    {
        bool temp0 = false;         //показывается на 4м диоде или нет
        for (int i=0; i<=n; ++i)
        {
            double x = a +i*h;
            if (!temp0 && (abs(fU(x,-1,1) - fU(x,0,1)) < 0.0001) && (fU(x,-1,1) > 0))   //заходим на 4 диод
                temp0 = true;
            if (temp0 && (abs(fU(x,-1,1) - fU(x,1,1)) < 0.0001))                        //уходим с него
                temp0 = false;
            if (temp0)
                std::cout << x << " " << 0 << std::endl;
            else
                std::cout << x << " " << fU(x,-1,1) << std::endl;
        }
    }*/
    std::cout << "\n\n\n";
    return;
}

void GraphUMax (double a, double b)
{
    int n = 400,                                //n- кол-во шагов. нужен достаточно большой
        var1,var2;                              //просто нужно для работы, значение не используется
    double h = (b-a)/n,
           Udd,                                 //максимальное значение
           min=300;
    for (int i=0; i<=n; ++i)
    {
        double x = a +i*h;
        Max6(fU(x,0,1),fU(x,-1,1),fU(x,1,1),-fU(x,0,1),-fU(x,-1,1),-fU(x,1,1), Udd,var1,var2);
        std::cout << x << " " << Udd << std::endl;
        if (min > Udd)
            min = Udd;
    }
    std::cout << "min = " << min;
    return;
}

inline double fUdsin (double t, int a, int b)
{
    return fU(t,t1,t2)*sin(2*PI*kglob*f*t);
}

inline double fUdcos (double t, int a, int b)
{
    return fU(t,t1,t2)*cos(2*PI*kglob*f*t);
}

inline double ABk (bool AB)
{
    double (*func)(double a, int b, int c);     //указатель на функцию
    if (AB)
        func = fUdsin;
    else
        func = fUdcos;
    return (2/T)*MSimpsonToch(-T/2,T/2,func,0.0001);
}

inline double Mk ()
{                                                   //pow(x,2) = x^2
    return sqrt(pow(ABk(false),2) + pow(ABk(true),2));
}

inline double fi ()
{
    return atan(ABk(false)/ABk(true));
}

double fUds (double t)
{
    double sumk = 0, temp = 10, tempPr = 0;
    int k = 1;
    while ((temp - tempPr) > 0.0001)
    {
       kglob = k;
       tempPr = temp;
       temp = Mk()*sin(2*PI*kglob*f*t + fi());
       sumk += temp;
       ++k;
    }
    return A0/2 + sumk;
}

int main()
{
    freopen("Ud.txt", "w", stdout);
    Ua = fU (0,0,1);
    Ub = fU (0,-1,1);
    Uc = fU (0,1,1);
    Ud = -Ua;
    Ue = -Ub;
    Uf = -Uc;
    double (*func)(double a, int b, int c);         //указатель на функцию
    func = fU;
    Max6 (Ua,Ub,Uc,Ud,Ue,Uf, Udd,t1,t2);
    std::cout << "ud = " << Udd << std::endl;
    std::cout << "Ud = " << MSimpsonToch (-PI,PI,func,0.0001)/(2*PI) << std::endl;
    fclose (stdout);
    freopen ("Graph1_6garm.txt", "w", stdout);
    GraphU (0,50,0,1,false);            //все 6 гармоник
    GraphU (0,50,-1,1,false);
    GraphU (0,50,1,1,false);
    GraphU (0,50,0,-1,false);
    GraphU (0,50,-1,-1,false);
    GraphU (0,50,1,-1,false);
    fclose (stdout);
    freopen ("Graph_Max.txt", "w", stdout);
    GraphUMax (0,50);
    fclose (stdout);
    freopen ("Graph2_6garm.txt", "w", stdout);
    GraphU (0,50,0,1,false);            //без 1 гармоники
    GraphU (0,50,-1,1,true);
    GraphU (0,50,1,1,false);
    GraphU (0,50,0,-1,false);
    GraphU (0,50,-1,-1,false);
    GraphU (0,50,1,-1,false);
    fclose (stdout);
    freopen ("Uds_Furie.txt", "w", stdout);
    for (int t=0; t<50; t++)            //ряд Фурье
        std::cout << fUds (t) << std::endl;
    fclose (stdout);
    freopen ("Ud_Uds.txt", "w", stdout);
    t1 = -1;
    t2 = 1;
    for (int t=0; t<50; t++)            //разница функций Ud - Uds
        std::cout << fUds (t) - fU(t,-1,-1) << std::endl;
    fclose (stdout);
    return 0;
}
