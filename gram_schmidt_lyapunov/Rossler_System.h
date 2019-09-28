#include<cmath>
#include<cstdlib>

const double a=0.15;
const double b=0.20;
const double c=10.0;

double f0(double x[], double t)
{
    return -(x[1]+x[2]);
}
double f1(double x[], double t)
{
    return x[0]+a*x[1];
}
double f2(double x[], double t)
{
    return b+x[2]*(x[0]-c);
}
double f3(double x[], double t)
{
    return -(x[4]+x[5]);
}
double f4(double x[], double t)
{
    return x[3]+a*x[4];
}
double f5(double x[], double t)
{
    return x[2]*x[3]+(x[0]-c)*x[5];
}