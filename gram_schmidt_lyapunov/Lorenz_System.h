#include<cmath>
#include<cstdlib>

const double sigma=16;
const double rho=45.92;
const double beta=4;

double f0(double x[], double t)
{
    return sigma*(x[1]-x[0]);
}

double f1(double x[], double t)
{
    return x[0]*(rho-x[2])-x[1];
}

double f2(double x[], double t)
{
    return x[0]*x[1]-(beta)*(x[2]);
}

double f3(double x[], double t)
{
    return sigma*(x[4]-x[3]);
}

double f4(double x[], double t)
{
    return (rho-x[2])*x[3]-x[0]*x[5]-x[4];
}

double f5(double x[], double t)
{
    return x[3]*x[1]+x[0]*x[4]-(beta)*(x[5]);
}
