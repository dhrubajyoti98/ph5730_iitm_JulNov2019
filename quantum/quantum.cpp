#include<complex>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<cstdio>
#include<omp.h>

using namespace std;

#define sigma0 0.05
#define x0 0.25

long double V(int j, double eps)
{
    long double x=eps*j;
    long double V0=-2*pow(50.0*M_PI,2);
    if(x>=0.484 && x<=0.516)
        return V0;
    else 
        return 0;
}

std::complex<long double> initial_wave(int j, double k0, long double eps)
{
    long double x=eps*j;
    double t1=(1.0/(2*pow(sigma0,2)))*pow(x-x0,2);
    //printf("%lf %lf\n",x,t1);
    std::complex<long double> w0(exp(-t1)*cos(x*k0),exp(-t1)*sin(x*k0));
    return w0;
}

std::complex<long double> E_func(int j,long double lambda, double eps)
{
    return std::complex<long double> (2L+pow(eps,2)*V(j,eps),-lambda);
}

std::complex<long double> omega_func(int j,long double lambda, std::complex<long double> wave[], double eps)
{
    return -wave[j+1]+ (std::complex<long double> (2L+pow(eps,2)*V(j,eps),lambda))*wave[j] - wave[j-1];
}

void copy(std::complex<long double> a[], std::complex<long double> b[],int J)
{
    #pragma omp parallel for
    for(int i=0;i<J+1;i++)
    {
        b[i]=a[i];
    }
}

int main(int argc, char* argv[])
{
    int J=(int)atof(argv[1]);
    long double eps=(long double)atof(argv[2]);

    long double delta=2*pow(eps,2);
    long double lambda=2*pow(eps,2)/delta;
    long double k0=atof(argv[3])*M_PI;

    std::complex<long double> e[J+1];
    std::complex<long double> f[J+1];

    std::complex<long double> omega[J+1];

    std::complex<long double> wave[J+1];

    //double t=0;

    e[1]=E_func(1,lambda,eps);
    for(int j=2;j<J+1;j++)
    {
        e[j]=E_func(j,lambda,eps)-(1.0L/e[j-1]);
    }


    for(int j=0;j<J+1;j++)
    {
        if(j==0)
        {
            wave[j]=std::complex<long double> (0,0);
        }
        else if(j==J)
        {
            wave[j]=std::complex<long double> (0,0);
        }
        else
        {
            wave[j]=initial_wave(j,k0,eps);
        }
    }
    int c=0;
    while(c<=2500)
    {
        /*if((t/5.0)-(int)(t/5.0)<0.001)
        {*/   
            char filename[50];
            sprintf(filename,"t%d.dat",c);
            FILE *fp=fopen(filename,"w");
            for(int i=0;i<J+1;i++)
            {
                fprintf(fp,"%Lf %Lf %Lf\n",eps*i,pow(std::abs(wave[i]),2),V(i,eps)/fabs(V(J/2,eps)));
            }
            fclose(fp);
            
        //} 

        printf("#Solving for n=%d\n",c);

        #pragma omp parallel for
        for(int j=1;j<J+1;j++)
        {
            omega[j]=omega_func(j,lambda,wave,eps);
        }

        f[1]=omega[1];
        for(int j=2;j<J;j++)
        {
            f[j]=omega[j]+(f[j-1]/e[j-1]);
        }

        std::complex<long double> wave_next[J+1];
        for(int j=J;j>=0;j--)
        {
            if(j==J)
            {
                wave_next[j]=std::complex<long double>(0,0);
            }
            else if(j==0)
            {
                wave_next[j]=std::complex<long double>(0,0);
            }
            else if(j==J-1)
            {
                wave_next[j]=-f[j-1]/e[j-1];
            }
            else
            {
                wave_next[j]=(wave_next[j+1]-f[j])/e[j];
            }
        }

        copy(wave_next,wave,J);

        //t=t+delta;
        c++;
    }

    return 0;
}