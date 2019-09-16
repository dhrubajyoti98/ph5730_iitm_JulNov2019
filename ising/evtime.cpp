#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include "ising.cpp"

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>


int main(int argc, char *argv[])
{
    double t=atof(argv[1]);
    gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r,time(0));

    Ising model1(1.0/t,r,1,0);
    Ising model2(1.0/t,r,0,0);

    int k=1;
    while(k<=500000)
    {
        #pragma omp parallel
        {
            int ID=omp_get_thread_num();

            if(ID==1)
            {
                model1.getNextConfig();
            }
            if(ID==2)
            {
                model2.getNextConfig();
            }
        }

        printf("%d %lf %lf\n",k,model1.calculateEnergy(),model2.calculateEnergy());
        k++;
    }
    return 0;
}