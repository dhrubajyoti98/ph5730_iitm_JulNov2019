#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include "ising.cpp"

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include<omp.h>


int main(int argc, char* argv[])
{
    gsl_rng *h=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(h,time(0));
    int NTHREADS;

    int steps=600000;

    double n1= 1.0/(steps);

    double TIN=1.5, TFIN=3.5;
    
    #pragma omp parallel 
    {
        int ID=omp_get_thread_num();
        int NT=omp_get_num_threads();

        #pragma omp master
        {
            NTHREADS=NT;
        }

        
        
        char fname[10];
        sprintf(fname,"f%d.dat",ID);
        FILE *f1=fopen(fname,"w");

        double H=(TFIN-TIN)/NT;
        double tstart=TIN+H*ID;
        double tend=TIN+H*(ID+1);

        #pragma omp critical
        {
            printf("THREAD ID=%d IS CALCULATING FOR T=%lf TO T=%lf\n",ID,tstart,tend);
        }

        double t=tstart;
        
        while(t<=tend)
        {
            Ising model(1.0/t,h,1,0);
            
            double E=0,M1=0,M2=0,mag=0;

            int k=1;
            while(k<=steps)
            {
                model.getNextConfig();
                k++;
            }
            int v=1;
            while(v<=steps)
            {
                model.getNextConfig();
                E=E+model.calculateEnergy();
                mag=model.sum();
                M1=M1+mag/pow(N,2);
                M2=M2+(mag*mag)/pow(N,4);
                v++;
            }
            fprintf(f1,"%lf %lf %lf %lf\n",t,n1*E,n1*M1,(1.0/t)*(n1*M2-n1*n1*M1*M1));

            #pragma omp critical
            {
                printf("THREAD ID=%d , T=%lf FINISHED.\n",ID,t);
            }

            t=t+0.001;            
        }
        fclose(f1);
        
    }

    FILE *fp=fopen("plot.plt","w");
    fprintf(fp,"plot ");
    for(int i=0;i<NTHREADS;i++)
    {
        fprintf(fp," \"f%d.dat\" u 1:b w lp,",i);
    }
    fclose(fp);
    
    return 0;
}