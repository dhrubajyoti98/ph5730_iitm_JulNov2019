#include<cstdio>
#include<cmath>
#include<cstdlib>
#include "predict.cpp"

int main(int argc, char *argv[])
{
    int SAMPLES=100;
    Prediction p[SAMPLES];
    int k=4;
    int m=2;

    int N=2;
    while(N<=5)
    {
        #pragma omp parallel for
        for(int i=0;i<SAMPLES;i++)
        {
            p[i].init((int)pow(10,N),30,k,m);
            p[i].populate(drand48()/10.0,drand48()/10.0);
            p[i].workout();
        }
       
        double SD=p[0].getSD();

        double value=0;
        for(int j=0;j<SAMPLES;j++)
        {
            value=value+p[j].ERROR[2];
        }
        value=sqrt(value/SAMPLES);
        printf("%d %lf\n",(int)pow(10,N),value/SD);

        N=N+1;
    }
    return 0;
}