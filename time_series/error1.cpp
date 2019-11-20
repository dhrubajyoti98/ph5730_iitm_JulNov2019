#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "predict.cpp"

int main(int argc, char* argv[])
{
    srand48(time(NULL));

    int SAMPLES=200; // 100 sample time-series.
    int PREDICT_LENGTH=15; //Predict length = 15.

    Prediction p[SAMPLES]; 
    
    int N=pow(10,(int)atof(argv[1]));
    int k=(int)atof(argv[2]);
    int m=(int)atof(argv[3]);
    int decide=(int)atof(argv[4]); // decide=0 is logistic, decide=1 is henon.

    #pragma omp parallel for
    for(int i=0;i<SAMPLES;i++)
    {
        p[i].init(N,PREDICT_LENGTH,k,m); 
        if(decide==0)
        {
            p[i].populate(drand48()); //logistic
        }
        else
        {
            p[i].populate(drand48()/10.0,drand48()/10.0); //henon
        }
        p[i].workout();
    }

    double SD=p[0].getSD();

    for(int T=1;T<=PREDICT_LENGTH;T++)
    {
        double value=0;
        for(int j=0;j<SAMPLES;j++)
        {
            value=value+p[j].ERROR[T-1];
        }
        value=sqrt(value/SAMPLES);
        printf("%d %lf\n",T,value/SD);
    }

    return 0;
}