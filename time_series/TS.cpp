#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include "predict.cpp"

int main(int argc, char* argv[])
{
    srand48(time(NULL));
    Prediction p;
    p.init(1000000,35,2,2);
    p.populate(drand48(),drand48());
    p.workout();

    for(int i=0;i<35;i++)
    {
        printf("%d %lf %lf %lf\n",i,p.TEST[i],p.PREDICT[i],p.ERROR[i]);
    }

    return 0;
}