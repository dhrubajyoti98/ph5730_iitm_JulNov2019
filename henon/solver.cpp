#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<omp.h>
#include "N_ODE.cpp"

double f0(double x[], double t)
{
    return x[1];
}

double f1(double x[], double t)
{
    return -x[0]-2*x[0]*x[2];
}

double f2(double x[], double t)
{
    return x[3];
}

double f3(double x[], double t)
{
    return -x[2]-(pow(x[0],2)-pow(x[2],2));
}

double U(double x, double y)
{
    return 0.5*(pow(x,2)+pow(y,2))+y*pow(x,2)-(1/3.0)*pow(y,3);
}

int main(int argc, char *argv[])
{
    srand48(time(NULL));
    int F=1;
    char filename[50];
    double E=atof(argv[1]);
    double TE=0;


    double YMIN=-0.5;
    double YMAX=1.1;

    double Y=YMIN;

    while(Y<=YMAX)
    {
        double PY=-0.6;
        while(PY<=0.7)
        {
            
            int N=4;
            double H=0.001,T=2000;

            N_ODE sys(N,H);  
            sprintf(filename,"data%d.dat",F);
            FILE *fp=fopen(filename,"w");
            
            double X=drand48();
            double PX=sqrt(2*(E-U(X,Y))-pow(PY,2));
    

            double WORKSPACE[N][4];

            double x[4]={X,PX,Y,PY};

            double (*f[4])(double[],double);
            f[0]=&f0;
            f[1]=&f1;
            f[2]=&f2;
            f[3]=&f3;

            while(sys.t<=T)
            {
                if(x[1]>0 && fabs(x[0])<0.0001)
                {
                    TE=U(x[0],x[2])+0.5*(pow(x[1],2)+pow(x[3],2));
                    fprintf(fp,"%lf %lf %lf %lf\n",x[2],x[3],sys.t,TE);
                }
                sys.solveRK4(x,x,f,WORKSPACE,N);
            }
            fclose(fp);

            PY=PY+0.1;
            F=F+1;
        }
        Y=Y+0.1;
    }
    

    FILE *plot=fopen("plot.plt","w");
    fprintf(plot,"set term pngcairo\n");
    fprintf(plot,"set output \"plot%lf.png\"\n",E);
    fprintf(plot,"set xlabel \"y\"");
    fprintf(plot,"set ylabel \"y-dot\"");
    fprintf(plot,"unset key\n");
    fprintf(plot,"plot ");
    for(int i=1;i<F;i++)
    {
        fprintf(plot,"\"data%d.dat\" u 1:2 w points ps 0.2, ",i);
    }
    fclose(plot);

    printf("Calculation Done for E=%lf\n",E);


    return 0;

}