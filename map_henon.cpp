#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<omp.h>

/*
Program to plot the x-y phase space of the Henon-like 
area-preserving 2D map.
*/

double f(double a, double x, double y)
{
    return x+a*(y-pow(y,3));
}

double g(double a, double x, double y)
{
    return y-a*(x-pow(x,3));
}

int main(int argc, char* argv[])
{
    double a=atof(argv[1]);
    double del=atof(argv[2]);
    int N=100000;
    int F=1;

    double XMIN=-1.0;
    double XMAX=1.0;

    double delX=del, delY=del;

    #pragma omp parallel
    {
        int ID=omp_get_thread_num();
        int NT=omp_get_num_threads();

        double H=(XMAX-XMIN)/NT;

        double X_START=XMIN+H*ID;
        double X_END=XMIN+H*(ID+1);

        double YMIN=-1.0;
        double YMAX=1.0;

        double x=X_START;
        while(x<=X_END)
        {
            double y=YMIN;
            while(y<=YMAX)
            {
                char filename[50];
                sprintf(filename,"data%d.dat",F);
                FILE *fp=fopen(filename,"w");

                double x0=x,y0=y;
                for(int i=0;i<=N;i++)
                {
                    x0=f(a,x0,y0);
                    y0=g(a,x0,y0);
                    fprintf(fp,"%lf %lf\n",x0,y0);
                }
                y=y+delY;
                fclose(fp);
                F++;
            }
            x=x+delX;
        }
        
    }

    FILE *plot=fopen("plot.plt","w");
    fprintf(plot,"set term pngcairo\n");
    fprintf(plot,"set output \"plot_%lf.png\"\n",a);
    fprintf(plot,"set xlabel \"x\"\n");
    fprintf(plot,"set ylabel \"y\"\n");
    fprintf(plot,"set xrange [-1:1]\n");
    fprintf(plot,"set yrange [-1:1]\n");
    fprintf(plot,"unset key\n");
    fprintf(plot,"plot ");
    for(int i=1;i<F;i++)
    {
        fprintf(plot,"\"data%d.dat\" u 1:2 w points ps 0.2, ",i);
    }
    fclose(plot);

    return 0;
}