#include<cstdio>
#include<cstdlib>
#include<cmath>

/*
Program to generate the equipotential lines for the given potential
*/

double U(double x, double y)
{
    return 0.5*(pow(x,2)+pow(y,2))+y*pow(x,2)-(1/3.0)*pow(y,3);
}

int main(int argc, char* argv[])
{
    double E=atof(argv[1]);
    double h=0.0001;
    double xmin=-0.9,xmax=0.9;
    double x=xmin;

    char filename[50];
    sprintf(filename,"file%0.5lf.dat",E);
    FILE *fp=fopen(filename,"w");

    while(x<=xmax)
    {
        double ymin=-0.6, ymax=1.0;
        double y=ymin;
        while(y<=ymax)
        {
            if(fabs(U(x,y)-E)<=0.0001)
            {
                fprintf(fp,"%lf %lf\n",x,y);
            }
            y=y+h;
        }
        x=x+h;
    }
    return 0;
}