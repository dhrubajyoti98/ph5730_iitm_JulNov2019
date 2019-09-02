#include<omp.h>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "N_ODE.cpp"

using namespace std;

class Points
{
    public:
    double y,py,x,px;

    double distance2(Points p)
    {
        return pow(y-p.y,2)+pow(py-p.py,2)+pow(px-p.px,2)+pow(x-p.x,2);
    }
};


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

void solve_trajectory(double E, Points p[], gsl_rng* r, int STEPS, double y, double py, double x, double px)
{
    double traj[4]={x,px,y,py};

    N_ODE sys(4,0.01);

    double (*f[4])(double[],double);
    f[0]=&f0;
    f[1]=&f1;
    f[2]=&f2;
    f[3]=&f3;

    double WORKSPACE[4][4];

    for(int i=0;i<STEPS;i++)
    {
        p[i].x=traj[0];
        p[i].px=traj[1];
        p[i].y=traj[2];
        p[i].py=traj[3];  
         
        sys.solveRK4(traj,traj,f,WORKSPACE,1);
    }
}

double calculate_dn(double E, int STEPS, gsl_rng* r, double del)
{
    Points p1[STEPS];
    Points p2[STEPS];

    double y,py,x,px;

    do
    {
        y=gsl_ran_flat(r,-0.4,0.4);
        py=gsl_ran_flat(r,-0.4,0.4);
        x=gsl_ran_flat(r,-0.4,0.4);
        px=sqrt(2*(E-U(x,y))-pow(py,2));
    } while (isnan(px));

    omp_set_num_threads(2);
    
    #pragma omp parallel 
    {
        int ID=omp_get_thread_num();
        if(ID==0)
            solve_trajectory(E,p1,r,STEPS,y,py,x,px);
        else
            solve_trajectory(E,p2,r,STEPS,y+del,py,x,px);
    }

    return p1[STEPS-1].distance2(p2[STEPS-1]); 
}

int main(int argc, char *argv[])
{
    double EMIN=atof(argv[1]);
    double EMAX=atof(argv[2]);

    double d0=pow(10,-7),LYA;

    int STEPS=pow(10,3);
    int IC=pow(10,2);

    gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r,time(0));

    FILE *fp=fopen(argv[3],"w");

    double E=EMIN;

    while(E<=EMAX)
    {
        double dn=0;
        printf("E=%lf started.\n",E);

        for(int i=1;i<=IC;i++)
        {
            dn=dn+calculate_dn(E,STEPS,r,d0);
        }

        dn=dn/IC;

        LYA=(1.0/STEPS)*log(sqrt(dn)/d0);
        fprintf(fp,"%lf %lf\n",E,LYA);
        printf("E=%lf done.\n",E);
        E=E+atof(argv[4]);
    }
    
    return 0;
}
