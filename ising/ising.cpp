#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<ctime>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include<omp.h>

#define N 64 /*N x N lattice*/

class Ising
{
    private:
    int lattice[N][N];
    int a,b;
    double BETA;
    gsl_rng *r;
    int START_TYPE; //START_TYPE=0 is hot-start, START_TYPE=1 is cold start.
    int UPDATE_TYPE; //UPDATE_TYPE=0 is Serial Site Update, UPDATE_TYPE=1 is Random Site Update

    void update_ab()
    {
        if(UPDATE_TYPE==0)
        {
            b=b+1;
            if(b==N)
            {
                b=0;
                a=a+1;
            }
            if(a==N)
            {
                a=0;
            }
        }
        else if(UPDATE_TYPE==1)
        {
            a=(int)(gsl_ran_flat(r,0,N));
            b=(int)(gsl_ran_flat(r,0,N));
        }
    }

    public:

    double calculateEnergy()
    {
        double energy=0;
        #pragma omp parallel for reduction(+:energy)
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int s1=lattice[i][j];
                int s2=lattice[i][(j+1)%N]+lattice[(i+1)%N][j];
                energy=energy+s1*s2;
            }
        }
        return (-1)*energy;
    }

    Ising(double B, gsl_rng *h, int start, int update)
    {
        r=h;
        BETA=B;

        a=0;
        b=0;

        START_TYPE=start;
        UPDATE_TYPE=update;

        double rnd;
        #pragma omp parallel for
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                if(START_TYPE==0) // HOT-START
                {
                    rnd=gsl_rng_uniform(r);
                    if(rnd<0.5)
                        lattice[i][j]=1;
                    else
                        lattice[i][j]=-1;
                }
                else if(START_TYPE==1) // COLD-START
                {
                    lattice[i][j]=1;
                }
            }
        }
    }

    void getNextConfig()
    {
        double energy_now=calculateEnergy(), energy_next, del_energy;
        do
        {
            update_ab();
            lattice[a][b]=(-1)*lattice[a][b];
            energy_next=calculateEnergy();
            del_energy=energy_next-energy_now;
            if(del_energy<0)
                break;
            else if(gsl_rng_uniform(r)<exp(-BETA*del_energy))
                break;
            else
            {
                lattice[a][b]=(-1)*lattice[a][b];
            }
        } while (1);
    }

    double sum()
    {
        double sum=0;
        #pragma omp parallel for reduction(+:sum)
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
                sum=sum+lattice[i][j];
        }
        return sum;
    }

    void printLattice()
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
                printf("%d  ",lattice[i][j]);
            printf("\n");
        }
        printf("\n");
    }
};