#include<cstdio>
#include<cmath>
#include<cstdlib>
#include "Rossler_System.h"
#include "N_ODE.cpp"

#define SIZE 3

class my_vector
{
    public:
    double VECTOR[SIZE];

    my_vector()
    {

    }

    void init(double a[])
    {
        for(int i=0;i<SIZE;i++)
        {
            this->VECTOR[i]=a[i];
        }
    }

    void init(my_vector a)
    {
        for(int i=0;i<SIZE;i++)
        {
            this->VECTOR[i]=a.VECTOR[i];
        }
    }

    void get_vector(double a[])
    {
        for(int i=0;i<SIZE;i++)
        {
            a[i]=this->VECTOR[i];
        }
    }

    double get(int i)
    {
        return this->VECTOR[i];
    }

    double inner_product(my_vector b)
    {
        double val=0;
        for(int i=0;i<SIZE;i++)
        {
            val=val+(b.VECTOR[i])*(this->VECTOR[i]);
        }
        return val;
    }

    double modulus()
    {
        double val=0;
        for(int i=0;i<SIZE;i++)
        {
            val=val+pow(VECTOR[i],2);
        }
        return sqrt(val);
    }

    my_vector subtract(my_vector b)
    {
        double h[3];
        for(int i=0;i<SIZE;i++)
        {
            h[i]=this->VECTOR[i]-b.VECTOR[i];
        }
        my_vector temp;

        temp.init(h);

        return temp;
    }

    my_vector add(my_vector b)
    {
        double h[3];
        for(int i=0;i<SIZE;i++)
        {
            h[i]=this->VECTOR[i]+b.VECTOR[i];
        }
        my_vector temp;

        temp.init(h);

        return temp;
    }

    my_vector divide(double value)
    {
        my_vector temp;
        for(int i=0;i<SIZE;i++)
        {
            temp.VECTOR[i]=this->VECTOR[i]/value;
        }
        return temp;
    }

    my_vector multiply(double value)
    {
        my_vector temp;
        for(int i=0;i<SIZE;i++)
        {
            temp.VECTOR[i]=this->VECTOR[i]*value;
        }
        return temp;
    }

};



int main(int argc, char* argv[])
{
    int N=2*SIZE;
    double H=atof(argv[1]);
    int T=(int)atof(argv[2]);

    int count=0;

    FILE *fp=fopen("data.dat","w");

    N_ODE sys(N,H);

    double x1[6]={1.0,2.0,3.0,1,0,0};
    double x2[6]={1.0,2.0,3.0,0,1,0};
    double x3[6]={1.0,2.0,3.0,0,0,1};

    double WORKSPACE1[N][4];
    double WORKSPACE2[N][4];
    double WORKSPACE3[N][4];

    double (*f[6])(double[],double);
    f[0]=&f0;
    f[1]=&f1;
    f[2]=&f2;
    f[3]=&f3;
    f[4]=&f4;
    f[5]=&f5;

    double lya[3]={0,0,0};
    double mod=0, inner1=0, inner2=0;
    
    while(sys.t<T)
    {
        fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",sys.t,x2[0],x2[1],x2[2],x2[3],x2[4],x2[5]);

        my_vector d[SIZE];

        #pragma omp parallel
        {
            int NT=omp_get_thread_num();
            if(NT==1)
            {
                sys.solveRK4(x1,x1,f,WORKSPACE1,1);
                double t1[3]={x1[3],x1[4],x1[5]};
                d[0].init(t1);
            }
            else if(NT==2)
            {
                sys.solveRK4(x2,x2,f,WORKSPACE2,1);
                double t2[3]={x2[3],x2[4],x2[5]};
                d[1].init(t2);
            }
            else if(NT==3)
            {
                sys.solveRK4(x3,x3,f,WORKSPACE3,1);
                double t3[3]={x3[3],x3[4],x3[5]};
                d[2].init(t3);
            }
        }

        my_vector v[SIZE], u[SIZE];

        /******GRAM SCHMIDT ORTHONORMALIZATION*******/

        v[0].init(d[0]);
        mod=v[0].modulus();
        lya[0]=lya[0]+log2(mod);
        u[0].init(v[0].divide(mod));

        inner1=d[1].inner_product(u[0]);
        v[1].init(d[1].subtract(u[0].multiply(inner1)));
        mod=v[1].modulus();
        lya[1]=lya[1]+log2(mod);
        u[1].init(v[1].divide(mod));

        inner1=d[2].inner_product(u[0]);
        inner2=d[2].inner_product(u[1]);
        my_vector temp1, temp2, temp3;
        temp1.init(u[0].multiply(inner1));
        temp2.init(u[1].multiply(inner2));

        temp3.init(temp1.add(temp2));

        v[2].init(d[2].subtract(temp3));
        mod=v[2].modulus();
        lya[2]=lya[2]+log2(mod);
        u[2].init(v[2].divide(mod));

        /******GRAM SCHMIDT ORTHONORMALIZATION*******/

        for(int i=0;i<SIZE;i++)
        {
            x1[i+SIZE]=u[0].VECTOR[i];
            x2[i+SIZE]=u[1].VECTOR[i];
            x3[i+SIZE]=u[2].VECTOR[i];
        }
        count++;

    }
    printf("Lyapunov1=%lf Lyapunov2=%lf Lyapunov3=%lf\n",lya[0]/(count*H), lya[1]/(count*H), lya[2]/(count*H));
    fclose(fp);

    return 0;
}