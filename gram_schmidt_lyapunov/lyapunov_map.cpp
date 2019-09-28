#include<cstdio>
#include<cmath>
#include<cstdlib>

#define SIZE 2

const double AA=1.4;
const double BB=0.3;

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

/**********MAP AND ITS PERTUBATION EVOLUTION EQUATIONS*********/

double f0(double x, double y, double dx, double dy)
{
    return 1-AA*pow(x,2)+y;
}
double f1(double x, double y, double dx, double dy)
{
    return BB*x;
}
double f2(double x, double y, double dx, double dy)
{
    return -2*AA*x*dx+dy;
}
double f3(double x, double y, double dx, double dy)
{
    return BB*dx;
}
/**********MAP AND ITS PERTUBATION EVOLUTION EQUATIONS*********/

void update(double x[])
{
    double tx,ty,tdx,tdy;
    tx=f0(x[0],x[1],x[2],x[3]);
    ty=f1(x[0],x[1],x[2],x[3]);
    tdx=f2(x[0],x[1],x[2],x[3]);
    tdy=f3(x[0],x[1],x[2],x[3]);

    x[0]=tx; x[1]=ty; x[2]=tdx; x[3]=tdy;
}


int main(int argc, char* argv[])
{
    int N=2*SIZE;

    int T=(int)atof(argv[1]);

    int count=0;

    double x1[4]={0.3,0.3,1,0};

    double x2[4]={0.3,0.3,0,1};

    double tx,ty,tdx,tdy;


    double lya[SIZE]={0,0};
    double mod=0, inner1=0;
    
    while(count<=T)
    {

        update(x1);
        update(x2);

        my_vector d[SIZE];
        
        double t1[SIZE]={x1[2],x1[3]};
        d[0].init(t1);
           
        double t2[SIZE]={x2[2],x2[3]};
        d[1].init(t2);
           
    
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

        /******GRAM SCHMIDT ORTHONORMALIZATION*******/

        for(int i=0;i<SIZE;i++)
        {
            x1[i+SIZE]=u[0].VECTOR[i];
            x2[i+SIZE]=u[1].VECTOR[i];
        }
        count++;

    }
    printf("Lyapunov1=%lf Lyapunov2=%lf\n",lya[0]/(count), lya[1]/(count));
    return 0;
}