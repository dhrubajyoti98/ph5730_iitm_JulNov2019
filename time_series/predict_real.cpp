#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<ctime>

#include<omp.h>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>
#include<gsl/gsl_fit.h>

using namespace std;

class Prediction
{
    private:
    vector<double> DATA;
    int N, PL, K, M;

    public:
    vector<double> TEST, PREDICT, ERROR;

    void init(int n, int pl, int k, int m)
    {
        N=n;
        PL=pl;
        K=k;
        M=m;

        DATA.resize(N);
        PREDICT.resize(PL);
        TEST.resize(PL);
        ERROR.resize(PL);
    }

    void load_data_from_file(FILE *fp, vector<double>& datafile)
    {
        char line[1024];
        double x;
        while(fgets(line,sizeof(line),fp)!=NULL)
	    {
		    sscanf(line,"%lf",&x);
		    datafile.push_back(x);
	    }
    }

    void populate() //Populate the DATA;
    {
        vector<double> datafile;
        FILE *fp=fopen("temp.data","r");
        load_data_from_file(fp,datafile);

        for(int i=0;i<N+PL;i++)
        {
            if(i<N)
            {
                DATA[i]=datafile[i];
            }
            else
            {
                TEST[i-N]=datafile[i];  
            }
        }
    }

    double getSD()
    {
        double mean=0,sd=0;
        for(int i=0;i<N;i++)
        {
            mean=mean+DATA[i];
        }
        mean=mean/N;
        for(int i=0;i<N;i++)
        {
            sd=sd+pow(DATA[i]-mean,2);
        }
        return sqrt(sd/N);
    }

    double getDistance(int i, int j)
    {
        double dist=0;
        for(int k=0;k<M;k++)
        {
            dist=dist+pow(DATA[i-k]-DATA[j-k],2);
        }
        return sqrt(dist);
    }

    gsl_vector* convertToGSL_Vector(vector<double> v, int size)
    {
        gsl_vector *r=gsl_vector_calloc(size);
        for(int i=0;i<size;i++)
        {
            gsl_vector_set(r,i,v[i]);
        }
        return r;
    }

    double predict(int T)
    {
        vector<double> distance;
        vector<double> index;

        double base=DATA[N-1];
        for(int j=N-T-1;j>=M-1;j--)
        {
            distance.push_back(getDistance(N-1,j));
            index.push_back(j);
        }
        int size=distance.size();

        gsl_vector *dist=convertToGSL_Vector(distance,size);
        gsl_vector *ind=convertToGSL_Vector(index,size);

        gsl_sort_vector2(dist,ind);

        double t[K];
        double tT[K];

        for(int j=0;j<K;j++)
        {
            t[j]=DATA[(int)gsl_vector_get(ind,j)];
            tT[j]=DATA[(int)(gsl_vector_get(ind,j))+T];
        }

        if(K==1)
        {
            return tT[0];
        }

        int temp;
        double val=0, val_err=0;
        double c0,c1,cov00,cov01,cov11;
        double sumsq;


        temp=gsl_fit_linear(t,1,tT,1,K,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
        temp=gsl_fit_linear_est(base, c0, c1, cov00, cov01, cov11, &val, &val_err);

        return val;
    }

    void workout()
    {
        for(int T=1;T<=PL;T++)
        {
            PREDICT[T-1]=predict(T);
            ERROR[T-1]=pow(PREDICT[T-1]-TEST[T-1],2);
        }
    }
};

int main(int argc, char* argv[])
{
    Prediction p;
    int N=(int)atof(argv[1]);
    int PL=(int)atof(argv[2]);
    int K=(int)atof(argv[3]);
    int M=(int)atof(argv[4]);
    p.init(N,PL,K,M);
    p.populate();
    p.workout();
    for(int i=0;i<PL;i++)
    {
        printf("%d %lf %lf\n",i,p.TEST[i],p.PREDICT[i]);
    }
    return 0;
}