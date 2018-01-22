/*--------------------------------------------------------------*\
|  Title: Shape of Mobius Band                                   |
|  Methods: Euler integration and Gradient descent               |
|  References: The Shape of a Mobius Band by L.Mahadevan         |
|              A Treatise on Mathematical theory of elasticity   |
|              by A.E.H Love                                     |
|  Author: Srikanth sarma                                        |
|          of house Stallion                                     |
|          First of the name                                     |
|          The king of the Cosmos and                            |
|          The Guardian of the Galaxy                            |
|          A sword against the Darkness and                      |
|          The sheild that guards the realms of men              |
|          A watchfull Protector and                             |
|          The Dark Knight                                       |
\*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double pi,ds;

void Solve(double**,double*,double*,int,double,double); //function takes initial values including guesses as input and computes the grid point values of all variables.
double absmax(double*);                         //Maximum absolte value in the array.
double mod(double);                            //absolute value.
void Shape(double,double,double*,double,int,double,double);

main()
{
    int N;
    double a,b,tolerance,FDiff,GDStep;
    a=0.5;
    b=0.5;
    tolerance=0.001;
    N=1000;
    FDiff=0.0001;
    GDStep=0.00001;

    double g[6];
    g[0]=-1;                    //initial guess values.
    g[1]=0;
    g[2]=0.5;
    g[3]=0;
    g[4]=0;
    g[5]=0.5;

    while(a<=1.5)
    {
        Shape(a,b,g,tolerance,N,GDStep,FDiff);
        a+=0.1;                     
    }
}

void Shape(double a,double b,double*g,double tol,int n,double c,double h)
{
    int i,j,k;
    double gt[6],dg[6],f[6],ft[6],J[6][6],grad[6];
    double *x[10];

    char S[20];
    sprintf(S,"MobiusBandGDfine%f.txt",a);
    FILE*fp;
    fp=fopen(S,"w");                              //file pointer to output final results.

    pi=3.141593;

    ds=2*pi/n;
    for(i=0;i<10;i++)
    {
        x[i]=(double*)malloc(n*sizeof(double));    //allocating memory space for x.
    }

    for(i=4;i<10;i++)                            //Known initial values
    {
        x[i][0]=0;
    }
    x[5][0]=pi/2;

    while(1)
    {
        Solve(x,g,f,n,a,b);                            //perform mesh computation.
        printf("%f\n",absmax(f));
        if(absmax(f)<tol)                     //break when absolute maxximum error is lower than the threshold.
        break;
        for(i=0;i<6;i++)
        {
            for(j=0;j<6;j++)                 //for numerical differentiation,we are using (f(x+h)-f(x))/h.
            {
                gt[j]=g[j];
            }
            gt[i]+=h;
            Solve(x,gt,ft,n,a,b);

            grad[i]=0;
            for(j=0;j<6;j++)
            {
                J[i][j]=(ft[j]-f[j])/h;        //collecting jacobian terms.
                grad[i]+=J[i][j]*f[j];
            }
        }
        for(i=0;i<6;i++)                       //Gradient descent: move against the gradient to reach the root.
        {
            g[i]-=c*grad[i];
        }
    }

    for(i=0;i<10;i++)                          //Output to a file.
    {
        for(j=0;j<n;j++)
        {
            fprintf(fp,"%f\t",x[i][j]);
        }
        fprintf(fp,"\n");
    }
    return;
}

void Solve(double**y,double*g,double*ef,int n,double a,double b)
{
    int i,j;
    double v[10],ac[2];                   //define 1st derivatives and 2nd derivatives.
                                         //In the following we can consider the analogy of 1st deivative as velocity and 2nd as acceleration

    for(i=0;i<4;i++)
    {
        y[i][0]=g[i];
    }
    v[0]=g[4];
    v[1]=g[5];

    for(i=1;i<n;i++)                     //Loop to cover all grid points.compute accelerations and veloities at present grid point by using the differential equations.
    {
        ac[0]=(-(1-b)*((a-b)*y[0][i-1]*y[1][i-1]*y[1][i-1]+y[2][i-1]*v[1])+(a-1)*y[0][i-1]*y[2][i-1]*y[2][i-1]+b*y[2][i-1]*v[1]+y[3][i-1]*y[0][i-1])/a;
        ac[1]=((1-a)*((a-b)*y[1][i-1]*y[0][i-1]*y[0][i-1]+y[2][i-1]*v[0])+(b-1)*y[1][i-1]*y[2][i-1]*y[2][i-1]-a*y[2][i-1]*v[0]+y[3][i-1]*y[1][i-1])/b;
        v[2]=(a-b)*y[0][i-1]*y[1][i-1];
        v[3]=-a*y[0][i-1]*v[0]-b*y[1][i-1]*v[1]-(a-b)*y[0][i-1]*y[1][i-1]*y[2][i-1];
        v[4]=(-y[0][i-1]*cos(y[6][i-1])+y[1][i-1]*sin(y[6][i-1]))/sin(y[5][i-1]);
        v[5]=y[0][i-1]*sin(y[6][i-1])+y[1][i-1]*cos(y[6][i-1]);
        v[6]=y[2][i-1]+(y[0][i-1]*cos(y[6][i-1])-y[1][i-1]*sin(y[6][i-1]))*cos(y[5][i-1])/sin(y[5][i-1]);
        v[7]=sin(y[5][i-1])*cos(y[4][i-1]);
        v[8]=sin(y[5][i-1])*sin(y[4][i-1]);
        v[9]=cos(y[5][i-1]);

        for(j=0;j<10;j++)               //moving to next grid point. 'ds' is like time step.
        {
            y[j][i]=y[j][i-1]+v[j]*ds;
            if(j<2)
            {
                y[j][i]+=0.5*ac[j]*ds*ds;
                v[j]+=ac[j]*ds;
            }
        }
    }

    for(i=0;i<6;i++)                  //calculate difference between the expected and calculated final grid value .
    {
        ef[i]=y[i+4][n-1];
    }
    ef[0]-=2*pi;
    ef[1]-=pi/2;
    ef[2]-=pi;
    return;
}

double absmax(double*y)
{
    double m=0;
    int i;

    for(i=0;i<6;i++)
    {
        if(mod(y[i])>m)
        m=mod(y[i]);
    }
    return m;
}

double mod(double y)
{
    if(y>0)
    return y;
    else
    return -y;
}
