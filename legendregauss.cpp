#include "legendregauss.h"
#include <iostream>
#include <cmath>
#include <cfloat>


nodesAndWeightsResult LegendreGaussNodesAndWeights(int N){
    double* xi = new double [N];
    double* w = new double [N];

    if (N == 0){
        xi[0] = 0.0;
        w[0] = 2.0;
    }
    else if (N == 1){
        xi[0] = -sqrt(1.0/3.0);
        w[0] = 1.0;
        xi[1] = -xi[0];
        w[1] = w[0];
    }
    else{
        double h = 0.001;
        int i = 0;
        double root; 

        for(double x=-1.0;x<=1.0;x=x+h){
            //set the accuracy to approx. 10^-15 but there is also a limit on maxSteps. (Modify these acc. to your needs)
            root=Bisection(N+1,x,x+h,0.00000000001,100000); 
            if(root!=999){ 
                xi[i]=root;
                i++;
            }
        };

        for(int i=0;i<N+1;i++){
            w[i] = (2/(double(1.0 - pow(xi[i], 2)) * pow(legendreFunction(N+1,xi[i]).LdashN,2)));
            // std::cout << w[i] << std::endl;
        };
    }

    nodesAndWeightsResult Result = {xi, w};
    return Result;
}


double Bisection(int N ,double a,double b, double tol, double max_iter){
double c;
if((legendreFunction(N,a).Ln*legendreFunction(N,b).Ln)<=0){  
    int iter=1;
    /*Bisection Method begins that tabulates the various values at each iteration*/
    do{
    c=(a+b)/2;
    if(legendreFunction(N,a).Ln*legendreFunction(N,c).Ln>0){
    a=c;
    }
    else if(legendreFunction(N,a).Ln*legendreFunction(N,c).Ln<0){
    b=c;
    }
    else if(legendreFunction(N,c).Ln==0){
        return c;
    }
    iter++;
        
    }while(abs(a-b)>=tol&&iter<=max_iter);
    return c;
}
else{
    return 999;
}
}

legendreResults legendreFunction(int N, double x){
    double LN;
    double LdashN;
    if( N == 0 ){
        LN= 0;
        LdashN = 1;
        legendreResults Result = {LN, LdashN};
        return Result;
    }
    else if( N == 1 ){
        LN = x;
        LdashN = 1;
        legendreResults Result = {LN, LdashN};
        return Result;
    }
    else{
        double LN2 = 1;
        double LN1 = x;
        double LdashN2 = 0;
        double LdashN1 = 1; 

        double i = 2;
        for(int k = 2; k < N; k++){
            LN = (((2 * i - 1) / i) * x * LN1) - (((i - 1) / k) * LN2);
            LdashN = LdashN2 + (((2 * k) - 1) * LN1);

            LN2 = LN1;
            LN1 = LN; 

            LdashN2 = LdashN1;
            LdashN1 = LdashN;
            i++;
        };
    };
    legendreResults Result = {LN, LdashN};
    return Result;
}

int AlmostEqual(double a, double b){
    int AlmEql;
    if(a == 0 || b == 0){
        if(abs(a-b) <= 2 * DBL_EPSILON){ AlmEql = 1; }
        else{AlmEql = 0;};
        }
    else{
        if( (abs(a-b) <= abs(a) * DBL_EPSILON) && (abs(a-b) <= abs(b) * DBL_EPSILON) ){ AlmEql = 1; }
        else{AlmEql = 0; };
        }
    return AlmEql;
}
    