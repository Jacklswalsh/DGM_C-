#include "output_file.h"
#include "Discontinuous_Galerkin_Method.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <string>

// Compile with: g++ -std=c++11 Discontinuous_Galerkin_Method.cpp


class NDGM{

    private:
        double* Wb;

    public:
        int N;
        nodesAndWeightsResult GLnodes_weights;
        double c = 1.0;
        double* lj1;
        double* ljn1;
        double** Dhij;
        double** D;
        double* xij;

        NDGM(int N){
            // initialise solution & Legendre Gauss Nodes and Weights
            this->N = N;
            this->xij = new double [N]; // Solution Array
            GLnodes_weights = LegendreGaussNodesAndWeights(N);
            Wb = barycentricWeights(N, GLnodes_weights.x);
            this->lj1 = lagrangeInterpolatingPolynomials(double(1.0), GLnodes_weights.x, Wb);
            this->ljn1 = lagrangeInterpolatingPolynomials(double(-1.0), GLnodes_weights.x, Wb);
            // std::cout << this->lj1[3] << this->ljn1[3] << std::endl;
            
            this->D = polynomialDerivativeMatrix(GLnodes_weights.x, N);
            this->Dhij = AllocateMatrixMemory(N,N);
            for(int j=0; j<N; j++){
                // std::cout << GLnodes_weights.w[j] <<std::endl;
                for (int i=0; i<N; i++){
                    this->Dhij[i][j] = -this->D[j][i] * (GLnodes_weights.w[j]/GLnodes_weights.w[i]);
                    // std::cout << "i: " << i << "   j: " << j << "   Dhij: " << Dhij[i][j] <<std::endl;
                };
            };
        }
    
        double* MxVDerivative(int N, double** Dij, double* fj){
            double t = 0;
            double* Infd = new double[N]; 
            for(int i=0; i<N; i++){
                t = 0; 
                for(int j=0; j<N; j++){
                    t = t + Dij[i][j] * fj[j];
                }; 
                Infd[i] = t;
            };
            return Infd;
        }

        double* DGDerivative(double XiL, double XiR){
            double* xijd = MxVDerivative(this->N, this->Dhij, this->xij);
            for(int j=0; j<this->N; j++){
                xijd[j] = xijd[j] + ((XiR * this->lj1[j] - XiL * this->ljn1[j])/this->GLnodes_weights.w[j]);
            };
            // std::cout << xijd[3] << std::endl;
            return xijd;
        }

        double interpolateToBoundary(double* lj){
            double interpolatedValue = 0.0;
            for(int j=0; j<this->N; j++){
                interpolatedValue = interpolatedValue + lj[j] * this->xij[j];
            };
            std::cout << lj[3] << "     Xij: " << xij[3] << std::endl;
            return interpolatedValue;
        }

        double g(double t){
            return 0.0;
        }

        double* DGTimeDerivative(double t){
            double XiL;
            double XiR;
            if(this->c > 0){
                XiL = g(t); 
                XiR = interpolateToBoundary(this->lj1);
                // std::cout << XiL << "     XiR: " << this->lj1[3] << std::endl;
            }
            else{
                XiR = g(t);
                XiL = interpolateToBoundary(this->ljn1);
            }
            double* xijdt = DGDerivative(XiL, XiR);
            for(int i=0; i<this->N; i++){
                xijdt[i] = -(this->c) * xijdt[i]; // PP
            }
            return xijdt;
        }

        double** polynomialDerivativeMatrix(double* xj, int N){
            double** D = AllocateMatrixMemory(N, N);
            double* wj = barycentricWeights(N, xj);
            for(int i=0; i<N; i++){
                D[i][i] = 0.0;
                for(int j=0; j<N; j++){
                    if( i != j){
                        D[i][j] = (wj[j]/wj[i]) * 1.0/(xj[i] - xj[j]);
                        D[i][i] = D[i][i] - D[i][j];
                    };
                };
            };
            return D;
        }


        double* lagrangeInterpolatingPolynomials(double x, double* xj, double* Wb){
            int N = this->N;
            double* lj = new double[N];

            bool xMatchesNode = 0;
            for(int j=0; j<N; j++){
                lj[j] = 0;
                if(AlmostEqual(x, xj[j]) == 1){
                    lj[j] = 1;
                    xMatchesNode = 1;
                };
            };
            if(xMatchesNode == 1){
                return lj;
            };
            double t = 0; 
            double s = 0;

            for(int j=0; j<N; j++){
                t = Wb[j]/(x-xj[j]);
                lj[j] = t;
                s = s + t;
            };
            for(int j=0; j<N; j++){
                lj[j] = lj[j] / s;
            };
            return lj;
        }


        double* barycentricWeights(int N, double* xj){

            double* W = new double[N];
            for(int j=0; j<N; j++){
                W[j] = 1.0;};
            // Clever looping trick to include evrything but the xn-xn term
            for(int j=1; j<N; j++){
                for(int k=0; k<j; k++){
                    W[k] = W[k] * (xj[k] - xj[j]);
                    W[j] = W[j] * (xj[j] - xj[k]);
                };
            };
            for(int j=0; j<N; j++){
                W[j] = 1/W[j];
            };
            return W;
        }

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
    };      


void printMatrix(double** M, int N){
    for(int i=0;i<N;i++){
    std::cout << "[ ";
        for(int j=0; j<N; j++){
            std::cout << M[i][j] << ", ";
        };
    std::cout << " ]" << std::endl;
    };
}


NDGM DGStepByRK3(double tn, double dt, NDGM dg){
    double am [3] = {0.0, -5.0/9.0, -153.0/128.0};
    double bm [3] = {0.0, 1.0/3.0, 3.0/4.0};
    double gm [3] = {1.0/3.0, 15.0/16.0, 8.0/15.0};
    double Gj[dg.N];
    double* xijd;
    double t;

    for(int m=0; m<3; m++){
        t = tn + bm[m] * dt;
        xijd = dg.DGTimeDerivative(t);
        // std::cout << xijd[3] << std::endl;

        for(int j=0; j<dg.N; j++){
            Gj[j] = am[m] * Gj[j] + xijd[j];
            dg.xij[j] = dg.xij[j] + gm[m] * dt * Gj[j];
        };
    };
    return dg;
}

integratorResults Integrator(double Nt, double Nout, double T, NDGM dg){
    double dt = T/Nt;
    double tn = 0;
    std::cout << "\n---------------------- \n\n    NEW CALCULATION        \n\n----------------------\n";

    // Defining the Initial Function
    double sigma = 0.2;
    for(int i=0; i<dg.N; i++){
        dg.xij[i] = exp(-log(2)*pow((dg.GLnodes_weights.x[i]+0.5), 2)/(pow(sigma,2)));
    };

    for(int n=0; n<Nt; n++){
        // std::cout << dg.xij[3] << std::endl;
        dg = DGStepByRK3(tn, dt, dg);
        tn = (n+1) * dt;
    };

    integratorResults Result = {dg.xij, dg.GLnodes_weights.x};
    return Result;
}

int main(){

    int N = 30;
    NDGM DGM(N);
    NDGM DGM2(N);
    NDGM DGM3(N);
    double Nt = 1000;
    double Nout = N;
    double T=0;
    // double T2 = 1;

    integratorResults Result = Integrator(Nt, Nout, T, DGM);
    integratorResults Result2 = Integrator(Nt*4, Nout, double(0.5), DGM2);
    integratorResults Result3 = Integrator(Nt*4, Nout, double(1.0), DGM3);

    // for (int i=0; i<5; i++){
    // std::cout << "X:  " << Result.Lgx[i] << std::endl;
    // std::cout << "xij:  " << Result.xij[i] << std::endl;
    // };

    writeToFile(Result, Result2, Result3, N);

    // delete[] DGM.xij;
    // delete[] DGM.lj1;
    // delete[] DGM.ljn1;
    // FreeMatrixMemory(N, DGM.D);
    // FreeMatrixMemory(N, DGM.Dhij);

    // delete[] DGM2.xij;
    // delete[] DGM2.lj1;
    // delete[] DGM2.ljn1;
    // FreeMatrixMemory(N, DGM2.D);
    // FreeMatrixMemory(N, DGM2.Dhij);

    // delete[] DGM3.xij;
    // delete[] DGM3.lj1;
    // delete[] DGM3.ljn1;
    // FreeMatrixMemory(N, DGM3.D);
    // FreeMatrixMemory(N, DGM3.Dhij);

return 1;}



