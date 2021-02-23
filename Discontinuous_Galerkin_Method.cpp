#include "output_file.h"
#include "legendregauss.h"
#include "Discontinuous_Galerkin_Method.h"
#include "mesh_generator.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <string>

// Compile with: g++ -std=c++11 Discontinuous_Galerkin_Method.cpp
// g++ -std=c++11 Discontinuous_Galerkin_Method.cpp -c -o Discontinuous_Galerkin_Method.o
// g++ -std=c++11 legendregauss.cpp -c -o legendregauss.o
// g++ -std=c++11 output_file.cpp -c -o output_file.o
// g++ -std=c++11 Discontinuous_Galerkin_Method.o legendregauss.o output_file.o -o DGM_C
// ./DGM_C   ---> outout 
// python3 CPlotter.py
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

    MeshGeneration Mesh;
    Mesh.mesh_generation(4, -8.0, 8.0);
    Mesh.print_mesh(0);

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



