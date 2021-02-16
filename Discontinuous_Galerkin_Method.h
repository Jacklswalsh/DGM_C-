#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H
struct integratorResults { 
    double* xij;
    double* Lgx;
};
void writeToFile(integratorResults Result, integratorResults Result2, integratorResults Result3, int N);
#endif



// Function Templates
nodesAndWeightsResult LegendreGaussNodesAndWeights(int N);
double Bisection(int N ,double a,double b, double tol, double max_iter);
legendreResults legendreFunction(int N, double x);
int AlmostEqual(double a, double b);

// Matrix Memory Allocations
double** AllocateMatrixMemory(int rows, int cols){
    double** A; 
    A = new double* [rows];
    for(int i = 0; i<rows; ++i){
        A[i] = new double [cols];
    } 
    return A;
}

void FreeMatrixMemory(int rows, double** Matrix){
    for(int i=0; i < rows; i++){
        delete[] Matrix[i];
        }
    delete[] Matrix;
    return;
}