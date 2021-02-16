#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H
struct integratorResults { 
    double* xij;
    double* Lgx;
};
void writeToFile(integratorResults Result, integratorResults Result2, integratorResults Result3, int N);
#endif


struct legendreResults { 
    double Ln; 
    double LdashN; 
};

struct nodesAndWeightsResult { 
    double* x; 
    double* w; 
};

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