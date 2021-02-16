#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H
struct integratorResults { 
    double* xij;
    double* Lgx;
};
void writeToFile(integratorResults Result, integratorResults Result2, integratorResults Result3, int N);
#endif