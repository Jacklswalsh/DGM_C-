#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include "output_file.h"


void writeToFile(integratorResults Result, integratorResults Result2, integratorResults Result3, int Nout){
    std::ofstream myfile;
    std::string filename;
    std::cout<<"Please enter the filename: ___________.csv:    ";
    std::cin >> filename;
    myfile.open(filename + ".csv", std::ofstream::out | std::ofstream::trunc);
    myfile.close();
    myfile.open(filename + ".csv");
    for(int i=0; i<Nout; i++){
        myfile << Result.Lgx[i] << "," << Result.xij[i] << "," << Result2.xij[i] << "," << Result3.xij[i] << std::endl;
    };
    myfile.close();
    return;
}