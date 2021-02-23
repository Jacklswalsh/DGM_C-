#ifndef MESHGENERATOR_H
#define MESHGENERATOR_H
#include <iostream>
#include <list>
#include <vector>
#include <iterator>

class MeshGeneration{

    private:
        int elements;
        double start;
        double end;
        double dx;

    public:
        std::list< std::vector<double> > mesh;
        MeshGeneration();
        void mesh_generation(int elements, double start, double end);
        void print_mesh(int i);

};

#endif