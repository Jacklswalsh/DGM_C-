// #include <cstdio>
#include "mesh_generator.h"
#include <iostream>
#include <list>
#include <vector>
#include <iterator>

    MeshGeneration::MeshGeneration(){
        this->mesh = mesh;
    };
    
    void MeshGeneration::mesh_generation(int elements, double start, double end){
        dx = (end - start)/(elements);
        for(double i = start; i < end; i += dx){
            std::vector<double> x = {i, i+dx};
            this->mesh.push_back(x);
            // std::cout << x[1] << std::endl;
        };
    };

    void MeshGeneration::print_mesh(int i){
        for (auto v : this->mesh) std::cout << v[i] << "\n"; 
        return;
    };

// int main(){

//     // MeshGeneration Mesh;
//     // Mesh.mesh_generation(4, -8.0, 8.0);
//     // Mesh.print_mesh(0);

//     return 1;
// };