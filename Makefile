CXX = g++
CPPFLAGS = -std=c++11



all: DGM_C

DGM_C: Discontinuous_Galerkin_Method.o legendregauss.o output_file.o mesh_generator.o
	${CXX} ${CPPFLAGS} Discontinuous_Galerkin_Method.o legendregauss.o output_file.o mesh_generator.o -o DGM_C
Discontinous_Galerkin_Method.o: Discontinuous_Galerkin_Method.cpp
	${CXX} ${CPPFLAGS} Discontinuous_Galerkin_Method.cpp -c -o Discontinuous_Galerkin_Method.o    
legendregauss.o: legendregauss.cpp
	${CXX} ${CPPFLAGS} legendregauss.cpp -c -o legendregauss.o
output_file.o: output_file.cpp
	${CXX} ${CPPFLAGS} output_file.cpp -c -o output_file.o
mesh_generator.o: mesh_generator.cpp
	${CXX} ${CPPFLAGS} mesh_generator.cpp -c -o mesh_generator.o
