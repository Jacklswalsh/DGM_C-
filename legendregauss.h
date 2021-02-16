#if !defined(LEGENDREGAUSS_H)
#define LEGENDREGAUSS_H

struct legendreResults { 
    double Ln; 
    double LdashN; 
};

struct nodesAndWeightsResult { 
    double* x; 
    double* w; 
};

#endif // LEGENDREGAUSS_H

nodesAndWeightsResult LegendreGaussNodesAndWeights(int N);
double Bisection(int N ,double a,double b, double tol, double max_iter);
legendreResults legendreFunction(int N, double x);
int AlmostEqual(double a, double b);