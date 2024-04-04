//
// Created by ed242 on 7/21/2022.
//

#include <iostream>
#include <string>
#include <cmath>

#include "mie_calc.h"

#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

using namespace std;

// r1, r2, lambda, rm1, im1, rm2, im2
// Output: extinction, scattering
int main(int argc, char ** argv) {
    if (argc <= 7) {
        cerr << "Not enough arguments provided" << endl;
        return -1;
    }

    double r1, r2, lambda, alpha1, alpha2, rm1, im1, rm2, im2;
    double k_coef[3];

    r1 = stod(argv[1]);
    r2 = stod(argv[2]);
    lambda = stod(argv[3]);
    rm1 = stod(argv[4]);
    im1 = stod(argv[5]);
    rm2 = stod(argv[6]);
    im2 = stod(argv[7]);

    alpha1 = 2.0*M_PI*r1/lambda;
    alpha2 = 2.0*M_PI*r2/lambda;
    coated_mie(alpha1, rm1, im1, alpha2, rm2, im2, k_coef);
    cout << k_coef[0]*M_PI*r2*r2 << " " << k_coef[1]*M_PI*r2*r2 << endl;

    return 0;
}