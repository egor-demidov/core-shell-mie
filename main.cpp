//
// Created by ed242 on 7/21/2022.
//

#include <iostream>
#include <string>
#include <cmath>

#include "mie_calc.h"

using namespace std;

// r1, r2, lambda, rm1, im1, rm2, im2
// Output: extinction, scattering
int main(int argc, char ** argv) {
    if (argc != 8) {
        cerr << "ERROR: exactly 7 floating-point arguments must be provided\n";
        cerr << "(1) radius of inner sphere\n";
        cerr << "(2) radius of outer sphere\n";
        cerr << "(3) incident light wavelength\n";
        cerr << "(4) real part of the refractive index for the inner sphere\n";
        cerr << "(5) imaginary part of the refractive index for the inner sphere\n";
        cerr << "(6) real part of the refractive index for the inner sphere\n";
        cerr << "(7) imaginary part of the refractive index for the inner sphere\n";
        return EXIT_FAILURE;
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

    cout << "C_ext\tC_sca\tC_abs\n";
    cout << k_coef[0]*M_PI*r2*r2 << "\t" << k_coef[1]*M_PI*r2*r2 << "\t" << k_coef[2]*M_PI*r2*r2 << endl;

    return 0;
}