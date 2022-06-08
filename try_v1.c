//This is file is for trying the FFTW library and seeing the normization constant

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

// Macros for real and imaginary parts
#define REAL 0
#define IMAG 1




int nx=100;
double x[nx];
double y[nx];
fftw_complex kx[nx/2], sumkx[nx/2], yc[nx]; // This is equivalent to: double x[n][2]
double amp[nx];

main()
{

    for(int i=0; i<sizeof(x);i++)
    {
        y[i] = 3 * sin( x[i] );
    }

    fftw_plan plan;

    plan = fftw_plan_dft(nx, &y[0], &kx[0],FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_execute(plan);

    g = fopen("pky.dat", "w");
    for (j = 1; j < ny / 2; j++){
        
        fprintf(g, "%.6e  %.6e\n", 1/100*j*2*pi, kx[j]);
    }
    fclose(g);

}