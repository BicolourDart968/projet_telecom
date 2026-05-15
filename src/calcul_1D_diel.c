#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define MU_0 12.566e-7

int calcul_1D_diel(int m, double **E, double **B, double (*source)(double), int step_time, double length,
double dt, double eps_0, double w, double (*eps_r)(double, double), double (*sigma)(double, double), double dx) {

    double k;

    *E = malloc(sizeof(double) * m);
    *B = malloc(sizeof(double) * m);
    double *E_bgauche = malloc(sizeof(double) * step_time);
    double *E_bdroit = malloc(sizeof(double) * step_time);
    if (*B == NULL || *E == NULL || E_bdroit == NULL || E_bgauche == NULL) {
        printf("ERREUR lors de l'allocation des tableaux des champs\n");
        return 1;
    }

    for (int i = 0; i < m; i++) {
        (*E)[i] = 0;
        (*B)[i] = 0;
    }

    for(int q = 0; q < step_time; q++) {

        double esp_r_val = eps_r(dx, length);
        //Coeffs utilises pour les calculs  
        k = sigma(dx, length) * dt / (2 * eps_0 * esp_r_val ); 

        (*E)[(int)round(m/2)] +=  dt/(esp_r_val*eps_0)*source(w*q*dt);

        E_bgauche[q] = (*E)[2];
        if (q >= 2) {
            (*E)[1] = E_bgauche[q-2];
            (*E)[m-1] = E_bdroit[q-2];
        }

        for(int j = 1; j < m; j++) {
            double eps_r_val = eps_r(dx*j, length);
            k = sigma(dx*j, length) * dt / (2.0 * eps_0 * eps_r_val);
            double cE = 1 / (2 * eps_r_val);
            (*E)[j] = 1/(1+k) * ((*E)[j]*(1-k) + cE*((*B)[j] - (*B)[j-1]));
        }

        E_bdroit[q] = (*E)[m-1];

        for(int j = 0; j < m-1; j++) {
            (*B)[j] += 0.5*((*E)[j+1] - (*E)[j]);
        }
    }
    return 0;
}