#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "plot_2D.h"

// Alloue et initialise à zéro un tableau 2D, retourne NULL en cas d'échec.
static double **alloc_field(int m) {
    double **f = malloc(sizeof(double *) * m);
    if (!f) return NULL;

    for (int i = 0; i < m; i++) {
        f[i] = calloc(m, sizeof(double));
        if (!f[i]) {
            for (int k = 0; k < i; k++) free(f[k]);
            free(f);
            return NULL;
        }
    }
    return f;
}

int calcul_2D_libre(int m, double ***E, double ***Bx, double ***By, double ***E_phas, double (*source)(double), int step_time, double dt, double eps_0, double w, double dx, double (*mur_conduct)(double, double, double), double length) {

    double A = 1.0;       //amplitude de la source
    int frame = 20;

    /* --- Allocation des champs --- */
    *E  = alloc_field(m);
    *Bx = alloc_field(m);
    *By = alloc_field(m);
    *E_phas = alloc_field(m);
    if (!*E || !*Bx || !*By || !*E_phas) {
        printf("ERREUR allocation champs\n");
        return 1;
    }

    /* --- Ouverture gnuplot --- */
    FILE *gp = gp_open();
    if (!gp) return 1;
    gp_setup_image(gp, m, -0.02, 0.02, "Champ Ez");

    int q_regime = round(step_time / 2);

    //Dans les boucles suivantes, on ne parcourt pas le bord afin de ne pas toucher aux CL

    for(int q = 0; q < step_time; q++) {
        for(int i = 1; i < m-1; i++) {
            for(int j = 1; j < m-1; j++) 
                (*E)[i][j] += 0.5*((*By)[i][j] - (*By)[i-1][j] - ((*Bx)[i][j] - (*Bx)[i][j-1]));
        }

        (*E)[(int)round(m/2)][(int)round(m/2)] -= dt/eps_0*A*source(w*dt*q);  //On impose la source au milieu

        if (q > q_regime) {
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    double value = fabs((*E)[i][j]);
                    if (value > (*E_phas)[i][j])
                        (*E_phas)[i][j] = value;
                }
            }
        }

        for(int i = 0; i < m-1; i++) {
            for(int j = 0; j < m-1; j++) {
                (*Bx)[i][j] -= 0.5*((*E)[i][j+1]-(*E)[i][j]);
                (*By)[i][j] += 0.5*((*E)[i+1][j]-(*E)[i][j]);
            }
        }

        if(q % frame == 0) {
            gp_plot_field(gp, *E, m);
        }
    }

    gp_close(gp);

    return 0;
} 