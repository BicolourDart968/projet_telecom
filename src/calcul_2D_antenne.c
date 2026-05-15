#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fonctions.h"
#include "plot_2D.h"

#define PI 3.14159265359
#define EPS_0 8.854e-12
#define MU_0 12.566e-7

// Réutilise la même allocation que calcul_2D_diel
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

int calcul_2D_antenne(int m, double ***E, double ***Bx, double ***By,
                      source src, int xr, int yr, int step_time,
                      double dt, double eps_0, double dx,
                      double (*eps_r_2D)(double, double, double),
                      double (*sigma_2D)(double, double, double),
                      double length, double lambda, double he, double Ra)
{
    double k = 0;
    int frame = 50;

    // Allocation 
    *E  = alloc_field(m);
    *Bx = alloc_field(m);
    *By = alloc_field(m);
    if (!*E || !*Bx || !*By) {
        fprintf(stderr, "ERREUR allocation champs\n");
        return 1;
    }
    double *P_sim  = calloc(step_time, sizeof(double));
    double *P_theo = calloc(step_time, sizeof(double));
    double *P_ref  = calloc(step_time, sizeof(double));
    double *Voc    = calloc(step_time, sizeof(double));
    double *Voc_theo = calloc(step_time, sizeof(double));
    double *Ez_ref = calloc(step_time, sizeof(double));
    if (!P_sim || !P_theo || !P_ref || !Voc || !Voc_theo || !Ez_ref) {
        fprintf(stderr, "ERREUR allocation puissance ou Voc\n");
        free(P_sim); free(P_theo); free(P_ref); free(Voc);
        return 1;
    }

    // Gnuplot 
    FILE *gp1 = gp_open();
    FILE *gp2 = gp_open();
    if (!gp1 || !gp2) {
        if (gp1) gp_close(gp1);
        if (gp2) gp_close(gp2);
        free(P_sim); free(P_theo); free(P_ref); free(Voc);
        return 1;
    }

    /* Paramètres de référence */
    double d = sqrt(pow((xr - src.x)*dx, 2)
                  + pow((yr - (src.y_start + src.y_end)/2)*dx, 2));
    int d0_cells = (int)(2 * lambda / dx);
    int x_ref = src.x + d0_cells;
    int y_ref = (src.y_start + src.y_end) / 2;
    double d0 = d0_cells * dx;
    int delay = (int)((d - d0) / (1.0/sqrt(EPS_0*MU_0) * dt));
    int q_start = (int)((d0 / (1.0/sqrt(MU_0*EPS_0) * dt)) * 1.2);
    double ratio_theo = sqrt(d0 / d);

    // Boucle temporelle 
    for (int q = 0; q < step_time; q++) {

        for (int i = 1; i < m - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                k = sigma_2D(dx*i, dx*j, length) * dt
                    / (2 * eps_0 * eps_r_2D(dx*i, dx*j, length));
                (*E)[i][j] = 1.0/(1.0+k) * ((*E)[i][j] * (1.0-k) + 1.0/(2.0 * eps_r_2D(i*dx, j*dx, length))* ((*By)[i][j] - (*By)[i-1][j] - ((*Bx)[i][j] - (*Bx)[i][j-1])));
            }
        }

        for (int j = src.y_start; j < src.y_end - 1; j++)
            (*E)[src.x][j] += src.A * src.forme(src.w * dt * q);

        for (int i = 0; i < m - 1; i++) {
            for (int j = 0; j < m - 1; j++) {
                (*Bx)[i][j] -= 0.5 * ((*E)[i][j+1] - (*E)[i][j]);
                (*By)[i][j] += 0.5 * ((*E)[i+1][j] - (*E)[i][j]);
            }
        }

        double signal_src = (*E)[x_ref][y_ref];
        double signal_rec = (*E)[xr][yr];
        P_ref[q]  = he * he * signal_src * signal_src / (8.0*Ra);
        P_sim[q]  = he * he * signal_rec * signal_rec / (8.0*Ra);
        P_theo[q] = (q - delay >= 0)
                   ? P_ref[q - delay] * ratio_theo * ratio_theo
                   : 0.0;

        // Calcul de la tension reçue par l'antenne
        Ez_ref[q] = signal_src;  // = (*E)[x_ref][y_ref]

        Voc[q] = he * fabs(signal_rec);
        Voc_theo[q] = (q - delay >= 0)
                    ? he * fabs(Ez_ref[q - delay]) * ratio_theo
                    : 0.0;
        if (q % frame == 0 && q > delay) {
            gp_plot_power(gp1, P_sim, P_theo, q, dt, q_start);
            gp_plot_voc(gp2, Voc, Voc_theo, q, dt, q_start);
        }
    }

    free(P_sim);
    free(P_theo);
    free(P_ref);
    free(Voc);
    gp_close(gp1);
    gp_close(gp2);
    return 0;
}