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

/**
 * Plot spécifique antenne : courbe puissance simulée vs théorique.
 * On garde ce helper ici car il n'est utilisé que par cette simulation.
 */
static void gp_plot_power(FILE *gp, double *P_sim, double *P_theo,
                          int q, double dt, int q_start)
{
    if (!gp) return;

    fprintf(gp, "set xlabel 't (s)'; set ylabel 'P (V^2/m^2)'\n");
    fprintf(gp, "set title 'Puissance recue'\n");
    fprintf(gp, "set style data lines\n");
    fprintf(gp, "set xrange [%g:*]\n", dt * q_start);
    fprintf(gp, "set yrange [0:*]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '-' w lines lc rgb 'red' lw 2 t 'Simulee',"
                " '-' w lines lc rgb 'blue' lw 2 t 'Theorique',"
                " '-' w lines lc rgb 'green' lw 2 t 'Difference'\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, P_sim[i]);
    fprintf(gp, "e\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, P_theo[i]);
    fprintf(gp, "e\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, fabs(P_sim[i] - P_theo[i]));
    fprintf(gp, "e\n");

    fflush(gp);
}

int calcul_2D_antenne(int m, double ***E, double ***Bx, double ***By,
                      source src, int xr, int yr, int step_time,
                      double dt, double eps_0, double dx,
                      double (*eps_r_2D)(double, double, double),
                      double (*sigma_2D)(double, double, double),
                      double length, double lambda)
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
    if (!P_sim || !P_theo || !P_ref) {
        fprintf(stderr, "ERREUR allocation puissance\n");
        free(P_sim); free(P_theo); free(P_ref);
        return 1;
    }

    // Gnuplot 
    FILE *gp = gp_open();
    if (!gp) { free(P_sim); free(P_theo); free(P_ref); return 1; }

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
                (*E)[i][j] = 1.0/(1.0+k) * ((*E)[i][j] * (1.0-k) + 1.0/(2.0 * eps_r_2D(i*dx, j*dx, length))* ((*By)[i][j] - (*By)[i-1][j] - ((*Bx)[i][j] - (*Bx)[i][j-1]))
                );
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
        P_ref[q]  = signal_src * signal_src;
        P_sim[q]  = signal_rec * signal_rec;
        P_theo[q] = (q - delay >= 0)
                   ? P_ref[q - delay] * ratio_theo * ratio_theo
                   : 0.0;

        if (q % frame == 0 && q > delay)
            gp_plot_power(gp, P_sim, P_theo, q, dt, q_start);
    }

    free(P_sim);
    free(P_theo);
    free(P_ref);
    gp_close(gp);
    return 0;
}