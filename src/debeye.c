#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

static FILE* plot_init_debye(int m) {
    FILE *pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
        printf("ERREUR: impossible d'ouvrir gnuplot\n");
        return NULL;
    }
    fprintf(pipe, "set xrange [0:%d]\n", m - 1);
    fprintf(pipe, "set yrange [-1.5:1.5]\n");
    fprintf(pipe, "set xlabel 'Position'\n");
    fprintf(pipe, "set ylabel 'E'\n");
    return pipe;
}

static void plot_frame_debye(FILE *pipe, double *E, int m, int q) {
    fprintf(pipe, "set title 'Debye - step %d'\n", q);
    fprintf(pipe, "plot '-' with lines lc rgb 'blue' title 'E'\n");
    for (int i = 0; i < m; i++) {
        fprintf(pipe, "%d %e\n", i, E[i]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
}

int debeye(int m, double **E, double **B, double **E_phas, int step_time, double dt, double eps_0, double w) {

    double A = 50.0;       //amplitude de la source

    double *P = malloc(sizeof(double)*m);   //Tableau de polarisation
    double *E_bgauche = malloc(sizeof(double) * step_time);
    double *E_bdroit = malloc(sizeof(double) * step_time);

    *E = malloc(sizeof(double) * m);
    *B = malloc(sizeof(double) * m);
    *E_phas = malloc(sizeof(double) * m);
    if (*B == NULL || *E == NULL || *E_phas == NULL || P == NULL || E_bdroit == NULL || E_bgauche == NULL) {
        printf("ERREUR lors de l'allocation des tableaux des champs\n");
        return 1;
    }

    FILE *pipe = plot_init_debye(m);
    if(pipe == NULL) {
        printf("Erreur lors de l'ouverture du pipe\n");
        return 1;
    }

    for (int i = 0; i < m; i++) {
        (*E)[i] = 0;
        (*B)[i] = 0;
        (*E_phas)[i] = 0.0;
        P[i] = 0.0;
    }

    int q_regime = round(step_time / 2);

    // Parametres Debye
    double tau = 1.0 / w;
    double chi_s = 1.25;       // eps_r_statique = 2.25 (verre), chi = eps_r - 1

    // Zone dielectrique : tiers central du domaine
    int diel_start = 0;
    int diel_end   = m;

    // Source : pulse gaussien court (broadband)
    int i_src = m / 8;
    double sigma_t = 1.0 * dt;
    double t0 = 1 * sigma_t;

    for(int q = 0; q < step_time; q++) {
        double t = q*dt;

        E_bgauche[q] = (*E)[2];
        if (q >= 2) {
            (*E)[1] = E_bgauche[q-2];
            (*E)[m-1] = E_bdroit[q-2];
        }

        //boucle spatiale pour E
        for(int i = 2; i < m-1; i++) {
            if(i >= diel_start && i < diel_end) {
                double dP = 1/tau*(chi_s*eps_0*(*E)[i] - P[i]);
                (*E)[i] += 0.5*((*B)[i] - (*B)[i-1]) - dt/eps_0*dP;
                P[i] += dt*dP;
            }

            else
                (*E)[i] += 0.5*((*B)[i] - (*B)[i-1]);
        }

        E_bdroit[q] = (*E)[m-1];
        if (q > q_regime) {
            for (int i = 0; i < m; i++) {
                double value = fabs((*E)[i]);
                if (value > (*E_phas)[i])
                    (*E_phas)[i] = value;
            }
        }

        double env = exp(-(t - t0) * (t - t0) / (2.0 * sigma_t * sigma_t));
        (*E)[i_src] += A * env * sin(w*t);

        //boucle spatiale pour B
        for(int i = 0; i < m-1; i++) {
            (*B)[i] += 0.5*((*E)[i+1] - (*E)[i]);
        }
        if (pipe && q % 5 == 0)      // affiche 1 frame sur 5
            plot_frame_debye(pipe, *E, m, q);
    }

    pclose(pipe);
    return 0;
}