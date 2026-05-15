/*Fonction qui initialise la grille en 1D
Elle prend en paramètre :
- la taille du domaine m et crée un tableau de taille 2m
- la source : une valeur true correspond à une sinusoïdale, false à une gaussienne
-le pointeur vers la grille*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

static FILE* plot_init_1D_libre(int m) {
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

static void plot_frame_1D_libre(FILE *pipe, double *E, int m, int q) {
    fprintf(pipe, "set title '1D Free space - step %d'\n", q);
    fprintf(pipe, "plot '-' with lines lc rgb 'blue' title 'E'\n");
    for (int i = 0; i < m; i++) {
        fprintf(pipe, "%d %e\n", i, E[i]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
}

int calcul_1D_libre(int m, double **E, double **B, double (*source)(double), int step_time, double dt, double eps_0, double w, double A) {

    *E = malloc(sizeof(double) * m);
    *B = malloc(sizeof(double) * m);
    double *E_bgauche = malloc(sizeof(double) * step_time);
    double *E_bdroit = malloc(sizeof(double) * step_time);
    if (*B == NULL || *E == NULL || E_bdroit == NULL || E_bgauche == NULL) {
        printf("ERREUR lors de l'allocation des tableaux des champs\n");
        return 1;
    }

    FILE *pipe = plot_init_1D_libre(m);
    if(pipe == NULL) {
        printf("Erreur lors de l'ouverture du pipe\n");
        return 1;
    }

    for (int i = 0; i < m; i++) {
        (*E)[i] = 0;
        (*B)[i] = 0;
    }

    for(int i = 0; i < step_time; i++) {
        E_bgauche[i] = (*E)[2];
        if (i >= 2) {
            (*E)[1] = E_bgauche[i-2];
            (*E)[m-1] = E_bdroit[i-2];
        }

        (*E)[1] += 0.5*((*B)[1] - (*B)[0]) - dt/eps_0*A*(*source)(w*i*dt);  // On place la source à côté du bord gauche

        //boucle spatiale pour E
        for(int j = 2; j < m-1; j++) {
            (*E)[j] += 0.5*((*B)[j] - (*B)[j-1]);
        }

        E_bdroit[i] = (*E)[m-2];

        //boucle spatiale pour B
        for(int j = 0; j < m-1; j++) {
            (*B)[j] += 0.5*((*E)[j+1] - (*E)[j]);
        }

        if (pipe && i % 5 == 0)      // affiche 1 frame sur 5
            plot_frame_1D_libre(pipe, *E, m, i);
    }

    pclose(pipe);
    return 0;
}