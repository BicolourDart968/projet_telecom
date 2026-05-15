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


int calcul_1D_libre(int m, double **E, double **B, double (*source)(double), int step_time, double dt, double eps_0, double w) {

    double A = 1.0;       //amplitude de la source

    *E = malloc(sizeof(double) * m);
    *B = malloc(sizeof(double) * m);
    if (*B == NULL || *E == NULL) {
        printf("ERREUR lors de l'allocation des tableaux des champs\n");
        return 1;
    }

    for (int i = 0; i < m; i++) {
        (*E)[i] = 0;
        (*B)[i] = 0;
    }

    for(int i = 0; i < step_time; i++) {
        (*E)[1] += 0.5*((*B)[1] - (*B)[0]) - dt/eps_0*(*source)(w*i*dt);  // On place la source à côté du bord gauche

        //boucle spatiale pour E
        for(int j = 2; j < m-1; j++) {
            (*E)[j] += 0.5*((*B)[j] - (*B)[j-1]);
        }

        //boucle spatiale pour B
        for(int j = 0; j < m-1; j++) {
            (*B)[j] += 0.5*((*E)[j+1] - (*E)[j]);
        }
    }

    return 0;
}