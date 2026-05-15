#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fonctions.h"
#include "plot_2D.h"

// Alloue et initialise à zéro un tableau 2D, etourne NULL en cas d'échec.
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


//Calcule la norme du vecteur de Poynting en un point (qualitatif donc pas de facteur 1/µ)
static double **build_poynting(double **E, double **Bx, double **By, int m) {
    double **S = alloc_field(m);
    if (!S) return NULL;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            S[i][j] = fabs(E[i][j]) * sqrt(Bx[i][j] * Bx[i][j] + By[i][j] * By[i][j]);
    return S;
}

//Librer tableaux 2D
static void free_field(double **f, int m) {
    if (!f) return;
    for (int i = 0; i < m; i++) 
        free(f[i]);
    free(f);
}

int calcul_2D_diel(int m, double ***E, double ***Bx, double ***By,
                   source src, int step_time, double dt, double eps_0,
                   double w, double dx,
                   double (*eps_r_2D)(double, double, double),
                   double (*sigma_2D)(double, double, double),
                   double length, int mask_type, int plot_mode)
{
    double k = 0;
    int frame = 20;

    // Allocation des champs 
    *E  = alloc_field(m);
    *Bx = alloc_field(m);
    *By = alloc_field(m);
    if (!*E || !*Bx || !*By) {
        fprintf(stderr, "ERREUR allocation champs\n");
        return 1;
    }

    // Ouverture gnuplot : 1 ou 2 fenêtres selon le mode 
    FILE *gp1 = gp_open();
    if (!gp1) return 1;

    if (plot_mode == 0) {
        gp_setup_image(gp1, m, -2.0, 2.0, "Champ Ez");
    } else if (plot_mode == 1) {
        gp_setup_image(gp1, m, -0.75, 0.75, "S = |E| * sqrt(Bx²+By²)");
    }

    for (int q = 0; q < step_time; q++) {

        // Mise à jour E 
        for (int i = 1; i < m - 1; i++) {
            for (int j = 1; j < m - 1; j++) {
                k = sigma_2D(dx*i, dx*j, length) * dt / (2 * eps_0 * eps_r_2D(dx*i, dx*j, length));
                (*E)[i][j] = 1.0/(1.0+k) * ( (*E)[i][j] * (1.0-k) + 1.0/(2.0*eps_r_2D(i*dx, j*dx, length))* ((*By)[i][j] - (*By)[i-1][j] - ((*Bx)[i][j] - (*Bx)[i][j-1])));
                if (mask_type == 1)
                    (*E)[i][j] *= young_slits(dx*i, dx*j, length, dx);
                else if (mask_type == 2)
                    (*E)[i][j] *= single_slit(dx*i, dx*j, length, dx);
                else if (mask_type == 3)
                    (*E)[i][j] *= fabry_perrot_mask(dx*i, dx*j, length, dx);
            }
        }

        // Source 
        for (int j = src.y_start; j < src.y_end - 1; j++)
            (*E)[src.x][j] += src.A * src.forme(src.w * dt * q + src.k * j * dx);

        // Mise à jour B 
        for (int i = 0; i < m - 1; i++) {
            for (int j = 0; j < m - 1; j++) {
                (*Bx)[i][j] -= 0.5 * ((*E)[i][j+1] - (*E)[i][j]);
                (*By)[i][j] += 0.5 * ((*E)[i+1][j] - (*E)[i][j]);
            }
        }

        // Affichage live 
        if (q % frame == 0) { //On toutes les x frame
            if (plot_mode == 0) {
                gp_plot_field(gp1, *E, m);

            } else if (plot_mode == 1) {
                double **S = build_poynting(*E, *Bx, *By, m);
                if (S) { 
                    gp_plot_field(gp1, S, m); 
                    free_field(S, m); 
                }
            }
        }
    }

    // Plot final diffracté : intensité bord droit dans une 2e fenêtre 
    if (plot_mode == 0 && mask_type <= 2) {
        FILE *gp_intensity = gp_open();
        if (gp_intensity) {
            gp_setup_curve(gp_intensity, "Intensite bord droit", "y", "E^2");
            int edge = m - 5;
            double *intensity = malloc(sizeof(double) * m);
            if (intensity) {
                for (int j = 0; j < m; j++)
                    intensity[j] = (*E)[edge][j] * (*E)[edge][j];
                gp_plot_curve(gp_intensity, intensity, m, "red");
                free(intensity);
            }
            gp_close(gp_intensity);
        }
    }

    gp_close(gp1);
    return 0;
}