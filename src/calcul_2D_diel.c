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


int calcul_2D_diel(int m, double ***E, double ***Bx, double ***By, double ***E_phas,
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
    *E_phas = alloc_field(m);
    double **E_re = alloc_field(m);
    double **E_im = alloc_field(m);
    if (!*E || !*Bx || !*By || !*E_phas) {
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

    int q_regime = round(3*step_time / 4);
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

        if (q > q_regime) {
            double phase = src.w * dt * q;
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    E_re[i][j] += (*E)[i][j] * cos(phase);
                    E_im[i][j] += (*E)[i][j] * sin(phase);
                }
            }
        }
         
        // Mise à jour B 
        for (int i = 0; i < m - 1; i++) {
            for (int j = 0; j < m - 1; j++) {
                (*Bx)[i][j] -= 0.5 * ((*E)[i][j+1] - (*E)[i][j]);
                (*By)[i][j] += 0.5 * ((*E)[i+1][j] - (*E)[i][j]);
            }
        }

        // Affichage live 
        /*if (q % frame == 0) { //On toutes les x frame
            if (plot_mode == 0) {
                gp_plot_field(gp1, *E, m);

            } else if (plot_mode == 1) {
                double **S = build_poynting(*E, *Bx, *By, m);
                if (S) { 
                    gp_plot_field(gp1, S, m); 
                    free_field(S, m); 
                }
            }
        }*/
    }

    int n_accum = step_time - q_regime;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            (*E_phas)[i][j] = sqrt(E_re[i][j]*E_re[i][j] + E_im[i][j]*E_im[i][j]) / n_accum;

    // Plot final diffracté : intensité bord droit dans une 2e fenêtre 
    //Généré par IA
    if (plot_mode == 0 && mask_type == 2) {
 
        double c_loc = 1.0 / sqrt(eps_0 * MU_0);  // si c n'est pas dans le scope
        double slit_width  = dx*40;         // largeur fente (cf single_slit)
        double slit_center = length / 2.0;          // centre de la fente en y
        double wall_x      = length / 2.0;          // position du mur en x
 
        // Ligne d'observation : à 3/4 du domaine en x
        // (assez loin du mur pour être en champ lointain,
        //  assez loin du bord droit pour éviter les réflexions)
        int obs_col = (int)( m*0.9);
        double dist = obs_col * dx - wall_x;        // distance mur → observation

        double lambda_sim = PI * c_loc / src.w;  // au lieu de 2*PI*c_loc/src.w
        printf("=== DIAGNOSTIC ===\n");
        printf("src.w   = %g rad/step\n", src.w);
        printf("dx      = %g m\n", dx);
        printf("dt      = %g s\n", dt);
        printf("lambda  = %g m\n", lambda_sim);
        printf("slit_width = %g m\n", slit_width);
        printf("a/lambda   = %g\n", slit_width / lambda_sim);
        printf("dist       = %g m\n", dist);
        printf("Fresnel number N = a²/(lambda*dist) = %g\n",
        slit_width*slit_width / (lambda_sim * dist));
 
        double *E_sim  = malloc(sizeof(double) * m);
        double *E_theo = malloc(sizeof(double) * m);
        if (E_sim && E_theo) {
 
            // 1) Champ simulé : enveloppe (phaseur) sur la colonne d'observation
            double max_sim = 0.0;
            for (int j = 0; j < m; j++) {
                E_sim[j] = (*E_phas)[obs_col][j];   // <-- E_phas, pas E
                if (E_sim[j] > max_sim) max_sim = E_sim[j];
            }
 
            // 2) Théorie : sinc évalué au même point
            //    Pour chaque point y, on calcule theta_s = atan2(y - y_centre, dist)
            //    puis |sinc(pi * a / lambda * sin(theta_s))|
            double max_theo = 0.0;
            for (int j = 0; j < m; j++) {
                double y  = j * dx;
                double dy = y - slit_center;
                
                // 1. Calcul de la vraie distance fente -> point (x_obs, y)
                double rho = sqrt(dist * dist + dy * dy); 
                
                double theta_s = atan2(dy, dist);
                double arg = PI * slit_width / lambda_sim * sin(theta_s);
                double sinc_val = (fabs(arg) < 1e-12) ? 1.0 : sin(arg) / arg;
                
                // 2. Facteur d'obliquité de l'optique physique
                double obliquity = cos(theta_s); 
                
                // 3. Le champ décroît en cos(theta)/sqrt(rho)
                E_theo[j] = fabs(sinc_val);
                
                if (E_theo[j] > max_theo) max_theo = E_theo[j];
            }
 
            // 3) Normalisation
            if (max_sim > 0.0)
                for (int j = 0; j < m; j++) E_sim[j] /= max_sim;
            if(max_sim > 0.0) 
                for (int j = 0; j < m; j++) E_theo[j] /= max_theo;
 
            // 4) Plot gnuplot superposé
            FILE *gp_diff = gp_open();
            if (gp_diff) {
                fprintf(gp_diff, "set title 'Diffraction simple fente — FDTD vs sinc'\n");
                fprintf(gp_diff, "set xlabel 'y (m)'\n");
                fprintf(gp_diff, "set ylabel '|E| normalise'\n");
                fprintf(gp_diff, "set grid\n");
                fprintf(gp_diff,
                    "plot '-' with lines lc rgb 'red'  lw 2 title 'FDTD (phaseur)', "
                    "     '-' with lines lc rgb 'blue' lw 2 title 'sinc theorique'\n");
 
                for (int j = 0; j < m; j++)
                    fprintf(gp_diff, "%g %g\n", j * dx, E_sim[j]);
                fprintf(gp_diff, "e\n");
 
                for (int j = 0; j < m; j++)
                    fprintf(gp_diff, "%g %g\n", j * dx, E_theo[j]);
                fprintf(gp_diff, "e\n");
 
                fflush(gp_diff);
                gp_close(gp_diff);
            }
 
            free(E_sim);
            free(E_theo);
        }
    }
    // Pour les autres mask_type (Young etc.), on garde l'ancien plot
    else if (plot_mode == 0 && mask_type <= 1) {
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