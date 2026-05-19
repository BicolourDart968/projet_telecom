#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "calcul_1D_diel.h"
#include "calcul_1D_libre.h"
#include "plot_1D.h"
#include "fonctions.h"
#include "plot_2D.h"
#include "calcul_2D_libre.h"
#include "calcul_2D_diel.h"
#include "calcul_2D_antenne.h"
#include "simulations_1D.h"
#include "estimate_refraction_angle.h"
#include "debeye.h"
#include "emission.h"

#define PI 3.14159265359
#define EPS_0 8.854e-12
#define MU_0 12.566e-7



int main() {
    double **E, **Bx, **By, **E_phas_2D, *E_1D, *B_1D, *E_phas_1D, dx, dt, lambda, c, time, length;
    int m, step_time;

    c = 1 / sqrt(MU_0 * EPS_0);
    lambda = 0.01;
    length = 20 * lambda;
    dx = lambda / 20;
    dt = dx / (2 * c);
    time = length / (2*c);

    m = round(length / dx);
    step_time = round(time / dt);

    int choix = 0;
    char buf[32];
    char message[256];
    printf("%d\n", step_time);
    while (1) {
        printf("\n1. Fentes de Young\n"
               "2. Poynting - fentes de Young\n"
               "3. Poynting - prisme\n"
               "4. Prisme\n"
               "5. Antenne\n"
               "6. Decroissance exponentielle\n"
               "7. Eau de mer - profondeur de peau\n"
               "8. 1D espace libre - source sinusoidale\n"
               "9. 1D espace libre - source gaussienne\n"
               "10. Guide d\'onde\n"
               "11. 1D libre - paquet d\'onde gaussien modulate\u00e9\n"
               "12. Simple fente\n"
               "13. Lentille\n"
               "14. Fabry Perot\n"
               "15. Envoi d'informations\n"
               "16. Phaseur 2D\n"
               "0. Quitter\n> ");
        fgets(buf, sizeof(buf), stdin);

        if (sscanf(buf, "%d", &choix) != 1) {
            printf("Entree invalide\n");
            continue;
        }

        switch (choix) {
            case 0:
                return 0;

            case 1: {
                //Young
                source src = {.A = 1.0, .y_start = 1, .y_end = m - 1,
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_0, sigma_0, length, 1, 0))
                    printf("Erreur simulation\n");
                break;
            }

            case 2: {
                //Poynting Young
                source src = {.A = 1.0, .y_start = 1, .y_end = m - 1,
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_0, sigma_0, length, 1, 1))
                    printf("Erreur simulation\n");
                break;
            }

            case 3: {
                //Poynting refraction
                source src = {.A = 1.0, .y_start = 1, .y_end = m - 1,
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_prisme, sigma_0, length, 0, 1))
                    printf("Erreur simulation\n");
                else
                    estimate_refraction_angle(E, m, length, dx);
                break;
            }

            case 4: {
                //Refraction
                source src = {.A = 1.0, .y_start = 1, .y_end = m - 1,
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_prisme, sigma_0, length, 0, 0))
                    printf("Erreur simulation\n");
                else
                    estimate_refraction_angle(E, m, length, dx);
                break;
            }

            case 5: {
                //Antenne
                source src = {.A = 5.0,
                              .y_start = round((m - 5) / 2),
                              .y_end   = round((m + 5) / 2),
                              .x = round(m / 2), .forme = sin,
                              .w = 2 * PI * c / lambda};
                double he = lambda / PI;  // hauteur equivalente pour dipole lambda/2
                double Ra = 73.0;         // resistance d'antenne pour dipole lambda/2
                if (calcul_2D_antenne(m, &E, &Bx, &By, &E_phas_2D, src,
                                      round(m / 4), round(m / 4),
                                      step_time, dt, EPS_0, dx,
                                      eps_r_prisme, sigma_0, length, lambda, he, Ra))
                    printf("Erreur simulation\n");
                break;
            }

            case 6: {
                //Decroissance exp
                double w = 2 * PI * c / lambda;
                if (calcul_1D_diel(m, &E_1D, &B_1D, &E_phas_1D, sin, step_time,
                                   length, dt, EPS_0, w,
                                   eps_r_milieu, sigma_milieu, dx, 1.0))
                    printf("Erreur simulation 1D\n");
                else
                    plot_1D(E_1D, m, dx);
                break;
            }

            case 7: {
                //Profondeur de peau
                if (run_skin_depth_validation(c))
                    printf("Erreur simulation eau de mer\n");
                break;
            }

            case 8: {
                //1D libre sin
                double w = 2 * PI * c / lambda;
                if (calcul_1D_libre(m, &E_1D, &B_1D, &E_phas_1D, sin, step_time,
                                   dt, EPS_0, w, 10.0))
                    printf("Erreur simulation 1D libre sinus\n");
                plot_1D(E_phas_1D, m, dx);
                break;
            }

            case 9: {
                //1D libre gauss
                double w = 2 * PI * c / lambda;
                if (calcul_1D_libre(m, &E_1D, &B_1D, &E_phas_1D, gaussienne, step_time,
                                    dt, EPS_0, w, 10.0))
                    printf("Erreur simulation 1D libre gaussienne\n");
                break;
            }

            case 10: {
                //Guide d'onde
                double theta = 45.0 * PI / 180.0;
                source src = {.A = 1.0, .y_start = round(m*0.4), .y_end = round(m*0.6),
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda,
                              .k = (PI / lambda) * sin(theta)};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_guide, sigma_0, length, 0, 0))
                    printf("Erreur simulation\n");
                break;
            }

            case 11: {
                // Prisme dispersif Debye
                double w = PI * c /lambda;
                if(debeye(m, &E_1D, &B_1D, &E_phas_1D, step_time, dt, EPS_0, w))
                    printf("Erreur simulation prisme Debye\n");
                
                break;
            }

            case 12: {
                // Simple fente
                source src = {.A = 1.0, .y_start = 1, .y_end = m - 1,
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_0, sigma_0, length, 2, 0))
                    printf("Erreur simulation\n");
                break;
            }

            case 13 :{
                source src = {.A = 1.0, .y_start = round(m*0.3), .y_end = round(m*0.7),
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_lentille, sigma_0, length, 0, 0))
                    printf("Erreur simulation\n");
                break;
            }

            case 14 :{
                source src = {.A = 1.0, .y_start = 1, .y_end = m-1,
                              .x = round(m / 8), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_diel(m, &E, &Bx, &By, &E_phas_2D, src, step_time,
                   dt, EPS_0, src.w, dx,
                   eps_r_0, sigma_0, length, 3, 0))
                    printf("Erreur simulation\n");
                break;
            }

            case 15: {
                printf("Entrez le message a envoyer : ");
                if (!fgets(message, sizeof(message), stdin)) {
                    printf("Erreur lecture message\n");
                    break;
                }
                size_t len = strlen(message);
                if (len > 0 && message[len - 1] == '\n')
                    message[len - 1] = '\0';

                //Adaptation des paramètres pour la simulation
                double case_length = length*10;
                double case_time = case_length / c;
                int case_m = round(case_length / dx);
                int case_step_time = round(case_time / dt);
                int xr = round(3.0 * case_m / 5.0);
                int yr = round(3.0 * case_m / 5.0);

                if (emission(case_m, &E, &Bx, &By, xr, yr, case_step_time,
                             dt, EPS_0, dx, case_length, lambda, message))
                    printf("Erreur simulation emission\n");
                break;
            }

            case 16 :{
               source src = {.A = 1.0, .y_start = round(m/2), .y_end = round(m/2),
                              .x = round(m / 2), .forme = sin,
                              .w = PI * c / lambda};
                if (calcul_2D_libre(m, &E, &Bx, &By, &E_phas_2D, sin, step_time,
                                   dt, EPS_0, src.w, dx,
                                   eps_r_0, length)) {
                    printf("Erreur simulation\n");
                } else {
                    FILE *gp2d = gp_open();
                    if (gp2d) {
                        double maxval = 0.0;
                        for (int i = 0; i < m; i++) {
                            for (int j = 0; j < m; j++) {
                                if (E_phas_2D[i][j] > maxval)
                                    maxval = E_phas_2D[i][j];
                            }
                        }
                        gp_setup_image(gp2d, m, 0.0, maxval > 0 ? maxval : 1e-6,
                                       "Phaseur 2D en espace libre");
                        gp_plot_field(gp2d, E_phas_2D, m);
                        gp_close(gp2d);
                    }

                    int center = round(m / 2);
                    int margin = 3;
                    int start = center + margin;
                    int n = m - start;
                    if (n > 0) {
                        double *r = malloc(sizeof(double) * n);
                        double *Eline = malloc(sizeof(double) * n);
                        double *Etheo = malloc(sizeof(double) * n);
                        if (r && Eline && Etheo) {
                            for (int i = 0; i < n; i++) {
                                int ix = start + i;
                                r[i] = dx * (ix - center);
                                Eline[i] = E_phas_2D[ix][center];
                            }
                            
                            // Calcul de la constante A basé sur le premier point valide
                            // Cela évite que les réflexions aux frontières (fin du domaine) ne faussent A
                            double A_const = Eline[0] * sqrt(r[0]);
                            
                            for (int i = 0; i < n; i++) {
                                Etheo[i] = A_const / sqrt(r[i]);
                            }

                            FILE *gpcurve = gp_open();
                            if (gpcurve) {
                                gp_setup_curve(gpcurve,
                                               "Comparaison E_phas centre vs A/sqrt(r)",
                                               "distance r (m)", "E_phas");
                                
                                // Injection de la valeur de A_const directement dans le titre de la légende Gnuplot
                                fprintf(gpcurve,
                                        "plot '-' with lines lw 2 lc rgb 'blue' title 'E_phas centre',"
                                        " '-' with lines lw 2 lc rgb 'red' title '%.4f / sqrt(r)'\n", A_const);
                                
                                for (int i = 0; i < n; i++)
                                    fprintf(gpcurve, "%g %g\n", r[i], Eline[i]);
                                fprintf(gpcurve, "e\n");
                                for (int i = 0; i < n; i++)
                                    fprintf(gpcurve, "%g %g\n", r[i], Etheo[i]);
                                fprintf(gpcurve, "e\n");
                                fflush(gpcurve);
                                gp_close(gpcurve);
                            }
                        }
                        free(r);
                        free(Eline);
                        free(Etheo);
                    }
                }
                break;
            }
            default:
                printf("Choix invalide\n");
        }
    }

    return 0;
}