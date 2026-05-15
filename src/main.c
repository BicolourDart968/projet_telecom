#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "calcul_1D_diel.h"
#include "calcul_1D_libre.h"
#include "plot_1D.h"
#include "fonctions.h"
#include "plot_2D.h"
#include "calcul_2D_diel.h"
#include "calcul_2D_antenne.h"
#include "simulations_1D.h"
#include "estimate_refraction_angle.h"
#include "debeye.h"

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
    time = length / c;

    m = round(length / dx);
    step_time = round(time / dt);

    int choix = 0;
    char buf[32];

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

            default:
                printf("Choix invalide\n");
        }
    }

    return 0;
}