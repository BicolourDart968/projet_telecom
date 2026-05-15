//Simualtion generee par IA

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "calcul_1D_diel.h"
#include "fonctions.h"
#include "simulations_1D.h"

#define PI 3.14159265359
#define EPS_0 8.854e-12
#define MU_0 12.566e-7

// Proprietes de l'eau de mer
static const double EAU_MER_EPS_R = 80.0;
static const double EAU_MER_SIGMA = 4.0;


// Calcule l'attenuation theorique alpha dans un milieu conducteur-dielectrique.
static double compute_alpha(double omega, double eps_r, double sigma) {
    double x = sigma / (omega * EPS_0 * eps_r);
    double sqrt_term = sqrt(1.0 + x * x);
    return omega * sqrt(MU_0 * EPS_0 * eps_r * 0.5) * sqrt(sqrt_term - 1.0);
}

// Regression linéaire sur log(|E[j]|) vs x = j*dx pour estimer alpha.
// En regime permanent : |E(x)| ~ A * exp(-alpha * x)
// => log(|E|) = log(A) - alpha * x   (droite de pente -alpha)
//
// j_start, j_end : indices de la plage spatiale à fitter (exclut la source et les bords).
// Retourne 0 si OK, 1 si pas assez de points exploitables.
//Code IA
static int fit_log_decay(const double *E, int m, double dx,
                         int j_start, int j_end,
                         double *alpha_out, double *A_out)
{
    double sum_x  = 0.0, sum_y  = 0.0;
    double sum_xx = 0.0, sum_xy = 0.0;
    int n = 0;

    // Ne garder que les maxima locaux de |E|
    for (int j = j_start + 1; j < j_end && j < m - 1; j++) {
        double val  = fabs(E[j]);
        double left = fabs(E[j - 1]);
        double right= fabs(E[j + 1]);

        if (val <= left || val <= right) continue;  // pas un max local
        if (val < 1e-20) continue;

        double x_j = j * dx;
        double y_j = log(val);

        sum_x  += x_j;
        sum_y  += y_j;
        sum_xx += x_j * x_j;
        sum_xy += x_j * y_j;
        n++;
    }

    if (n < 2) return 1;

    double denom = n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-30) return 1;

    double slope     = (n * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / n;

    *alpha_out = -slope;
    *A_out     = exp(intercept);
    return 0;
}

// Lance la validation de la profondeur de peau pour l'eau de mer à plusieurs frequences.
int run_skin_depth_validation(double c) {
    //Declarartion des tableaux de stockage des valeurs
    double freqs[4] = {1e3, 1e5, 1e7, 1e9};
    double delta_theo[4];
    double delta_sim[4];
    double freq_valid[4];
    int valid_count = 0;

    printf("\nFréquence (Hz) | tan(delta) | delta_theo (m) | delta_sim (m) | erreur (%%)\n");
    printf("---------------------------------------------------------------------\n");

    for (int i = 0; i < 4; i++) {
        double f = freqs[i];
        double lambda_f = c / f;
        double omega = 2.0 * PI * f;

        // Théorie
        double alpha_th = compute_alpha(omega, EAU_MER_EPS_R, EAU_MER_SIGMA);
        double delta_th = alpha_th > 0.0 ? 1.0 / alpha_th : INFINITY;

        // Discrétisation adaptee à fréquence spécifique : on recalcule les paramètres de la discretisation en fct de w 
        // dx doit resoudre à la fois lambda et la profondeur de peau
        double dx_f = fmin(lambda_f / 20.0, delta_th / 10.0);   
        double dt_f = dx_f / (2.0 * c);
        double length_f = 15.0 * delta_th;
        int m_f = (int)round(length_f / dx_f);

        if (m_f > 1000000) {
            printf("%g Hz : m_f=%d > 1e6, fréquence ignorée\n", f, m_f);
            continue;
        }
        if (m_f < 10) m_f = 10;

        int steps_per_period = (int)(1.0 / (f * dt_f)) + 1;
        int step_time_f = 20 * steps_per_period;

        if (step_time_f > 16000000) {
            printf("f=%.0e Hz ignorée : %d pas requis\n", f, step_time_f);
            continue;
        }

        // Simulation 1D
        double *E_1D = NULL;
        double *B_1D = NULL;
        double *E_1D_phas = NULL;
        if (calcul_1D_diel(m_f, &E_1D, &B_1D, &E_1D_phas, sin, step_time_f, length_f,
                           dt_f, EPS_0, omega, eps_r_eau_mer, sigma_eau_mer, dx_f, 1.0)) {
            printf("Erreur simulation 1D pour %g Hz\n", f);
            free(E_1D);
            free(B_1D);
            continue;
        }

        // Estimation de l'attenuation simulee
        int j_end_f = m_f - 2; // On cherche sur presque tout le domaine

        double alpha_sim_f = 0.0;
        double A_fit_f = 1.0;
        double *absE_f = malloc(sizeof(double) * m_f);
        for (int j = 0; j < m_f; j++) absE_f[j] = fabs(E_1D[j]);

       int j_src = m_f / 2;  // adapter selon l'implémentation réelle
        int j_fit_start = j_src + 5;  // marge pour éviter le champ proche
        if (fit_log_decay(absE_f, m_f, dx_f, j_fit_start, j_end_f, &alpha_sim_f, &A_fit_f)) {
            printf("Impossible d'estimer alpha pour %g Hz\n", f);
        }
        free(absE_f);

        double delta_sim_f = alpha_sim_f > 0.0 ? 1.0 / alpha_sim_f : INFINITY;
        double tan_d = EAU_MER_SIGMA / (omega * EPS_0 * EAU_MER_EPS_R);
        double error_pct = delta_th > 0.0
                         ? fabs(delta_sim_f - delta_th) / delta_th * 100.0
                         : INFINITY;

        delta_theo[valid_count] = delta_th;
        delta_sim[valid_count] = delta_sim_f;
        freq_valid[valid_count] = f;
        valid_count++;

        printf("%-10g | %.3g | %-10.4g | %-10.4g | %.3g %%\n",
               f, tan_d, delta_th, delta_sim_f, error_pct);

        free(E_1D);
        free(B_1D);
    }

    // Tracé gnuplot
    if (valid_count > 0) {
        FILE *gp = popen("gnuplot -persistent", "w");
        if (gp == NULL) {
            printf("Erreur ouverture gnuplot\n");
            return 1;
        }

        fprintf(gp, "set logscale xy\n");
        fprintf(gp, "set xlabel 'Fréquence (Hz)'\n");
        fprintf(gp, "set ylabel 'Profondeur de peau δ (m)'\n");
        fprintf(gp, "set title 'Profondeur de peau eau de mer — validation FDTD vs théorie'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set style data lines\n");
        fprintf(gp, "plot '-' with lines lc rgb 'blue' lw 2 title 'δ théorique', "
                     "'-' with points pt 7 lc rgb 'red' title 'δ simulé'\n");

        for (int i = 0; i < valid_count; i++)
            fprintf(gp, "%g %g\n", freq_valid[i], delta_theo[i]);
        fprintf(gp, "e\n");

        for (int i = 0; i < valid_count; i++)
            fprintf(gp, "%g %g\n", freq_valid[i], delta_sim[i]);
        fprintf(gp, "e\n");

        fflush(gp);
        pclose(gp);
    }

    return 0;
}