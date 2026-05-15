//Code généré par IA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estimate_refraction_angle.h"

#define PI 3.14159265359

void estimate_refraction_angle(double **E, int m, double length, double dx) {
    int x_start = (int)ceil((length / 2.0) / dx);
    if (x_start < 1) x_start = 1;
    if (x_start >= m - 1) x_start = m - 2;

    int count = 0;
    double sum_y = 0.0;
    double sum_x = 0.0;
    double sum_yy = 0.0;
    double sum_xy = 0.0;

    for (int y = 1; y < m - 1; y++) {
        double x_zero = NAN;
        // Chercher le dernier zéro dans le prisme, pas le premier
        for (int x = m - 2; x >= x_start; x--) {
            double v1 = E[x][y];
            double v2 = E[x + 1][y];
            if (v1 * v2 < 0.0) {
                x_zero = x + fabs(v1) / (fabs(v1) + fabs(v2));
                break;
            }
        }
        if (!isnan(x_zero)) {
            // Vérifier que le point est bien dans le prisme
            double px = x_zero * dx;
            double py = (double)y * dx;
            if (!(py < px - length/2.0 && px > length/2.0 && py < length/2.0))
                continue;
            double yi = (double)y;
            sum_y += yi;
            sum_x += x_zero;
            sum_yy += yi * yi;
            sum_xy += x_zero * yi;
            count++;
        }
    }

    if (count < 2) {
        printf("Impossible d\'estimer l\'angle de réfraction : pas assez de zéros détectés.\n");
        return;
    }

    double denom = count * sum_yy - sum_y * sum_y;
    if (fabs(denom) < 1e-12) {
        printf("Impossible d\'estimer l\'angle de réfraction : régression instable.\n");
        return;
    }

    double a = (count * sum_xy - sum_x * sum_y) / denom;
    double theta_t = atan(a);
    double theta_t_deg = theta_t * 180.0 / PI;

    double sin_theta_i = sqrt(2.0) / 2.0;
    double sin_theta_t = sin_theta_i / sqrt(5.0);
    if (sin_theta_t > 1.0) sin_theta_t = 1.0;
    if (sin_theta_t < -1.0) sin_theta_t = -1.0;
    double theta_t_theo = asin(sin_theta_t);
    double theta_t_theo_deg = theta_t_theo * 180.0 / PI;

    printf("Angle de réfraction simulé : %.2f deg\n", theta_t_deg);
    printf("Angle de réfraction théorique (eps_r=5, theta_i=45deg) : %.2f deg\n", theta_t_theo_deg);
}