#include <math.h>
#include <stdlib.h>

#define LAMBDA 0.01
#define DX LAMBDA/20

double sigma(double x, double length) {
    /*if (x < 0 || x > length)
        return 0;
    else if (x >= 0 && x <= length/3
    )
        return 0;
    else if (x >= length/3 && x <= 2*length/3)
        return 1e-1;
    else
        return 0;*/
    return 1.3e-1;
}

double esp_r(double x, double length) {
    if (x < 0 || x > length)
        return 1;
    else if (x >= 0 && x <= length/3)
        return 1;
    else if (x >= length/3 && x <= 2*length/3)
        return 1;
    else
        return 1;
}

double eps_r_prisme(double x, double y, double length) {
    if (y < x - length/2 && x > length/2 && y < length/2)
        return 5.0;
    return 1.0;
}

double eps_r_guide(double x, double y, double length) {
    double border = length * 0.25;
    if (y < border || y > length - border)
        return 2.0;
    return 5.0;
}

double eps_r_0(double x, double y, double length) {
    return 1.0;
}

double sigma_0(double x, double y, double length) {
    return 0.0;
}

double eps_r_milieu(double x, double length) {
    (void) length;
    return 4.0;
}

double sigma_milieu(double x, double length) {
    (void) length;
    return 0.5;
}

double eps_r_eau_mer(double x, double length) {
    (void) length;
    return 80.0;
}

double sigma_eau_mer(double x, double length) {
    (void) length;
    return 4.0;
}

double constante(double t){return 1.0;}

double gaussienne(double t) {
    double sigma = 1.0;  // largeur de la pulsation
    return exp(-t*t / (2*sigma*sigma));
}

static double paquet_w = 1.0;
static double paquet_t0 = 0.0;
static double paquet_sigma_t = 1.0;

void set_paquet_params(double w, double t0, double sigma_t) {
    paquet_w = w;
    paquet_t0 = t0;
    paquet_sigma_t = sigma_t;
}

double paquet_gauss_modul(double omega_t) {
    double t = omega_t / paquet_w;
    double dt = t - paquet_t0;
    return exp(-dt*dt / (2.0 * paquet_sigma_t * paquet_sigma_t)) * sin(paquet_w * t);
}

double eps_r_libre_1D(double x, double length) {
    (void) x; (void) length;
    return 1.0;
}

double sigma_libre_1D(double x, double length) {
    (void) x; (void) length;
    return 0.0;
}

double young_slits(double x, double y, double length, double dx) {
    //On place le centre du mur au centre du domaine
    double center_x = length / 2.0;
    double thickness = 5.0 * dx;
    //position des fentes : un peu au-dessus et en dessous du milieu du mur
    double slit_center1 = length / 2.0 - length / 16.0;   
    double slit_center2 = length / 2.0 + length / 16.0;
    double slit_width = length / 20.0;

    if (fabs(x - center_x) <= thickness / 2.0) {
        // in the barrier
        if (fabs(y - slit_center1) <= slit_width / 2.0 || fabs(y - slit_center2) <= slit_width / 2.0) {
            return 1.0;  // slit
        } else {
            return 0.0;  // barrier
        }
    } else {
        return 1.0;
    }
}

double single_slit(double x, double y, double length, double dx) {
    double center_x = length / 2.0;
    double thickness = 5.0 * dx;
    double slit_center = length / 2.0;
    double slit_width = length / 20.0;

    if (fabs(x - center_x) <= thickness / 2.0) {
        if (fabs(y - slit_center) <= slit_width / 2.0) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    return 1.0;
}

double eps_r_lentille(double x, double y, double length)
{
    double n      = 2.5;
    double R      = 0.25;              // rayon de courbure [m]
    double delta0 = 0.04;             // épaisseur au centre [m]
    double x_c    = length / 3.0;     // lentille au tiers du domaine
    double y_c    = length / 2.0;     // centrée verticalement

    double dy = y - y_c;
    double delta_y = delta0 - (dy * dy) / R;   // épaisseur locale

    if (delta_y <= 0.0) return 1.0;

    if (fabs(x - x_c) <= delta_y / 2.0)
        return n * n;                

    return 1.0;
}

double eps_r_fabry_perot(double x, double y, double length)
{
    (void) y;

    double eps_mir = 10000000.0;                  /* n = 3, R ≈ 25 % par interface */
    double thick   = 2.0 * DX;            /* épaisseur miroir ≈ 2 cellules  */

    /* cavité : L = 5 lambda (= 10 * lambda/2 → résonance ordre m=10) */
    double L_cav   = 5.0 * LAMBDA;

    /* centre de la cavité au milieu du domaine */
    double center  = length / 2.0;
    double mir1    = center - L_cav / 2.0; /* miroir gauche */
    double mir2    = center + L_cav / 2.0; /* miroir droit  */

    if (fabs(x - mir1) <= thick / 2.0 || fabs(x - mir2) <= thick / 2.0)
        return eps_mir;

    return 1.0;
}

double fabry_perrot_mask(double x, double y, double length, double dx)
{
    (void) y;
    double thick = 3.0 * dx;
    double L_cav = 4.0 * LAMBDA;
    double center = length / 2.0;
    double mir1 = center - L_cav / 2.0;
    double mir2 = center + L_cav / 2.0;

    if (fabs(x - mir1) <= thick / 2.0 || fabs(x - mir2) <= thick / 2.0)
        return 0.25;
    return 1.0;
}