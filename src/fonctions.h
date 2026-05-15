#ifndef FONCTIONS_H
#define FONCTIONS_H

#define PI 3.14159265359
#define EPS_0 8.854e-12
#define MU_0 12.566e-7

double sigma(double x, double length);
double esp_r(double x, double length);
double mur_conduct(double x, double y, double length);
double sigma_absorbant(double x, double y, double length);
double eps_r_prisme(double x, double y, double length);
double eps_r_guide(double x, double y, double length);
double double_slit(double x, double y, double length);
double sigma_0(double x, double y, double length);
double eps_r_0(double x, double y, double length);
double eps_r_milieu(double x, double length);
double sigma_milieu(double x, double length);
double eps_r_eau_mer(double x, double length);
double sigma_eau_mer(double x, double length);
double constante(double t);
double gaussienne(double t);
void set_paquet_params(double w, double t0, double sigma_t);
double paquet_gauss_modul(double omega_t);
double eps_r_libre_1D(double x, double length);
double sigma_libre_1D(double x, double length);
double young_slits(double x, double y, double length, double dx);
double single_slit(double x, double y, double length, double dx);
double eps_r_lentille(double x, double y, double length);
double eps_r_fabry_perot(double x, double y, double length);
double fabry_perrot_mask(double x, double y, double length, double dx);

typedef struct{         //Structure du type source
    double A;
    int y_start, y_end;
    int x;
    double (*forme)(double);    //Type de fonction (sin, gaussienne etc)
    double w;                   // Initialiser à 1 si pas sinusoïde
    double k;                   // phase gradient le long de y pour source inclinée
}source;

#endif
