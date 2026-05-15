#include "fonctions.h"

int calcul_2D_diel(int m, double ***E, double ***Bx, double ***By, double ***E_phas, source src, int step_time, double dt, double eps_0, double w, double dx, double (*eps_r_2D)(double, double, double), double (*sigma_2D)(double, double, double), double length, int mask_type, int plot_mode);