#include "fonctions.h"

int calcul_2D_antenne(int m, double ***E, double ***Bx, double ***By, double ***E_phas, source src, int xr, int yr, int step_time, double dt, double eps_0, double dx, double (*eps_r_2D)(double, double, double), double (*sigma_2D)(double, double, double), double length, double lambda, double he, double Ra);