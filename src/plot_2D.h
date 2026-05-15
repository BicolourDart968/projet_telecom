#ifndef PLOT_2D_H
#define PLOT_2D_H

#include <stdio.h>

/**
 * Module de plot 2D unifié pour gnuplot.
 *
 * Toutes les fonctions utilisent "plot '-' matrix with image" (2D pur).
 * Pas de multiplot — pour afficher deux champs, on ouvre deux fenêtres.
 *
 * Workflow simple :
 *   FILE *gp = gp_open();
 *   gp_setup_image(gp, m, -2.0, 2.0, "Ez");
 *   for (q ...) { gp_plot_field(gp, E, m); }
 *   gp_close(gp);
 *
 * Workflow deux fenêtres (Poynting) :
 *   FILE *gp1 = gp_open();  FILE *gp2 = gp_open();
 *   gp_setup_image(gp1, m, -cb, cb, "S_x");
 *   gp_setup_image(gp2, m, -cb, cb, "S_y");
 *   for (q ...) {
 *       gp_plot_field(gp1, Sx, m);
 *       gp_plot_field(gp2, Sy, m);
 *   }
 *   gp_close(gp1); gp_close(gp2);
 */

/* --- Cycle de vie --- */
FILE *gp_open(void);
void  gp_close(FILE *gp);

/* --- Configuration (appeler une fois avant la boucle) --- */
void gp_setup_image(FILE *gp, int m, double cb_min, double cb_max,
                    const char *title);

void gp_setup_curve(FILE *gp, const char *title,
                    const char *xlabel, const char *ylabel);

/* --- Envoi de données (appeler à chaque frame) --- */
void gp_plot_field(FILE *gp, double **field, int m);

void gp_plot_curve(FILE *gp, double *data, int n, const char *color);

#endif /* PLOT_2D_H */