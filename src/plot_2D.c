//Code genere par IA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "plot_2D.h"

/* ================================================================
 *  Helper interne
 * ================================================================ */

static void write_matrix(FILE *gp, double **field, int m) {
    for (int j = m - 1; j >= 0; j--) {
        for (int i = 0; i < m; i++)
            fprintf(gp, "%g ", field[i][j]);
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
}

/* ================================================================
 *  Cycle de vie
 * ================================================================ */

static int gp_window_id = 0;

FILE *gp_open(void) {
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        fprintf(stderr, "ERREUR : impossible d'ouvrir gnuplot\n");
        return NULL;
    }
    fprintf(gp, "set terminal qt %d\n", gp_window_id++);
    fflush(gp);
    return gp;
}

void gp_close(FILE *gp) {
    if (gp) pclose(gp);
}

/* ================================================================
 *  Configuration (une seule fois par fenetre)
 * ================================================================ */

void gp_setup_image(FILE *gp, int m, double cb_min, double cb_max,
                    const char *title) {
    if (!gp) return;
    fprintf(gp, "set title '%s'\n", title ? title : "");
    fprintf(gp, "set xrange [0:%d]\n", m);
    fprintf(gp, "set yrange [0:%d]\n", m);
    fprintf(gp, "set cbrange [%g:%g]\n", cb_min, cb_max);
    fprintf(gp, "set xlabel 'x'; set ylabel 'y'\n");
    fprintf(gp, "set palette defined (0 '#50c878', 1 '#8a2be2')\n");
    fflush(gp);
}

void gp_setup_curve(FILE *gp, const char *title,
                    const char *xlabel, const char *ylabel) {
    if (!gp) return;
    fprintf(gp, "set title '%s'\n", title ? title : "");
    fprintf(gp, "set xlabel '%s'; set ylabel '%s'\n",
            xlabel ? xlabel : "", ylabel ? ylabel : "");
    fprintf(gp, "set grid\n");
    fflush(gp);
}

/* ================================================================
 *  Envoi de donnees (par frame)
 * ================================================================ */

void gp_plot_field(FILE *gp, double **field, int m) {
    if (!gp || !field) return;
    fprintf(gp, "plot '-' matrix with image notitle\n");
    write_matrix(gp, field, m);
    fflush(gp);
}

void gp_plot_curve(FILE *gp, double *data, int n, const char *color) {
    if (!gp || !data) return;
    fprintf(gp, "plot '-' with lines lc rgb '%s' notitle\n",
            color ? color : "red");
    for (int i = 0; i < n; i++)
        fprintf(gp, "%d %g\n", i, data[i]);
    fprintf(gp, "e\n");
    fflush(gp);
}

void gp_plot_xy_curve(FILE *gp, double *x, double *y, int n,
                      const char *color, const char *label) {
    if (!gp || !x || !y) return;
    fprintf(gp, "plot '-' with lines lc rgb '%s' title '%s'\n",
            color ? color : "red", label ? label : "data");
    for (int i = 0; i < n; i++)
        fprintf(gp, "%g %g\n", x[i], y[i]);
    fprintf(gp, "e\n");
    fflush(gp);
}

void gp_plot_power(FILE *gp, double *P_sim, double *P_theo,
                          int q, double dt, int q_start)
{
    if (!gp) return;

    fprintf(gp, "set xlabel 't (s)'; set ylabel 'P (W)'\n");
    fprintf(gp, "set title 'Puissance recue'\n");
    fprintf(gp, "set style data lines\n");
    fprintf(gp, "set xrange [%g:*]\n", dt * q_start);
    fprintf(gp, "set yrange [0:*]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '-' w lines lc rgb 'red' lw 2 t 'Simulee',"
                " '-' w lines lc rgb 'blue' lw 2 t 'Theorique',"
                " '-' w lines lc rgb 'green' lw 2 t 'Difference'\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, P_sim[i]);
    fprintf(gp, "e\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, P_theo[i]);
    fprintf(gp, "e\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, fabs(P_sim[i] - P_theo[i]));
    fprintf(gp, "e\n");

    fflush(gp);
}

/**
 * Plot de la tension reçue par l'antenne
 */
void gp_plot_voc(FILE *gp, double *Voc, double *Voc_theo,
                        int q, double dt, int q_start)
{
    if (!gp) return;

    fprintf(gp, "set xlabel 't (s)'; set ylabel 'Voc (V)'\n");
    fprintf(gp, "set title 'Tension reçue'\n");
    fprintf(gp, "set style data lines\n");
    fprintf(gp, "set xrange [%g:*]\n", dt * q_start);
    fprintf(gp, "set yrange [0:*]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '-' w lines lc rgb 'red' lw 2 t 'Simulee',"
                " '-' w lines lc rgb 'blue' lw 2 t 'Theorique',"
                " '-' w lines lc rgb 'green' lw 2 t 'Difference'\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, Voc[i]);
    fprintf(gp, "e\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, Voc_theo[i]);
    fprintf(gp, "e\n");

    for (int i = q_start; i <= q; i++)
        fprintf(gp, "%g %g\n", dt * i, fabs(Voc[i] - Voc_theo[i]));
    fprintf(gp, "e\n");

    fflush(gp);
}
