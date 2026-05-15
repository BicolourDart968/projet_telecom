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