#include <stdio.h>
#include <stdlib.h>

int plot_1D(double *field, int m, double dx) {
    FILE *fout = popen("gnuplot -persistant", "w");
    if(fout == NULL) {
        printf("ERREUR lors de l'allocation du fichier\n");
        return 1;
    }

    fprintf(fout, "set xlabel 'x'; set ylabel 'y'\n");
    fprintf(fout, "plot '-' with lines\n");

    for(int i = 0; i < m-1; i++)
        fprintf(fout, "%g %g\n", dx*i, field[i]);

    fprintf(fout, "e\n"); //fin du dataset
    fflush(fout);
    pclose(fout);

    return 0;
}