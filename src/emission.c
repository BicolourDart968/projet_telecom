#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "fonctions.h"
#include "plot_2D.h"

#define PI 3.14159265359
#define EPS_0 8.854e-12
#define MU_0 12.566e-7

// Réutilise la même allocation que calcul_2D_diel
static double **alloc_field(int m) {
    double **f = malloc(sizeof(double *) * m);
    if (!f) return NULL;
    for (int i = 0; i < m; i++) {
        f[i] = calloc(m, sizeof(double));
        if (!f[i]) {
            for (int k = 0; k < i; k++) free(f[k]);
            free(f);
            return NULL;
        }
    }
    return f;
}
 
// Fonction qui convertit un message en bits
// Retourne le nombre total de bits
// Code IA
int message_to_bits(char *message, int **bits) {
    int nb_bits = 0;
    int len = strlen(message);
    (*bits) = calloc(len*8 + 16, sizeof(int)); //Taille du message + 2 bits d'init/fin
    if(!bits)
        return(1);

        // Préambule : 10101010
        for(int i = 0; i < 8; i++)
            (*bits)[nb_bits++] = (i % 2 == 0) ? 1 : 0;

    for (int i = 0; i < len; i++) {
        char c = message[i];
        for (int b = 7; b >= 0; b--) {
            (*bits)[nb_bits++] = (c >> b) & 1;
        }
    }

    // Séquence de fin : 10101010
    for(int i = 0; i < 8; i++)
        (*bits)[nb_bits++] = (i % 2 == 0) ? 1 : 0;

    return 0;
}

// Convertit un tableau de bits en message texte
// Code IA
void bits_to_message(int *bits, int nb_bits, char *message) {
    int nb_chars = nb_bits / 8;
    
    for (int i = 0; i < nb_chars; i++) {
        char c = 0;
        for (int b = 0; b < 8; b++) {
            c |= bits[i * 8 + b] << (7 - b);
        }
        message[i] = c;
    }
    message[nb_chars] = '\0';
}

int emission(int m, double ***E, double ***Bx, double ***By,
                      int xr, int yr, int step_time,
                      double dt, double eps_0, double dx,
                      double length, double lambda, char *message)
{
    double w0 = 2*PI * 1/sqrt(EPS_0*MU_0) / lambda;
    int frame = 50;
    int x_src = round(m/2);
    int y_src = round(m/2);
    int *bits;
    if(message_to_bits(message, &bits)) {
        printf("Erreur lors de l'allcation du tableau deds bits\n");
        return 1;
    }

    // Allocation 
    *E  = alloc_field(m);
    *Bx = alloc_field(m);
    *By = alloc_field(m);
    if (!*E || !*Bx || !*By) {
        fprintf(stderr, "ERREUR allocation champs\n");
        return 1;
    }
    
    double *Voc = calloc(step_time, sizeof(double));
    if (!Voc) {
        fprintf(stderr, "ERREUR allocation puissance\n");
        return 1;
    }

    // Gnuplot 
    FILE *gp = gp_open();
    if (!gp) {
        gp_close(gp);
        free(Voc);
        return 1;
    }

    /* Paramètres de référence */
    double d = sqrt(pow((xr - x_src)*dx, 2)
                  + pow((yr - y_src)*dx, 2));
    int d0_cells = (int)(2 * lambda / dx);
    double d0 = d0_cells * dx;

    int N_period = 120;       //Pas de temps entre chaque bits
    int n_bits = 0;     //Compteur de bits
    int len = strlen(message)*8+16;    //Taille message + 16bits d'init/fin
    int bit = bits[n_bits];

    //Réception
    double rms_acc = 0.0;                  // accumulateur RMS
    int rms_count = 0;                     // compteur dans la fenêtre
    int delay_steps = (int)(d / (1.0/sqrt(EPS_0*MU_0) * dt));
    int rx_start = delay_steps;            // premier échantillon utile au récepteur
    int rx_phase = 0;                      // 0=préambule, 1=données, 2=fin
    int *rx_bits = calloc(len * 8 + 64, sizeof(int));
    int rx_nb_bits = 0;
    double rx_envelope_max = 0.0;          // pour calibrer le seuil
    char rx_message[256] = {0};

    // Boucle temporelle 
    int progress_step = step_time / 100;
    if (progress_step < 1) progress_step = 1;
    for (int q = 0; q < step_time; q++) {
        double t = dt*q;
        for(int i = 1; i < m-1; i++) {
            for(int j = 1; j < m-1; j++) 
                (*E)[i][j] += 0.5*((*By)[i][j] - (*By)[i-1][j] - ((*Bx)[i][j] - (*Bx)[i][j-1]));
        }
        
        if(q % N_period == 0 && n_bits < len)
            bit = bits[n_bits++];
        else if(q % N_period == 0 && n_bits >= len)
            bit = 0;

        (*E)[x_src][y_src] += sin(w0*t)*bit;

        for(int i = 0; i < m-1; i++) {
            for(int j = 0; j < m-1; j++) {
                (*Bx)[i][j] -= 0.5*((*E)[i][j+1]-(*E)[i][j]);
                (*By)[i][j] += 0.5*((*E)[i+1][j]-(*E)[i][j]);
            }
        }

        // Démodulation récepteur
        //Code IA
        if (q >= rx_start) {
            double ez = (*E)[xr][yr];
            rms_acc += ez * ez;
            rms_count++;

            // Évaluer l'enveloppe tous les N pas de temps (= 1 bit)
            if (rms_count >= N_period) {
                double envelope = sqrt(rms_acc / rms_count);  //Moyenne quadratique
                
                rms_acc = 0.0;
                rms_count = 0;

                // Calibrer le seuil sur les premiers bits
                if (envelope > rx_envelope_max)
                    rx_envelope_max = envelope;

                double seuil = 0.4 * rx_envelope_max;
                int rx_bit = (envelope > seuil) ? 1 : 0;
                printf("Bit reçu : %d\n", rx_bit);
                if (rx_phase == 0) {
                    // Attente du préambule 10101010
                    static int rx_pattern_count = 0;
                    int expected = (rx_pattern_count % 2 == 0) ? 1 : 0;
                    if (rx_bit == expected)
                        rx_pattern_count++;
                    else
                        rx_pattern_count = (rx_bit == 1) ? 1 : 0;

                    if (rx_pattern_count >= 8) {
                        rx_phase = 1;
                        rx_nb_bits = 0;
                        rx_pattern_count = 0;
                        printf("Préambule détecté, réception en cours...\n");
                    }
                }
                else if (rx_phase == 1) {
                    // Réception des données
                    rx_bits[rx_nb_bits++] = rx_bit;

                    // Vérifier si les 8 derniers bits sont tous à 1 (séquence de fin)
                    if (rx_nb_bits >= 8) {
                        int match = 1;
                        for (int b = 0; b < 8; b++) {
                            int expected = (b % 2 == 0) ? 1 : 0;
                            if (rx_bits[rx_nb_bits - 8 + b] != expected) {
                                match = 0;
                                break;
                            }
                        }
                        if (match) {
                            rx_nb_bits -= 8;
                            rx_phase = 2;
                            bits_to_message(rx_bits, rx_nb_bits, rx_message);
                            printf("Message reçu : %s\n", rx_message);
                            return 0;
                        }
                    }
                }
            }

            Voc[q] = (*E)[xr][yr];

            if (q % frame == 0) {  
                fprintf(gp, "set title 'Tension reçue au récepteur'\n");
                fprintf(gp, "set xlabel 'Pas de temps'\n");
                fprintf(gp, "set ylabel 'Voc (V)'\n");
                fprintf(gp, "set xrange [%d:%d]\n", rx_start, step_time);
                fprintf(gp, "plot '-' with lines title 'Voc'\n");
                for (int i = rx_start; i <= q; i++)
                    fprintf(gp, "%d %e\n", i, Voc[i]);
                fprintf(gp, "e\n");
                fflush(gp);
            }
        }
    }
    printf("\n");

    free(Voc);
    gp_close(gp);
    return 0;
}