CC = gcc
CFLAGS = -Wall -Wextra -g -std=gnu11 -Wno-unused-parameter
LDFLAGS = -lm

TARGET = main
SRC = main.c calcul_1D_libre.c plot_1D.c calcul_1D_diel.c fonctions.c calcul_2D_libre.c plot_2D.c calcul_2D_diel.c calcul_2D_antenne.c simulations_1D.c estimate_refraction_angle.c debeye.c
OBJ = $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.c calcul_1D.h plot_1D.h calcul_1D_diel.h fonctions.h calcul_2D_libre.h plot_2D.h calcul_2D_diel.h calcul_2D_antenne.h simulations_1D.h estimate_refraction_angle.h debeye.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
