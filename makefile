all:
	gcc molecular.c -o mol -lm -Wall
debug:
	gcc molecular.c -o mol -lm -Wall -g
