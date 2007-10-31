BEGIN {go = 0;}

(NF == 6 && go > 0) {
	for (i = 1; i <= 6; i++) p[k++] = $i;
	go++; if (go == 3) {go = 0;}
}

/eta,q/ {p[k++] = $2;}

/second partial in parameter space/ {go++; k = 0;}

END {
	for (k = 0; k < 12; k++) printf ("%5.0f", p[k]);
	printf ("%10.5f\n", p[12]);
}
