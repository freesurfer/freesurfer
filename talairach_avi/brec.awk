BEGIN {
	depth = 0;
	n = 0;
}

($1 == "rec") {n++;}
($1 == "endrec") {depth--;}

(NF > 0) {
##########################
# find first printing char
##########################
	k = match($0, /["-}]/);
	text = substr($0, k);

##################################
# test for presence of leading "-"
##################################
	lneg = 0;
	if (index(text, "-") == 1) lneg++;

#################################
# convert embedded tabs to spaces
#################################
	l = length(text);
	temp = ""
	for (k = 1; k <= l; k++) {
		c = substr(text, k, 1);
		if (c == "\t") {
			m = length(temp) % 8;
			for (i = 0; i < 8 - m; i++) temp = temp" ";
		}
		else {
			temp = temp"&";
			gsub (/&/, c, temp);
		}
	}
	text = temp;

########
# output
########
	if (dlimit == 0 || n <= dlimit) {
		printf ("%d", n);
		for (i = 0; i < 6*depth - lneg; i++) printf(" ");
		printf ("%s\n", substr (text, 1, 108 - (6*depth)));
	}
}

($1 == "rec") {depth++;}
($1 == "endrec") {n--;}
