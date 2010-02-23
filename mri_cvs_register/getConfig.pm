1;

sub getConfig (%) {
	my %args = @_;
	my %CONFIG = ();

	open(CONFIG, "<$args{'conf_file'}")
		or print STDERR "No conf file in $FILE_LOCATION/conf $!";

	while (<CONFIG>) {
		next if /^\#/;
		chomp;
		next unless $_; # skip empty lines
		s/\=/___EQ___/; # replace first equal sign
		my ($key, $val) = split /___EQ___/;
		$key =~ s/[\t ]//g;
		# trim leading & trailing spaces
		$val = join " ", grep { $_ } split /[\t ]/, $val;
		$CONFIG{$key} = $val;
	}

	close CONFIG;

	\%CONFIG;

}

