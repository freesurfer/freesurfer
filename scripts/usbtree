#!/usr/bin/perl
#
# Reads /proc/bus/usb/devices and selectively lists and/or
# interprets it.
# Options:
#	"-l" lists vendor/product ID and revision.
#	[filename] reads USB product data from filename instead of
#	  the default location (/proc/bus/usb/devices).
#
# Originally written by Randy Dunlap. 

$PROGNAME = $0;
$DEVFILENAME = "/proc/bus/usb/devices";
$showconfig = "yes";
$listids = 0;

if ($#ARGV > -1)
{
	if ("$ARGV[0]" eq "-l")
	{
		$listids = 1;
	}
	else
	{
		$DEVFILENAME = $ARGV[0];
	}
}

if (!open(DEVNUM, "<$DEVFILENAME"))
{
	print "$PROGNAME: cannot open '$DEVFILENAME'\n";
	exit 1;
}

while ($line = <DEVNUM>)	# read a text line from DEVNUM
{
	# skip all lines except those we recognize:
	if (($line !~ "^C:")		# Configuration: one is active
	    && ($line !~ "^D:")	# Device:
	    && ($line !~ "^I:")	# Interface: protocol group
	    && ($line !~ "^S:") # String: used with root hub
	    && ($line !~ "^P:")	# Ids
	    && ($line !~ "^T:")	# Topology: starts each device
	    )
	{
		next;	# to the next line
	}

	chomp $line;		# remove line endings

	# First convert '=' signs to spaces.
	$line =~ tr/=/ /;

	# and convert all '(' and ')' to spaces.
	$line =~ tr/(/ /;
	$line =~ tr/)/ /;

	# split the line at spaces.
	@fields = split / +/, $line;

	# T:  Bus=01 Lev=01 Prnt=01 Port=03 Cnt=01 Dev#=  3 Spd=1.5 MxCh= 0
	if ($line =~ "^T:")
	{
		# split yields: $bus, $level, $parent, $port, $count, $devnum, $speed, $maxchild.

		$bus    = @fields [2];
		$level  = @fields [4];
		$parent = @fields [6];		# parent devnum
		$port   = @fields [8] + 1;	# make $port 1-based
		$count  = @fields [10];
		$devnum = @fields [12];
		$speed  = @fields [14];
		$maxchild = @fields [16];
		$devclass = "?";
		$intclass = "?";
		$driver   = "?";
		$ifnum    = "?";
		$showclass = "?";	# derived from $devclass or $intclass
		$lastif = "?";			# show only first altsetting
		$HCtype = "?";
		$showconfig = "no";
		$nconfig = "0";
		next;
	} # end T: line

	# only show the _active_ configuration
	# C:* #Ifs= 1 Cfg#= 1 Atr=a0 MxPwr=100mA
	elsif ($line =~ "^C:") {
	    if ($line =~ "^C:\\*") {
		$showconfig = @fields[4];
	    } else {
		$showconfig = "no";
	    }
	    next;
	}

	# D:  Ver= 1.00 Cls=00(>ifc ) Sub=00 Prot=00 MxPS= 8 #Cfgs=  1
	elsif ($line =~ "^D:")
	{ # for D: line
		$devclass = @fields [5];
		$nconfig = @fields [13];
		next;
	}

	# P:  Vendor=0aa7 ProdID=0304 Rev= 0.52
	elsif ($line =~ "^P:")
	{
		$ids = $line;
		$ids =~ s/P: +//;
		$ids =~ s/ +/ /g;
		$ids =~ s/Vendor /Vendor=/g;
		$ids =~ s/ ProdID /, ProdID=/g;
		$ids =~ s/ Rev /, Rev=/g;
		next;
	}
	# in case this is a root hub, look at the device strings.
	#  - S:  Manufacturer:Linux 2.6.5 ehci_hcd	[all 2.6]
	#  - S:  Product=USB UHCI Root Hub 		[all 2.4, 2.2]
	#  - S:  Product=OPTi Inc 82C861		[2.6/PCI_NAMES]
	elsif ( $line =~ "^S:" )
	{ # for S: line
		if ($level == 00 && $line =~ "hcd")
		{
		    $HCtype = @fields [4];
		}
		elsif ($level == 00 && $line =~ "HCI" && $HCtype eq "?")
		{
		    $HCtype = @fields [3];
		}
		
		if ($line =~ "Product")
		{
			$product = $line;
			$product =~ s/Product //;
			$product =~ s/S: +//;
			$product =~ s/ +/ /g;
			$product =~ s/ $//g;
		}
		next;
	}

	# the rest of this code:
	#  - only shows interface descriptors
	#  - for the active configuration
	#  - for the first (prefer: active!) altsetting
	elsif (!($line =~ "^I:")
		|| "$showconfig" eq "no") {
	    next;
	}

	
	# I:  If#= 0 Alt= 0 #EPs= 1 Cls=03(HID  ) Sub=01 Prot=02 Driver=hid
	$intclass = @fields [9];
	$ifnum    = @fields [2];
	$driver   = @fields [15];

	if (($devclass eq ">ifc") || ($devclass eq "unk."))
	{	# then use InterfaceClass, not DeviceClass
		$showclass = $intclass;
	}
	else
	{	# use DeviceClass
		$showclass = $devclass;
	}

	if ($level == 0)
	{
	    # substitute real driver name
	    if ($HCtype =~ "UHCI-alt")
	    {
		$HC = "uhci";
	    }
	    elsif ($HCtype =~ "UHCI")
	    {
		$HC = "usb-uhci";
	    }
	    elsif ($HCtype =~ "OHCI")
	    {
		$HC = "usb-ohci";
	    }
	    else
	    {
		$HC = $HCtype;
	    }

		print sprintf ("/: Bus $bus.Port $port: Dev $devnum, Class=root_hub, Drv=%s/%sp, %sM\n",
			 $HC, $maxchild, $speed );
	}
	elsif ($lastif ne $ifnum)
	{
		$temp = $level;
		while ($temp >= 1)
		{
			print "    ";
			$temp--;
		}

		if ($nconfig ne "1") {
		    $temp = " Cfg $showconfig/$nconfig";
		} else {
		    $temp = "";
		}

		print sprintf ("|_ Port $port: Dev $devnum$temp, If $ifnum, Prod=$product, Class=$showclass, Drv=$driver%s, %sM",
			($maxchild == 0) ? "" : ("/" . $maxchild . "p"),
			$speed);

		if ($listids == 1)
		{
			print ", $ids";
		}
		print "\n";
		$lastif = $ifnum;
	}
	$product="";
} # end while DEVNUM

close (DEVNUM);

# END.
