#!/usr/bin/env perl

use Set::IntSpan  qw(grep_set map_set);
#use Getopt::Long;
use Getopt::Std;

$[ = 1;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#Define -f option as an string
#&GetOptions("f=s")|| &usage();
getopts('f:')|| &usage();
$fields = $opt_f if $opt_f;

die "Invalid range $fields" unless valid Set::IntSpan $fields;
$set = new   Set::IntSpan $fields;
@elements = elements $set;
print "#Printed elements", @elements;

while (<>) {
    if ( /Size of step / ) { 
                          @Fld = split;
                          $size_of_step = $Fld[5];
	print "#Size of Step ",$size_of_step;
    }

    if ( /^========/ ) {
                        @Fld = split;
                        $timestep = $Fld[3] * $size_of_step;
                        $/ = '--------'; # change the record separator
                        $_ = <>;
                        @Fld = split;
                        print $timestep,@Fld[@elements];
                        $/ = "\n";
    }
}

sub usage {
    my $msg  = shift;
    print STDERR "$0: $msg\n" if $msg;
    die <<EOF;
usage: $0 [-f n,m-p,q] file

        -f      Fields to be output.
EOF
} 
