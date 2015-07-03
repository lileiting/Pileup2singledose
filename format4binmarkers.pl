#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;

sub usage{
    print <<USAGE;

  $FindBin::Script <in>

    Convert format of results from pileup2singledose.pl 
    for BinMarkers

    Examples
        perl $FindBin::Script input.txt > output.txt

USAGE
    exit;
}

sub convert_format{
    my $fh = shift;
    while(<$fh>){
        s/\.\./-/g;
        my @F = split /\t/;
        print join "\t", @F[0,3..$#F];
    }
}

sub main{
    usage unless @ARGV == 1;
    my ($infile, $outfile) = @ARGV;
    open my $fh, "<", $infile or die;
    convert_format($fh);
    close $fh;
}

main() unless caller;
