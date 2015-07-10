#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

# URL: https://github.com/lileiting/Pileup2singledose

sub usage{
    print <<USAGE;

perl filtergenotypes.pl [OPTIONS]

  [-i,--in] <GENOTYPES.MATRIX.TXT>
    Results from pileup2singledose.pl

  -t,--threshold DECIMAL 
    Allowed rate of missing data
    Default: 0.5

  -o,--out OUTPUT.txt
    Output file name
    Default: STDOUT

  -h,--help
    Print help

USAGE
    exit;
}

sub read_commands{
    usage unless @ARGV;
    my $infile = '';
    my $outfile = '';
    my $threshold = 0.5;
    my $help;
    GetOptions(
        "i|in=s" => \$infile,
        "o|out=s" => \$outfile,
        "t|threshold=f" => \$threshold,
        "h|help" => \$help
    );
    usage if $help;
    $infile = shift @ARGV if !$infile and @ARGV > 0;
    usage unless $infile;
    open my $in_fh, "<", $infile or die "$infile: $!";
    usage unless $threshold >= 0 and $threshold <= 1;
    my $out_fh = \*STDOUT;
    if($outfile){
        open my $fh, ">", $outfile or die "$outfile: $!";
        $out_fh = $fh;
    }
    return {
        in_fh => $in_fh,
        out_fh => $out_fh,
        threshold => $threshold
    };
}

sub missing_rate{
    my @genotypes = @_[3..$#_];
    my %count = (lm=>0, ll=>0, nn=>0, np => 0, q/-/ => 0);
    map{$count{$_}++}@genotypes;
    my $total_genotypes = $count{lm} + $count{ll} + 
                          $count{nn} + $count{np} + 
                          $count{q/../};
    die "Total genotypes equal to zero ?" if $total_genotypes == 0;
    return $count{q/../}/$total_genotypes;
}

sub message{local $\ = "\n"; print STDERR @_}

sub main{
    my $para = read_commands;
    my $in_fh = $para->{in_fh};
    my $out_fh = $para->{out_fh};
    my $threshold = $para->{threshold};
    my $num_of_lines = 0;
    my $num_of_markers = 0;
    my $num_of_filtered_markers = 0;
    while(<$in_fh>){
        $num_of_lines++;
        next if /^\s*#/ or /^\s*$/;
        chomp;
        my @F = split /\s+/;
        next unless $F[2] eq q/<genotype>/;
        $num_of_markers++;
        next if missing_rate(@F) > $threshold;
        $num_of_filtered_markers++;
        print $out_fh qq/$_\n/;
    }
    message "Threshold: $threshold";
    message "Number of lines in input file: $num_of_lines";
    message "Number of markers in input file: $num_of_markers";
    message "Number of final filtered markers: $num_of_filtered_markers";
    exit;
}

main() unless caller;

