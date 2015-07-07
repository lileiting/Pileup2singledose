#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;

sub usage{
    print <<USAGE;

  $FindBin::Script <Matrix>

    Convert a, h, b genotypes to lmxll or nnxnp
    based on the first genotype. Assume the 
    first genotype as maternal parent. Then if
    first genotype is h, then convert to lmxll,
    and if the first genotype is a, then convert
    to nnxnp.

USAGE
    exit;
}

sub convert{
    my %hash = ( 
       h => {h => 'lm', a => 'll', '-' => '-'},
       a => {h => 'np', a => 'nn', '-' => '-'}
    );

    my @new;
    push @new, $_[0];
    for my $i (1..$#_){
        my $new_gt = $hash{$_[1]}->{$_[$i]} // $_[$i];
        push @new, $new_gt;
    }
    return @new;
}

sub main{
    usage unless @ARGV;
    my $infile = shift @ARGV;
    open my $fh, "<", $infile or die;
    while(<$fh>){
        chomp;
        my @F = split /\s+/;
        if ($. == 1 and $F[1] ne q/h/ and $F[1] ne q/a/){
            print "$_\n";
            next;
        }
        print join("\t", convert(@F))."\n";
    }
    close $fh;
}

main unless caller;
