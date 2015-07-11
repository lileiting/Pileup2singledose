#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;

sub usage{
    print <<USAGE;

  $FindBin::Script <Matrix>

    Convert a, h, b genotypes to lmxll or nnxnp
    based on the first two genotypes (parent). 
    Assume the first genotype as maternal 
    parent and second genotype as paternal 
    parent. Then if first genotype is h, convert 
    to lmxll, and if the first genotype is a, 
    convert to nnxnp.

    Use this script after filter_genotypes.pl
    and use the output from filter_genotypes.pl
    as input for this script.

USAGE
    exit;
}

sub which_type{
    my ($a, $b) = @_;
    if($a eq 'h'){
        return 'lmxll';
    }elsif($b eq 'h'){
        return 'nnxnp';
    }elsif($a eq '..' and $b eq 'a'){
        return 'lmxll';
    }elsif($a eq 'a' and $b eq '..'){
        return 'nnxnp';
    }else{
        return 0;
    }
}

sub convert{
    my %hash = ( 
       lmxll => {h => 'lm', a => 'll'},
       nnxnp => {h => 'np', a => 'nn'}
    );
    my $type = which_type(@_[3,4]);

    my @new;
    push @new, @_[0..2];
    for my $i (3..$#_){
        my $new_gt = $hash{$type}->{$_[$i]} // $_[$i];
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
        if ($. == 1 and $F[3] ne q/h/ and $F[3] ne q/a/ and $F[3] ne q/../){
            print "$_\n";
            next;
        }
        print join("\t", convert(@F))."\n";
    }
    close $fh;
}

main unless caller;
