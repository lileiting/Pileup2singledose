#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;

sub usage{
    print <<"usage";

    $FindBin::Script <table_process_markers.txt>

usage
    exit;
}

sub count_bases{
    my @bases = @_;
    my %stat = (A => 0, C => 0, G => 0, T => 0);
    for my $base (@bases){
        next if $base eq '*';
        while($base =~ /([ATGC])(\d+)/g){
            $stat{$1} += $2;
        }
    }
    return %stat;
}

sub format_hash{
    my %hash = @_;
    my @keys = sort{$hash{$b} <=> $hash{$a}}(keys %hash);
    my @values = @hash{@keys};
    my $ratio = sprintf "%.5f", $values[1] / $values[0];
    return (join("|", @keys), @values, $ratio);
}

sub is_single_dose{
    my $ratio = shift;
    return ($ratio >= 1 / 30 && $ratio <= 1/6) ? 1 : 0;
}

sub process_data{
    my $fh = shift;
    while(my $genotype = <$fh>){
        next unless $genotype =~ /position/;
        next unless $genotype =~ /SNP/;
        chomp $genotype;
        my @F = split /\t/, $genotype; 

        # $F[4] first is maternal parent
        my ($order, $n1, $n2, $n3, $n4, $ratio) =  format_hash(count_bases(@F[5..$#F]));
        next unless is_single_dose($ratio);
        next unless $n2 >= 4 * 30;
        ($order, $n1, $n2, $n3, $n4) =  format_hash(count_bases(@F[4..$#F]));
        print join("\t", $F[0], 
                         join("|", @F[1..3]),
                         $order, $n1, $n2, $n3, $n4, ($n1 + $n2) / 60           
              ), "\n";
    }
}

sub main{
    usage unless @ARGV;
    my $file = shift @ARGV;
    open my $in_fh, "<", $file or die "$file:$!";
    process_data($in_fh);
    close $in_fh;
}

main unless caller;
