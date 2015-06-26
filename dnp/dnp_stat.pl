#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use FindBin;

sub usage{
    print <<USAGE;

  $FindBin::Script CSV [OPTINOS]

     CSV file should be obtained from clean_data_for_dnp_stat.pl

     -o,--output FILE
     -h,--help

USAGE
    exit;
}

sub get_variation{
    my ($type, $str) = @_;
    my $re = $type eq 'dnp' ? '  ' : ' ';
    my $n = 0;
    $n++ while $str =~ /\|\|$re\|\|/g;
    return $n;
}

sub get_snp{ get_variation( 'snp', @_ )}
sub get_dnp{ get_variation( 'dnp', @_ )}

sub analyze_variation{
    my ($type, $str, $size) = @_;
    my @stat;
    my @char = split //, $str;
    my $valid_length = int(@char / $size) * $size; 
    for(my $i = 0; $i < $valid_length; $i += $size){
        my $substr = join("", @char[$i .. $i+ $size - 1]);
        my $num_of_dnp = $type eq 'dnp' ? get_dnp($substr) : get_snp($substr);
        push @stat, $num_of_dnp;
    }
    return @stat;
}

sub analyze_snp { analyze_variation ('snp', @_)}
sub analyze_dnp { analyze_variation ('dnp', @_)}

sub analyze_file{
    my ($infile, $out_fh) = @_;
    
    open my $in_fh, $infile or die "$infile:$!";
    my ($file,$name, $str) = ('', '', '');
    while(<$in_fh>){
        chomp;
        next if /^\s*#/ or /^\s*$/;
        ($file, $name, $str) = split /,/;
        die unless $str;
        for my $size (100, 250, 500,1000){
            my @snp = analyze_snp($str, $size);
            my @dnp = analyze_dnp($str, $size);
            for(0..$#dnp){
                print $out_fh join("\t", $file,$name, 
                                         $size, 
                                         $_+1, 
                                         $dnp[$_], 
                                         $snp[$_]).
                                   "\n";
            }
        }
    }
}

sub main{
    GetOptions(
        "output=s" => \my $outfile,
        "help"  => \my $help
    );
    usage if $help or @ARGV == 0;

    my $out_fh = \*STDOUT;
    open $out_fh, ">", $outfile or die "$outfile: $!" if $outfile;
    print $out_fh join("\t", qw/File       Alignment     Region_size 
                                Region_pos Number_of_DNP Number_of_SNP
                               /)."\n";

    my $file  = shift @ARGV;
    analyze_file($file, $out_fh);
}

main() unless caller;
