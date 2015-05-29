#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

sub usage{
    print <<USAGE;

perl get_matrix.pl [OPTIONS]

    -m,--male <male.pileup>
    -f,--female <female.pileup>
    -p,--progeny <progeny1.pileup>
    -p,--progeny <progeny2.pileup>,<progeny3.pileup>
    
    -o,--output <pileup.matrix.txt>
    
    -h,--help
    
USAGE
    exit;
}

sub read_commands{
    usage unless @ARGV;
    my $male = '';
    my $female = '';
    my @progenies = ();
    my $output = '';
    my $help;
    GetOptions(
        "m|male=s" => \$male,
        "f|female=s" => \$female,
        "p|progeny=s" => \@progenies,
        "o|output=s" => \$output,
        "h|help" => \$help    
    );
    usage if $help;
    @progenies = split(/,/, join(',', @progenies));
    
    my $num_infile = 0;
    $num_infile++ if $male;
    $num_infile++ if $female;
    $num_infile += scalar(@progenies);
    usage unless $num_infile > 0;
    
    my $out_fh = \*STDOUT;
    if($output){
        open my $fh, ">", $output or die "$output: $!";
        $out_fh = $fh;
    }
    return {
        male => $male, 
        female => $female, 
        prog_ref =>\@progenies, 
        out_fh => $out_fh};
}

sub main{
    my $para = read_commands;
    my $male = $para->{male};
    my $female = $para->{female};
    my $prog_ref = $para->{prog_ref};
    my $out_fh = $para->{out_fh};
    my %super_hash;
    
    my @files;
    for my $file ($female, $male, @$prog_ref){
        next unless $file and -e $file;
        push @files, $file;
        warn "Loading pileup data from $file ...\n";
        my $num_line = 0;
        my $num_valid = 0;
        open my $fh, "<", $file or die "$file: $!";
        while(<$fh>){
            $num_line++;
            next unless /^(\S+)\s+(\d+)\s+([ATGCatgcNn])\s+(\d+)\s+(\S+)\s+(\S+)/;
            $num_valid++;
            my($chr, $pos, $nt, $num, $bases, $qual) = ($1, $2, $3, $4, $5, $6);
            $super_hash{"$chr-$pos"}->{$file} = $bases;
        }
        close $fh;
        warn "$num_line lines, $num_valid valid lines\n";
    }
    
    my @positions = sort{$a cmp $b}(keys %super_hash);

    # print title
    print $out_fh join("\t", "#position", @files), "\n";
    
    # print bases
    for my $position (@positions){
        my @line = ($position);
        for my $file (@files){
            if($super_hash{$position}->{$file}){
                push @line, $super_hash{$position}->{$file};
            }else{
                push @line, "*";
            }
        }
        print $out_fh join("\t", @line), "\n";
    }
}

main;