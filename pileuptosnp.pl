#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

############################################################
# Usage
############################################################

sub usage{
    print <<USAGE;

perl pileuptosnp.pl <genome.pileup.matrix.txt> [OPTIONS]

  -t,--threshold NUM
    threshold

  -h,--help

USAGE
    exit;
}

############################################################
# Read commands
############################################################

sub read_commands{
    my %para = (threshold => 3);

    GetOptions ("help" => \$para{help},
                "threshold" => \$para{threshold}
               );

    usage if $para{help} or @ARGV == 0;
    $para{infile} = shift @ARGV;
    return \%para;
}

############################################################
# Process pileup
############################################################

sub parse_reads{
    my @reads = @_;
    my %nt=();
    my $read;
    while(@reads){
        $read = shift @reads;
        while($read =~ /([.ATGCatgc][+\-](\d+))|([.ATGCatgc](?![+\-]))/g){
            my $match = $&;
            my $n = $2;
            if($match =~ /[+\-]/){
                $read =~ /[A-Za-z]{$n}/g;
                $match .= $&;
                $nt{$match}++;
            }else{
                $nt{$match}++;
            }
        }
    }
    return %nt ? %nt : ("*" => "");
}

sub judge{
    my ($read,$nt,$cut_off)=@_;
    my %count = &parse_reads($read);
    my $nt_count=0;
    if(exists $count{$nt}){
        $nt_count = $count{$nt};
    }
    return $nt_count >= $cut_off ? 1 : 0;
}

sub decide_genotype{
    my ($type, $bases, $main, $mutant, $para) = @_;
    my $threshold = $para->{threshold};
    die unless $type eq 'lmxll' or $type eq 'nnxnp';
    if(&judge($bases,$mutant,$threshold) && &judge($bases,$main,$threshold)){
        return $type eq 'lmxll' ? 'lm' : 'np';
    }elsif(!&judge($bases,$mutant,1) && &judge($bases,$main,3)){
        return $type eq 'lmxll' ? 'll' : 'nn';
    }else{
        return "..";
    }
}

sub variation_type {
    my($main, $mutant) = @_;
    return ($main =~ /[+\-]/ or $mutant =~ /[+\-]/) ? "Indel" : "SNP";
}

sub is_single_dose{
    my($main, $mutant);
    die "Number of mutant allele is 0!" unless $mutant;
    return ($main/$mutant <= 30 and $main/$mutant >= 6) ? 1 : 0;
}

sub process_pileup{
    my $para = shift;
    my($infile, $threshold) = ($para->{infile}, $para->{threshold});
    open my $in_fh, "<", $infile or die;
    
    while(my $line = <$in_fh>){
        my ($id,$female,@progeny)=split(/ /,$line);
        die "Data format incorrect: $line" 
            unless $id and $female and @progeny;
        
        map{$_ =~ s/,/./g; $_ =~ tr/atgc/ATGC/}($female,@progeny);
        my %hash=&parse_reads(@progeny);
        next unless keys %hash > 1;
        my ($main,$mutant,@other_nt) = sort {$hash{$b} <=> $hash{$a}} (keys %hash);
        next unless is_single_dose($hash{$main}, $hash{$mutant});
        my $type = variation_type($main, $mutant);
        
        # SNP
        print "$id\t<$type>\t<position>";
        map{@_ = &parse_reads($_);print "\t",@_}($female,@progeny);
        print "\n";
        print "$id\t<$type>\t";
        # lmxll
        if(&judge($female,$main,$threshold) && &judge($female,$mutant,$threshold)){
            print "<lmxll>\tlm";
            foreach(@progeny){
                my $genotype = decide_genotype('lmxll', $_, $main, $mutant, $para);
                print "\t$genotype";
            }
        }elsif(&judge($female,$main,$threshold) && !&judge($female,$mutant,1)){
            # nnxnp
            print "<nnxnp>\tnn";
            foreach(@progeny){
                my $genotype = decide_genotype('nnxnp', $_, $main, $mutant, $para);
                print "\t$genotype";
            }
        }else{
            #other
            print "<Unknown>";
        }
        print "\n";
    }
    close $in_fh;
}

############################################################
# Main
############################################################

sub main{
    my $para = read_commands;
    process_pileup($para);
}

main();
