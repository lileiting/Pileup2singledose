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

  -i,--min NUM
  -x,--max NUM
    the range for main / mutant 
    Default: min is 6, max is 30

  -t,--threshold NUM
    threshold, default: 3

  -h,--help

USAGE
    exit;
}

############################################################
# Read commands
############################################################

sub read_commands{
    my %para = (
        min => 6,
        max => 30,
        lines => 0,
        not_single_dose => 0,
        is_single_dose => 0,
        threshold => 3
    );

    GetOptions (
        "i|min=i" => \$para{min},
        "x|max=i" => \$para{max},
        "threshold=i" => \$para{threshold},
        "help" => \$para{help}
    );

    usage if $para{help} or @ARGV == 0;
    $para{infile} = shift @ARGV;
    die "Input file not exist" unless -e $para{infile};
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
    }elsif(!&judge($bases,$mutant,1) && &judge($bases,$main,$threshold)){
        return $type eq 'lmxll' ? 'll' : 'nn';
    }else{
        return "..";
    }
}

sub decide_h_a_b{
    my ($bases, $main, $mutant, $para) = @_;
    my $threshold = $para->{threshold};
    if(&judge($bases,$mutant,$threshold) and &judge($bases,$main,$threshold)){
        return 'h';
    }elsif(!&judge($bases,$mutant,1) and &judge($bases,$main,$threshold)){
        return 'a';
    }elsif(&judge($bases,$mutant,$threshold) and !&judge($bases,$main,1)){
        return 'b';
    }else{die}
}

sub variation_type {
    my($main, $mutant) = @_;
    return ($main =~ /[+\-]/ or $mutant =~ /[+\-]/) ? "Indel" : "SNP";
}

sub is_single_dose{
    my($main, $mutant, $para) = @_;
    die "Number of mutant allele is 0!" unless $mutant;
    return ($main/$mutant <= $para->{max} and $main/$mutant >= $para->{min}) ? 1 : 0;
}

sub print_hash{
    my %hash = @_;
    for my $key (keys %hash){
        my $value = $hash{$key};
        print "Key: $key; Value: $value\n";
    }
}



sub process_pileup{
    my $para = shift;
    my($infile, $threshold) = ($para->{infile}, $para->{threshold});
    open my $in_fh, "<", $infile or die "$infile: $!";
    
    while(my $line = <$in_fh>){
        next if $line =~ /^\s*#/ or $line =~ /^\s*$/;
        $para->{lines}++;
        #print $line;
        #next if $line =~ /^\[/;
        my ($id,$female,@progeny)=split(/\s+/,$line);
        #die "Data format incorrect: $line" 
        #    unless $id and $female and @progeny;
        
        map{$_ =~ s/,/./g; $_ =~ tr/atgc/ATGC/}($female,@progeny);
        my %hash=&parse_reads(@progeny);
        next unless keys %hash > 1;
        my ($main,$mutant,@other_nt) = sort {$hash{$b} <=> $hash{$a}} (keys %hash);
        #print_hash(%hash);
        #print "Main: $main $hash{$main}; Mutant: $mutant $hash{$mutant}\n";
        unless(is_single_dose($hash{$main}, $hash{$mutant}, $para)){
            $para->{not_single_dose}++;
            next;
        }else{
            $para->{is_single_dose}++;
        }
        my $type = variation_type($main, $mutant);
        
        # SNP
        print "$id\t<$type>\t<bases>";
        map{@_ = &parse_reads($_);print "\t",@_}($female,@progeny);
        print "\n";
        print "$id\t<$type>\t<genotype>";
        # lmxll
        #if(&judge($female,$main,$threshold) && &judge($female,$mutant,$threshold)){
        #    print "<lmxll>\tlm";
        #    foreach(@progeny){
        #        my $genotype = decide_genotype('lmxll', $_, $main, $mutant, $para);
        #        print "\t$genotype";
        #    }
        #}elsif(&judge($female,$main,$threshold) && !&judge($female,$mutant,1)){
        #    # nnxnp
        #    print "<nnxnp>\tnn";
        #    foreach(@progeny){
        #        my $genotype = decide_genotype('nnxnp', $_, $main, $mutant, $para);
        #        print "\t$genotype";
        #    }
        #}else{
        #    #other
        #    print "<Unknown>";
        #}
        for my $progeny_bases ($female, @progeny){
            my $genotype = decide_h_a_b($progeny_bases, $main, $mutant, $para);
            print "\t$genotype";
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
    
    ###
    printf "Number of lines in file: %d\n", 
        $para->{lines};
    printf "Number of positions exist mutants: %d\n", 
        $para->{not_single_dose} + $para->{is_single_dose};
    printf "Number of positions is single dose: %d\n", 
        $para->{is_single_dose};   
}

main();
