#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

sub usage{
    print <<USAGE;

perl pileuptosnp.pl <genome.pileup.matrix.txt> [OPTIONS]

  -t,--threshold NUM
    threshold

  -h,--help

USAGE
    exit;
}

sub main{
    my $help;
    my $threshold = 3;

    GetOptions ("help" => \$help,
                "threshold" => \$threshold
               );

    usage if $help or @ARGV == 0;
    my $infile = shift @ARGV;

    process_pileup($infile, $threshold);
}

sub process_pileup{
    my($infile, $threshold) = @_;
    open my $in_fh, "<", $infile or die;
    
    while(my $line = <$in_fh>){
        my ($id,$female,@progeny)=split(/ /,$line);
        map{$_ =~ s/,/./g; $_ =~ tr/atgc/ATGC/}($female,@progeny);
        my %hash=&parse_reads(@progeny);
        unless((keys %hash) > 1){next;}
        my ($main,$mutant,@other_nt) = sort {$hash{$b} <=> $hash{$a}} (keys %hash);
        #<--single dose
        unless($hash{$main}/$hash{$mutant} <= 30 && $hash{$main}/$hash{$mutant} >= 6){next;}
        #single dose -->
        my $type;
        if($main =~ /[+\-]/ || $mutant =~ /[+\-]/){$type = "Indel";}
        else{$type = "SNP";}
        # SNP
        print "$id\t<$type>\t<position>";
        map{@_ = &parse_reads($_);print "\t",@_}($female,@progeny);
        print "\n";
        print "$id\t<$type>\t";
        # lmxll
        if(&judge($female,$main,$threshold) && &judge($female,$mutant,$threshold)){
            print "<lmxll>\tlm";
            foreach(@progeny){
                if(&judge($_,$mutant,$threshold) && &judge($_,$main,$threshold)){print "\tlm";}
                elsif(!&judge($_,$mutant,1) && &judge($_,$main,3)){print "\tll";}
                else{print "\t..";}
            }
        }
        # nnxnp
        elsif(&judge($female,$main,$threshold) && !&judge($female,$mutant,1)){
            print "<nnxnp>\tnn";
            foreach(@progeny){
                if(&judge($_,$mutant,$threshold)  && &judge($_,$main,$threshold)){print "\tnp";}
                elsif(!&judge($_,$mutant,1) && &judge($_,$main,3)){print "\tnn";}
                else{print "\t..";}
            }
        }
        #other
        else{print "<Unknown>";}
        print "\n";
    }
    close $in_fh;
}

sub judge{
    my ($read,$nt,$cut_off)=@_;
    my %count = &parse_reads($read);
    my $nt_count=0;
    if(exists $count{$nt}){
        $nt_count = $count{$nt};
    }
    if($nt_count >= $cut_off){1;}
    else{0;}
}
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
            }else{$nt{$match}++;}
        }
    }
    if(%nt){%nt;}
    else{("*","");}
}

main();
