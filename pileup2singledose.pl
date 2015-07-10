#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin;

# URL: https://github.com/lileiting/Pileup2singledose

############################################################
# Usage
############################################################

sub usage{
    print <<USAGE;

  $FindBin::Script <genome.pileup.matrix.txt> [OPTIONS]

    -i,--min NUM
    -x,--max NUM
      the range for main / mutant 
      Default: min is 6, max is 30

    -t,--threshold NUM
      threshold, default: 3

    -v,--verbose

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
        threshold => 3,
        verbose => 0
    );

    GetOptions (
        "i|min=i" => \$para{min},
        "x|max=i" => \$para{max},
        "threshold=i" => \$para{threshold},
        "verbose" => \$para{verbose},
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
    for my $read (@reads){
        while($read =~ /(([.ATGCatgc][+\-](\d+))|([.ATGCatgc](?![+\-])))/g){
            my $match = $1;
            my $n = $3;
            if($match =~ /[+\-]/){
                $read =~ /([A-Za-z]{$n})/g;
                $match .= $1;
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

sub decide_segregation_type{
    my ($female, $male, $main, $mutant, $para) = @_;
    my $threshold = $para->{threshold};

    my $verbose = $para->{verbose};
    print "# $.: Female: $female; Male: $male; Main: $main; Mutant: $mutant\n"
        if $verbose;

    if(&judge($female,$mutant,$threshold) and &judge($female,$main,$threshold)
#        and not &judge($male,$mutant,1) && &judge($male,$main,$threshold)
    ){
        return 'lmxll';
    }elsif(&judge($male,$mutant,$threshold) and &judge($male,$main,$threshold)
#        and not &judge($female,$mutant,1) && &judge($female,$main,$threshold)
    ){
        return 'nnxnp';
    }else{
        return "unexpected";
    }
}

sub decide_genotype{
    my ($type, $bases, $main, $mutant, $para) = @_;
    my $threshold = $para->{threshold};
    return 'Undef' unless $type eq 'lmxll' or $type eq 'nnxnp';
    if(&judge($bases,$mutant,$threshold) && &judge($bases,$main,$threshold)){
        return $type eq 'lmxll' ? 'lm' : 'np';
    }elsif(!&judge($bases,$mutant,1) && &judge($bases,$main,$threshold)){
        return $type eq 'lmxll' ? 'll' : 'nn';
    }else{
        return "-";
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
    }else{return ".."}
}

sub variation_type {
    my($main, $mutant) = @_;
    return ($main =~ /[+\-]/ or $mutant =~ /[+\-]/) ? "Indel" : "SNP";
}

sub is_single_dose{
    my($main, $mutant, $para) = @_;
    my $threshold = $para->{threshold};
    die "Number of mutant allele is 0!" unless $mutant;
    return (
        $main/$mutant <= $para->{max} and 
        $main/$mutant >= $para->{min} and 
        $mutant >= $threshold
    ) ? 1 : 0;
}

sub print_hash{
    my %hash = @_;
    for my $key (keys %hash){
        my $value = $hash{$key};
        print "Key: $key; Value: $value\n";
    }
}

sub print_hash_str{
    my %hash = @_;
    if(keys %hash == 0){
        return 'empty';
    }elsif(keys %hash == 1){
        my ($key, $value) = each %hash;
        return "$key|$value"
    }else{
        my @array;
        for my $key (keys %hash){
            my $value = $hash{$key};
            push @array, "$key|$value";
        }
        return join(";", @array);
    }
}

# Assume the first column is maternal parent, second column
# is paternal parent.
sub process_pileup{
    my $para = shift;
    my($infile, $threshold, $verbose) = @{$para}{qw/infile threshold verbose/};
    open my $in_fh, "<", $infile or die "$infile: $!";
    
    while(my $line = <$in_fh>){
        # Assume first line as title and print it directly
        print $line and next if $. == 1;

        # Ignore comments or blank lines
        print $line and next if $line =~ /^\s*#/ or $line =~ /^\s*$/;

        # Count valid positions
        $para->{lines}++;

        my ($id,$female,$male, @progeny)=split(/\s+/,$line);
        die "Data format incorrect: $line" 
            unless $id and $female and @progeny;
        
        map{$_ =~ s/,/./g; $_ =~ tr/atgc/ATGC/}($female,$male,@progeny);
        
        my %hash=&parse_reads($female, $male, @progeny);
        
        my $tmp_str = print_hash_str(%hash);
        #print "# $.: No mutant ($tmp_str)\n" and 
        next if keys %hash <= 1;
        
        my ($main,$mutant,@other_nt) = 
            sort {$hash{$b} <=> $hash{$a}} (keys %hash);
            
        #print_hash(%hash);
        #print "# $.: Main: $main $hash{$main}; Mutant: $mutant $hash{$mutant}\n";
        unless(is_single_dose($hash{$main}, $hash{$mutant}, $para)){
            $para->{not_single_dose}++;
            print "# $.: Not single dose ".print_hash_str(%hash)."\n" 
                if $verbose;
            next;
        }else{
            $para->{is_single_dose}++;
            print "# $.: Is single dose ".print_hash_str(%hash)."\n" 
                if $verbose;
        }
        my $type = variation_type($main, $mutant);
       
        my $seg_type = decide_segregation_type($female, $male, $main, $mutant, $para);
        print "# $.: Seg type: $seg_type\n" if $verbose;
        if ($seg_type eq "unexpected"){
            $para->{seg_type_unexpected}++;
            next;
        }

        # SNP
        my $base_info = "$id\t<$type>\t<bases>";
        map{my %tmp = &parse_reads($_);$base_info .= "\t".print_hash_str(%tmp)}(
            $female,$male, @progeny);
        my $gt_info = "$id\t<$type>\t<genotype>";

        for my $progeny_bases ($female, $male, @progeny){
            my $genotype = decide_genotype($seg_type, $progeny_bases, $main, $mutant, $para);
            $gt_info .= "\t$genotype";
        }
        print "$base_info\n$gt_info\n";
    }
    close $in_fh;	
}


############################################################
# Main
############################################################

sub print_stat{
    my $para = shift;
    printf "Number of lines in file: %d\n",
        $para->{lines};
    printf "Number of positions exist mutants: %d\n",
        $para->{not_single_dose} + $para->{is_single_dose};
    printf "Number of positions is single dose: %d\n",
        $para->{is_single_dose};
    print "Segeration type could not be decided: ",
        $para->{seg_type_unexpected},"\n";
}

sub main{
    my $para = read_commands;
    process_pileup($para);
    print_stat($para);
}

main() unless caller;
