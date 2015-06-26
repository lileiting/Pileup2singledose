#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;

sub usage{
    print <<USAGE;

$FindBin::Script FILE [FILE ...]

USAGE
    exit;
}

sub clean_data{
    my $file = shift;
    open my $in_fh, "<", $file or die "$file: $!";
    my ($name, $str) = ('', '');
    my $status;
    while(<$in_fh>){
        s/[\t\r\n]//g;
        ($name) = ($_ =~ /^(\S+)/) and next if $. == 1;
        next if /^\s*#/ or /^\s*$/;
        my $length = length($_);
        if($length >= 60 + 13){
             $status = $length - 60;
             $str .= substr($_, $status, 60);
        }elsif($status){
             $str .= substr($_, $status, $length - $status);
        }else{die "CAUTION: Status missing! File: $file, line $.\n"}     
    }
    print "$file,$name,$str\n";
    close $in_fh;
}

sub main{
    usage if @ARGV == 0;
    my @files = @ARGV;
    for my $file (@files){
        clean_data($file);
    }
}

main() unless caller;
