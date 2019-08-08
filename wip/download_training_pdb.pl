#!/usr/bin/env perl

use strict;
use warnings;

my @input_list = `cat $ARGV[0]`;

my %hash;

foreach my $line (@input_list){
    $hash{substr($line, 0, 4)} = 1;
}

foreach my $key (keys %hash){
    system "wget -q https://files.rcsb.org/download/$key.pdb.gz";
    unless (-e "$key.pdb.gz" and -f "$key.pdb.gz"){
        print "FAIL: $key\n";
    }
}
