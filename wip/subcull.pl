#!/usr/bin/env perl

use strict;
use warnings;

my @input_list = `cat $ARGV[0]`;

foreach my $line (@input_list){
    my @A = split(/\s+/, $line);
    if ($A[1] >= 80 and $A[1] <= 250){
        print $line;
    }
}

