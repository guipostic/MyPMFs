#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $input_file_path;
my $rmsd_files_list;
my $rmsd_threshold = 0;
my $help = ($#ARGV < 0) ? 1 : 0;
my $write_correlation_files;
my $forbidden = "";

GetOptions(
    'input=s' => \$input_file_path,
    'rmsd=s' => \$rmsd_files_list,
    'write' => \$write_correlation_files,
    'threshold=s' => \$rmsd_threshold,
    'xclude_native' => \$forbidden,
    'help' => \$help
);
helper() if ($help);

if (!$input_file_path){
    die "Error: missing -i arg\n";
}
if (!$rmsd_files_list){
    die "Error: missing -r arg\n";
}

if ($forbidden){
    $forbidden = 'native.pdb';
}

my @input_file = `cat $input_file_path`;
my @rmsd_files = `cat $rmsd_files_list`;

# ANALYZE DECOYS RMSD
my %hashdecoy_rmsd;
foreach my $file (@rmsd_files){
    chomp $file;
    my $pdb = "";
    if ($file =~ m/\/home\/\w+\/GITHUB\/3DRobot_set\/(\w+)\/list.txt/){
        $pdb = $1;
    }
    my @rmsd_file = `cat $file`;
    shift @rmsd_file;
    foreach my $line (@rmsd_file){
	chomp $line;
	my ($decoy, $rmsd) = split(/\s+/, $line);
        $hashdecoy_rmsd{"$pdb-$decoy"} = $rmsd;
    }
}

# ANALYZE DECOYS ENERGIES
my %hash5char_decoys;
my %hash5char_energy;

my $pdb5char = "";
my $energy_found = 1;

my $total_energy = 0;
my $total_models = 0;

foreach my $line (@input_file){
    chomp $line;
    if ($line =~ m/^Now processing .+(\d\w\w\w\w)\/(\w+\.pdb)/){
        if ($energy_found){
            $pdb5char = $1;
            my $decoy = $2;
            if (exists $hash5char_decoys{$pdb5char}){
                push(@{$hash5char_decoys{$pdb5char}}, $decoy);
            }
            else{
                $hash5char_decoys{$pdb5char} = [];
                push(@{$hash5char_decoys{$pdb5char}}, $decoy);
    	    }
            $energy_found = 0;
        }
        else{
            die "Error at line: $line\n";
        }
    }
    if ($line =~ m/^Pseudo-energy = (.+)/){
        if ($pdb5char){
            my $energy = $1;
            if (exists $hash5char_energy{$pdb5char}){
                push(@{$hash5char_energy{$pdb5char}}, $energy);
            }
            else{
                $hash5char_energy{$pdb5char} = [];
                push(@{$hash5char_energy{$pdb5char}}, $energy);
    	    }
            $pdb5char = "";
            $energy_found = 1;
            $total_energy+=$energy;
            $total_models++;
        }
        else{
            die "Error at line: $line\n";
        }
    }
}

my $avg_energy = $total_energy/$total_models;
print "AVG energy = $avg_energy\n";

my $correct_native = 0;
my $incorrect_native = 0;

my $total_accuracy = 0;
my $overall_delta = 0;

# FOR EACH PROTEIN
foreach my $key (sort keys %hash5char_decoys){
    my @decoys = @{$hash5char_decoys{$key}};
    my @energies = @{$hash5char_energy{$key}};

    if ($write_correlation_files){
        open(OUT, '>', "$key.cor") or die "Error: cannot write $key.cor\n\t$!\n";
        for (my $i=0; $i<=$#decoys; $i++){
            my $rmsd = $hashdecoy_rmsd{"$key-$decoys[$i]"};
            print OUT "$decoys[$i] $energies[$i] $rmsd\n";
        }
        close OUT;
    }

    # CHECK RANKING OF NATIVES
    my $index = -1;
    my $lowest_energy = 1000000;
    for (my $i=0; $i<=$#decoys; $i++){
        if ($energies[$i] < $lowest_energy){
            $lowest_energy = $energies[$i];
            $index = $i;
        }
    }

    my $lowest_rmsd = $hashdecoy_rmsd{"$key-$decoys[$index]"};
    #print "$key lowest : $decoys[$index] (E=$energies[$index]; RMSD=$lowest_rmsd)\n";

    if ($lowest_rmsd == 0){
        $correct_native++;
    }
    else{
        $incorrect_native++;
    }

    # COMPUTE ACCURACY
    my $correct = 0;
    my $incorrect = 0;
    my $total_delta = 0;

    for (my $i=0; $i<=$#decoys; $i++){
        for (my $j=$i+1; $j<=$#decoys; $j++){
            if ($j<=$#decoys){
                my $rmsdi = $hashdecoy_rmsd{"$key-$decoys[$i]"};
                my $rmsdj = $hashdecoy_rmsd{"$key-$decoys[$j]"};
                my $delta = $rmsdi - $rmsdj;
                if ($delta < 0){
                    $delta*=-1;
                }
		if ($delta > $rmsd_threshold and $decoys[$i] ne $forbidden and $decoys[$j] ne $forbidden){
                    if ($rmsdi < $rmsdj and $energies[$i] < $energies[$j]){
                        $correct++;
                    }
                    elsif ($rmsdi > $rmsdj and $energies[$i] > $energies[$j]){
                        $correct++;
                    }
                    else{
                        $incorrect++;
                        $total_delta+=$delta;
                    }
                }
            }
        }
    }
    my $accuracy = ($correct/($correct+$incorrect))*100;
    #print "$key: $accuracy%\n";

    my $avg_delta = $total_delta/$incorrect;
    #print "$key: $avg_delta A\n";

    $total_accuracy+=$accuracy;
    $overall_delta+=$avg_delta;
}

my $total_native = $correct_native+$incorrect_native;
my $accuracy_native = ($correct_native/($total_native))*100;
print "Native found = $accuracy_native% ($correct_native found; $incorrect_native not found)\n";

my $avg_accuracy = $total_accuracy/$total_native;
print "Accuracy = $avg_accuracy%\n";

my $delta_rmsd = $overall_delta/$total_native;
print "DeltaRMSD = $delta_rmsd\n";

#sprintf("%.2f", 1.555);


sub helper {
    print
      "\n",
      "\tOptions:\n",
      "\t--input            TSV file; column 0: UniProt AC; column 1: UniProt ID;\n",
      "\t--pdbhtml          The PDBhtml database (CSV file)\n",
      "\t--exclude_natives  For the input proteins that have at least one PDB structure: do not search into the HHsearch results\n",
      "\t--template_shared  Allow N input proteins to share the same template (default = 1)\n",
      "\t--results_hh       TSV file representing the HHsearch results\n",
      "\t--smr              Use the SWISS-MODEL REPOSITORY INDEX file, instead of the HHsearch results\n",
      "\t--only_pdb         Do not search into the HHsearch results\n",
      "\t--no_filter        Do not filter the redundancy within the complexes\n",
      "\t--max_title        Maximum character length for the title (default = 300)\n",
      "\t--help             This help\n",
      "\n";
    exit 1;
}

