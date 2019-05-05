#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;

my $version = "1.0";
my $date_of_version = "15-02-2017";


########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72-4-8_Get_stats_on_superVCF_files__Make_synthTable__feat2.pl

SYNOPSIS

perl PP72-4-8_Get_stats_on_superVCF_files__Make_synthTable__feat2.pl       myfile.superVCF


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message


DESCRIPTION
Generate the final table containing breadth of coverage infos and number of poly-allelic sites. Required to calculate density of PA sites
Analyze the .superVCF file


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72 xx-xx-xx-xx
Project ID connection: P25 - Density / Script based on version 2

ENDUSAGE
########################################################################################################

#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage


my @interm;
my @interm2;
my $ID;

my %record;

my %report;

my $Nb_total_pos = 0;
my $coverage;
my $Max_cov = 0;

my $Nb_total_pos_cov10mini = 0;
my $coverage10mini;

#......................................................................................................#


if(@ARGV != 1){
  print "version: $version\n";
  print "date: $date_of_version\n";
  print "$usage\n"; 
  exit(0);
}

if($printVersion){
    print "version $version\n";
    print "date: $date_of_version\n";    
    exit(0);
}

if($help){
  print "version: $version\n";
  print "date: $date_of_version\n";  
  print $usage;
  exit(0);
}



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                             Files\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

#####  Input file  ##########################

my ($superVCF_file) = @ARGV;

open (IN1, $superVCF_file) || die "Can't open $superVCF_file: $!";


#####  sample name  ##########################  

@interm = split(/\./, $superVCF_file);
$ID = $interm[0];

print "Sample name:", $ID, "\n";

#####  Output file name  ########################## 
my $nameOUT1 = "Table_report_" . $ID . ".tmpTB1";
open (OUT1, ">$nameOUT1")  || die "Can't open $nameOUT1: $!";
print "Output TB1 = $nameOUT1\n";

my $nameOUT2 = "Table_report_" . $ID . ".tmpTB2";
open (OUT2, ">$nameOUT2")  || die "Can't open $nameOUT2: $!";
print "Output TB2 = $nameOUT2\n";



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                             Parsing .superVCF file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";


foreach my $line (<IN1>) {
    chomp($line);
    
    @interm = split(/\t/, $line);
    
    # discard comment line 
    if ($line =~ /^\s*#/) {
	    next;			
	    
    # discard blank line
    } elsif ($line =~ /^\s*$/) {
	    next;	
	    
    # record 
    } else {

	#		scaff	  position	rad
	$record{$interm[1]}{$interm[3]}{$interm[2]}{seqID} = $interm[4];
	if ($interm[6] >0) {
	    $record{$interm[1]}{$interm[3]}{$interm[2]}{rep} = 1;
	} else {
	    $record{$interm[1]}{$interm[3]}{$interm[2]}{rep} = 0;
	}
	$record{$interm[1]}{$interm[3]}{$interm[2]}{cod} = $interm[7];
	$record{$interm[1]}{$interm[3]}{$interm[2]}{ref} = $interm[8];
	$record{$interm[1]}{$interm[3]}{$interm[2]}{type} = $interm[9];
	$record{$interm[1]}{$interm[3]}{$interm[2]}{alleles} = $interm[10];
	$record{$interm[1]}{$interm[3]}{$interm[2]}{occ} = $interm[11];
	$record{$interm[1]}{$interm[3]}{$interm[2]}{covgen} = $interm[12];
	
	#warn Dumper(\%record); getc();
    }
}


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                               Generate and write general table .tmpTB1:\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";


%report = ();

foreach my $scaffold (nsort  (keys %record)) {
    foreach my $position (sort {$a <=> $b; } (keys %{$record{$scaffold}})) {
	foreach my $RAD (keys %{$record{$scaffold}{$position}}) {
	    
	    #General coverage
	    $Nb_total_pos ++;
	    $coverage += $record{$scaffold}{$position}{$RAD}{covgen};
	    if ($Max_cov < $record{$scaffold}{$position}{$RAD}{covgen}) {
		$Max_cov = $record{$scaffold}{$position}{$RAD}{covgen};
	    }
	    
	    if ($record{$scaffold}{$position}{$RAD}{covgen} >= 10) {
		$Nb_total_pos_cov10mini ++;
		$coverage10mini += $record{$scaffold}{$position}{$RAD}{covgen};
	    }
	    

	    #repeats
	    if ($record{$scaffold}{$position}{$RAD}{rep} == 0 ) {
		$report{rep}{0}{cov} += $record{$scaffold}{$position}{$RAD}{covgen};
		$report{rep}{0}{pos} ++;
		
		if ($record{$scaffold}{$position}{$RAD}{covgen} >= 10) {
		    $report{rep}{0}{cov10mini} += $record{$scaffold}{$position}{$RAD}{covgen};
		    $report{rep}{0}{pos10mini} ++;
		}
		
	    } else {
		$report{rep}{1}{cov} += $record{$scaffold}{$position}{$RAD}{covgen};
		$report{rep}{1}{pos} ++;
		
		if ($record{$scaffold}{$position}{$RAD}{covgen} >= 10) {
		    $report{rep}{1}{cov10mini} += $record{$scaffold}{$position}{$RAD}{covgen};
		    $report{rep}{1}{pos10mini} ++;
		}
	    }

	    #coding
	    if ($record{$scaffold}{$position}{$RAD}{cod} == 0 ) {
		$report{cod}{0}{cov} += $record{$scaffold}{$position}{$RAD}{covgen};
		$report{cod}{0}{pos} ++;
		
		if ($record{$scaffold}{$position}{$RAD}{covgen} >= 10) {
		    $report{cod}{0}{cov10mini} += $record{$scaffold}{$position}{$RAD}{covgen};
		    $report{cod}{0}{pos10mini} ++;
		}
		
	    } else {
		$report{cod}{1}{cov} += $record{$scaffold}{$position}{$RAD}{covgen};
		$report{cod}{1}{pos} ++;
		
		if ($record{$scaffold}{$position}{$RAD}{covgen} >= 10) {
		    $report{cod}{1}{cov10mini} += $record{$scaffold}{$position}{$RAD}{covgen};
		    $report{cod}{1}{pos10mini} ++;
		}
	    }	    
	    
	    
	    #variants
	    if ($record{$scaffold}{$position}{$RAD}{type} ne "" ) {
		$report{site_with_variants} ++ ;
	    }	    
	    
	    
	    
	}
    }
}


#write Table
print OUT1 "\tMean_coverage\tMean_coverage_10mini\tMax_coverage\tTotal_number_common_positions\tNumber_of_covered_positions(cov>=10)\t%sites_covered\tMean_coverage_in_repeats\tMean_coverage_in_non_repeats\tMean_coverage_in_coding\tMean_coverage_in_non_coding\tTotal_number_of_variable_sites\n";

print OUT1 $ID,"\t";
print OUT1 sprintf("%.2f", $coverage/$Nb_total_pos ), "\t";
print OUT1 sprintf("%.2f", $coverage10mini/$Nb_total_pos_cov10mini ), "\t";
print OUT1 $Max_cov , "\t";
print OUT1 $Nb_total_pos, "\t", $Nb_total_pos_cov10mini , "\t", sprintf("%.2f", $Nb_total_pos_cov10mini/$Nb_total_pos*100) , "%\t";
print OUT1 sprintf("%.2f", $report{rep}{1}{cov10mini}/$report{rep}{1}{pos10mini}), "\t", sprintf("%.2f", $report{rep}{0}{cov10mini}/$report{rep}{0}{pos10mini}), "\t";
print OUT1 sprintf("%.2f", $report{cod}{1}{cov10mini}/$report{cod}{1}{pos10mini}), "\t", sprintf("%.2f", $report{cod}{0}{cov10mini}/$report{cod}{0}{pos10mini}), "\t";
print OUT1 $report{site_with_variants},"\n";



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                       Generate and write polyallelic sites and coverage table .tmpTB2:\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";



%report = ();

my $total_position;
my $total_NBpolyallelesAll;
my $total_NBpolyallelesSNP;
my $total_NBpolyallelesInDel;

foreach my $scaffold (nsort  (keys %record)) {
    foreach my $position (sort {$a <=> $b; } (keys %{$record{$scaffold}})) {
	foreach my $RAD (keys %{$record{$scaffold}{$position}}) {
	    
	    $report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{position} ++;
	    $total_position ++;
	    

	    $report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{coverage_glob} += $record{$scaffold}{$position}{$RAD}{covgen};
	    $report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{coverage_glob_nb_pos} ++;
		

	    if ($record{$scaffold}{$position}{$RAD}{covgen} >= 10) {
    
		$report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{coverage10mini} += $record{$scaffold}{$position}{$RAD}{covgen};
		$report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{coverage10mini_nb_pos} ++;
	    }
	    
	    if ($record{$scaffold}{$position}{$RAD}{type} =~ m/,/) {
		$report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{NBpolyallelesAll} += 1;
		$total_NBpolyallelesAll ++;
		
		if (($record{$scaffold}{$position}{$RAD}{type} !~ m/I/) and ($record{$scaffold}{$position}{$RAD}{type} !~ m/D/)) {
		    $report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{NBpolyallelesSNP} += 1;
		    $total_NBpolyallelesSNP ++;		    
		} 
		if ($record{$scaffold}{$position}{$RAD}{type} !~ m/S/) {									# elsif ou if ???
		    $report{$record{$scaffold}{$position}{$RAD}{rep}}{$record{$scaffold}{$position}{$RAD}{cod}}{NBpolyallelesInDel} += 1;
		    $total_NBpolyallelesInDel ++;	    
		} 
	    }
	}
    }
}

	


print OUT2 "Isolate\tRepeats\tCoding\tObserved_Pos\t%Observed_Pos\tMean_Cov\tMean_cov(>=10)\tNB_of_all_polyalleles\tNB_of_polyallelic_SNPs\tNB_of_polyallelic_InDels\n";

foreach my $Repeat (nsort keys %report) {
    foreach my $coding (nsort keys %{$report{$Repeat}}) {	
	    print OUT2 $ID,"\t";
	    print OUT2 "$Repeat\t$coding\t", $report{$Repeat}{$coding}{position},"\t", sprintf("%.0f", $report{$Repeat}{$coding}{position}/$total_position*100),"%\t";
	    print OUT2 sprintf("%.1f", $report{$Repeat}{$coding}{coverage_glob}/$report{$Repeat}{$coding}{coverage_glob_nb_pos}),"\t"; 
	    print OUT2 sprintf("%.1f", $report{$Repeat}{$coding}{coverage10mini}/$report{$Repeat}{$coding}{coverage10mini_nb_pos}),"\t"; 
	    print OUT2 $report{$Repeat}{$coding}{NBpolyallelesAll},"\t";
	    print OUT2 $report{$Repeat}{$coding}{NBpolyallelesSNP},"\t";
	    print OUT2 $report{$Repeat}{$coding}{NBpolyallelesInDel},"\n";	    
	
    }
    
}


print OUT2 "\t\t\t\t", $total_position ,"\t\t\t\t", $total_NBpolyallelesAll,"\t", $total_NBpolyallelesSNP,"\t",  $total_NBpolyallelesInDel,"\n";



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                                 Completed!!\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";




__END__


