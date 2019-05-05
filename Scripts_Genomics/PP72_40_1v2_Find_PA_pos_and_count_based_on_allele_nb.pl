#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;


my $version = "1.0";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_04_02_1v1_Find_PA_pos_and_make_file_for_Venn.pl

SYNOPSIS

perl PP72_04_02_1v1_Find_PA_pos_and_make_file_for_Venn.pl      ***.subVCF     type   


OPTIONS
'type' indicates the options and it can be:
  all                 > all positions are included
  norep               > repeat positions are removed
  norepcod            > only non repeat and coding positions are reported
  norep-noInDels      > repeat positions and indels-containing positions are removed
  norepcod-noInDels   > repeat positions, non-coding positions and indels-containing positions are removed
  all-but-noInDels    > all positions are included but indels-containing positions are removed
  

DESCRIPTION
Analyze .subVCF files and report the poly-allelic sites under the format "scaffoldxxxx_position"
Designed for venn diagram generation


contact:
frederic.masclaux@unil.ch
University of Lausanne

DATE
14/10/2016

ENDUSAGE
########################################################################################################


my $numArgs = $#ARGV + 1;

if($numArgs != 2){
  print "version: $version\n";
  print "$usage\n"; 
  exit(0);
}

#..............Variables...............................................................................#

my @int1;
my @int2;
my $Filename;
my $SampleName;

my $count;
my $count_1allele = 0;
my $count_2alleles = 0;
my $count_3alleles = 0;
my $count_4alleles = 0;
my $count_morethan4alleles = 0;

my @interm;
my %record;

my $scaffold; 
my $position;
my $RAD;
my $REFallele;
my @TypeAllele;
my @ListAllele;


#......................................................................................................#


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                            Manage the files\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

my ($FileIN, $type) = @ARGV;

if ($type ne "all" and $type ne "norepcod" and $type ne "norep" and $type ne "norep-noInDels"
    and $type ne "norepcod-noInDels" and $type ne "all-but-noInDels") {
  die "type is not correct: $type\n$usage\n";
}

#input
open (IN1, $FileIN) || die "Can't open $FileIN: $!";
print "open file: $FileIN\n";
print "\n";

#sample name
@int1 = split(/\./,$FileIN);
$Filename = $int1[0];
@int2 = split(/\-/,$FileIN);
$SampleName = $int2[0] . "-" . $int2[1];

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);
my $NameOUT1 = $Filename . "_" . $type . "_" . $date . ".Count-PA-types.txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";
print "Output file: $NameOUT1\n";

print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print">>>>>  Parsing the file $FileIN\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";


foreach my $line (<IN1>) {
  chomp($line);
  
  if 	($line =~ /^\s*#/) {
    next;
      
      
  } elsif ($line =~ /^\s*$/) {
    next;
      
  } else {

    @interm = split(/\t/,$line);
   
    #Scaffold and position
    $scaffold = $interm[1];
    $RAD = $interm[2];
    $position = $interm[3];
    $REFallele = $interm[8];
    @TypeAllele = split(/\,/, $interm[9]);
    @ListAllele = split(/\,/, $interm[10]);

    
    if ($type eq "norep") {
      if ($interm[6] == 0 and $interm[9] =~ m/,/) {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        } elsif (@TypeAllele == 2) {
          $count_2alleles ++;
        } elsif (@TypeAllele == 3) {
          $count_3alleles ++;
        } elsif (@TypeAllele == 4) {
          $count_4alleles ++;
        } elsif (@TypeAllele > 4) {
          $count_morethan4alleles ++;
        }
      } else {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        }          
      }
    } elsif ($type eq "norepcod") {
      if ($interm[6] == 0 and $interm[7] == 1 and $interm[9] =~ m/,/) {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        } elsif (@TypeAllele == 2) {
          $count_2alleles ++;
        } elsif (@TypeAllele == 3) {
          $count_3alleles ++;
        } elsif (@TypeAllele == 4) {
          $count_4alleles ++;
        } elsif (@TypeAllele > 4) {
          $count_morethan4alleles ++;
        }
      } else {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        }          
      }
    } elsif ($type eq "norep-noInDels") {
      if ($interm[6] == 0 and $interm[9] !~ m/I/ and $interm[9] !~ m/D/ and $interm[9] =~ m/,/) {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        } elsif (@TypeAllele == 2) {
          $count_2alleles ++;
        } elsif (@TypeAllele == 3) {
          $count_3alleles ++;
        } elsif (@TypeAllele == 4) {
          $count_4alleles ++;
        } elsif (@TypeAllele > 4) {
          $count_morethan4alleles ++;
        }
      } else {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        }          
      }
    } elsif ($type eq "norepcod-noInDels") {
      if ($interm[6] == 0 and $interm[7] == 1 and $interm[9] !~ m/I/ and $interm[9] !~ m/D/ and $interm[9] =~ m/,/) {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        } elsif (@TypeAllele == 2) {
          $count_2alleles ++;
        } elsif (@TypeAllele == 3) {
          $count_3alleles ++;
        } elsif (@TypeAllele == 4) {
          $count_4alleles ++;
        } elsif (@TypeAllele > 4) {
          $count_morethan4alleles ++;
        }
      } else {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        }          
      }
    } elsif ($type eq "all-but-noInDels") {
      if ($interm[9] !~ m/I/ and $interm[9] !~ m/D/ and $interm[9] =~ m/,/) {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        } elsif (@TypeAllele == 2) {
          $count_2alleles ++;
        } elsif (@TypeAllele == 3) {
          $count_3alleles ++;
        } elsif (@TypeAllele == 4) {
          $count_4alleles ++;
        } elsif (@TypeAllele > 4) {
          $count_morethan4alleles ++;
        }
      } else {
        if (@TypeAllele == 1) {
          $count_1allele ++;
        }          
      }
    } else {
      print $scaffold , "_", $position, "\n";
    }

  }  
}

print OUT1 "\t", "1_allele\t2_alleles\t3_alleles\t4_alleles\t>4_alleles" ,"\n";

print OUT1 $SampleName ,"\t", $count_1allele, "\t", $count_2alleles,  "\t", $count_3alleles, "\t", $count_4alleles, "\t", $count_morethan4alleles,"\n";


__END__