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
my $Filename;

my $count;

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

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);
my $NameOUT1 = $Filename . "_" . $type . "_" . $date . ".ForVenn.txt";
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

    
    if ($interm[9] =~ m/,/) {
      
      if ($type eq "norep") {
        if ($interm[6] == 0) {
          print OUT1 $scaffold , "_", $position, "\n";
          print $line, "\n";
        }
      } elsif ($type eq "norepcod") {
        if ($interm[6] == 0 and $interm[7] == 1) {
          print OUT1 $scaffold , "_", $position, "\n";
        }        
      } elsif ($type eq "norep-noInDels") {
        if ($interm[6] == 0 and $interm[9] !~ m/I/ and $interm[9] !~ m/D/) {
          print OUT1 $scaffold , "_", $position, "\n";
        }
      } elsif ($type eq "norepcod-noInDels") {
        if ($interm[6] == 0 and $interm[7] == 1 and $interm[9] !~ m/I/ and $interm[9] !~ m/D/) {
          print OUT1 $scaffold , "_", $position, "\n";
        }
      } elsif ($type eq "all-but-noInDels") {
        if ($interm[9] !~ m/I/ and $interm[9] !~ m/D/) {
          print OUT1 $scaffold , "_", $position, "\n";
        }        
      } else {
        print OUT1 $scaffold , "_", $position, "\n";
      }
    }
  }
    
}


__END__