#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use Getopt::Long;

my $version = "1.0";
my $date = "10-03-2017";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_04_16v1_Filter_LUMPY_results.pl

SYNOPSIS

perl PP72_04_16v1_Filter_LUMPY_results.pl  --file=Lumpy.vcf  --GenomeFile=myGenome.fa    --Size=10000     --SU_thresh=10  --PE_thresh=10  --SR_thresh=10


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-f / --file                   Output file from LUMPY
-g / --GenomeFile             Genome .fasta file
-S / --Size                   Size threshold (bp) for the scaffolds
     --SU_thresh              Threshold for the number of pieces of evidence supporting the variant
     --PE_thresh              Threshold for the number of paired-end reads supporting the variant
     --SR_thresh              Threshold for the number of split reads of evidence supporting the variant


DESCRIPTION
Filter LUMPY output files (like DNA2-PE-1_on_DNA2.vcf)
Filters are based on scaffold size, SU, PE and SR


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $count_filt;                             # counter for retained scaffolds
my $count_tot;                              # counter for total nb scaffolds
my $cursor;                                 # cursor to record status 0 or 1
my @interm;                                 # array to store split line
my @interm2;                                # array to store split elements
my $file;                                   # input file from LUMPY
my $GenomeFile;                             # Genome file
my $name;                                   # scaffold name
my $NameOUT1;                               # Name of the output file
my $PE;                                     # Number of paired-end reads supporting the variant
my $PE_thresh;                              # Threshold for the number of paired-end reads supporting the variant
my %record;                                 # main hash to record data
my %report_counting;                        # hash to record counts
my $scaffold_name;                          # Current scaffold name
my $seq = "";                               # to store sequence
my $Size;                                   # Size threshold to keep scaffold
my $SR;                                     # Number of split reads of evidence supporting the variant
my $SR_thresh;                              # Threshold for the number of split reads of evidence supporting the variant
my $SU;                                     # Number of pieces of evidence supporting the variant
my $SU_thresh;                              # Threshold for the number of pieces of evidence supporting the variant
my $temp;








#......................................................................................................#


if(@ARGV==0 or @ARGV!=6){
  print "version: $version\n";
  print "date: $date\n";
  print "$usage\n"; 
  exit(0);
}


GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion,
            'file=s'                       => \$file,
            'GenomeFile=s'                 => \$GenomeFile,
            'Size=s'                       => \$Size,
            'SU_thresh=s'                  => \$SU_thresh,
            'PE_thresh=s'                  => \$PE_thresh,
            'SR_thresh=s'                  => \$SR_thresh);


if($printVersion){
    print "version $version\n";
    print "date: $date\n";    
    exit(0);
}

if($help){
  print "version: $version\n";
  print "date: $date\n";  
  print $usage;
  exit(0);
}





warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                               Files\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";


print "Genome file = ",$GenomeFile,"\n";
open (IN1, $GenomeFile) || die "Can't open $GenomeFile: $!";

print "LUMPY file = ",$file,"\n";
open (IN2, $file) || die "Can't open $file: $!";



#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

$temp = "_filtered_" . $date_file . ".vcf";
$NameOUT1 = $file;
$NameOUT1 =~ s/\.vcf/$temp/;
print "Output file = ", $NameOUT1, "\n";

open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                       Parsing genome .fasta file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$seq ="";

foreach my $line (<IN1>) {
    chomp($line);
	
	if 	($line =~ /^\s*#/) {
		next;			
	} elsif ($line =~ /^\s*$/) {
		next;	
	} elsif ($line =~ /^>/) {

		if (length($seq) > 0){
      $record{$name} = length($seq);
		}
		$seq="";
		$name=$line;
		$name =~ s/\>//;
		$count_tot ++;
		
		next;
	} else {
		$seq .= $line;
	}
}

if (length($seq) > 0){
    $record{$name} = length($seq);
}


print "Total nb of sequences: $count_tot\n";
print "**************************************","\n";


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                              Parsing and filtering output file from LUMPY\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $line (<IN2>) {
    chomp($line);
	
	if 	($line =~ /^\s*#/) {
    print OUT1 $line,"\n";		
	} elsif ($line =~ /^\s*$/) {
		next;	
	} else {
    
    @interm = split(/\t/, $line);
    $scaffold_name = $interm[0];
    
    @interm2 = split(/\;/, $interm[7]);
    foreach my $t (@interm2) {
      if ($t =~/SU=/) {
        $SU = [split(/\=/, $t)] -> [1];
      } elsif ($t =~/PE=/) {
        $PE = [split(/\=/, $t)] -> [1];
      } elsif   ($t =~/SR=/) {
        $SR = [split(/\=/, $t)] -> [1];
      }
    }
    
    # Filtering steps
    $cursor = "report";
    if ($record{$scaffold_name} < $Size ) {
      $cursor = "NO_report";
      $report_counting{Size} ++;
    } elsif ($SU <= $SU_thresh) {
      $cursor = "NO_report";
      $report_counting{NBevents} ++;      
    } elsif ($PE <= $PE_thresh or $SR <= $PE_thresh) {
      $cursor = "NO_report";
      $report_counting{OneIsZero} ++      
    }
    
    if ($cursor eq "report") {
      print OUT1 $line,"\n";
    }
  }
}

foreach my $t (nsort keys %report_counting) {
  print $t," : ", $report_counting{$t}, "\n";
}

__END__






