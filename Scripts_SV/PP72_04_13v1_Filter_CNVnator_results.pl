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

PP72_04_13v1_Filter_CNVnator_results.pl

SYNOPSIS

perl PP72_04_13v1_Filter_CNVnator_results.pl  --file=Output_CNVnator   --GenomeFile=myGenome.fa    --Size=10000


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-f / --file                   Output file from CNVnator
-g / --GenomeFile             Genome .fasta file
-S / --Size                   Size threshold (bp) for the scaffolds



DESCRIPTION
Filter CNVnator output files (like output_CNVnator_Nu6-WG1_on_Rhiir2.only2mapped.txt)
Filters are based on scaffold size, p-value n¡1 (< 0.05), q0 (< 0.5)


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
my $file;                                   # input file = Output file from CNVnator
my $GenomeFile;                             # Genome file
my $name;                                   # scaffold name
my $NameOUT1;                               # Name of the output file
my $pval1;                                  # CNVnator first p-value
my $q0;                                     # CNVnator q0 value
my $RD;                                     # CNVnator RD value
my %record;                                 # main hash to record data
my $scaffold_name;                          # Current scaffold name
my $seq = "";                               # to store sequence
my $Size;                                   # Size threshold to keep scaffold
my $temp;








#......................................................................................................#


if(@ARGV==0 or @ARGV!=3){
  print "version: $version\n";
  print "date: $date\n";
  print "$usage\n"; 
  exit(0);
}


GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion,
            'file=s'                       => \$file,
            'GenomeFile=s'                 => \$GenomeFile,
            'Size=s'                       => \$Size);


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

print "CNVnator file = ",$file,"\n";
open (IN2, $file) || die "Can't open $file: $!";



#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

$temp = "_filtered_" . $date_file . "_.txt";
$NameOUT1 = $file;
$NameOUT1 =~ s/\.txt/$temp/;
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
warn"                              Parsing and filtering output file from CNVnator\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $line (<IN2>) {
    chomp($line);
	
	if 	($line =~ /^\s*#/) {
		next;			
	} elsif ($line =~ /^\s*$/) {
		next;	
	} else {
    
    @interm = split(/\t/, $line);
    $scaffold_name = [split(/\:/, $interm[1])] -> [0];
    $pval1 = $interm[4];
    $q0 = $interm[8];
    $RD = $interm[3];
    
    # Filtering steps
    $cursor = "report";
    if ($record{$scaffold_name} < $Size ) {
      $cursor = "NO_report";
    } elsif ($pval1 >= 0.05) {
      $cursor = "NO_report";
    } elsif ($q0 >= 0.5) {
      $cursor = "NO_report";
    } elsif ($RD == 0) {
      $cursor = "NO_report";
    }
   
    
    if ($cursor eq "report") {
      print OUT1 $line,"\n";
    }
  }
}

__END__






