#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Sortkeys  = 1;
use Sort::Naturally;
use Getopt::Long;

my $version = "1.0";
my $date = "05-05-2017";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_30_8v1_Make_a_list_of_common_transcripts_1st_column.pl

SYNOPSIS

perl PP72_30_8v1_Make_a_list_of_common_transcripts_1st_column.pl  --pattern=Detail --fasta=transcripts.fa


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-p / --pattern                Pattern to recognize the results files
-f / --fasta                  Transcript file


DESCRIPTION
Read automatically all files with the pattern in the folder
Extract contents of the first column
Report all elements found in more than one file in an output file
Report all unique elements in a second output file
Generate additional files with the selected transcripts

contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72_30

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $count;                                  # count
my $cursor;                                 # cursor
my $fasta;                                  # Fasta file of the transcript sequences
my @interm;                                 # array to store split line
my $IsolateName;                            # Name of the isolate
my $pattern;                                # pattern to recognize the results files
my %record;                                 # hash to record data to be treated
my %record_per_file;                        # hash to record data to be treated
my $name;                                   # scaffold name
my $seq = "";                               # to store sequence
my %seqtot;                                 # hash to record sequences

#......................................................................................................#


if(@ARGV==0){
  print "version: $version\n";
  print "date: $date\n";
  print "$usage\n"; 
  exit(0);
}


GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion,
            'pattern=s'                    => \$pattern,
            'fasta=s'                      => \$fasta);

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
warn"                                             Files\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

#####  Input files  ##########################

print "Pattern: $pattern\n";
my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { /^$pattern/ } readdir(DIR);
closedir DIR;

$count = 0;
warn "Files\n";
for my $i (@files) {
    warn $i,"\n";
    $count ++;
}
warn "\n";
if ($count == 0) {
    die "no *$pattern* file found ------\n\n";
} else {
    print ">> $count files\n";    
}

warn $fasta,"\n";
open (IN2, $fasta) || die "Can't open $fasta: $!";

$IsolateName = [split(/\./, $fasta )] -> [0];

#####  Output files  ##########################

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

my $NameOUT1 = "Common_elements__" . $IsolateName . "__" . $date_file . ".txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";
my $NameOUT2 = "Unique_elements__" . $IsolateName . "__" . $date_file . ".txt";
open (OUT2, ">$NameOUT2")  || die "Can't open $NameOUT2: $!";

my $NameOUT3 = "Common_transcripts__" . $IsolateName . "__" . $date_file . ".fasta";
open (OUT3, ">$NameOUT3")  || die "Can't open $NameOUT3: $!";
my $NameOUT4 = "Unique_transcripts__" . $IsolateName . "__" . $date_file . ".fasta";
open (OUT4, ">$NameOUT4")  || die "Can't open $NameOUT4: $!";


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                         Parsing files\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$cursor = 0;

foreach my $f (nsort @files) {
  open (IN1, "<$f");
  warn " treatment of $f\n";
  $cursor ++;
  %record_per_file = ();
  
  foreach my $line (<IN1>) {
    chomp($line);
    
    if 	($line =~ /^\s*#/) {
      next;
    } elsif ($line =~ /^\s*$/) {
      next;
    } else {
      @interm = split(/\t/,$line);
      $record_per_file{$interm[0]} ++;
    }
  }
  foreach my $elem ( keys %record_per_file) {
    $record{$elem} ++;
  }
}


#warn Dumper(\%record); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                           Parsing .fasta file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $line (<IN2>) {
    chomp($line);
	if 	($line =~ /^\s*#/) {
		next;			
	} elsif ($line =~ /^\s*$/) {
		next;	
	} elsif ($line =~ /^>/) {

		if (length($seq) >0){
			$seqtot{$name} = $seq;
		}
		$seq="";
		$line = [split(" ", $line)] -> [0];
		$name=substr($line,1);
		$count ++;
		next;

	} else {
		$seq .= $line;
	}

}

$seqtot{$name} = $seq;

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                        Report final data\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $elem (nsort keys %record) {
  if ($record{$elem} == 1) {
    print OUT2 $elem, "\n";
    print OUT4 ">", $elem, "\n", $seqtot{$elem}, "\n";
  } else {
    print OUT1 $elem, "\n";
    print OUT3 ">", $elem, "\n", $seqtot{$elem}, "\n";
  }
  
}



__END__

  
  



