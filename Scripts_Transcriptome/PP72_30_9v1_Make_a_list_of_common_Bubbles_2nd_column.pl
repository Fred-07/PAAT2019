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

PP72_30_9v1_Make_a_list_of_common_Bubbles_2nd_column.pl

SYNOPSIS

perl PP72_30_9v1_Make_a_list_of_common_Bubbles_2nd_column.pl    --resultfile=Detail     --fasta=transcripts.fa     --list=Common_elements.txt


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-r / --resultfile             Detailed_final_results_ file
-f / --fasta                  bubble file correpsonding to the Detailed_final_results_ file
-c / --common                 single column file containing the elements (transcript names) to select in the fasta file


DESCRIPTION
Based on the list of transcripts, extract the bubbles that match the transcripts
Report the sequences of the matching bubbles in a new file.
Output filename is deduced from the --list file (the first field between '_' characters)

contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72_30

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $basename;                               # to name the output file
my $count;                                  # count
my $cursor;                                 # cursor
my $fasta;                                  # Fasta file of the transcript sequences
my @interm;                                 # array to store split line
my $IsolateName;                            # Name of the isolate
my $list;                                   # List of elements (transcript names) to select in the fasta file
my $name;                                   # scaffold name
my %record;                                 # hash to record data to be treated
my %record_bubbles;                         # hash to record bubbles
my $resultfile;                             # Detailed_final_results_ file

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
            'resultfile=s'                 => \$resultfile,
            'fasta=s'                      => \$fasta,
            'list=s'                       => \$list);

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

warn $list,"\n";
open (IN1, $list) || die "Can't open $list: $!";

warn $resultfile,"\n";
open (IN2, $resultfile) || die "Can't open $resultfile: $!";

warn $fasta,"\n";
open (IN3, $fasta) || die "Can't open $fasta: $!";


$IsolateName = [split(/\_/, $fasta )] -> [2];
print "Isolate name = $IsolateName\n";

@interm = split(/\_/, $list );
$basename = $interm[0] . "_bubbles__";

#####  Output files  ##########################

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

my $NameOUT3 = $basename . $IsolateName . "__" . $date_file . ".fasta";
open (OUT3, ">$NameOUT3")  || die "Can't open $NameOUT3: $!";



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                             Parsing list file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $line (<IN1>) {
  chomp($line);
  
  if 	($line =~ /^\s*#/) {
    next;
  } elsif ($line =~ /^\s*$/) {
    next;
  } else {
    $record{$line} = "";
  }
}

#warn Dumper(\%record); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                 Parsing Detailed_final_results_ file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$cursor = 0;

 
foreach my $line (<IN2>) {
  chomp($line);
  
  if 	($line =~ /^\s*#/) {
    next;
  } elsif ($line =~ /^\s*$/) {
    next;
  } else {
    @interm = split(/\t/,$line);
    if (exists $record{$interm[0]}) {
      $name = $interm[1];
      $name =~ s/\|/;/g;
      $record_bubbles{$name} = $interm[1];
    }
  }
}



#warn Dumper(\%record_bubbles); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                           Parsing .fasta file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $line (<IN3>) {
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
		@interm = split(/\|/, $line);
		$name=substr($interm[0],1) . ";" . $interm[1] . ";" . $interm[2] . ";" . $interm[3];
		$count ++;
		next;

	} else {
		$seq .= $line;
	}

}

$seqtot{$name} = $seq;

#warn Dumper(\%seqtot); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                        Report final data\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $bubble (nsort keys %record_bubbles) {
  foreach my $long_name (keys  %seqtot) {
    @interm = split(/\;/, $long_name);
    $name= $interm[0] . ";" . $interm[1] . ";" . $interm[2];
    #print "$bubble  -  $long_name  -  $name\n";
    if ($bubble eq $name) {
      #print "1 $bubble $name\n";
      #print "2 ", $record_bubbles{$bubble} , "\n";
      #print "3 ", $seqtot{$long_name}, "\n";
      print OUT3 ">", $long_name , "\n", $seqtot{$long_name}, "\n";
      #print "$bubble -> $long_name\n";
    }
  }
  #print "-" x 30, "\n\n";
}




__END__

  
  



