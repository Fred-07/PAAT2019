#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use Getopt::Long;
use List::Util qw[min max];

my $version = "1.0";
my $date_of_version = "13-02-2017";


########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72-XXX_Get_breadth_coverage_for_WG_data.pl

SYNOPSIS

perl PP72-XXX_Get_breadth_coverage_for_WG_data.pl   --cluster=Novoalignmulticopy.gff    --repeat=Repeats.gff     --coding=AugustusPrediction.gff
																										--genome=ReferenceGenome.fa					--coverage=Sample.coverage



OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
--cluster                     'Cluster' file. Can be gff3 or xxxx
--repeat                      'Repeats' file generated with RepeatMasker. .gff3 format.
--coding                      'Gene prediction' file generated with Augustus. .gff3 format
--genome											Reference genome in fasta format
--coverage										Coverage for the sample as calculated with samtools depth (Sample.coverage)

DESCRIPTION
Calculate the breadth of coverage for whole genome data. Required to calculate density of PA sites

contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72 xx-xx-xx-xx

ENDUSAGE
########################################################################################################

#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $count;
my $fastaheader;
my $fileCluster;
my $fileRepeat;
my $fileCoding;
my $fileGenome;
my $fileCoverage;
my @int1;
my @int2;
my %record_fasta;
my %record_pos;
my $repeat;
my %report;
my $SampleName;
my $seq;



#old
my $RepScore;


my $counter;                                # counter
my $cursor_type;



my (@interm, @interm2,@interm3, @interm4);
   
my $scaffold; 
my $position;
my $type;
my $coverage;
my @ALLallelesID;
my $value;

my @test_type;
my $REFallele;
my @ALTalleleList;

my %AlleleDist;


my @test_ref;
my $test_cond;

my @ordered_isolates;
my $str_ordered_isolates;

#......................................................................................................#


if(@ARGV != 5){
  print "version: $version\n";
  print "date: $date_of_version\n";
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion,
            'cluster=s'                    => \$fileCluster,
            'repeat=s'                     => \$fileRepeat,
            'coding=s'                     => \$fileCoding,
						'genome=s'                     => \$fileGenome,
						'coverage=s'                   => \$fileCoverage);

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

#input
open (IN2, $fileCluster) || die "Can't open $fileCluster: $!";
open (IN3, $fileRepeat) || die "Can't open $fileRepeat: $!";
open (IN4, $fileCoding) || die "Can't open $fileCoding: $!";
open (FASTA, $fileGenome) || die "Can't open $fileGenome: $!";
open (IN6, $fileCoverage) || die "Can't open $fileCoverage: $!";

print "file for Cluster: $fileCluster\n";
print "file for Repeat: $fileRepeat\n";
print "file for Coding: $fileCoding\n";
print "file for Genome: $fileGenome\n";
print "file for Coverage: $fileCoverage\n";


#sample name
@int1 = split(/\./,$fileCoverage);
$SampleName = $int1[0];

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);
my $NameOUT1 = "Positions_per_feature_CRC_".$SampleName."_" . $date . "_PRV3.tmpTB3";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";

print "Output file: $NameOUT1\n";

print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                             Parsing Genome file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";


#####  Parsing Genome file  ##########################

$seq = "";

foreach my $line (<FASTA>) {
    chomp($line);
    
    if ($line =~ /^\s*#/) {
	    next;
    } elsif ($line =~ /^\s*$/) {
	    next;	
    } elsif ($line =~ /^>/) {
			if (length($seq) > 0){
				$record_fasta{$fastaheader} = length($seq);
			}
			$seq="";
			$fastaheader = substr($line,1);
			$count ++;
			next;
		} else {
			$seq .= $line;
		}
}
$record_fasta{$fastaheader} = length($seq);

print "Number of sequences: $count\n\n"; 
#print Dumper(\%record); getc(); 


foreach my $header (keys %record_fasta) {
	for (my $cursor=1; $cursor <= $record_fasta{$header}; $cursor ++) {
		$record_pos{$header}{$cursor}{cluster} = 0;
		$record_pos{$header}{$cursor}{repeat} = 0;
		$record_pos{$header}{$cursor}{coding} = 0;
		$record_pos{$header}{$cursor}{coverage} = 0;
	}
}

#print Dumper(\%record); getc();



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                             Parsing 'cluster' file (or 'MultiCopElem' file)\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $line (<IN2>) {
    chomp($line);
    
    if 	($line =~ /^\s*#/) { 			# discard blank line
		next;			
		
    } else {
        @interm = split(/\t/, $line);
        #print">>$interm[0] $interm[1] $interm[2] $interm[3] $interm[4]\n";#getc();
        $RepScore  = [split(/MulticopyMeanScore=/, $interm[8])] -> [1];
        if ($RepScore > 0 and $RepScore < 1) {
            $RepScore = 1;
        }
        for (my $i=$interm[3]; $i <= $interm[4]; $i +=1) {
            #print ">>>>$i\n"; getc();
						$record_pos{$interm[0]}{$i}{cluster} += $RepScore;
						$record_pos{$interm[0]}{$i}{repeat} = 2;
        }
    }
}

#warn Dumper(\%record); getc();

print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                         Parsing 'repeats' file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $line (<IN3>) {
    chomp($line);
    
    if 	($line =~ /^\s*#/) { 			# discard blank line
		next;			
		
    } else {
        @interm = split(/\t/, $line);
        #print">>$interm[0] $interm[1] $interm[2] $interm[3] $interm[4]\n";#getc();
        for (my $i=$interm[3]; $i <= $interm[4]; $i +=1) {
            #print ">>>>$i\n"; getc();
            $record_pos{$interm[0]}{$i}{repeat} += 1;

        }
    }

}


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                         Parsing 'coding' file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $line (<IN4>) {
    chomp($line);
    
    if 	($line =~ /^\s*#/) { 			# discard blank line
		next;			
		
    } else {
        @interm = split(/\t/, $line);
        #print">>$interm[0] $interm[1] $interm[2] $interm[3] $interm[4]\n";#getc();
        
        if ($interm[2] eq "CDS") {
            for (my $i=$interm[3]; $i <= $interm[4]; $i +=1) {
                #print ">>>>$i\n"; getc();
                $record_pos{$interm[0]}{$i}{coding} = 1;

            }
        }
    }

}
#warn Dumper(\%record); getc();


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                         Parsing 'coverage' file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $line (<IN6>) {
    chomp($line);
    
    if 	($line =~ /^\s*#/) { 			# discard blank line
		next;			
		
    } else {
        @interm = split(/\t/, $line);
        $record_pos{$interm[0]}{$interm[1]}{coverage} = $interm[2];
    }

}
#warn Dumper(\%record_pos); getc();


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                              Computation\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $scaffold (  keys %record_pos) {
    foreach my $position (keys %{$record_pos{$scaffold}}) {

			if ($record_pos{$scaffold}{$position}{coverage} >= 10) {
					if ($record_pos{$scaffold}{$position}{repeat} > 0) {
						$repeat = 1;
					} else {
						$repeat = 0;
					}
			
					$report{$repeat}{$record_pos{$scaffold}{$position}{coding}} ++;		
			}
	
    }
}


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                             Write file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";


print OUT1 "Isolate\tRepeats\tCoding\tCovered_positions\n";

foreach my $Repeat (nsort keys %report) {
    foreach my $coding (nsort keys %{$report{$Repeat}}) {
			print OUT1 $SampleName,"\t";
			print OUT1 "$Repeat\t$coding\t", $report{$Repeat}{$coding}, "\n";
	    
    }
}

