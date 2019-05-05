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

PP72_04_12v1_ReadFasta_filterSize_separateScaffolds.pl

SYNOPSIS

perl PP72_04_12v1_ReadFasta_filterSize_separateScaffolds.pl  --file=GenomeFasta.fa  --Size=10000  --out_dir FASTAs


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-f / --file                   Input fasta file
-S / --Size                   Size threshold (bp) for the scaffolds
-o / --out_dir                Output directory


DESCRIPTION
Read a .fasta file and keep the scaffolds bigger than the threshold.
Output a single .fasta file for each retained scaffold (ex: scaffold1.fa)
Fasta files are formatted to 80 characters/line


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
my $file;                                   # input file
my @interm;                                 # array to store split line
my $name;                                   # scaffold name
my $out_dir;                                # Folder to output individual scaffold files
my %record;                                 # main hash to record data
my $seq = "";                               # to store sequence
my $Size;                                   # Size threshold to keep scaffold


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
            'Size=s'                       => \$Size,
            'out_dir=s'                    => \$out_dir);


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
warn"                                             Genome file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

print "Input file = ",$file,"\n";
open (IN1, $file) || die "Can't open $file: $!";

if (-e $out_dir and -d $out_dir) {
    print "Directory $out_dir is found\n";
} else {
    print "Directory not found\n =>Making directory $out_dir\n";
    mkdir $out_dir;
}

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                       Parsing genome .fasta file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $line (<IN1>) {
    chomp($line);
	
	if 	($line =~ /^\s*#/) {
		next;			
	} elsif ($line =~ /^\s*$/) {
		next;	
	} elsif ($line =~ /^>/) {

		if (length($seq) >= $Size){
			open(OUT1, ">$out_dir/$name.fa");
            $count_filt ++;
            $seq =~ s/([^\n]{80}|[^\n]{1,79}$)/$1\n/g;
            print OUT1 ">", $name, "\n";
            print OUT1 $seq, "\n";
            close(OUT1);
            
		}
		$seq="";
		$name=$line;
		$name =~ s/\>//;
		$count_tot ++;
		
		next;

	
	# keep line, add to sequence string
	} else {
		$seq .= $line;
	}
}

if (length($seq) >= $Size){
    open(OUT1, ">$out_dir/$name.fa");
    $count_filt ++;
    $seq =~ s/([^\n]{80}|[^\n]{1,79}$)/$1\n/g;
    print OUT1 ">", $name, "\n";
    print OUT1 $seq, "\n";
    close(OUT1);
}



print "nb of retained sequences: $count_filt\n";
print "Total nb of sequences: $count_tot\n";
print "**************************************","\n";



__END__

warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Parsing List of SNP file and print results   \n\n";

#####....................................................................................................#####
##############################################################################################################



my $line1 = <IN2>;


my $nameOUT1 = "Target_" . $target . ".fasta";
open (OUT1, ">$nameOUT1")  || die "Can't open $nameOUT1: $!";

my $nameOUT2 = "Primer3_Target_" . $target . ".fasta";
open (OUT2, ">$nameOUT2")  || die "Can't open $nameOUT2: $!";

$count = 0;


foreach my $line (<IN2>) {
    chomp($line);
	    
    @interm = split(/\t/, $line);

	
    # discard comment line 
    if 	($line =~ /^\s*#/) {
	next;			
	    
    # discard blank line
    } elsif ($line =~ /^\s*$/) {
	next;
	
    } else {
		
		foreach my $i ( @interm ) {
			warn $i,"\t";
		}
		warn "\n";
		
		$test = 1;
		$name = $interm[0];
		
		#extremity 5'
		if ($interm[2] < 500) {
			$start = 0;
			$PositionSNP = $interm[2]-1;
		} else {
			$start = $interm[2] -1 -500;
			$PositionSNP = 500; 
		}
		
		#extremity 3'
		#print $interm[2]," ", $interm[1], " ", $seqtot{$interm[1]}, "\n"; getc();
		if ($interm[2] > length($seqtot{$interm[1]}) -1 -500) {
			$end = 500 + length($seqtot{$interm[1]}) -1;
		} else {
			$end = 1000 -1;
		}
		
		$TempSeq = substr($seqtot{$interm[1]}, $start, $end);
		
	
		$newBase = substr ($TempSeq, $PositionSNP, 1);
		warn "newBase= $newBase\n";
		$LeftSeq = substr ($TempSeq, 0, $PositionSNP);
		$RightSeq = substr ($TempSeq, $PositionSNP+1);
		warn "TempSeq= $TempSeq\n";
		warn "FinalSe=", $LeftSeq, $newBase, $RightSeq, "\n";
		warn ">> PositionSNP= $PositionSNP\n";
		warn "SNP in TempSeq= ", substr($TempSeq, $PositionSNP,1),"\n";
		warn "around SNP in TempSeq= ", substr($TempSeq, $PositionSNP-5,11),"\n";
	
		$test = 1;
	
		if ($test == 1) {
			$count ++;
			print OUT1 ">seq",$count,"_", $interm[0],"_on_", $interm[1], "\n";
			print OUT1 $LeftSeq, $newBase, $RightSeq, "\n\n";
			print OUT2 ">seq",$count,"_", $interm[0],"_on_", $interm[1], "\n";
			print OUT2 $LeftSeq, "[", substr($newBase,0,1), "]", $RightSeq, "\n\n";
		}
		warn "1 car:", substr($seqtot{$interm[1]}, $interm[2] -1 , 1),"\n";
		warn "10   :", substr($seqtot{$interm[1]}, $interm[2] -1 -5, 11),"\n";
		
		
		#print"\n\n"; getc();
    }

}

__END__





