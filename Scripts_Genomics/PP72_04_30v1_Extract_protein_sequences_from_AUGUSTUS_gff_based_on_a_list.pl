#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;


warn "\n$0\n"; 				# nom du fichier
warn "\n","*********************","\n";


my $numArgs = $#ARGV + 1;

unless( $numArgs == 2) {
    
    print "\n";    
    print "===============================================================================\n";
    print "Usage: perl $0 Predicted_genes_by_Augustus.gff  list.txt ";
    print "\n";
    print "===============================================================================\n";
    print "   Author: Frederic Masclaux\n";
    print "   Date: 2018-07-06\n";
    print "   Version: 1\n";
    print "   Name: $0\n";
    print "   Purpose: \n";
    print "   Extract protein sequences from .gff file created by Augustus based on a provided list and prepare a fasta file\n";
    print "   Work with --Augustus-- files\n";      
    print "   Project ID connection: P72 (/P51)\n";
    print "   Comments :\n";
    print "===============================================================================\n";
    exit;
}


warn "\n","*********************","\n";

#####  Opening files  ##########################

my ($fileIN1, $fileIN2) = @ARGV;

open (IN1, $fileIN1) || die "Can't open $fileIN1: $!";
open (IN2, $fileIN2) || die "Can't open $fileIN2: $!";



#####  Output file  ##########################

$fileIN2 =~ s{.*/}{};
my $Filename = [split(/\./, $fileIN2)] -> [0] . "." . [split(/\./, $fileIN2)] -> [1];

my $nameOUT1 = $Filename . ".fasta";
open (OUT1, ">$nameOUT1")  || die "Can't open $nameOUT1: $!";



#.....Variables.....#

my @interm;
my @int1;
my @int2;

my $rec1 = 0;
my $rec2 = 0;

my $Locusname;
my $Transcriptname;
my $ProtSeq;

my %record;

my $genomelocation;



#....................#


warn "\n","*********************","\n";

#####  Parsing list file  ##########################

warn "Parsing list file\n\n";

foreach my $line (<IN2>) {
    chomp($line);

	if 	($line =~ /^\s*#/) {
		next;			

	} elsif ($line =~ /^\s*$/) {
		next;	

	} else {
        
        $record{$line} ++;
    }
}
        
foreach my $t (nsort keys %record) {
    print $t, " ", $record{$t}, "\n";
}


warn "\n","*********************","\n";

#####  Parsing GFF file  ##########################

warn "Parsing GFF file\n\n";

foreach my $line (<IN1>) {
    chomp($line);

    # discard blank line
    if ($line =~ /^\s*$/) {
        next;	
		
    } elsif ($line =~ m/^#\ start\ gene/) {
        if ($rec1 >0) {
            warn "====> $Locusname $Transcriptname $ProtSeq";
            die "rec1 was not reinitialized";
            
        }
        
        @interm = split(/\ /, $line);
        $Locusname = $interm[3];
        $rec1 = 1;
        warn "====> test1 : $Locusname","\n";
	
    } elsif ($rec1 == 1) {
        @interm = split(/\t/, $line);
        warn $interm[2];
        if ($interm[2] eq "transcript") {
            $genomelocation = $interm[0] . "_" . $interm[3] . "_" . $interm[4];
            @int1 = split(/;/, $interm[8]);
            foreach my $t (@int1){
                if (uc($t) =~ m/ID/) {
                    @int2 = split(/\=/, $t);
                    $Transcriptname = $int2[1];
                    $record{$Transcriptname} += 0;
                    $rec1 = 0;
                }
            }
    
        }
	
    } elsif ($line =~ /^#\ protein\ sequence =/) {
        $rec2 = 1;
        @interm = split(/\[/, $line);
        if ($interm[1] =~ m/\]$/) {
            $interm[1] = substr($interm[1], 0, -1);
            $rec2 = 0;
        }
        $ProtSeq = $interm[1];
    } elsif ($rec2 == 1) {
        @interm = split(/\ /, $line);
        
        if ($interm[1] =~ m/\]$/) {
            $interm[1] = substr($interm[1], 0, -1);
            $rec2 = 0;
        }
        $ProtSeq .= $interm[1];
        
    } elsif ($line =~ /^#\ end\ gene/) {
        if ($ProtSeq eq "") {
            die "empty prot"
	} elsif ($rec2 == 1) {
	    die "protein never ended"
	} elsif ($Transcriptname eq "" or $Locusname ="") {
	    die "names are missing", $Locusname, " ", $Transcriptname;

	} elsif ($record{$Transcriptname} > 0 ) {
	    print OUT1 ">", $Transcriptname, " ", $genomelocation, "\n";
	    print OUT1 $ProtSeq, "\n";
	}
	warn ">", $Transcriptname, "\n";
	$rec1 = 0;
	$rec2 = 0;
	$Locusname = "";
	$Transcriptname = "";
	$ProtSeq ="";
    }
    
}

__END__
	
	
	