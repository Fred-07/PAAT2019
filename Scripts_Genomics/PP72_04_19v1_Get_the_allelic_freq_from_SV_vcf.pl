#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use Getopt::Long;

my $version = "1.0";
my $date = "13-07-2017";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_04_19v1_Get_the_allelic_freq_from_SV_vcf.pl

SYNOPSIS

perl PP72_04_19v1_Get_the_allelic_freq_from_SV_vcf.pl


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message




DESCRIPTION
Analyze the field SVLEN= in .vcf created with CNVnator or LUMPY.
Determine the allelic frequency of the deletion


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $AOfieldpos;                             # position of the field AO in the string
my $AOfreq;                                 # number of observed AO
my $bin_low;
my $bin_high;
my @BinSizes;                               # Bin sizes for counting
my $count;                                  # counter
my @int1;                                   # temporary array
my @interm;                                 # array to store split line
my @interm2;                                # array to store split element (2nd level)
my @List_of_samples;                        # Contains all saple names
my $NameOUT1;                               # Name of the output file
my $NameOUT2;                               # Name of the output file
my $Rcomp_Name;                             # Sample name compatible with R 
my %record;                                 # main hash to record data
my %record_bis;
my $ROfieldpos;                             # position of the field RO in the string
my $ROfreq;                                 # number of observed RO
my $SampleName;
my $Size;                                   # Size of the SV
my $SVfreq;

#......................................................................................................#



if(@ARGV!=0){
  print "version: $version\n";
  print "date: $date\n";
  print "$usage\n"; 
  exit(0);
}


GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion);


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



my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { /.vcf$/ } readdir(DIR);
closedir DIR;

$count = 0;
warn "VCF files\n";
for my $i (@files) {
    warn $i,"\n";
    $count ++;
}
warn "\n";
if ($count == 0) {
    die "no .subVCF file found ------\n\n";
} else {
    print ">> $count files\n";    
}



#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

$NameOUT1 = "List_allelic_freq_for_Deletion_" . $date_file . ".txt";
print "Output file = ", $NameOUT1, "\n";

open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                       Parsing .vcf file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $f (nsort @files) {
    open (IN1, "<$f");
    warn " treatment of $f\n";
    

    #sample name and isolate name
    @int1 = split(/\_/,$f);
    $SampleName = $int1[0];
    #$SampleName = [split(/\./,$f)] -> [0];
    warn "Sample: $SampleName\n";
    #$sampleNB{$IsolateName} ++;
    push(@List_of_samples, $SampleName);
    
    foreach my $line (<IN1>) {
        chomp($line);
        if 	($line =~ /^\s*#/) {
            next;
        } elsif ($line =~ /^\s*$/) {
            next;
        } else {
            @interm = split(/\t/,$line);
            
            if ($interm[7] =~ m/SVTYPE=DEL/) {
              @interm2 = split(/\:/, $interm[8]);
              $ROfieldpos = 0;
              $AOfieldpos = 0;
              $count = 0;
              foreach my $elem (@interm2) {
                if ($elem eq "RO") {
                  $ROfieldpos = $count;
                } elsif ($elem eq "AO") {
                  $AOfieldpos = $count;
                }
                $count ++;
              }
              @interm2 = split(/\:/, $interm[9]);
              $ROfreq = $interm2[$ROfieldpos];
              $AOfreq = $interm2[$AOfieldpos];
              $SVfreq =  sprintf("%.2f", $AOfreq/($AOfreq+$ROfreq));
              $record{$SampleName}{$SVfreq} += 0;
            }
        }
    }
}

#warn Dumper(\%record); getc();


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                       Report the results\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $SampleName (nsort keys %record) {
  #$Rcomp_Name = [split(/\-/, $SampleName)] -> [0] . "-on-" . [split(/\_/, $SampleName)] -> [2];
  $Rcomp_Name = [split(/\-/, $SampleName)] -> [0];
  print OUT1 $Rcomp_Name;
  foreach my $freq (nsort keys %{$record{$SampleName}}) {
      print OUT1 "\t", $freq;
  }
  print OUT1 "\n";
}

  

__END__






