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

PP72_04_17v1_Get_the_variant_size_dist_from_SV_vcf.pl

SYNOPSIS

perl PP72_04_17v1_Get_the_variant_size_dist_from_SV_vcf.pl 


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message




DESCRIPTION
Analyze the field SVLEN= in .vcf created with CNVnator or LUMPY.
Classify the length of the deletion


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

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
my $SampleName;
my $Size;                                   # Size of the SV

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
warn"                                               Bin sizes\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

@BinSizes = ("10-100", "101-1000","1001-10000","10001-100000","100001-1000000","1000001-+");

#@BinSizes = ("10-100", "101-500","501-1000","1001-5000","5001-10000","10001-+");

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

$NameOUT1 = "Table_Deletion_size_" . $date_file . ".txt";
print "Output file = ", $NameOUT1, "\n";

open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";

$NameOUT2 = "List_Deletion_size_" . $date_file . ".txt";
print "Output file = ", $NameOUT2, "\n";

open (OUT2, ">$NameOUT2")  || die "Can't open $NameOUT2: $!";



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
              @interm2 = split(/;/, $interm[7]);
              $Size = "";
              foreach my $elem (@interm2) {
                if ($elem =~m/SVLEN=/) {
                  $Size = abs([split(/=/, $elem)] -> [1]);
                }
              }
              #print $Size,"\n";
              if ($Size ne "") {
                foreach my $bin (@BinSizes) {
                  $bin_low = "";
                  $bin_high = "";
                  ($bin_low, $bin_high) = split(/\-/, $bin);
                  if ($bin_high eq "+") {
                    $bin_high = 1000000000000000;
                  }
                  if ($Size >= $bin_low and $Size <= $bin_high) {
                    $record{$SampleName}{$bin} ++;
                  }
                }
                $record_bis{$SampleName}{$Size} ++;
                
              }
              
            }
        }
    }
}

foreach my $SampleName (@List_of_samples) {
  foreach my $bin (@BinSizes) {
    $record{$SampleName}{$bin} += 0;
    $record_bis{$SampleName}{0} ++;
  }
}
 


#warn Dumper(\%record); getc();


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                       Report the results\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

print OUT1 "bin(bp)";
foreach my $SampleName (nsort keys %record) {
  print OUT1 "\t", $SampleName;
}
print OUT1 "\n";

foreach my $bin (@BinSizes) {
  print OUT1 $bin;
  foreach my $SampleName (nsort keys %record) {
    print OUT1 "\t", $record{$SampleName}{$bin};
  }
  print OUT1 "\n";
}
  


foreach my $SampleName (nsort keys %record_bis) {
  #$Rcomp_Name = [split(/\-/, $SampleName)] -> [0] . "-on-" . [split(/\_/, $SampleName)] -> [2];
  #$Rcomp_Name = [split(/\-/, $SampleName)] -> [0];
  $Rcomp_Name = [split(/\_/, $SampleName)] -> [0];
  print OUT2 $Rcomp_Name;
  foreach my $size (nsort keys %{$record_bis{$SampleName}}) {
    if ($size != 0) {
      print OUT2 "\t", $size;
    }
    
  }
  print OUT2 "\n";
}

  

__END__






