#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Sortkeys  = 1;
use Sort::Naturally;
use Getopt::Long;

my $version = "1.0";
my $date = "10-03-2017";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_30_2_Filter_and_assemble_data_from_kissplice.pl

SYNOPSIS

perl PP72_30_2_Filter_and_assemble_data_from_kissplice.pl  --file1=mainOutput***.tsv


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-f1 / --file1                 Input file from KisSplice


DESCRIPTION
Filter the file mainOutput***.tsv to remove bubbles that match several sequences or bubbles that match several comonents
Report the final result


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72_30

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $counter;                                # counter
my $Error_control;                          # cursor to record if an error was found
my $file1;                                  # input file blastx.outfmt6.w_pct_hit_length 
my $file2;                                  # input file mainOutput***.tsv
my $gene;                                   # gene name
my $group;                                  # group name -> g1 g2...
my $header1;
my $header2;
my $highest_value;                          # to make the sorting
my @interm;                                 # array to store split line
my @interm2;                                # secondary array to store split line
my $isoform;                                # isoform name -> i1 i2...
my @isoform_list;                           # tab to contain the isoforms which requires additional filtering
my %record;                                 # hash to record data to be treated
my %record_counter;                         # hash to count
my %record_blastx;                          # hash to record all blastx.outfmt6.w_pct_hit_length data
my %record_kissplice;                       # hash to record all kisplice data
my $temp;                                   # temporary variable
my $temp2;                                  # temporary variable
my @temp;                                   # temporary array
my %temp_hash;                              # to store very temporary hash (counting)


#old
my %counterGlob;														# hash to count
my $dt1;														        # date/hour 1
my $dt2;														        # date/hour 2
my $dt_format;														  # format of the date
my $genome;                                 # reference genome file
my $header;                                 # fasta header

my $limit = 500;                            # size to exclude at the extrimities of the scaffolds
my @list_of_SampleNames;                    # 
my $outputfileprefix;                       # Prefix for the output file

my %record_coord;                           # hash to record data about genome coordinates
my %record_fasta;                           # hash to record genome sequences
my $ref;                                    # Reference coverage file
my $sample_name;                            # Contains the name of the sample from the coverage file
my $scaffold;                               # Store the current scaffold
my $seq = "";                               # to store sequence
my $temp_output;                            # temporary output = decompressed cov file
my $threshold_scaff_size = 10000;           # minimal size of the analyzed scaffolds
my $window;                                 # Size of the window to walk on the scaffolds



#......................................................................................................#


if(@ARGV==0){
  print "version: $version\n";
  print "date: $date\n";
  print "$usage\n"; 
  exit(0);
}


GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion,
            'file1=s'                      => \$file1,
            'file2=s'                      => \$file2);

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

print "Input file 1 = ",$file1,"\n";
open (IN1, $file1) || die "Can't open $file1: $!";


#####  Output files  ##########################

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

$temp = $file1;
$temp =~ s/.tsv//;
$temp2 = $temp;
$temp2 =~ s/mainOutput//;
my $NameOUT1 = $temp . "_filtered-Ks_" . $temp2 ."_" . $date_file . ".txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";

my $NameOUT3 = "Final_results" . "_filtered-Ks_" . $temp2 ."_" . $date_file . ".txt";
open (OUT3, ">$NameOUT3")  || die "Can't open $NameOUT3: $!";




warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                           Record data from the mainOutput***.tsv file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";
$header1 ="";

foreach my $line (<IN1>) {
  chomp($line);
  
  if ($header1 eq "") {
    $header1 = $line;
  } elsif ($line =~ /^\s*#/) {
      next;
  } elsif ($line =~ /^\s*$/) {
      next;	
  } else {
      @interm = split(/\t/, $line);
      
      $record_kissplice{$interm[0]} = $line;
      
  }
}

#warn Dumper(\%record_kissplice); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                   Filter data from mainOutput***.tsv file\n";
warn"                                              Save the file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

# 2 groups of the same cluster cannot match to the same element in the protein database.
# If they do match, it means that they are the same gene or they belong tothe same gene family.
# In this case, they are removed.

%record =();
$counter = 1; # give an unique identifier for fullname

foreach my $rec_gene (keys %record_kissplice) {
  @interm = split(/\t/, $record_kissplice{$rec_gene});
  if ($interm[11] eq "True" or $interm[12] eq "True") {             # 10 = several components = several "genes"    /   12 = sequencing error
    print "deleted:  " , $interm[1], " ", $interm[0]," ",  $interm[10] ," ",$interm[12],"\n";
    delete $record_kissplice{$rec_gene};
  }

}
#warn Dumper(\%record_kissplice); getc();
print OUT1 $header1,"\n";
foreach my $rec_gene (keys %record_kissplice) {
  print OUT1 $record_kissplice{$rec_gene}, "\n";
}
close(OUT1);


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                      Report final statistics\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";


print "Filename:\t", $file1, "\n";
print OUT3 "=====" x 30, "\n";
print OUT3 "Filename:\t", $file1, "\n";



# Number of genes with variants
%record_counter = ();
foreach my $rec_gene (keys %record_kissplice) {
  $record_counter{total}{$rec_gene} ++;
  
  @interm = split(/\t/, $record_kissplice{$rec_gene});
  #record genes with non-syn variants
  if ($interm[9] =~ m/\|/) {                  # this loop searching for "|" makes the script compatible with the type 0b of KisSplice
    if ($interm[9] =~ m/True/) {              # it simplifies the fields "Is_not_synonymous" to only one characteristics, either "TRUE" or "FALSE". If both T and F are found, T is kept
      $interm[9] = "True";                    # 
    } elsif ($interm[9] =~ m/False/) {
      $interm[9] = "False";
    } else {
      $interm[9] = "N/A";
    }
  }
  if ($interm[2] =~ m/\|/) {                  # this loop searching for "|" makes the script compatible with the type 0b of KisSplice
    if ($interm[2] =~ m/True/) {              # it simplifies the fields "Is_in_CDS" to only one characteristics, either "TRUE" or "FALSE".
      $interm[2] = "True";                    # 
    } else {
      $interm[2] = "False";
    }
  }  
  if ($interm[2] eq "True" and $interm[9] eq "True") {
    $record_counter{nonsyn}{$rec_gene} ++
  }
  
}

$counter = 0;
foreach my $t  (keys %{$record_counter{total}}) {
  $counter ++;
}
print "Number of genes with variants:\t", $counter, "\n";
print OUT3 "Number of genes with variants:\t", $counter, "\n";

$counter = 0;
foreach my $t  (keys %{$record_counter{nonsyn}}) {
  $counter ++;
}
print "Number of genes containing non-synonymous variants:\t", $counter, "\n";
print "\n";
print OUT3 "Number of genes containing non-synonymous variants:\t", $counter, "\n";
print OUT3 "\n";

# Proportions of each type of variant: non-coding, syn, non-syn
%record_counter = ();
$record_counter{ "True" }{ "True" } = 0;
$record_counter{ "True" }{ "False" } = 0;
$record_counter{ "False" }{ "N/A" } = 0;

foreach my $rec_gene (keys %record_kissplice) {
  @interm = split(/\t/, $record_kissplice{$rec_gene});
  if ($interm[9] =~ m/\|/) {                  # this loop searching for "|" makes the script compatible with the type 0b of KisSplice
    if ($interm[9] =~ m/True/) {              # it simplifies the fields "Is_not_synonymous" to only one characteristics, either "TRUE" or "FALSE". If both T and F are found, T is kept
      $interm[9] = "True";                    # 
    } elsif ($interm[9] =~ m/False/) {
      $interm[9] = "False";
    } else {
      $interm[9] = "N/A";
    }
  }
  if ($interm[2] =~ m/\|/) {                  # this loop searching for "|" makes the script compatible with the type 0b of KisSplice
    if ($interm[2] =~ m/True/) {              # it simplifies the fields "Is_in_CDS" to only one characteristics, either "TRUE" or "FALSE".
      $interm[2] = "True";                    # 
    } else {
      $interm[2] = "False";
    }
  }  
  
  $record_counter{ $interm[2] }{ $interm[9] } ++;
}


print "Number of non-synonymous variants:\t",  $record_counter{ "True" }{ "True" }, "\n";
print "Number of synonymous variants:\t",  $record_counter{ "True" }{ "False" }, "\n";
print "Number of non-coding variants:\t",  $record_counter{ "False" }{ "N/A" }, "\n";

print OUT3 "Number of non-synonymous variants:\t",  $record_counter{ "True" }{ "True" }, "\n";
print OUT3 "Number of synonymous variants:\t",  $record_counter{ "True" }{ "False" }, "\n";
print OUT3 "Number of non-coding variants:\t",  $record_counter{ "False" }{ "N/A" }, "\n";
print OUT3 "=====" x 30, "\n";
#warn Dumper(\%record_kissplice); getc();
  
  



