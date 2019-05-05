#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Sortkeys  = 1;
use Sort::Naturally;
use Getopt::Long;

my $version = "1.0";
my $date = "14-04-2017";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_30_4_Summarize_onlyKs_results_in_a_single_table.pl

SYNOPSIS

perl PP72_30_4_Summarize_onlyKs_results_in_a_single_table.pl  --pattern=Final__Ks


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-p / --pattern                Pattern to recognize the results files


DESCRIPTION
Read automatically all result files in the folder and create a summary file


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
my $pattern;                                # pattern to recognize the results files



#OLD
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
            'pattern=s'                    => \$pattern);

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

print "Pattrn: $pattern\n";
my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { /$pattern/ } readdir(DIR);
closedir DIR;

$count = 0;
warn "VCF files\n";
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


#####  Output files  ##########################

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

my $NameOUT1 = "Table_result_summary__" . $pattern . "__" . $date_file . ".txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                         Parsing files\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";


foreach my $f (nsort @files) {
  open (IN1, "<$f");
  warn " treatment of $f\n";
   
  foreach my $line (<IN1>) {
    chomp($line);
    
    if 	($line =~ /^\s*#/) {
      next;
        
    } elsif ($line =~ /^\s*$/) {
      next;
        
    } elsif ($line =~ /Filename/) {
      @interm = split(/\t/,$line);
      $temp = $interm[1];
      $temp =~ s/.tsv//;
      $temp =~ s/mainOutput//;
      $record{$temp}{filename} = $interm[1];
             
    } elsif ($line =~ /Number\ of\ genes\ with\ variants/) {
      @interm = split(/\t/,$line);
      $record{$temp}{geneswithvariants} = $interm[1];      
    } elsif ($line =~ /Number\ of\ genes\ containing\ non-synonymous\ variants/) {
      @interm = split(/\t/,$line);
      $record{$temp}{geneswithnon_syn} = $interm[1];
    } elsif ($line =~ /Number\ of\ non-synonymous\ variants/) {
      @interm = split(/\t/,$line);
      $record{$temp}{non_synvariants} = $interm[1];
    } elsif ($line =~ /Number\ of\ synonymous\ variants/) {
      @interm = split(/\t/,$line);
      $record{$temp}{synvariants} = $interm[1];
    } elsif ($line =~ /Number\ of\ non-coding\ variants/) {
      @interm = split(/\t/,$line);
      $record{$temp}{noncodingvariants} = $interm[1];
    }
  }
}

#warn Dumper(\%record); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                      Report final statistics\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $sample (nsort keys %record) {
  print OUT1 "\t", $sample;
}
print OUT1 "\n";

print OUT1 "Number of genes with variants:";
foreach my $sample (nsort keys %record) {
  print OUT1 "\t", $record{$sample}{geneswithvariants};
}
print OUT1 "\n";

print OUT1 "Number of genes with non-synonymous variants:";
foreach my $sample (nsort keys %record) {
  print OUT1 "\t", $record{$sample}{geneswithnon_syn};
}
print OUT1 "\n";


print OUT1 "Number of non-synonymous variants:";
foreach my $sample (nsort keys %record) {
  print OUT1 "\t", $record{$sample}{non_synvariants};
}
print OUT1 "\n";

print OUT1 "Number of synonymous variants";
foreach my $sample (nsort keys %record) {
  print OUT1 "\t", $record{$sample}{synvariants};
}
print OUT1 "\n";

print OUT1 "Number of non-coding variants";
foreach my $sample (nsort keys %record) {
  print OUT1 "\t", $record{$sample}{noncodingvariants};
}
print OUT1 "\n";


__END__



  
  



