#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use Getopt::Long;

my $version = "1.0";
my $date = "24-07-2018";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_04_21v1_Analyze_sma3_v2_output.pl

SYNOPSIS

perl PP72_04_21v1_Analyze_sma3_v2_output.pl


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message




DESCRIPTION
Generate a synthetic table from all the '*goslim_go_summary.tsv' files generated with sma3_v2.pl
Read automatically all files in the folders


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $counter = 0;
my $cursor;
my @List_of_obs_first_sec;
my @List_of_sections;
my $section;




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



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                               Files\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";



my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { /goslim_go_summary.tsv$/ } readdir(DIR);
closedir DIR;

$count = 0;
for my $i (@files) {
    print $i,"\n";
    $count ++;
}
print "\n";
if ($count == 0) {
    die "no file found ------\n\n";
} else {
    print ">> $count files\n";    
}
print "\n";


#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

@interm = split(/\./,$files[0]);

$NameOUT1 = "Global_table__" . $interm[1] . "__" . $date_file . ".txt";
print "Output file = ", $NameOUT1, "\n";

open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";



print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                       Parsing .vcf file\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $f (nsort @files) {
  open (IN1, "<$f");
  print " treatment of $f\n";
  

  #sample name and isolate name
  @int1 = split(/_on_/,$f);
  $int1[0] =~s/PE-/WG/g;
  if ($int1[0] eq "DNA1-WG1") {
    $int1[0] = "DAOM-DNA1";
  }
  if ($int1[0] eq "DNA2-WG1") {
    $int1[0] = "DAOM-DNA2";
  }    
  $SampleName = $int1[0];

  print "Sample: $SampleName\n";
  push(@List_of_samples, $SampleName);
  
  $counter = 0;
  foreach my $line (<IN1>) {
    chomp($line);
    
    if ($line =~ /^\s*$/) {
      $cursor = 0;
    } elsif ($line =~ /^\s*#/) {
      $section = $line;
      $section =~ s/\"//g;
      $record{$section}{title}{$SampleName} ++;
      $cursor = 1;
    } else {
      if ($line =~ m/^GO:/) {
        $line =~ s/\t/ /;
      }
      @interm = split(/\t/,$line);

      $record{$section}{data}{$interm[0]}{$SampleName} = $interm[1];

    }
  }
}

#fill empty cells
foreach my $sec (keys %record) {
  foreach my $obs (keys %{$record{$sec}{data}}) {
    foreach my $sample (@List_of_samples) {
      $record{$sec}{data}{$obs}{$sample} += 0;
    }
  }
}


#print Dumper(\%record); getc();


print"##############################################################################################################\n";
print"#####....................................................................................................#####\n";
print"                                       Report the results\n";
print"#####....................................................................................................#####\n";
print"##############################################################################################################\n";

foreach my $sample (nsort @List_of_samples) {
  print OUT1 "\t", $sample;
}
print OUT1 "\n";

@List_of_sections = ("#Annotation summary","##GO Slim", "#Category Biological process", "#Category Cellular component", "#Category Molecular function",
                     "#UniProt Pathways", "##UniProt Keyword categories", "#Category Biological_process", "#Category Cellular_component","#Category Molecular function");

@List_of_obs_first_sec = ("Number of query sequences:", "With Annotations","With GENENAME","With DESCRIPTION","With ENZYME","With GO","With KEYWORD","With PATHWAY",
                          "Annotator a1","Annotator a2","Annotator a3","Annotator a12","Annotator a123");

foreach my $sec (@List_of_sections ) {
  print OUT1 $sec, "\n";
  if ($sec !~ m/^##/ and $sec !~ m/^#Annotation/) {
    foreach my $obs (keys %{$record{$sec}{data}}) {
      print OUT1 $obs;
      foreach my $sample (nsort @List_of_samples) {
        print OUT1 "\t", $record{$sec}{data}{$obs}{$sample};
      }
      print OUT1 "\n";
    }
  } elsif ($sec !~ m/^##/ and $sec =~ m/^#Annotation/) {
    foreach my $obs (@List_of_obs_first_sec) {
      print OUT1 $obs;
      foreach my $sample (nsort @List_of_samples) {
        print OUT1 "\t", $record{$sec}{data}{$obs}{$sample};
      }
      print OUT1 "\n";
    }    
  }
}


__END__






