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

PP72_30_1_Filter_blastx.outfmt6.w_pct_hit_length.pl

SYNOPSIS

perl PP72_30_1_Filter_blastx.outfmt6.w_pct_hit_length.pl  --file1=blastx.outfmt6.w_pct_hit_length   --file2=mainOutput***.tsv


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-f1 / --file1                 Input file from blastx
-f2 / --file2                 Input file from KisSplice
-v / --verbose                To report the deleted elements


DESCRIPTION
Filter the file blastx.outfmt6.w_pct_hit_length to keep only the best isoform of each gene, remove the gene with multiple groups
Filter the file mainOutput***.tsv to remove bubbles that match 2 sequences

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
my $verbose;

my @arrayC1;                                # Array to store Read_counts_variant_1
my @arrayC2;                                # Array to store Read_counts_variant_2
my $count_discard_n1;                       # To count nb of discarded elements
my $count_discard_n2;                       # To count nb of discarded elements
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
my $sumC1;                                  # variable to store sum of Read_counts_variant_1
my $sumC2;                                  # variable to store sum of Read_counts_variant_2
my $temp;                                   # temporary variable
my $temp2;                                  # temporary variable
my @temp;                                   # temporary array
my %temp_hash;                              # to store very temporary hash (counting)

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
            'file2=s'                      => \$file2,
            'verbose!'                     => \$verbose);

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

print "Input file 2 = ",$file2,"\n";
open (IN2, $file2) || die "Can't open $file2: $!";

#####  Output files  ##########################

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

my $NameOUT1 = $file1 . "_filtered-BXKS_" . $date_file . ".txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";

$temp = $file2;
$temp =~ s/.tsv//;
$temp2 = $temp;
$temp2 =~ s/mainOutput//;
my $NameOUT2 = $temp . "_filtered-BxKs_" . $temp2 ."_" . $date_file . ".txt";
open (OUT2, ">$NameOUT2")  || die "Can't open $NameOUT2: $!";

my $NameOUT3 = "Final_results" . "_filtered-BxKs_" . $temp2 ."_" . $date_file . ".txt";
open (OUT3, ">$NameOUT3")  || die "Can't open $NameOUT3: $!";

my $NameOUT4 = "Detailed_final_results" . "_filtered-BxKs_" . $temp2 ."_" . $date_file . ".txt";
open (OUT4, ">$NameOUT4")  || die "Can't open $NameOUT4: $!";

my $NameOUT5 = "List_of_gene_names_" . $temp2 ."_" . $date_file . ".txt";
open (OUT5, ">$NameOUT5")  || die "Can't open $NameOUT5: $!";


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                           Record data from the blastx.outfmt6.w_pct_hit_length file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";
$header1 = "";

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
      
      $record_blastx{$interm[0]} = $line;
      
  }
}

#warn Dumper(\%record_blastx); getc();



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                              Data from blastx.outfmt6.w_pct_hit_length file\n";
warn"                                          Analyze the groups\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

# 2 groups of the same cluster cannot match to the same element in the protein database.
# If they do match, it means that they are the same gene or they belong to the same gene family.
# In this case, they are removed.

%record =();
$counter = 1; # give an unique identifier for fullname

foreach my $rec_gene (keys %record_blastx) {
  @interm = split(/\t/, $record_blastx{$rec_gene});
  
  @interm2 = split(/_/, $interm[0]);
  $isoform = pop(@interm2);
  $group = pop(@interm2);
  
  $gene = join('_', @interm2);
  #print $gene," ", $isoform, "\n";
  
  $record{$gene}{db_hit}{$group}{$interm[1]} ++;
  $record{$gene}{status}{$group} = "keep";
  $record{$gene}{fullname}{$counter} = $interm[0];
  $counter ++;

}


# Discard a cluster if the groups of this cluster (TRINITY_DN1456_c0) have the same match in the database
foreach my $rec_gene (keys %record) {
  %temp_hash = ();
  foreach my $gr (keys %{$record{$rec_gene}{db_hit}}) {
    foreach my $db_hit (keys %{$record{$rec_gene}{db_hit}{$gr}}) {
      $temp_hash{ $db_hit } ++;
    }
  }
  foreach my $t (keys %temp_hash) {
    if ($temp_hash{ $t } > 1) {
      foreach my $gr (keys %{$record{$rec_gene}{status}}) {
        $record{$rec_gene}{status}{$gr} = "discard";
        if ($verbose) {print "discarded $rec_gene \n";}
      }
    }
  }
}

#warn Dumper(\%record); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                     Remove unwanted groups from blastx.outfmt6.w_pct_hit_length file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $rec_gene (keys %record) {
  foreach my $gr (keys %{$record{$rec_gene}{status}}) {
    if ($record{$rec_gene}{status}{$gr} eq "discard") {
      foreach my $t (keys %{$record{$rec_gene}{fullname}}) {
        delete $record_blastx{ $record{$rec_gene}{fullname}{$t} }
      }
    }
  }
}

#warn Dumper(\%record_blastx); getc();



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                              Data from blastx.outfmt6.w_pct_hit_length file\n";
warn"                                          Analyze the isoforms\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

%record =();

foreach my $rec_gene (keys %record_blastx) {
  @interm = split(/\t/, $record_blastx{$rec_gene});
  
  @interm2 = split(/_/, $interm[0]);
  $isoform = pop(@interm2);
  
  $gene = join('_', @interm2);
  #print $gene," ", $isoform, "\n";
  
  $record{$gene}{nb}{$isoform} ++;
  $record{$gene}{db_hit_len}{$isoform} = $interm[12];
  $record{$gene}{pct_hit_len_aligned}{$isoform} = $interm[13];
  $record{$gene}{evalue}{$isoform} = $interm[10];
  $record{$gene}{length}{$isoform} = $interm[3];
  $record{$gene}{db_hit}{$isoform} = $interm[1];
  
  $record{$gene}{status}{$isoform} = "keep";
  $record{$gene}{fullname}{$isoform} = $interm[0];
  
}


# Check if an isoform is found more than one time, which should not happen!!
$Error_control = 0;
foreach my $rec_gene (keys %record) {
  foreach my $iso (keys %{$record{$rec_gene}{nb}}) {
    if ($record{$rec_gene}{nb}{$iso} > 1) {
      print "The isoform ", $record{$rec_gene}{fullname}{$iso}, " has more than one record\n";
      $Error_control = 1;
    }
  }
}
if ($Error_control == 0){
  print ">> PASSED  |  No problem with record redundancy of isoforms\n";
} else {
  print ">> ERROR  |  Record redundancy of isoforms was detected\n";
}



# Discard a gene if the isoforms of this gene have different matches in the database
foreach my $rec_gene (keys %record) {
  %temp_hash = ();
  foreach my $iso (keys %{$record{$rec_gene}{db_hit}}) {
    $temp_hash{ $record{$rec_gene}{db_hit}{$iso} } ++;
  }
   $counter = 0;
  foreach my $t (keys %temp_hash) {
    $counter ++;
  }
  if ($counter > 1) {
    foreach my $iso (keys %{$record{$rec_gene}{status}}) {
      $record{$rec_gene}{status}{$iso} = "discard";
      if ($verbose) {print "discarded $rec_gene \n";}
    }
  }
}




# Define the best isoform
foreach my $rec_gene (keys %record) {
  $highest_value = 0;
  @isoform_list =();
  $isoform ="";
  
  # test for pct_hit_len_aligned
  foreach my $iso (keys %{$record{$rec_gene}{pct_hit_len_aligned}}) {
    if ($record{$rec_gene}{pct_hit_len_aligned}{$iso} > $highest_value) {
      $highest_value = $record{$rec_gene}{pct_hit_len_aligned}{$iso};
      $isoform = $iso;
      
    }
  }
  push @isoform_list, $isoform;
  
  foreach my $iso (keys %{$record{$rec_gene}{pct_hit_len_aligned}}) {
    if ($record{$rec_gene}{pct_hit_len_aligned}{$iso} == $highest_value and $iso ne $isoform) {
      $highest_value = 0;
      push @isoform_list, $iso;
      $isoform ="";
    }
  }
  
  
  if ($isoform ne "") {
    foreach my $iso (keys %{$record{$rec_gene}{status}}) {
      if ($isoform ne $iso) {
        $record{$rec_gene}{status}{$iso} = "discard";
      }
    }
  } else {
    $highest_value = 0;
    $isoform ="";
    
    # test for length
    foreach my $iso (@isoform_list) {
      if ($record{$rec_gene}{length}{$iso} > $highest_value) {
        $highest_value = $record{$rec_gene}{length}{$iso};
        $isoform = $iso;
        
      }
    }
    
    foreach my $iso (@isoform_list) {
      if ($record{$rec_gene}{length}{$iso} == $highest_value and $iso ne $isoform) {
        $highest_value = 0;
        $isoform ="";
      }
    }
    if ($isoform ne "") {
      foreach my $iso (keys %{$record{$rec_gene}{status}}) {
        if ($isoform ne $iso) {
          $record{$rec_gene}{status}{$iso} = "discard";
        }
      }
    } else {
      $highest_value = 0;
      $isoform ="";
      
      # test for evalue
      foreach my $iso (@isoform_list) {
        if ($record{$rec_gene}{evalue}{$iso} > $highest_value) {
          $highest_value = $record{$rec_gene}{evalue}{$iso};
          $isoform = $iso;
          
        }
      }
      
      foreach my $iso (@isoform_list) {
        if ($record{$rec_gene}{evalue}{$iso} == $highest_value and $iso ne $isoform) {
          $highest_value = 0;
          $isoform ="";
        }
      }
      if ($isoform ne "") {    
        foreach my $iso (keys %{$record{$rec_gene}{status}}) {
          if ($isoform ne $iso) {
            $record{$rec_gene}{status}{$iso} = "discard";
          }
        }
      } else {
        @temp = nsort @isoform_list;
        $isoform = shift(@temp);
        foreach my $iso (keys %{$record{$rec_gene}{status}}) {
          if ($isoform ne $iso) {
            $record{$rec_gene}{status}{$iso} = "discard";
          }
        }
      }
    }
  }
   
}

#warn Dumper(\%record); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                     Remove unwanted isoforms from blastx.outfmt6.w_pct_hit_length file\n";
warn"                                              Save the file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$count_discard_n1 = 0;
$count_discard_n2 = 0;

foreach my $rec_gene (keys %record) {
  $count_discard_n1 ++;
  foreach my $iso (keys %{$record{$rec_gene}{status}}) {
    if ($record{$rec_gene}{status}{$iso} eq "discard") {
      delete $record_blastx{ $record{$rec_gene}{fullname}{$iso} };
      $count_discard_n2 ++;
    }
  }
}

print "Number of genes discarded because their isoforms: ", $count_discard_n1, "\n";
print "Number of isoforms discarded: ", $count_discard_n2, "\n";
#warn Dumper(\%record_blastx); getc();
print OUT1 $header1,"\n";
foreach my $rec_gene (nsort keys %record_blastx) {
  print OUT1 $record_blastx{$rec_gene}, "\n";
}
close(OUT1);


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                              Record data from the mainOutput***.tsv file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$header2 ="";
%temp_hash = ();

foreach my $line (<IN2>) {
  chomp($line);
  
  if ($header2 eq "") {
    $header2 = $line             . "\tCountsAllele1\tCountsAllele2\tfreqAllele1";
  } elsif ($line =~ /^\s*#/) {
      next;
  } elsif ($line =~ /^\s*$/) {
      next;	
  } else {
      @interm = split(/\t/, $line);
      $temp_hash{$interm[0]} ++;
      $record_kissplice{$interm[0]}{ $temp_hash{$interm[0]} } = $line;
      
  }
}



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                         Filter based on the read coverage \n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$count_discard_n1 = 0;
$count_discard_n2 = 0;
%record =();

# get total number of records per gene
foreach my $rec_gene (keys %record_kissplice) {
  $record{discarded}{$rec_gene} = 0;
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    $record{total}{$rec_gene} ++;
  }
}


foreach my $rec_gene (keys %record_kissplice) {
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    $sumC1 = 0;
    $sumC2 = 0;  
    @interm = split(/\t/, $record_kissplice{$rec_gene}{$t});
    
    @arrayC1 = split(/\|/, $interm[14]);
    for (@arrayC1) {
      s/C\d+_//g;
    }
    for my $t (@arrayC1) {
      $sumC1 += $t;
    }   
    
    @arrayC2 = split(/\|/, $interm[15]);
    for (@arrayC2) {
      s/C\d+_//g;
    }
    for my $t (@arrayC2) {
      $sumC2 += $t;
    }
  
    if ($sumC1+$sumC2 < 6 or $sumC1 < 2 or $sumC2 < 2 ) {
      if ($verbose) {print "deleted:  " , $interm[0], " ", $interm[1]," ",  $interm[14] ," ",$interm[15],"\n";}
      delete $record_kissplice{$rec_gene}{$t};
      $record{discarded}{$rec_gene} ++;
      $count_discard_n1 ++;
    } else {
      $record_kissplice{$rec_gene}{$t} = $record_kissplice{$rec_gene}{$t} . "\t" . $sumC1 . "\t" . $sumC2  . "\t" . sprintf("%.2g", $sumC1/($sumC1+$sumC2));
    }
  }

}

print "Number of bubbles discarded because low coverage: ", $count_discard_n1, "\n";

# Remove genes for which there is no record left
$count_discard_n1 = 0;
foreach my $rec_gene (keys %{$record{total}}) {
  if ($record{total}{$rec_gene} == $record{discarded}{$rec_gene}) {
    delete $record_kissplice{$rec_gene};
    if ($verbose) {print "deleted:  " , $rec_gene, "\n";}
    $count_discard_n1 ++;
  }
}
print "Number of genes discarded because filtering coverage removed all records: ", $count_discard_n1, "\n";
#warn Dumper(\%record_kissplice);# getc();



warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                   Filter data from mainOutput***.tsv file\n";
warn"                                              Save the file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

%record =();

# get total number of records per gene
foreach my $rec_gene (keys %record_kissplice) {
  $record{discarded}{$rec_gene} = 0;
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    $record{total}{$rec_gene} ++;
  }
}

# filter and delete unwanted records
$count_discard_n1 = 0;
foreach my $rec_gene (keys %record_kissplice) {
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    @interm = split(/\t/, $record_kissplice{$rec_gene}{$t});
    if ($interm[10] eq "True" or $interm[12] eq "True") {             # 10 = several components = several "genes"    /   12 = sequencing error
      if ($verbose) {print "deleted:  " , $interm[1], " ", $interm[0]," ",  $interm[10] ," ",$interm[12],"\n";}
      delete $record_kissplice{$rec_gene}{$t};
      $record{discarded}{$rec_gene} ++;
      $count_discard_n1 ++ ;
    }
  }
}
print "Number of records discarded because the 'component' and the 'sequencing' filters : ", $count_discard_n1, "\n";

# Remove genes for which there is no record left
$count_discard_n1 = 0;
foreach my $rec_gene (keys %{$record{total}}) {
  if ($record{total}{$rec_gene} == $record{discarded}{$rec_gene}) {
    delete $record_kissplice{$rec_gene};
    if ($verbose) {print "deleted:  " , $rec_gene, "\n";}
    $count_discard_n1 ++;
  }
}
print "Number of genes discarded because filtering removed all records: ", $count_discard_n1, "\n";


#warn Dumper(\%record_kissplice); getc();
print OUT2 $header2,"\n";
foreach my $rec_gene (nsort keys %record_kissplice) {
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    print OUT2 $record_kissplice{$rec_gene}{$t}, "\n";
  }
}
close(OUT2);

#warn Dumper(\%record_kissplice); #getc();


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                       Find common elements between blastx data and KisSplice data\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

$counter = 0;
my $cdiscard = 0;
foreach my $rec_gene (keys %record_kissplice) {
  if (exists $record_blastx{$rec_gene}) {
    $counter ++;
  } else {
    delete $record_kissplice{$rec_gene};
    $cdiscard ++;
  }
}
print "Removed elements (and not genes)" ,"\n";
print $counter,"\n";
print "total:", $cdiscard+$counter,"\n";

#warn Dumper(\%record_kissplice); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                      Report final statistics\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";


print "Filename:\t", $file2, "\n";
print OUT3 "=====" x 30, "\n";
print OUT3 "Filename:\t", $file2, "\n";

# Number of genes with a blastx match after filtration
%record_counter = ();
foreach my $rec_gene (keys %record_blastx) {
  $record_counter{total}{$rec_gene} ++; 
}

$counter = 0;
foreach my $t  (keys %{$record_counter{total}}) {
  $counter ++;
}
print "Number of genes with a blastx match after filtration:\t", $counter, "\n";
print OUT3 "Number of genes with a blastx match after filtration:\t", $counter, "\n";
#----


# Number of genes with variants
%record_counter = ();
foreach my $rec_gene (keys %record_kissplice) {
  $record_counter{total}{$rec_gene} ++;

  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {  
    @interm = split(/\t/, $record_kissplice{$rec_gene}{$t});
    #record genes with non-syn variants
    if ($interm[2] eq "True" and $interm[9] eq "True") {
      $record_counter{nonsyn}{$rec_gene} ++;
    }
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
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    @interm = split(/\t/, $record_kissplice{$rec_gene}{$t});
  
    $record_counter{ $interm[2] }{ $interm[9] } ++;
  }
}

print "Number of non-synonymous variants:\t",  $record_counter{ "True" }{ "True" }, "\n";
print "Number of synonymous variants:\t",  $record_counter{ "True" }{ "False" }, "\n";
print "Number of non-coding variants:\t",  $record_counter{ "False" }{ "N/A" }, "\n";

print OUT3 "Number of non-synonymous variants:\t",  $record_counter{ "True" }{ "True" }, "\n";
print OUT3 "Number of synonymous variants:\t",  $record_counter{ "True" }{ "False" }, "\n";
print OUT3 "Number of non-coding variants:\t",  $record_counter{ "False" }{ "N/A" }, "\n";
print OUT3 "=====" x 30, "\n";
#warn Dumper(\%record_kissplice); getc();
  

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                      Report detailed final table\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

print OUT4 $header2, "\t", $header1, "\n";
foreach my $rec_gene (nsort keys %record_kissplice) {
  foreach my $t (keys %{$record_kissplice{$rec_gene}}) {
    print OUT4 $record_kissplice{$rec_gene}{$t}, "\t";
    print OUT4 $record_blastx{$rec_gene}, "\n";
  }
}
close(OUT4);


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                      Report list of genes\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

%temp_hash = ();
%record = ();

foreach my $rec_gene (nsort keys %record_kissplice) {
  @interm = split(/\t/, $record_blastx{$rec_gene});
  $temp_hash{$interm[0]} ++;
  if ($temp_hash{$interm[0]} == 1) {
    $record{$interm[14]} ++;
  }
}
foreach my $gene_func (nsort keys %record) {
  print OUT5 $gene_func, "\t", $record{$gene_func}, "\n";
}
close(OUT5);
