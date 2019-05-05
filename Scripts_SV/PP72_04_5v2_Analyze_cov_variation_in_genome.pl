#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;
use Getopt::Long;

use DateTime;
use DateTime::Format::MySQL;
use DateTime::Format::Duration;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

my $version = "1.0";
my $date = "04-02-2017";

########################################################################################################
my $usage = <<'ENDUSAGE';

NAME

PP72_04_5v1_Analyze_cov_variation_in_genome.pl

SYNOPSIS

perl PP72_04_5v1_Analyze_cov_variation_in_genome.pl  --window=500   --ref=A1-onA1-genome.coverage   --outputfileprefix=genomeA1


OPTIONS

-V / --version                Print version number 
-h / --help                   Print help message
-w / --window                 Size of the window (bp)
-r / --ref                    Reference coverage file (.coverage.gz, obtained with the reads used to generate the genome)
-g / --genome                 Reference genome file (.fasta)
-o / --outputfileprefix       Prefix for the output file

DESCRIPTION
Calculate coverage ratio between a reference sample and a test sample.
Automatically read all coverage files (.coverage.gz) in the folder


contact:
frederic.masclaux@unil.ch
University of Lausanne

PROJECT = P72

ENDUSAGE
########################################################################################################


#..............Variables...............................................................................#
my $printVersion;                           # print version number, if set
my $help;                                   # print usage

my $count;                                  # counter
my $counter1;                               # counter
my $counter2;                               # counter
my %counterGlob;														# hash to count
my $dt1;														        # date/hour 1
my $dt2;														        # date/hour 2
my $dt_format;														  # format of the date
my $genome;                                 # reference genome file
my $header;                                 # fasta header
my @interm;                                 # array to store split line
my $limit = 500;                            # size to exclude at the extrimities of the scaffolds
my @list_of_SampleNames;                    # 
my $outputfileprefix;                       # Prefix for the output file
my %record;                                 # main hash to record data
my %record_coord;                           # hash to record data about genome coordinates
my %record_fasta;                           # hash to record genome sequences
my $ref;                                    # Reference coverage file
my $sample_name;                            # Contains the name of the sample from the coverage file
my $scaffold;                               # Store the current scaffold
my $seq = "";                               # to store sequence
my $temp_output;                            # temporary output = decompressed cov file
my $threshold_scaff_size = 10000;           # minimal size of the analyzed scaffolds
my $window;                                 # Size of the window to walk on the scaffolds



#OLD

my @interm2;                                # array to store split line
my $isoform_name;                           # to store the isoform name (= i1)
my $length;                                 # length of the sequence (read from the header)
my %Name_Hash;                              # Hash to store the corresponding between short names and long names

my $read_cluster;                           # to store the read cluste name (=TRINITY_DN9981_c0)
my $ref_infogene;                           #
my $ref_isof;                               # reference name

my $size;                                   # to store length of sequence
my $switcher = 1;                           # use as a switcher to start and stop recording
my $tmp;





#......................................................................................................#


if(@ARGV==0){
  print "version: $version\n";
  print "date: $date\n";
  print "$usage\n"; 
  exit(0);
}


GetOptions( 'help!'                        => \$help,
            'version!'                     => \$printVersion,
            'window=s'                     => \$window,
            'ref=s'                        => \$ref,
            'genome=s'                     => \$genome,
            'outputfileprefix=s'           => \$outputfileprefix);


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

####  Search for files  ##########################

my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { /.coverage.gz$/ } readdir(DIR);
closedir DIR;

$count = 0;
$counter1 = 0;
warn ">Compressed coverage files\n";
for my $i (nsort @files) {
    warn "  ",$i,"\n";
    $count ++;
    # search for ref file
    if ($ref eq $i) {
      $counter1 ++;
    }
}
warn "\n";
if ($count == 0) {
    die "ERROR: no .coverage.gz file found   ------\n\n";
} elsif ($counter1 != 1) {
  die "ERROR: no .coverage.gz file is found as reference: $ref   ------\n\n";
} else {
    print ">> $count files\n";    
}


#####  Output files  ##########################

#output
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date_file =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);

my $NameOUT1 = "Ratio_Coverage_" . $outputfileprefix . "_" . $date_file . ".txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";



#####  Format for duration calculation  ########################## 
$dt_format = DateTime::Format::Duration->new(
    pattern => '%e days, %H hours, %M minutes, %S seconds'
);


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                       Analyze the genome .fasta file\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

print "Genome = ",$genome,"\n";
open (FASTA, $genome) || die "Can't open $genome: $!";

foreach my $line (<FASTA>) {
  chomp($line);
  
  if ($line =~ /^\s*#/) {
      next;
  } elsif ($line =~ /^\s*$/) {
      next;	
  } elsif ($line =~ /^>/) {
      if (length($seq) >0){
        $record_fasta{$header} = $seq;
      }
      $seq="";
      $header = substr($line,1);
      $count ++;
      next;
  } else {
      $seq .= $line;
  }
}
$record_fasta{$header} = $seq;

close(FASTA);

foreach my $h (keys %record_fasta) {
  $counter1 = 0;
  $record_coord{$h}{start} = $limit +1;
  $record_coord{$h}{end} = length($record_fasta{$h}) - $limit -1;
  $record_coord{$h}{size_tot} = length($record_fasta{$h});
	
	
	if ($record_coord{$h}{size_tot} >= $threshold_scaff_size) {
		for (my $ii = $record_coord{$h}{start}; $ii <= $record_coord{$h}{end}-$window; $ii = $ii + $window) {
			$counter1 ++;
			$record_coord{$h}{intervalles}{$counter1}{min} = $ii;
			$record_coord{$h}{intervalles}{$counter1}{max} = $ii + $window;
		}
	}

}

#warn Dumper(\%record_coord); getc();


warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                         Parse the .coverage files\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $file (nsort @files) {
    print ">Treatment of $file\n";
    $dt1 = DateTime->now(time_zone => 'floating', formatter => 'DateTime::Format::MySQL');
    print "  Start at $dt1\n";

    #Uncompress the file
    print "  Unzipping $file\n";
    $temp_output = [split(/\.coverage/,$file)] -> [0] . ".coverage";
    gunzip $file => $temp_output or die "gunzip failed: $GunzipError\n";
    
    print "  Temporary file: $temp_output\n";    
    open (IN1, $temp_output) || die "Can't open $temp_output: $!";
    
    #Get sample name
    $sample_name = [split(/\./,$file)] -> [0];
    print "  Sample name: $sample_name\n";
    print "  Recording data...\n";
    
    push @list_of_SampleNames, $sample_name;
    
		$counter1 = 0;
		$counter2 = 0;
		%counterGlob = ();   # ?? ? ? ?
    
    foreach my $line (<IN1>) {
      chomp($line);

      if 	($line =~ /^\s*#/) {
        next;
      } elsif ($line =~ /^\s*$/) {
        next;
      } else {      
        @interm = split(/\t/, $line);
        $scaffold = $interm[0];
        
        if ( (exists $record_coord{$scaffold}) and ($ref eq $file)) {
          foreach my $interv (keys %{$record_coord{$scaffold}{intervalles}}) {
            if ($interm[1] >= $record_coord{$scaffold}{intervalles}{$interv}{min} and $interm[1] <= $record_coord{$scaffold}{intervalles}{$interv}{max}) {
              $record_coord{$scaffold}{intervalles}{$interv}{ref}{sum} += $interm[2];
              $record_coord{$scaffold}{intervalles}{$interv}{ref}{nbpos} += 1;
            }
          }
        }
				if (exists $record_coord{$scaffold}) {
          foreach my $interv (keys %{$record_coord{$scaffold}{intervalles}}) {
            if ($interm[1] >= $record_coord{$scaffold}{intervalles}{$interv}{min} and $interm[1] <= $record_coord{$scaffold}{intervalles}{$interv}{max}) {
              $record_coord{$scaffold}{intervalles}{$interv}{samples}{$sample_name}{sum} += $interm[2];
              $record_coord{$scaffold}{intervalles}{$interv}{samples}{$sample_name}{nbpos} += 1;
            }
          }
				}
      }
    }
    unlink($temp_output);
    close(IN1);
    $dt2 = DateTime->now(time_zone => 'floating', formatter => 'DateTime::Format::MySQL');
    print "  End at $dt2\n";
    print "  Duration:", $dt_format->format_duration($dt2-$dt1), "\n";
}

#warn Dumper(\%record_coord); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                         Complete missing data\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $scaff (keys %record_coord) {
  foreach my $interv (keys %{$record_coord{$scaff}{intervalles}}) {
    foreach my $sample (@list_of_SampleNames) {
      $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{sum} .="";
      $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{nbpos} .="";
    }
    $record_coord{$scaff}{intervalles}{$interv}{ref}{sum} .="";
    $record_coord{$scaff}{intervalles}{$interv}{ref}{nbpos} .="";
  }
}
      

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                        Calculations mean coverage\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $scaff (keys %record_coord) {
  foreach my $interv (keys %{$record_coord{$scaff}{intervalles}}) {
    foreach my $sample (keys %{$record_coord{$scaff}{intervalles}{$interv}{samples}}) {
      if ($record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{sum} ne ""  and  $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{nbpos} ne "" ){
          $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{meanCov} = $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{sum} / $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{nbpos};
      } else {
        $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{meanCov} = "NA";
      }
    }
    #ref
    if ($record_coord{$scaff}{intervalles}{$interv}{ref}{sum} ne ""  and  $record_coord{$scaff}{intervalles}{$interv}{ref}{nbpos} ne "" ){
      $record_coord{$scaff}{intervalles}{$interv}{ref}{meanCov} = $record_coord{$scaff}{intervalles}{$interv}{ref}{sum} / $record_coord{$scaff}{intervalles}{$interv}{ref}{nbpos};
    } else {
      $record_coord{$scaff}{intervalles}{$interv}{ref}{meanCov} = "NA";
    }
  }
}

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                           Calculations ratio\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

foreach my $scaff (keys %record_coord) {
  foreach my $interv (keys %{$record_coord{$scaff}{intervalles}}) {
    foreach my $sample (keys %{$record_coord{$scaff}{intervalles}{$interv}{samples}}) {
      if ($record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{meanCov} ne "NA" and $record_coord{$scaff}{intervalles}{$interv}{ref}{meanCov} ne "NA") {
        $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{ratio} =
              sprintf("%.3f", $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{meanCov} / $record_coord{$scaff}{intervalles}{$interv}{ref}{meanCov});
      } else {
        $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{ratio} = "NA";
      }
    }
  }
}

#warn Dumper(\%record_coord); getc();

warn"##############################################################################################################\n";
warn"#####....................................................................................................#####\n";
warn"                                            Report Results\n";
warn"#####....................................................................................................#####\n";
warn"##############################################################################################################\n";

print OUT1 "Scaffold\tInterv\tPositionRange";

foreach my $isolate (nsort @list_of_SampleNames) {
    print OUT1 "\t", $isolate;
}
print OUT1 "\n";

foreach my $scaff (nsort keys %record_coord) {
  foreach my $interv (nsort keys %{$record_coord{$scaff}{intervalles}}) {
    print OUT1 $scaff,"\t", $interv, "\t", $record_coord{$scaff}{intervalles}{$interv}{min},"-", $record_coord{$scaff}{intervalles}{$interv}{max};
    foreach my $sample (nsort keys %{$record_coord{$scaff}{intervalles}{$interv}{samples}}) {
      print OUT1 "\t", $record_coord{$scaff}{intervalles}{$interv}{samples}{$sample}{ratio};
    }
    print OUT1 "\n";
  }
}

print OUT1 "#" . $ref . "\n";




__END__
