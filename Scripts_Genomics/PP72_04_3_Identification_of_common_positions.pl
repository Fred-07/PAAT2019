#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Sort::Naturally;

warn "\n$0\n"; 				# nom du fichier
warn "\n","*********************","\n";


my $numArgs = $#ARGV + 1;

unless( $numArgs == 1) {
    
    print "\n";    
    print "===============================================================================\n";
    print "Usage: perl $0  freedomdegree\n";
    print "\n";
    print "===============================================================================\n";
    print "   Author: Frederic Masclaux\n";
    print "   Date: 2016-09-26\n";
    print "   Version: 1\n";
    print "   Name: $0\n";
    print "   Purpose: \n";
    print "   Analyze .subVCF files to determine which positions have the same allelic composition (among the variable positions)\n";
    print "   Report nb of common positions between replicates\n";
    print "   Generate .CommonPos.subVCF files\n";
    print "   freedomdegree can be 0, -1 or -2. It corresponds to the number of samples that can be removed from the commons\n";
    print "   The script excludes all subVCF files containing the string 'CommonPos' in their name\n";    
    print "   Project ID connection: P72\n";
    print "   Comments :\n";
    print "===============================================================================\n";
    exit;
}


my ($freedomdegree) = @ARGV;


#.....Variables.....#
my $count;

my @int1;
my @int2;

my $ReplicateName;
my $IsolateName;

my %sampleNB;

my $header = "";

my $scaffold; 
my $position;
my $RAD;
my $REFallele;
my @TypeAllele;
my @ListAllele;

my @interm;
my %record;
my %file_content;

my $temp;

my $Allele_is_common;
my $Allele_tested;
my %report_results;

my %test_hash;
my $nb_rep;

my $allele;
my $Allele_count;

my $Result_PA_Valid;
my $NbMissingAllele_perSite;
my $NbMissingAllele_in1rep;
my $NbMissingAllele_in2reps;
my $NbMissingAllele_in3reps;

my $freedomdegreeName;



warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Input and Output files\n\n";

#####....................................................................................................#####
##############################################################################################################


my $some_dir = ".";
opendir(DIR, $some_dir) || die "can't opendir $some_dir: $!";
my @files = grep { !/CommonPos/ } grep { /.subVCF$/ } readdir(DIR);
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


#output file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date =  sprintf ("%04d%02d%02d",$year+1900,$mon+1,$mday);
my $NameOUT1 = "Table_stats_from_subVCF__Find_common_positions_among_replicates_" . $date . ".txt";
open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";


# Freedom degree
if ($freedomdegree eq "0" ) {
    $freedomdegreeName = "all-";
} elsif ($count - $freedomdegree > 0 and $freedomdegree < 0 and $freedomdegree > -3) {
    $freedomdegreeName = "partial-";
} else {
    die("problem with freedomdegree...")
}
print ">> " . $freedomdegreeName . ($count +  $freedomdegree) . "(". $freedomdegree . ")\n";


warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Parsing .subVCF files - 1st step - To know where to look\n\n";

#####....................................................................................................#####
##############################################################################################################


foreach my $f (nsort @files) {
    open (IN1, "<$f");
    warn " treatment of $f\n";
    

    #sample name and isolate name
    @int1 = split(/\./,$f);
    $ReplicateName = $int1[0];  
    @int1 = split(/\-/,$f);
    $IsolateName = $int1[0];
    warn "Isolate: $IsolateName  -  Sample: $ReplicateName\n";
    $sampleNB{$IsolateName} ++;
    
    foreach my $line (<IN1>) {
        chomp($line);
        
        if 	($line =~ /^\s*#/) {
            next;
            
        } elsif ($line =~ /^\s*$/) {
            next;
            
        } else {
	    
            @interm = split(/\t/,$line);
           
            #Scaffold and position
            $scaffold = $interm[1];
            $position = $interm[3];
           
            if ($interm[6] == 0 and $interm[9] ne "") {
                $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{toCHECK} += 1;
            }
        }
    }
}



warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Parsing .subVCF files - 2nd step - To record the data\n\n";

#####....................................................................................................#####
##############################################################################################################


foreach my $f (nsort @files) {
    open (IN1, "<$f");
    warn " treatment of $f\n";
    

    #sample name
    @int1 = split(/\./,$f);
    $ReplicateName = $int1[0];  
    @int1 = split(/\-/,$f);
    $IsolateName = $int1[0];
    warn "Isolate: $IsolateName  -  Sample: $ReplicateName\n";

    
    foreach my $line (<IN1>) {
        chomp($line);
        
        if 	($line =~ /^\s*#/) {
            if ($header eq "") {
                $header = $line . "\n";
            }
            next;
            
        } elsif ($line =~ /^\s*$/) {
            next;
            
        } else {
	    
            @interm = split(/\t/,$line);
            
            #Scaffold and position
            $scaffold = $interm[1];
            $RAD = $interm[2];
            $position = $interm[3];
            $REFallele = $interm[8];
            
            if (exists $record{$scaffold} and exists $record{$scaffold}{$position} and exists $record{$scaffold}{$position}{SAMPLES}{$IsolateName} and $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{toCHECK} > 0) {
                
                if ($interm[9] ne "") {

                    $record{$scaffold}{$position}{REF} = $REFallele;
        
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{TYPES} = $interm[9];
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{ALLELES} = $interm[10];
                    
                    $temp = 0;
                    if ($interm[10] =~ m/,/) {
                        $temp = 1; 
                    }
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{foundPA} += $temp;
                                        
                    $temp = 0;
                    if ($interm[9] !~ m/S/) {
                        $temp = 1;
                    }
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{foundnotS} += $temp;
                    
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{foundNA} += 0;
                    
                } elsif ($interm[9] eq "" and $interm[10] ne "na") {
                    $record{$scaffold}{$position}{REF} = $REFallele;
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{TYPES} = "R";
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{ALLELES} = $interm[10];
                    #$record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{TYPES} = "na";
                    #$record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{ALLELES} = "na";
                    #$record{$scaffold}{$position}{SAMPLES}{$IsolateName}{foundNA} += 1;
                    
                } elsif ($interm[9] eq "" and $interm[10] eq "na") {
                    $record{$scaffold}{$position}{REF} = $REFallele;
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{TYPES} = "na";
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{REP}{$ReplicateName}{ALLELES} = "na";
                    $record{$scaffold}{$position}{SAMPLES}{$IsolateName}{foundNA} += 1;
                }
                
                # Record file content
                $file_content{ToReport}{$scaffold}{$position}{$IsolateName} += 0;
                $file_content{DATA}{$IsolateName}{$ReplicateName}{$scaffold}{$position} = $line . "\n";
                #warn Dumper(\%file_content); getc();
                
            }
        }
    }
}

#warn Dumper(\%record); getc();

warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Analyze data...\n\n";

#####....................................................................................................#####
##############################################################################################################

foreach my $scaffold (  (keys %record)) {
    foreach my $position ( (keys %{$record{$scaffold}})) {
        foreach my $iso ( (keys %{$record{$scaffold}{$position}{SAMPLES}})) {
            
            # SNP consistency
            $Allele_is_common = "true";
            $Allele_tested = "";
            if (exists $record{$scaffold} and exists $record{$scaffold}{$position} and exists $record{$scaffold}{$position}{SAMPLES}{$iso}) {
                if ($record{$scaffold}{$position}{SAMPLES}{$iso}{foundPA} == 0) {
                    
                    if ($record{$scaffold}{$position}{SAMPLES}{$iso}{foundnotS} == 0) {
                        %test_hash = ();
                        $nb_rep = 0;
                        $Allele_count = 0;
                        foreach my $rep (nsort (keys %{$record{$scaffold}{$position}{SAMPLES}{$iso}{REP}})) {
                            $nb_rep ++;
                            $allele = $record{$scaffold}{$position}{SAMPLES}{$iso}{REP}{$rep}{ALLELES};
                            $test_hash{$allele} ++;
                        }
                        $test_hash{na} += 0;
                        #warn Dumper(\%test_hash); getc();
                        
                        foreach my $allele (keys %test_hash) {
                            if ($allele ne "na") {
                                $Allele_count ++;
                            }
                        }
                        
                        if ($Allele_count == 1 and $test_hash{na} == 0) {
                            $report_results{"01commonSNP_noNA"}{$iso} ++;
                            $file_content{ToReport}{$scaffold}{$position}{$iso} = 1;
                        } elsif ($Allele_count == 1 and $test_hash{na} == 1) {
                            $report_results{"02commonSNP_1NA"}{$iso} ++;
                        } elsif ($Allele_count == 1 and $test_hash{na} == 2) {
                            $report_results{"03commonSNP_2NAs"}{$iso} ++;
                        } elsif ($Allele_count == 1 and $test_hash{na} > 2) {
                            $report_results{"04commonSNP_fewNAs"}{$iso} ++;                            
                        } elsif ($Allele_count > 1 and $test_hash{na} == 0) {
                            $report_results{"05inconsistentSNPs_noNA"}{$iso} ++;
                        } elsif ($Allele_count > 1 and $test_hash{na} > 0) {
                            $report_results{"06inconsistentSNPs_fewNAs"}{$iso} ++;
                            print "Inconsistent at: ", $scaffold, " ",  $position, "\n";
                        } else {
                            $report_results{"07undeterminedSNPsite"}{$iso} ++;
                        }
                        $report_results{"08TotalnbSNPsites"}{$iso} ++;

                    }
                }
            }


            

            # Poly-allelic site consistency
            $Allele_is_common = "true";
            $Allele_tested = "";
            if (exists $record{$scaffold} and exists $record{$scaffold}{$position} and exists $record{$scaffold}{$position}{SAMPLES}{$iso}) {
                if ($record{$scaffold}{$position}{SAMPLES}{$iso}{foundPA} > 0) {
                    %test_hash = ();
                    $nb_rep = 0;
                    foreach my $rep (nsort (keys %{$record{$scaffold}{$position}{SAMPLES}{$iso}{REP}})) {
                        $nb_rep ++;
                        foreach my $allele (split(",", $record{$scaffold}{$position}{SAMPLES}{$iso}{REP}{$rep}{ALLELES})) {
                            $test_hash{$allele}{nb} ++;
                            $test_hash{$allele}{$nb_rep} ++;
                        }
                    }
                    $test_hash{na}{nb} += 0;
                     
                    #warn Dumper(\%test_hash); getc();
                    
                    #Check for same alleles found in all rep
                    $Allele_tested = "tested";
                    $Result_PA_Valid = "true";
                    if ($test_hash{na}{nb} == 0) {
                        foreach my $allele (keys %test_hash) {
                            if ($allele ne "na" and $test_hash{$allele}{nb} < $nb_rep) {
                                $Result_PA_Valid = "wrong";
                            }
                        }
                    } elsif ($test_hash{na}{nb} == 1) {
                        foreach my $allele (keys %test_hash) {
                            if ($allele ne "na" and $test_hash{$allele}{nb} < $nb_rep-1) {
                                $Result_PA_Valid = "wrong"; 
                            }
                        }
                    } elsif ($test_hash{na}{nb} == 2) {
                        foreach my $allele (keys %test_hash) {
                            if ($allele ne "na" and $test_hash{$allele}{nb} < $nb_rep-2) {
                                $Result_PA_Valid = "wrong";
                            }
                        }
                    }
                    
                    #Test for inconsistency characteristics
                    $NbMissingAllele_perSite = 0;
                    $NbMissingAllele_in1rep = 0;
                    $NbMissingAllele_in2reps = 0;
                    if ($test_hash{na}{nb} == 0 and $Result_PA_Valid eq "wrong") {
                        foreach my $allele (keys %test_hash) {
                            if ($allele ne "na" and $test_hash{$allele}{nb} < $nb_rep) {
                                $NbMissingAllele_perSite += 1;
                            }
                            if ($allele ne "na" and $test_hash{$allele}{nb} == $nb_rep-3) {
                                $NbMissingAllele_in3reps += 1;
                            }
                            if ($allele ne "na" and $test_hash{$allele}{nb} == $nb_rep-2) {
                                $NbMissingAllele_in2reps += 1;
                            }   
                            if ($allele ne "na" and $test_hash{$allele}{nb} == $nb_rep-1) {
                                $NbMissingAllele_in1rep += 1;
                            }
                         
                        }
                    }
                    
                    #Record results
                    if ($test_hash{na}{nb} == 0 and $Result_PA_Valid eq "true") {
                        $report_results{"10commonPAsite_noNA"}{$iso} ++;
                        $file_content{ToReport}{$scaffold}{$position}{$iso} = 1;
                    } elsif ($test_hash{na}{nb} == 1 and $Result_PA_Valid eq "true") {
                        $report_results{"11commonPAsite_Minus1rep_1NA"}{$iso} ++;
                        if ($freedomdegree == -1) {
                            $file_content{ToReport}{$scaffold}{$position}{$iso} = 1;
                        }
                    } elsif ($test_hash{na}{nb} == 2 and $Result_PA_Valid eq "true") {
                        $report_results{"12commonPAsite_Minus2rep_2NAs"}{$iso} ++;
                        if ($freedomdegree == -2) {
                            $file_content{ToReport}{$scaffold}{$position}{$iso} = 1;
                        }                        
                    } elsif ($NbMissingAllele_perSite == 1 and $NbMissingAllele_in1rep == 1) {
                        $report_results{"13JustOneAlleleMissing_in1rep_noNA"}{$iso} ++;
                    } elsif ($NbMissingAllele_perSite == 1 and $NbMissingAllele_in2reps == 1) {
                        $report_results{"14JustOneAlleleMissing_in2rep_noNA"}{$iso} ++;                        
                    } elsif ($NbMissingAllele_perSite == 1 and $NbMissingAllele_in3reps == 1) {
                        $report_results{"15JustOneAlleleMissing_in3rep_noNA"}{$iso} ++;
                    } elsif ($NbMissingAllele_perSite > 1) {
                        $report_results{"16inconsistentPAsites"}{$iso} ++;
                    } elsif ($test_hash{na}{nb} > 0) {
                        $report_results{"17PAsites_WithfewNAsites"}{$iso} ++;
                    } else {
                        $report_results{"18PAsites_undetermined"}{$iso} ++;
                        foreach my $allele (keys %test_hash) {
                            print "Undetermined PA site:", $iso, " > ",$allele, " : ", $test_hash{$allele}{nb}, "\n";
                        }
                        print "\n";
                    }
                    $report_results{"19TotalnbPAsites"}{$iso} ++;
                }
            }                           
        }
    }
}

# Complete missing values
foreach my $class (nsort keys %report_results) {
    foreach my $iso ( keys %sampleNB) {
        $report_results{$class}{$iso} += 0
    }   
}



#warn Dumper(\%report_results); getc();

warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Report Results...\n\n";

#####....................................................................................................#####
##############################################################################################################

#print isolate names
foreach my $class (nsort keys %report_results) {
    foreach my $iso (nsort keys %{$report_results{$class}}) {
        print OUT1 "\t", $iso;
    }
    last;
}
print OUT1 "\n";

#report values for each class
foreach my $class (nsort keys %report_results) {
    print OUT1 substr($class, 2 );
    foreach my $iso (nsort keys %{$report_results{$class}}) {
        print OUT1 "\t", $report_results{$class}{$iso};
    }
    print OUT1 "\n";
}

warn "\n","*********************","\n";
##############################################################################################################
#####....................................................................................................#####

warn "Create subset of .subVCF files based on consensus positions\n\n";

#####....................................................................................................#####
##############################################################################################################

foreach my $iso ( nsort (keys %{$file_content{DATA}})) {
    foreach my $rep (nsort (keys %{$file_content{DATA}{$iso}})) {
        $NameOUT1 = $rep . ".CommonPos." . $freedomdegreeName . ($count + $freedomdegree) . "(". $freedomdegree . ")" . ".subVCF";
        open (OUT1, ">$NameOUT1")  || die "Can't open $NameOUT1: $!";
        print OUT1 $header;
        foreach my $scaffold (  nsort (keys %{$file_content{DATA}{$iso}{$rep}})) {
            foreach my $position ( nsort (keys %{$file_content{DATA}{$iso}{$rep}{$scaffold}})) {
                if ($file_content{ToReport}{$scaffold}{$position}{$iso} != 0) {
                    print OUT1 $file_content{DATA}{$iso}{$rep}{$scaffold}{$position};
                }
            }
        }
        close(OUT1);
    }
}



__END__
