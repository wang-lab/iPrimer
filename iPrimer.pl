#!/usr/bin/perl -w

# This software tool is distributed under Apache 2.0 with Commons Clause license.
# Copyrights @ University of Illinois at Chicago

use strict;
use warnings;
use diagnostics;
use lib ".";

my $inputFile = "";
my $filterFile = "";      # optional file
my $primerFile = "primerOutput.xls";

my $hasFilterFile = "n";
my $relaxed_condition = "n";
my $globalCrossDegree = 15;
my $filterCrossDegree = 16;
my $maxBlastScore = 32;

my $minSeqLength = 100;
my $maxAmpliconLength = 400;

my $minPrimerLength = 19;
my $optPrimerLength = 22;
my $maxPrimerLength = 25;

my $minF1cLength = 12;
my $maxF1cLength = 17;
my $minB1cLength = 12;
my $maxB1cLength = 17;
my $maxFIPLength = 37;
my $maxBIPLength = 37;

# Gap between adjacent primers
# F3 - F2 - FL - F1 - B1 - BL - B2 - B3
my ($minGapF3F2, $maxGapF3F2) = (3, 60);
my ($minGapF2FL, $maxGapF2FL) = (1, 20);
my ($minGapF2F1, $maxGapF2F1) = (40, 60);  # loop size including F2
my ($minGapFLF1, $maxGapFLF1) = (3, 20);

my ($minGapF1B1, $maxGapF1B1) = (1, 50);

my ($minGapB1BL, $maxGapB1BL) = (3, 20);
my ($minGapB1B2, $maxGapB1B2) = (40, 60); # loop size including B2
my ($minGapBLB2, $maxGapBLB2) = (1, 20);
my ($minGapB2B3, $maxGapB2B3) = (3, 60);

my %ff_pair;
my %fr_pair;
my %rf_pair;
my %rr_pair;

my %primer_set;

my $minGCPercent = 35;
my $maxGCPercent = 65;

my $optTm = 62;
my $minTm; # determined by $optTm and $tmRange
my $maxTm; # determined by $optTm and $tmRange
my $tmRange = 4; # range for all primers, 2x deviation from $optTm
my $stemTm = 45;
my $maxEndTmData = 46;
my $maxEndTmFilter = 46;

my $polyAT_reject_length = 5;
my $polyGC_reject_length = 4;
my $primer_self_reject_length = 5; # primer self-anneal is more stringent to avoid secondary structure
my $primer_pair_reject_length = 7; # primer pair annealing
my $sequence_self_reject_length = 9;
my $dimer_reject_length = 4;

my $primerConc = 0.00000025;  # 0.25 uM primer concentration
my $saltConc = 0.15; # 150 mM salt concentration
my ($dH, $dS, $dG, $dHi, $dSi, $dGi) = tmNeighborParam();

# minimal delta G for the five 3' end bases
my $maxEndStability = -8;

my $successCount = 0; # number of sequences that have picked primers
my $sequenceCount = 0;  # report number of processed sequences on screen

my $totalCheckedCount = 0;
my $complexityFailedCount = 0;
my $varFailedCount = 0;
my $polyXFailedCount = 0;
my $tmHighFailedCount = 0;
my $tmLowFailedCount = 0;
my $gcFailedCount = 0;
my $endStabilityFailedCount = 0;
my $selfPrimerFailedCount = 0;
my $dimerFailedCount = 0;
my $selfSequenceFailedCount = 0;
my $globalFilterFailedCount = 0;
my $globalPlusStrandFailedCount = 0;
my $globalMinusStrandFailedCount = 0;
my $endTmDataFailedCount = 0;
my $endTmFilterFailedCount = 0;
my $blastGlobalFailedCount = 0;
my $blastSelfFailedCount = 0;

my @filter = ();
my %filterUnit = ();

my @dna = ();      # original sequences
my %nmerUnit = (); # 15mer perfect words
my @masked_dna = (); # to store masked dna sequences
my @dna_var = ();  # contains sequences with base variations

my $dust_dir = "NCBI/dust";
my $blastDir = "NCBI/blast";
my $tmpDirectory = "iprimer_tmp";
my $blastDbName = "$tmpDirectory/iprimer";

my $oldTime = time();
my $time;

for (my $i = 0 ; $i < $#ARGV ; $i++) {
    if ($ARGV[$i] eq '-f' ) {
         $inputFile = $ARGV[$i+1];
         print "User selected input file: $inputFile.\n";
    }
    elsif ($ARGV[$i] eq '-r' ) {
          $relaxed_condition = $ARGV[$i+1];
          print "User selected relaxed condition: $relaxed_condition.\n";
          if ($relaxed_condition ne 'y' and $relaxed_condition ne 'yes') {
               print qq{Allowed input for -r is 'y' or 'yes'. Please try again.\n};
               exit;
          }
    }
}

if (! -e $inputFile) {
    print "The input sequence file cannot be found. Please use the correct command line (e.g., perl iPrimer.pl -f test.fasta).\n";
    exit;
}

my $variationFile = $inputFile;  # optional file

print "***Starting iPrimer to design primers for $inputFile...***\n";

unlink $primerFile if -e $primerFile;

if ($relaxed_condition eq "y" or $relaxed_condition eq "yes") {
     $globalCrossDegree = 17;
     $filterCrossDegree = 18;
     $maxBlastScore = 36;
     $tmRange = 10;
     $polyAT_reject_length = 6;
     $minGCPercent = 30;
     $maxGCPercent = 70;
     $sequence_self_reject_length = 12;
     $primer_self_reject_length = 7;
     $primer_pair_reject_length = 8;
}

print "\niPrimer is now importing your data file...\n";

@dna = importFasta($inputFile);
printStatistics(@dna);

@dna_var = importFasta($variationFile);

if ($hasFilterFile eq 'y') {
     print "\niPrimer is now including $filterFile...\n";
     @filter = importFasta($filterFile);

     print "\niPrimer is now parsing all filter sequences...\n";
     %filterUnit = screenKey(10, \@filter);
}

print "iPrimer is now parsing all input sequences...\n";
%nmerUnit = screenKey(10, \@dna);

print "iPrimer is now masking low-complexity regions...\n";
@masked_dna = maskRepeat($inputFile);
printStatistics(@masked_dna);

print "iPrimer is now reading the sequence annotations...";
for (my $index = 0; $index <= $#dna; $index++) {
     my $acc = $1 if $dna[$index]{'id'} =~ /^(\w+)/;
     $dna[$index]{'acc'} = $acc;
}
print "done.\n";

print "\nBuild index for blast program......\n";
prepareBlast();

print "All file parsing and formatting finished. Used time: ", usedTime(), "\n";
print "\niPrimer is now picking primers with the selected Tm range and degree of cross-hybridization......\n\n";

my $skip_count = 0;
for (my $index = 0; $index <= $#dna; $index++) {
     if (length($dna[$index]{'seq'}) < $minSeqLength) {
         print "$dna[$index]{'acc'} skipped as the sequence length < $minSeqLength...\n";
         next;
     }
     my $gc_content = gcPercent($dna[$index]{'seq'});
     # Set Tm range and adjust optTm is necessary
     $optTm-- if $gc_content < 0.4;
     $optTm++ if $gc_content > 0.6;

     $minTm = $optTm - $tmRange / 2;
     $maxTm = $optTm + $tmRange / 2;

     print "GC content: $gc_content; optimal Tm: $optTm; Tm range: $minTm - $maxTm.\n";
     my $primer_number = pickPrimer($index);
     print "For sequence #", ($index+1), ", $primer_number primers were picked. Now picking candidate primer sets...\n";

     print "Preparing the hash tables...";
     prepareHashTable($index);
     print "used time: ", usedTime(), "\n";

     my $primer_set_number = pickPrimerSet($index);
     print "\nFor $dna[$index]{'acc'}, $primer_set_number primer sets were designed...";
     print "used time: ", usedTime(), "\n";

     exportPrimer($index) if $primer_set_number > 0;

} # DNA loop

# remove blast-related temp files
system("rm -r $tmpDirectory");

print "\nThe primers are listed in file $primerFile in tab dilimited format.\n";
print "\n***Program ended.***\n";


#########################################################################################
# start of all subroutines
#########################################################################################

# purpose: to select an primer for each gene that meet tm and cross-hybridization requirement.
# method: for a given primer, all possible 15-mers are evaluated by using two overlapping 10mers.
#         All 15-mers should be unique in the retained primer.
sub pickPrimer {

     my ($index) = @_;
     # take the shorter length of $keyLength and $filterKeyLength
     my $shorterKeyLength = ($globalCrossDegree <= $filterCrossDegree) ? $globalCrossDegree : $filterCrossDegree;

     my $primer;   # temp variable to store a primer for check-up
     my $rPrimer;  # primer complement
     my $primer_number = 0; # a single qualified primer index
     my $position;
     POSITION: for ($position = 0; $position <= length($dna[$index]{'seq'}) - $maxPrimerLength; $position++) {

               # avoid masked regions for primer design
               $primer = substr($masked_dna[$index]{'seq'}, $position, $optPrimerLength);

               # flag for forward and reverse primers, default to 1 (true)
               my $forwardOK = 1;
               my $reverseOK = 1;

               $totalCheckedCount++;

               # first apply less CPU-intensive filters. The most time-consuming tasks are performed later

               ###########################################################################################
               # check for non-atcg characters
               ###########################################################################################

               if ($primer =~ /[^atcg]/) {
                    $complexityFailedCount++;
                    next POSITION;
               }

               ##########################################################################################
               # test for tm range
               ##########################################################################################

               my $tm = tmNeighbor($primer);
               # out of Tm range and < $optTm
               if ($tm < $minTm) {
                    $primer = substr($dna[$index]{'seq'}, $position, $optPrimerLength + 1);
                    next POSITION if substr($primer, -1, 1) =~ /[^atcg]/;
                    if (tmNeighbor($primer) < $minTm) {
                        $primer = substr($dna[$index]{'seq'}, $position, $optPrimerLength + 2);
                        next POSITION if substr($primer, -1, 1) =~ /[^atcg]/;
                        if (tmNeighbor($primer) < $minTm) {
                            $primer = substr($dna[$index]{'seq'}, $position, $optPrimerLength + 3);
                            next POSITION if substr($primer, -1, 1) =~ /[^atcg]/;
                            if (tmNeighbor($primer) < $minTm) {
                                    $tmLowFailedCount++;
                                    next POSITION;
                            }
                        }
                    }
               }
               # out of Tm range and > $optTm
               elsif ($tm > $maxTm) {
                       $primer = substr($dna[$index]{'seq'}, $position, $optPrimerLength - 1);
                       if (tmNeighbor($primer) > $maxTm) {
                           $primer = substr($dna[$index]{'seq'}, $position, $optPrimerLength - 2);
                           if (tmNeighbor($primer) > $maxTm) {
                               $primer = substr($dna[$index]{'seq'}, $position, $optPrimerLength - 3);
                               if (tmNeighbor($primer) > $maxTm) {
                                       $tmHighFailedCount++;
                                       next POSITION;
                               }
                           }
                       }
               }

               ###########################################################################################
               # check for base variations
               ###########################################################################################

               my $primer_var = substr($dna_var[$index]{'seq'}, $position, length($primer));
               if ($primer ne $primer_var) {
                    print "$primer_var skipped due to genome variation...\n";
                    $varFailedCount++;
                    next POSITION;
               }

               ##########################################################################################
               # check for poly-X bases
               ##########################################################################################

               my $polyX_flag = checkPolyX($primer);
               if ($polyX_flag) {
                    $polyXFailedCount++;
                    next POSITION;
               }

               ##########################################################################################
               # check for %GC
               ##########################################################################################

               my $GC_flag = checkGC($primer);
               if ($GC_flag) {
                    $gcFailedCount++;
                    next POSITION;
               }

               ##########################################################################################
               # check the last 5 end-bases (3') for stability (delta G)
               ##########################################################################################

               $rPrimer = dnaComplement($primer);

               $forwardOK = 0 if endStability($primer) < $maxEndStability;   # forward primer
               $reverseOK = 0 if endStability($rPrimer) < $maxEndStability;  # reverse primer

               if (!$forwardOK && !$reverseOK) {
                    $endStabilityFailedCount++;
                    next POSITION;
               }

               ##########################################################################################
               # primer self-annealing filters: primer-primer, primer-sequence, and primer end dimer
               ##########################################################################################

               # position where this self-annealing stretch is found
               my $self_anneal_jump = selfAnnealTest($primer, $primer, $primer_self_reject_length);

               # primer-primer self annealing filter failed
               if ($self_anneal_jump >= 0) {
                    $position += $self_anneal_jump;
                    $selfPrimerFailedCount += $self_anneal_jump + 1;
                    $totalCheckedCount += $self_anneal_jump;
                    next POSITION;
               }

               # if did not find any self match for primer-primer, continue to test primer-sequence self annealing
               $self_anneal_jump = selfAnnealTest($primer, $dna[$index]{'seq'}, $sequence_self_reject_length);
               # primer-sequence self annealing filter failed
               if ($self_anneal_jump >= 0) {
                    $position += $self_anneal_jump;
                    $selfSequenceFailedCount += $self_anneal_jump + 1;
                    $totalCheckedCount += $self_anneal_jump;
                    next POSITION;
               }

               # test for primer end dimer
               # $self_anneal_jump should be 0 if a dimer is found, equivalent to >= 0 (confirmed)
               if ($forwardOK) {
                    my $endBase = substr($primer, -$dimer_reject_length, $dimer_reject_length);
                    $self_anneal_jump = selfAnnealTest($endBase, $primer, $dimer_reject_length);
                    $forwardOK = 0 if $self_anneal_jump == 0;
               }
               if ($reverseOK) {
                    my $endBase = substr($rPrimer, -$dimer_reject_length, $dimer_reject_length);
                    $self_anneal_jump = selfAnnealTest($endBase, $rPrimer, $dimer_reject_length);
                    $reverseOK = 0 if $self_anneal_jump == 0;
               }
               if (!$forwardOK && !$reverseOK) {
                    $dimerFailedCount++;
                    next POSITION;
               }

               ###############################################################################################
               # test for global n-mer match
               ###############################################################################################

               my $unitIndex = -1;
               # step through the 10-mer unitKeys within the primer(from the 3'-end first)
               for (my $keyIndex = length($primer) - $shorterKeyLength; $keyIndex >= 0 ; $keyIndex--) {
                    my $n_mer = substr($primer, $keyIndex, $globalCrossDegree); # 15-mer
                    # do not test <15mer, because this would lead to problems in the sub
                    $unitIndex = nmerTest($n_mer, $index, \%nmerUnit, \@dna) if length($n_mer) == $globalCrossDegree;
                    if ($unitIndex >= 0) {
                         if ($unitIndex == 0.1) {
                              $globalPlusStrandFailedCount += $keyIndex + 1;
                         }
                         else {
                               $globalMinusStrandFailedCount += $keyIndex + 1;
                         }
                         $totalCheckedCount += $keyIndex;

                         # the start position of the 15-mer
                         $position += $keyIndex;
                         next POSITION;
                    }

                    if ($hasFilterFile eq 'y') {
                         $n_mer = substr($primer, $keyIndex, $filterCrossDegree); # 13-mer
                         $unitIndex = nmerTest($n_mer, $index, \%filterUnit, \@filter) if length($n_mer) == $filterCrossDegree;
                         if ($unitIndex >= 0 ) {
                              # the end position of the 13-mer
                              $position += $keyIndex; # the start position of the 14mer
                              $globalFilterFailedCount += $keyIndex + 1;

                              $totalCheckedCount += $keyIndex;
                              next POSITION;
                         }
                    }

               } # end for keyIndex

#                ################################################################################################
#                # test for end Tm match
#                ################################################################################################

               if ($forwardOK) {
                    $unitIndex = endTmForward($primer, $index, $position, $maxEndTmData, $globalCrossDegree, \%nmerUnit, \@dna);
                    $forwardOK = 0 if $unitIndex >= 0;
               }
               if ($reverseOK) {
                    $unitIndex = endTmReverse($primer, $index, $position, $maxEndTmData, $globalCrossDegree, \%nmerUnit, \@dna);
                    $reverseOK = 0 if $unitIndex >= 0;
               }
               if (!$forwardOK && !$reverseOK) {
                    $endTmDataFailedCount++;
                    next POSITION;
               }

               if ($hasFilterFile eq 'y') {
                    if ($forwardOK) {
                         $unitIndex = endTmForward($primer, $index, $position, $maxEndTmFilter, $filterCrossDegree, \%filterUnit, \@filter);
                         $forwardOK = 0 if $unitIndex >= 0;
                    }
                    if ($reverseOK) {
                         $unitIndex = endTmReverse($primer, $index, $position, $maxEndTmFilter, $filterCrossDegree, \%filterUnit, \@filter);
                         $reverseOK = 0 if $unitIndex >= 0;
                    }
                    if (!$forwardOK && !$reverseOK) {
                         # the end position of the end n-mer
                         $endTmFilterFailedCount++;
                         next POSITION;
                    }
               }

               ################################################################################################
               # test for global Blast score
               ################################################################################################

               my $blastPosition = blast($primer, $index);
               if ($$blastPosition[0] >= 0) {  # 0 index for start position
                  $position += $$blastPosition[0] - 1;
                  $blastGlobalFailedCount += $$blastPosition[0];

                  $totalCheckedCount += $$blastPosition[0] - 1;
                  next POSITION;
               }

               ################################################################################################
               ### test for sequence self-annealing Blast score
               ################################################################################################

               my $annealPosition = anneal($primer, $index);
               if ($$annealPosition[0] >= 0) {
                  $position += $$annealPosition[0] - 1;
                  $blastSelfFailedCount += $$annealPosition[0];

                  $totalCheckedCount += $$annealPosition[0] - 1;
                  next POSITION;
               }

               #***********************************************************************************************
               # All screenings are ok, recording info and proceed to the next Position
               #***********************************************************************************************
               $dna[$index]{'primerSeq'}[$primer_number] = $primer;
               $dna[$index]{'rPrimerSeq'}[$primer_number] = dnaComplement($primer);
               $dna[$index]{'primerPos'}[$primer_number] = $position;
               $dna[$index]{'forward'}[$primer_number] = $forwardOK;
               $dna[$index]{'reverse'}[$primer_number] = $reverseOK;

               my $primer_stem = truncate_primer($primer, $stemTm);
               $dna[$index]{'stemSeq'}[$primer_number] = $primer_stem;
               $dna[$index]{'rStemSeq'}[$primer_number] = dnaComplement($primer_stem);
               $primer_number++;

     } # POSITION loop

     return $primer_number;

} # end of sub pickPrimer

# Prepare the hash table for compatible pairs
sub prepareHashTable {

    my ($index) = @_; # gene index

    # positions for all qualified primers for the current sequence
    my @position = @{$dna[$index]{'primerPos'}};
    # sequences for all qualified primers for the current sequence
    my @primerSeq = @{$dna[$index]{'primerSeq'}};
    my @rPrimerSeq = @{$dna[$index]{'rPrimerSeq'}};

    # flag array for forward primers
    my @forward = @{$dna[$index]{'forward'}};
    # flag array for reverse primers
    my @reverse = @{$dna[$index]{'reverse'}};

    for (my $i = 0; $i < $#position; $i++) {
        my $primer_i = $primerSeq[$i];
        my $primer_ir = $rPrimerSeq[$i];
        for (my $j = $i + 1; $j <= $#position; $j++) {
            my $primer_j = $primerSeq[$j];
            my $primer_jr = $rPrimerSeq[$j];
            next if $position[$j] - $position[$i] + length($primer_j) > $maxAmpliconLength;
            if ($forward[$i] == 1 and $forward[$j] == 1) {
                 my $dimer_flag = checkPrimerPair($primer_i, $primer_j);
                 $ff_pair{$i}{$j} = 1 unless $dimer_flag;
            }
            if ($forward[$i] == 1 and $reverse[$j] == 1) {
                 my $dimer_flag = checkPrimerPair($primer_i, $primer_jr);
                 $fr_pair{$i}{$j} = 1 unless $dimer_flag;
            }
            if ($reverse[$i] == 1 and $forward[$j] == 1) {
                 my $dimer_flag = checkPrimerPair($primer_ir, $primer_j);
                 $rf_pair{$i}{$j} = 1 unless $dimer_flag;
            }
            if ($reverse[$i] == 1 and $reverse[$j] == 1) {
                 my $dimer_flag = checkPrimerPair($primer_ir, $primer_jr);
                 $rr_pair{$i}{$j} = 1 unless $dimer_flag;
            }
         }
    }
}

# picker at least one pair of primers
sub pickPrimerSet {

    my ($index) = @_; # gene index

    # positions for all qualified primers for the current sequence
    my @position = @{$dna[$index]{'primerPos'}};
    # sequences for all qualified primers for the current sequence
    my @primerSeq = @{$dna[$index]{'primerSeq'}};
    my @rPrimerSeq = @{$dna[$index]{'rPrimerSeq'}};
    my @stemSeq = @{$dna[$index]{'stemSeq'}};
    my @rStemSeq = @{$dna[$index]{'rStemSeq'}};

    my %FIP;
    my %BIP;
    my %F2;
    my %B2;
    my %F2_picked;
    my %B2_picked;
    my %setB;
    my %FB_hash; # map F and B primer sets
    my $setB_count = 0;
    my $setF_count = 0;
    my $primer_set_count = 0;

    print "\nStart selecting candidate primer sets...\n";
    # Loop for B1 primer
    B1: for (my $B1_i = 4; $B1_i < scalar(@position) - 3; $B1_i++) {
        my $primer_B1c = $stemSeq[$B1_i];
        next if length($primer_B1c) < $minB1cLength or length($primer_B1c) > $maxB1cLength;

        # Loop for BL primer
        BL: foreach my $BL_i (sort {$a<=>$b} keys %fr_pair) {
            my $primer_BL = $primerSeq[$BL_i];
            my $gapB1BL = $position[$BL_i] - $position[$B1_i] - length($primer_B1c);
            next if $gapB1BL < $minGapB1BL;
            last if $gapB1BL > $maxGapB1BL;

            # Loop for B2 primer
            B2: foreach my $B2_i (sort {$a<=>$b} keys %{$fr_pair{$BL_i}}) {
                next unless exists $rr_pair{$B2_i};
                my $primer_B2 = $rPrimerSeq[$B2_i];
                my $gapBLB2 = $position[$B2_i] - $position[$BL_i] - length($primer_BL);
                next if $gapBLB2 < $minGapBLB2;
                last if $gapBLB2 > $maxGapBLB2;
                my $gapB1B2 = $position[$B2_i] - $position[$B1_i] - length($primer_B1c) + length($primer_B2); # including B2 length
                next if $gapB1B2 < $minGapB1B2;
                last if $gapB1B2 > $maxGapB1B2;

                my $B2_base = substr($primer_B2, 0, 1);
                my $ext_base = substr($dna[$index]{'seq'}, $position[$B1_i] + length($primer_B1c), 1);
                next if $B2_base eq $ext_base;

                my $primer_BIP = join("", ($primer_B1c, $primer_B2));
                next if length($primer_BIP) > $maxBIPLength;

                my $polyX_flag1 = checkPolyX($primer_BIP);
                next if $polyX_flag1;
                my $GC_flag1 = checkGC($primer_BIP);
                next if $GC_flag1;

                my $dimer_flag1 = checkPrimerSelf($primer_BIP);
                my $dimer_flag2 = checkPrimerPair($primer_BIP, $primer_BL);
                next if $dimer_flag1 or $dimer_flag2;

                $BIP{$primer_BIP} = 1;
                $B2{$B2_i} = $position[$B2_i];

                # Loop for B3 primer
                B3: foreach my $B3_i (sort {$a<=>$b} keys %{$rr_pair{$B2_i}}) {
                    next unless defined $fr_pair{$BL_i}{$B3_i};
                    my $primer_B3 = $rPrimerSeq[$B3_i];
                    my $gapB2B3 = $position[$B3_i] - $position[$B2_i] - length($primer_B2);
                    next if $gapB2B3 < $minGapB2B3;
                    last if $gapB2B3 > $maxGapB2B3;

                    my $dimer_flag3 = checkPrimerPair($primer_BIP, $primer_B3);
                    next if $dimer_flag3;

                    ##### Passed all position tests, record results here
                    #*********************************************************************

                    $setB{$B1_i}{$BL_i}{$B2_i}{$B3_i} = $primer_BIP;
                    $setB_count++;
                    #*********************************************************************

                } # Loop for B3
            } # Loopf for B2
        }  # Loop for BL
    }  # Loop for B1

    #print "Set B count: $setB_count; unique BIP primers: ", scalar(keys %BIP), "; unique B2 primers: ", scalar(keys %B2), " *** \n", join(", ", (sort {$a<=>$b} values %B2)), ".\n";
    print "Set B count: $setB_count; unique BIP primers: ", scalar(keys %BIP), "; unique B2 primers: ", scalar(keys %B2), ".\n";
    print "$setB_count primer B sets designed.\n";

    # Loop for F1 primer
    F1: for (my $F1_i = 3; $F1_i < scalar(@position) - 4; $F1_i++) {
        my $primer_F1c = $rStemSeq[$F1_i];
        next if length($primer_F1c) < $minF1cLength or length($primer_F1c) > $maxF1cLength;

        # Loop for B1 primer
        B1: for (my $B1_i = $F1_i + 1; $B1_i < scalar(@position) - 3; $B1_i++) {
            my $gapF1B1 = $position[$B1_i] - $position[$F1_i] - length($primer_F1c);
            next if $gapF1B1 < $minGapF1B1;
            last if $gapF1B1 > $maxGapF1B1;
            next if !exists $setB{$B1_i};
            $FB_hash{$F1_i}{$B1_i} = 1;
        }
    }

    # Loop for F3 primer
    F3: foreach my $F3_i (sort {$a<=>$b} keys %ff_pair) {
        next unless exists $fr_pair{$F3_i};
        my $primer_F3 = $primerSeq[$F3_i];

        # Loop for F2 primer
        F2: foreach my $F2_i (sort {$a<=>$b} keys %{$ff_pair{$F3_i}}) {
            # next if current primer is already designed
            next if exists $F2_picked{$F2_i};
            next unless exists $ff_pair{$F2_i} and exists $fr_pair{$F2_i};

            my $primer_F2 = $primerSeq[$F2_i];
            my $gapF3F2 = $position[$F2_i] - $position[$F3_i] - length($primer_F3);
            next if $gapF3F2 < $minGapF3F2;
            last if $gapF3F2 > $maxGapF3F2;

            # Loop for FL primer
            FL: foreach my $FL_i (sort {$a<=>$b} keys %{$fr_pair{$F2_i}}) {
                next unless defined $fr_pair{$F3_i}{$FL_i} and exists $rr_pair{$FL_i} and exists $rf_pair{$FL_i};
                my $primer_FL = $rPrimerSeq[$FL_i];
                my $gapF2FL = $position[$FL_i] - $position[$F2_i] - length($primer_F2);
                next if $gapF2FL < $minGapF2FL;
                last if $gapF2FL > $maxGapF2FL;

                # Loop for F1 primer
                F1: for (my $F1_i = $FL_i + 1; $F1_i < scalar(@position) - 4; $F1_i++) {
                    my $primer_F1 = $stemSeq[$F1_i];
                    my $gapFLF1 = $position[$F1_i] - $position[$FL_i] - length($primer_FL);
                    next if $gapFLF1 < $minGapFLF1;
                    last if $gapFLF1 > $maxGapFLF1;
                    my $gapF2F1 = $position[$F1_i] - $position[$F2_i]; # including F2 length
                    next if $gapF2F1 < $minGapF2F1;
                    last if $gapF2F1 > $maxGapF2F1;

                    next unless exists $FB_hash{$F1_i};

                    my $F2_base = substr($rPrimerSeq[$F2_i], -1, 1);
                    my $ext_base = substr($dna[$index]{'seq'}, $position[$F1_i] - 1, 1);
                    next if $F2_base eq $ext_base;

                    my $primer_F1c = $rStemSeq[$F1_i];
                    my $primer_FIP = join("", ($primer_F1c, $primer_F2));
                    next if length($primer_FIP) > $maxFIPLength;

                    my $polyX_flag2 = checkPolyX($primer_FIP);
                    next if $polyX_flag2;
                    my $GC_flag2 = checkGC($primer_FIP);
                    next if $GC_flag2;

                    my $dimer_flag4 = checkPrimerSelf($primer_FIP);
                    my $dimer_flag5 = checkPrimerPair($primer_FIP, $primer_FL);
                    my $dimer_flag6 = checkPrimerPair($primer_FIP, $primer_F3);
                    next if $dimer_flag4 or $dimer_flag5 or $dimer_flag6;

                    $FIP{$primer_FIP} = 1;
                    $F2{$F2_i} = $position[$F2_i];
                    $setF_count++;

                    foreach my $B1_i (sort {$a<=>$b} keys %{$FB_hash{$F1_i}}) {

                        # $setB{$B1_i}{$BL_i}{$B2_i}{$B3_i} = 1;
                        foreach my $BL_i (keys %{$setB{$B1_i}}) {
                            next unless defined $ff_pair{$F3_i}{$BL_i} and defined $ff_pair{$F2_i}{$BL_i} and defined $rf_pair{$FL_i}{$BL_i};

                            foreach my $B2_i (keys %{$setB{$B1_i}{$BL_i}}) {
                                # skip current primer is already designed
                                next if exists $B2_picked{$B2_i};
                                next unless defined $fr_pair{$F3_i}{$B2_i} and defined $fr_pair{$F2_i}{$B2_i} and defined $rr_pair{$FL_i}{$B2_i};

                                foreach my $B3_i (keys %{$setB{$B1_i}{$BL_i}{$B2_i}}) {
                                    next unless defined $fr_pair{$F3_i}{$B3_i} and defined $fr_pair{$F2_i}{$B3_i} and defined $rr_pair{$FL_i}{$B3_i};

                                    my $primer_BIP = $setB{$B1_i}{$BL_i}{$B2_i}{$B3_i};
                                    my $dimer_flag7 = checkPrimerPair($primer_BIP, $primer_F3);
                                    my $dimer_flag8 = checkPrimerPair($primer_BIP, $primer_FIP);
                                    my $dimer_flag9 = checkPrimerPair($primer_BIP, $primer_FL);
                                    next if $dimer_flag7 or $dimer_flag8 or $dimer_flag9;

                                    $F2_picked{$F2_i} = [$F3_i, $F2_i, $FL_i,$F1_i, $B1_i, $BL_i, $B2_i, $B3_i, $primer_FIP, $primer_BIP];
                                    $B2_picked{$B2_i} = $F2_i;
                                    push @{$dna[$index]{'primer_set'}}, $F2_picked{$F2_i};                                     #
                                    $primer_set_count++;
                                    next F2;
                                }
                            }
                        }
                    }
                }  # Loop for F1
            }  # Loop for FL
        }  # Loop for F2
    }   # Loop for F3

    print "Set F count: $setF_count; unique FIP primers: ", scalar(keys %FIP), "; unique F2 primers: ", scalar(keys %F2), ".\n";
    print "$setF_count primer F sets designed.\n";

    return $primer_set_count;
}


sub checkPrimerPair {
    my ($primer1, $primer2) = @_;

    # check the end bases of primer1 against the f2 primer sequence first
    my $endBase = substr($primer1, -$dimer_reject_length, $dimer_reject_length);
    my $dimer_pair = selfAnnealTest($endBase, $primer2, $dimer_reject_length);

    # $primer1 dimer OK, continue to check $primer2
    if ($dimer_pair < 0) {
         $endBase = substr($primer2, -$dimer_reject_length, $dimer_reject_length);
         $dimer_pair = selfAnnealTest($endBase, $primer1, $dimer_reject_length);

         # if f2 dimer OK, continue to check primer-primer annealing
         if ($dimer_pair < 0) {
              # position where this self-annealing stretch is found
              # check for self-annealing of 5' and 3' primer sequences
              $dimer_pair = selfAnnealTest($primer1, $primer2, $primer_pair_reject_length);

              if ($dimer_pair < 0) {
                   return 0;
              }
         }
    }
    return 1;
}

sub checkPrimerSelf {
    my ($primer) = @_;

    # check the end bases of primer against the primer sequence first
    my $endBase = substr($primer, -$dimer_reject_length, $dimer_reject_length);
    my $dimer = selfAnnealTest($endBase, $primer, $dimer_reject_length);

    # $primer dimer OK, continue to check $primer2
    if ($dimer < 0) {
         # check for self-annealing of 5' and 3' primer sequences
         my $dimer1 = selfAnnealTest($primer, $primer, $primer_self_reject_length);
         if ($dimer1 < 0) {
              return 0;
         }
    }
    return 1;
}

# self-annealing test for both primer-primer and primer-sequence
sub selfAnnealTest {

    my ($primer, $sequence, $reject_length) = @_;

    # reverse the primer instead of the whole sequence because primer is much shorter and it takes less time
    my $primer_complement = dnaComplement($primer);
    my (%forward_hash, @reverse_array);

    # fill in two hashes for primer and sequence/primer
    for (my $i = length($primer) - $reject_length; $i >= 0; $i--) {
         push @reverse_array, substr($primer_complement, $i, $reject_length);
    }
    for (my $j = 0; $j <= length($sequence) - $reject_length; $j++) {
         $forward_hash{substr($sequence, $j, $reject_length)} = 1;
    }
    # for every key in the primer hash, check in sequential order
    foreach (@reverse_array) {
          if (exists $forward_hash{$_}) {
              # $position -= $oligoLength - ($keyIndex + $keyLength) from OligoPicker 3' end
              my $jumpLen = length($primer) - index($primer_complement, $_) - $reject_length;
              return $jumpLen;
          }
    }
    return -1;
}

# check for poly-X bases
sub checkPolyX {
    my $primer = shift;
    my $flag = 0;
    if ($primer =~ /a{$polyAT_reject_length}/ || $primer =~ /c{$polyGC_reject_length}/ ||
       $primer =~ /g{$polyGC_reject_length}/ || $primer =~ /t{$polyAT_reject_length}/
       ) {
         $flag = 1;;
       }
       return $flag;
}


# check for GC content
sub checkGC {
    my $primer = shift;
    my $gcContent = 100 * gcPercent($primer);
    my $flag = 0;
    if ($gcContent > $maxGCPercent || $gcContent < $minGCPercent) {
         $flag = 1;;
       }
       return $flag;
}


# design stem primer (shortened based on $stemTm)
sub truncate_primer {
    my ($primer, $tm_threshold) = @_;
    my $tm = tmNeighbor($primer);
    my $last_primer = $primer;
    while ($tm > $tm_threshold) {
          $last_primer = $primer;
          $primer = substr($primer, 0, length($primer) - 1);
          $tm = tmNeighbor($primer);
    }
    return $last_primer;
}


sub maskRepeat {

    my $file = shift;
    my $maskedFile = "masked_" . $file;
    system("$dust_dir/dust $file > $maskedFile");

    my @masked_seq = importFasta($maskedFile);

    unlink $maskedFile if -e $maskedFile; # the masked fasta file

    print "Sequence file has been masked and imported. Number of DNA sequences: ", scalar(@masked_seq), "\n";

    return @masked_seq;
}


# N contiguous base match; n_mer is 15 here
sub nmerTest {

    my ($n_mer, $index, $unitRef, $seqRef) = @_;

    my $unitSeq = substr($n_mer, 0, 10);  # 10 bases, smaller than $keyLength
    my $unitIndex = nmerTestOneStrand($n_mer, $index, $unitSeq, $unitRef, $seqRef, 1);

    # if did not find any match for this 15mer, continue to test for its complementary strand
    if ($unitIndex < 0) {
         $n_mer = dnaComplement($n_mer);
         $unitSeq = substr($n_mer, 0, 10);

         # unitSeq is now the last 10 bases of the 15mer; the order is reversed
         # There should be no match if to the 15mer itself, so keyIndex is not important (ignored)
         $unitIndex = nmerTestOneStrand($n_mer, $index, $unitSeq, $unitRef, $seqRef, 2);
    }
    else {
          $unitIndex = 0.1; # plus strand failed, value used for the plus strand failed counter
    }
    return $unitIndex;
}

# internally used by nmerTest sub
sub nmerTestOneStrand {

    # position and keyIndex are used to find the exact n-mer location in the current seq
    my ($n_mer, $index, $unitSeq, $unitRef, $seqRef, $strand_orientation) = @_;

    my @geneList;   # array for unitSeq gene indexes and positions

    # almost always exists
    if (exists $$unitRef{$unitSeq}) {

         # create an array that contain indices for specific $keyUnits
         @geneList = split(/ /, $$unitRef{$unitSeq});    # contain gene indexes

         # step through each gene index in the geneList array
         foreach my $unitIndex (@geneList) {
              # skip if from the same gene and the strand orientation is the same;
              if (($strand_orientation == 1) and ($$seqRef[$unitIndex]{'id'} eq $dna[$index]{'id'})) {
                  next;
              }
              # Other sequences examined
              elsif ($$seqRef[$unitIndex]{'seq'} =~ /$n_mer/) {
                    return 1;
              }
         }

    } # end if

    return -1; # not found
}

# the end perfect match Tm can not exceed a limit
sub endTmForward {

    my ($primer, $index, $position, $maxEndTm, $crossDegree, $unitRef, $seqRef) = @_;

    my $endLength = endLength($primer, $maxEndTm);

    # no need to check if the endLength is already longer than crossDegree (tested before)
    return -1 if $endLength >= $crossDegree;

    my $n_mer = substr($primer, -$endLength, $endLength);
    my $unitIndex = nmerTest($n_mer, $index, $unitRef, $seqRef);

    return $unitIndex;

}

# the end perfect match Tm can not exceed a limit
sub endTmReverse {

    my ($primer, $index, $position, $maxEndTm, $crossDegree, $unitRef, $seqRef) = @_;

    my $rPrimer = dnaComplement($primer);
    my $endLength = endLength($rPrimer, $maxEndTm);

    # no need to check if the endLength is already longer than crossDegree (tested before)
    return -1 if $endLength >= $crossDegree;

    my $n_mer = substr($primer, 0, $endLength);
    my $unitIndex = nmerTest($n_mer, $index, $unitRef, $seqRef);

    return $unitIndex;

}

sub endLength {

    my ($primer, $maxEndTm) = @_;
    my $endTm;
    my $endLength;
    for ($endLength = 10; $endLength <= $minPrimerLength; $endLength++) {
          $endTm = tmNeighbor(substr($primer, -$endLength, $endLength));
          last if $endTm > $maxEndTm;
    }
    return $endLength;
}

# tm calculation using the Nearest Neighbor Method
sub tmNeighbor {

    my ($oligo) = @_;

    my $H = 0;
    my $S = 0;

    for (my $position = 0; $position < length($oligo)-1; $position++) {

         my $dimer = substr($oligo, $position, 2);
         $H += $$dH{$dimer};
         $S += $$dS{$dimer};
    }

    my $firstBase = substr($oligo, 0, 1);
    my $lastBase = substr($oligo, -1, 1);

    #initiation energy
    $H += $$dHi{$firstBase} + $$dHi{$lastBase};
    $S += $$dSi{$firstBase} + $$dSi{$lastBase};

    # salt correction
    $S = $S + 0.368 * (length($oligo) - 1) * log($saltConc);

    my $tm = 1000 * $H / ($S + 1.9872 * log($primerConc / 4)) - 273.15;

    return sprintf("%.1f", $tm);
}

# Check the stability of the last 5 3' end-bases by calculating the delta G
sub endStability {

    my ($oligo) = @_;

    my $G = 0;

    # count the last five bases
    for (my $position = 2; $position <= 6; $position++) {

         my $dimer = substr($oligo, -$position, 2);
         $G += $$dG{$dimer};
    }
    return $G;
}


sub tmNeighborParam {

    my %dH = ();
    my %dS = ();
    my %dG = ();

    my %dHi = ();
    my %dSi = ();
    my %dGi = ();

    ($dH{'aa'}, $dS{'aa'}, $dG{'aa'}) = (-7.9, -22.2, -1);
    ($dH{'tt'}, $dS{'tt'}, $dG{'tt'}) = ($dH{'aa'}, $dS{'aa'}, $dG{'aa'});
    ($dH{'at'}, $dS{'at'}, $dG{'at'}) = (-7.2,-20.4,-0.88);
    ($dH{'ta'}, $dS{'ta'}, $dG{'ta'}) = (-7.2, -21.3, -0.58);
    ($dH{'ca'}, $dS{'ca'}, $dG{'ca'}) = (-8.5, -22.7, -1.45);
    ($dH{'tg'}, $dS{'tg'}, $dG{'tg'}) = ($dH{'ca'}, $dS{'ca'}, $dG{'ca'});
    ($dH{'ct'}, $dS{'ct'}, $dG{'ct'}) = (-7.8, -21, -1.28);
    ($dH{'ag'}, $dS{'ag'}, $dG{'ag'}) = ($dH{'ct'}, $dS{'ct'}, $dG{'ct'});
    ($dH{'ga'}, $dS{'ga'}, $dG{'ga'}) = (-8.2, -22.2, -1.3);
    ($dH{'tc'}, $dS{'tc'}, $dG{'tc'}) = ($dH{'ga'}, $dS{'ga'}, $dG{'ga'});
    ($dH{'gt'}, $dS{'gt'}, $dG{'gt'}) = (-8.4, -22.4, -1.44);
    ($dH{'ac'}, $dS{'ac'}, $dG{'ac'}) = ($dH{'gt'}, $dS{'gt'}, $dG{'gt'});
    ($dH{'cg'}, $dS{'cg'}, $dG{'cg'}) = (-10.6, -27.2, -2.17);
    ($dH{'gc'}, $dS{'gc'}, $dG{'gc'}) = (-9.8, -24.4, -2.24);
    ($dH{'gg'}, $dS{'gg'}, $dG{'gg'}) = (-8, -19.9, -1.84);
    ($dH{'cc'}, $dS{'cc'}, $dG{'cc'}) = ($dH{'gg'}, $dS{'gg'}, $dG{'gg'});


    ($dHi{'a'}, $dHi{'g'}, $dSi{'a'}, $dSi{'g'}, $dGi{'a'}, $dGi{'g'}) = (2.3, 0.1, 4.1, -2.8, 1.03, 0.98);
    ($dHi{'t'}, $dHi{'c'}, $dSi{'t'}, $dSi{'c'}, $dGi{'t'}, $dGi{'c'}) = ($dHi{'a'}, $dHi{'g'}, $dSi{'a'}, $dSi{'g'}, $dGi{'a'}, $dGi{'g'});

    return (\%dH, \%dS, \%dG, \%dHi, \%dSi, \%dGi);

}

sub prepareBlast {

    mkdir $tmpDirectory;

    my $blastFile = $inputFile;
    $blastFile = "'$blastFile $filterFile'" if -e $filterFile;

    #     -i Input file for formatting
    #     -p Type of file: T - protein; F - nucleotide
    #     -o Parse options: T - true: parse seqId and create indexes; F - false: do not parse
    #     -n  Base name for BLAST files [String]  Optional

    system("$blastDir/formatdb -i $blastFile -p F -o F -n $blastDbName");
    print ("Blast preparation is done!\n");
}

# return the failed position (the end position of the blasted region) if any
sub blast {

    my ($oligo, $index) = @_;

    #>NC_045512.2:21563-25384 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
    my $acc = $dna[$index]{'acc'};
    my $fasta = ">$acc\n$oligo";
    #     -p Program Name: blastn, blastp, etc
    #     -d Database: must be first formatted with formatdb
    #     -i Query File: default - stdin
    #     -e Expectation value
    #     -o Blast output file: default - stdout
    #     -F Filter query sequence (F/T): DUST with blastn
    #     -m Output format: 8 for tab
    #     -I  Show GI's in deflines [T/F]  default = F
    #     -S  Query strands to search against database (for blast[nx], and tblastx).  3 is both, 1 is top, 2 is bottom [Integer] default = 3

    my @result;

    open (BLAST, "echo '$fasta' | $blastDir/blastall -p blastn -m 8 -F F -S 3 -e 1000 -d $blastDbName |");
    my @field = ();
    while (<BLAST>) {

            s/\s+$//;

            # the blast result format is as following:
            # gfp_coding      gfp_coding      100.00     20      0       0       1       20      1       20      8.5e-04   40.14
            @field = split "\t", $_;
            my $hit_acc = $field[1];
            $hit_acc = $1 if $hit_acc =~ /^(\w+)/;

            # Next if the hit is itself; may not exist if mapped to non-coding RNAs
            if (($hit_acc eq $acc) && ($field[2] eq '100.00') && ($field[3] == length($oligo))) {
                 next;
            }
            else {
                  if ($field[11] > $maxBlastScore) {
                      #print "$_", "\n";
                      return [$field[6], $field[7]];
                  }
                  last; # only need to check the first non-self score since the result is sorted by score!
            }
    }
    close(BLAST);
    return [-1, -1];
}

sub anneal {

    my ($primer, $index) = @_;
    my $sequence = $dna[$index]{'seq'};
    $sequence =~ s/(.{70})/$1\n/g; # to make fasta format

    my $oligoFile = "$tmpDirectory/query_oligo.$$";
    my $seqFile = "$tmpDirectory/query_seq.$$";

    open (OLIGO, ">$oligoFile") || die("Can not open file for writing!\n"); # for writing
    print OLIGO ">oligo\n$primer\n";
    close(OLIGO);

    open (SEQ, ">$seqFile") || die("Can not open file for writing!\n"); # for writing
    print SEQ ">sequence\n$sequence";
    close(SEQ);


    #   -i  First sequence [File In]
    #   -j  Second sequence [File In]
    #   -p  Program name: blastp, blastn, blastx, tblastn, tblastx. For blastx 1st sequence should be nucleotide, tblastn 2nd sequence nucleotide [String]
    #   -o  alignment output file [File Out]  default = stdout
    #   -F  Filter query sequence (DUST with blastn, SEG with others) [String] default = T
    #   -e  Expectation value (E) [Real] default = 10.0
    #   -S  Query strands to search against database (blastn only).  3 is both, 1 is top, 2 is bottom [Integer] default = 3
    #   -D  Output format: 0 - traditional, 1 - tabulated [Integer]  default = 0
    #   -W  Wordsize (zero invokes default behavior) [Integer] default = 0

    # tmpseq_0        tmpseq_1        9.09    11      10      0       36      45      2970    2980    0.054   22.30

    # `system command` will not return anything if there is no output (!exists); system("command") returns 0 if no output
    my @result = `$blastDir/bl2seq -i $oligoFile -j $seqFile -p blastn -S 2 -D 1 -F F -W 7 -e 10000`;

    if (exists $result[0]) {

        my @field = split "\t", $result[0];
        if ($field[11] > 2 * $sequence_self_reject_length) {
             #print $result[0], "\n";
             return [$field[6]-1, $field[7]-1]; # to fix a bug in bl2seq. 35-44 is wrongly displayed as 36-45
        }

    }
    return [-1, -1];
}


sub usedTime {

     my $usedTime;
     $usedTime = time() - $oldTime;
     $oldTime = time();
     return $usedTime;
}


sub printStatistics {

     my @sequence = @_;
     my $totalLength = 0;
     my $index;
     for ($index = 0; $index <= $#sequence; $index++) {
           $totalLength += length($sequence[$index]{'seq'});
     }
     my $averageLength = int($totalLength / $index);
     print scalar(@sequence) . " sequence(s) with average length of $averageLength.\n";

}


sub exportPrimer {

     my $index = shift;

     # positions for all qualified primers for the current sequence
     my @position = @{$dna[$index]{'primerPos'}};
     # sequences for all qualified primers for the current sequence
     my @primerSeq = @{$dna[$index]{'primerSeq'}};
     my @rPrimerSeq = @{$dna[$index]{'rPrimerSeq'}};
     my @stemSeq = @{$dna[$index]{'stemSeq'}};
     my @rStemSeq = @{$dna[$index]{'rStemSeq'}};

     open (OUT, ">>$primerFile") || die("Can not open file for writing!\n");
     print OUT "Accession\tLength\tPrimers\tF3\tF2\tFL\tF1c\tB1c\tBL\tB2\tB3\t";
     print OUT "F3\tlen\tTm\tF2\tlen\tTm\tFL\tlen\tTm\tF1c\tlen\tTm\tB1c\tlen\tTm\tBL\tlen\tTm\tB2\tlen\tTm\tB3\tlen\tTm\tFIP\tlen\tTm\tBIP\tlen\tTm\t\n";

     my ($dG, $foldingTm);
     foreach my $indx_ref (@{$dna[$index]{'primer_set'}}) {

             print OUT $dna[$index]{'acc'};
             print OUT "\t";
             print OUT length($dna[$index]{'seq'});
             print OUT "\t";
             print OUT scalar(@position);
             print OUT "\t";

             my ($F3_i, $F2_i, $FL_i,$F1_i, $B1_i, $BL_i, $B2_i, $B3_i, $primer_FIP, $primer_BIP) = @{$indx_ref};
             my $primer_F3 = $primerSeq[$F3_i];
             my $primer_F2 = $primerSeq[$F2_i];
             my $primer_FL = $rPrimerSeq[$FL_i];
             my $primer_F1c = $rStemSeq[$F1_i];
             my $primer_B1c = $stemSeq[$B1_i];
             my $primer_BL = $primerSeq[$BL_i];
             my $primer_B2 = $rPrimerSeq[$B2_i];
             my $primer_B3 = $rPrimerSeq[$B3_i];

             foreach my $pos ($F3_i, $F2_i, $FL_i,$F1_i, $B1_i, $BL_i, $B2_i, $B3_i) {
                     print OUT $position[$pos], "\t";
             }
             foreach my $seq ($primer_F3, $primer_F2, $primer_FL, $primer_F1c, $primer_B1c, $primer_BL, $primer_B2, $primer_B3) {
                     print OUT $seq, "\t", length($seq), "\t", tmNeighbor($seq), "\t";
             }
             print OUT $primer_FIP, "\t", length($primer_FIP), "\t", tmNeighbor($primer_F2), "\t";
             print OUT $primer_BIP, "\t", length($primer_BIP), "\t", tmNeighbor($primer_B2), "\t";
             print OUT "\n";
     }
     close(OUT);
}


sub fastaToTab {

     my ($fastaFile, $tabFile) = @_;
     my $id = "";
     my $dna= "";
     my $lastLine = "";

     open (IN, "$fastaFile") || die("Can not open $fastaFile file for reading in fastaToTab sub!\n");
     open (OUT, ">$tabFile") || die("Can not open $tabFile file for writing!\n");

     while (<IN>) {

          s/\s+$//;
          next if ($_ !~ /\S/);

          if ($_ =~ /^\>/) {
               $id = $_;
               $id =~ s/^\>//;
               if ($lastLine =~ /^\>/) {   # for some casual format with more than one def line per seq
                    $id .= $_;
               }
               else {
                    print OUT "\n" if ($dna ne ""); # not the first line of the file
                    print OUT "$id\t";
                    $id = "";
               }
          }
          else {
               $_ =~ s/\s//g;
               $dna = $_;
               print OUT $dna;
          }
          $lastLine = $_;
     }

     close(IN);
     close(OUT);
}

# fasta file import
sub importFasta {

    my ($fastaFile) = @_;

    # intermediate tmp file with PID to store formatted and un-duplicate sequence data
    my $tabFile = "$fastaFile $$.tab";

    # get a tab file format from the original fasta file
    fastaToTab($fastaFile, $tabFile);

    my @seq = importTabSeq($tabFile);

    unlink $tabFile if -e $tabFile;

    return @seq;
}

# will take two input param: tab file name and starting index position for @dna
# the second param is used to import a filter fasta file
# return the number of elements in the @dna array
sub importTabSeq {

     my ($tabFile) = @_;
     my @sequence = ();
     my $index = 0;

     open (IN, "$tabFile") || die("Cannot open $tabFile file for reading in importTab sub!\n");
     while (<IN>) {

          s/\s+$//;
          my ($id, $sequence) = split /\t/, $_;

          $sequence[$index]{'id'} = $id;
          # all subsequent sequence manipulations are on small letters a, t, c, and g
          $sequence =~ tr/A-Z/a-z/;
          $sequence[$index]{'seq'} = $sequence;

          $index++;

     }
     close(IN);

     return @sequence;
}

# get a list of 10-mer duplicate n.t. sequences, and push into a hash
# $keyUnit{'aaatcgagatcaa'} = "1, 45, 3,444, 3, 555"  # gene index1, position1, gene index2, position2......
sub screenKey {

     my ($keyLength, $seq) = @_;
     my $keySeq;
     my %nmer = (); # keyUnit used in the main sub
     for (my $index = 0; $index < scalar(@$seq); $index++) {
          for (my $position = 0; $position <= length($$seq[$index]{'seq'}) - $keyLength; $position++) {

               $keySeq = substr($$seq[$index]{'seq'}, $position, $keyLength);
               next if $keySeq =~ /[^atcg]/;

               if (exists $nmer{$keySeq}) {
                     $nmer{$keySeq} .= " $index";
               }
               else {
                     $nmer{$keySeq} = "$index";
               }
          }
     }

     # remove index redundancy
     foreach my $word(keys %nmer) {
             my @indx = split / /, $nmer{$word};
             my %ind_hash = ();
             foreach (@indx) {
                     $ind_hash{$_} = 1;
             }
             $nmer{$word} = join(" ", keys(%ind_hash));
     }
     return %nmer;

}

# calculate gc pecent range 0-1
sub gcPercent {

    my ($sequence) = @_;
    $sequence =~ tr/AUTCG/autcg/;
    my $length = length($sequence);
    my $gcCount= 0;

    for (my $j = 0; $j < $length; $j++) {
         if (substr($sequence, $j, 1) eq 'c' || substr($sequence, $j, 1) eq 'g') {
              $gcCount++;
         }
    }

    return sprintf("%.2f", $gcCount / $length);

}

# return the self-complementary strand of the input sequence
sub dnaComplement {

    my ($sequence) = @_;
    $sequence =~ tr/atcgATCG/tagcTAGC/;
    $sequence = reverse($sequence);
    return $sequence;
}






