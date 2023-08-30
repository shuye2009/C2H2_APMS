#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Dec 06, 2005
#Last modified on: Feb 21, 2018


use strict;
use warnings;
use Cwd;
#use Chart::Gnuplot;
my $dir = getcwd;

print "$dir\n";



my $score_files = $ARGV[0]; #comma seperated strings
my $score_columns = $ARGV[1]; #comma seperated integers
my $sizes = $ARGV[2]; #comma seperated strings
my $gi_col = $ARGV[3]; #single integer
my $baitfile = $ARGV[4]; #single string
my $preyfile = $ARGV[5]; #single string
my $dataset = $ARGV[6]; #single string

#`chdir $wd`;
#check if all input files exist
my @sizes = split(/,/, $sizes);
my @score_files = split(/,/, $score_files);
my @score_columns = split(/,/, $score_columns);
foreach(@score_files){
	print "$_\n";
	die "The score file, $_ does not exist\n" unless -e $_;
}

foreach(@score_columns){
	print "$_\n";
	die "The score column, $_ does not exist\n" unless defined $_;
}

my @bp_files = ($baitfile, $preyfile);
foreach my $f (@bp_files){
	print "$f\n";
	die "The file, $f does not exist\n" unless -e $f;
}


my %baits = ();
my %preys = ();

open(B, "<$baitfile") or die "Can not open bait file!";
open(P, "<$preyfile") or die "Can not open prey file!";

#die "debugging stop\n";
my $number_baits = 0;
while(<B>){         
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;
    my @line = split(/\t/, $temp);
    my $bait = $line[0];
    my @bait_ids = split(/\|/, $bait);
			
	my $bait_gi = $bait_ids[$gi_col];
	
    $baits{$bait_gi} = 1;
    $number_baits++;
}
while(<P>){
    my $temp = $_;
           
    $temp =~ s/\n//;
    $temp =~ s/\r//;
    my @line = split(/\t/, $temp);
        
    my $prey = $line[0];

	my @prey_ids = split(/\|/, $prey);
	my $prey_gi = $prey_ids[$gi_col];
	
    $preys{$prey_gi} = 1;
}
        
foreach my $size (@sizes){
   my $positive_reference_file = $dataset."_positive_reference_PPI_".$size.".tab";
   my $negative_reference_file = $dataset."_negative_reference_PPI_".$size.".tab";
   #my $reference_file = $wd.$dataset."_positive_reference_PPI_".$size.".tab"; #for FLY only 
   #my $negative_reference_file = $wd.$dataset."_negative_reference_PPI_".$size.".tab"; #for FLY only 
   open (R, "<$positive_reference_file") or die "can not open positive reference file";
   open (NR, "<$negative_reference_file") or die "can not open negative reference file";
   my %positive_reference_PPI = ();
   my %negative_reference_PPI = ();
   my $size_positive_ref = 0;
   my $size_negative_ref = 0;
   #my $header = <R>;
   
   while (<R>){
            
            chomp($_);
            my $temp = $_;
            $temp =~ s/\n//;
            $temp =~ s/\r//;
        
            my @tary = split(/\t/,$temp);
            my $bait = $tary[0];
            my $prey = $tary[1];
            
            
            my @bait_ids = split(/\|/, $bait);
			my @prey_ids = split(/\|/, $prey);
			my $bait_gi = $bait_ids[$gi_col];
			my $prey_gi = $prey_ids[$gi_col];  
        	$positive_reference_PPI{"$bait_gi\t$prey_gi"} = 1;
        	$positive_reference_PPI{"$prey_gi\t$bait_gi"} = 1;
            
           # print "$bait_gi\t$prey_gi\n";
            $size_positive_ref++;
    }
        
    while (<NR>){
            
        chomp($_);
        my $temp = $_;
        $temp =~ s/\n//;
        $temp =~ s/\r//;
        
        my @tary = split(/\t/,$temp);
        my $bait = $tary[0];
        my $prey = $tary[1];
            
        my @bait_ids = split(/\|/, $bait);
		my @prey_ids = split(/\|/, $prey);
		my $bait_gi = $bait_ids[$gi_col];
		my $prey_gi = $prey_ids[$gi_col];  
        $negative_reference_PPI{"$bait_gi\t$prey_gi"} = 1;
        $negative_reference_PPI{"$prey_gi\t$bait_gi"} = 1;
            #print "$bait_gi\t$prey_gi\n";
        $size_negative_ref++;
   }
    my $sets = scalar(@score_files);
	for (my $set=0; $set<$sets; $set++){
		
		my $infile = $score_files[$set];
		my $score_column = $score_columns[$set]; #
		
		
		my $prefix = $infile;
		$prefix =~ s/\.txt//g;
		$prefix =~ s/\.tab//g;
		$prefix =~ s/\.csv//g;
		$prefix =~ s/\.tsv//g;
		
	
		#my $outfile = $prefix."_ROC_".$size.".tab";
		my $outfile = $prefix."_GS_ROC_".$size.".tab";
		
		print "start\t$size\t$infile\n";
		print "$positive_reference_file\n$negative_reference_file\n";
		print "size positive $size_positive_ref\nsize negative $size_negative_ref\n";
		
		open(OUT, ">$outfile") or die "Can not open file!";
		
		open(IN, "<$infile") or die "Can not open $infile!";
		

		my %net = (); # interactions that match positive or negative references, not entire dataset
		my %net_all = (); #entire dataset
		
		
		my $max = -9999999999999;
		my $min = 999999999999;
		
		my %intersection = (); # intersection with positive reference set
		my $total_positive = 0;
		my $total_negative = 0;
		

		my $header = <IN>;
		while(<IN>){
		     my $temp = $_;
		   
		     $temp =~ s/\n//;
		     $temp =~ s/\r//;
		     my @line = split(/\t/, $temp);
		
		     my $orf1 = $line[0];
		     my $orf2 = $line[1];
		     my $score = $line[$score_column];
		     unless($orf1 && $orf2){next;}
		     if($orf1 eq $orf2){next;}
		     
		     
		    my @orf1_ids = split(/\|/, $orf1);
		    my @orf2_ids = split(/\|/, $orf2);
		    my $orf1_gi = $orf1_ids[$gi_col]; #$gi_col
		    my $orf2_gi = $orf2_ids[$gi_col];
		    
		  
		     
		  #  if($orf1_name eq "gfp"){next;}
		      
		      
			if($orf1_gi and $orf2_gi){
				
		    	$net_all{"$orf1_gi\t$orf2_gi"} = $score;
			     #if(exists $reference_PPI{"$orf1_gi\t$orf2_gi"}){ #for positive reference that only has gi as id
			     if(exists $positive_reference_PPI{"$orf1_gi\t$orf2_gi"}){
			         $intersection{"$orf1_gi\t$orf2_gi"} = 1;
			         $total_positive++;
			         $net{"$orf1_gi\t$orf2_gi"} = $score;
			         if($score > $max){$max = $score;}
			         if($score < $min){$min = $score;}
			     }elsif(exists $negative_reference_PPI{"$orf1_gi\t$orf2_gi"}){
			
			         $total_negative++;
			         $net{"$orf1_gi\t$orf2_gi"} = $score;
			         if($score > $max){$max = $score;}
			         if($score < $min){$min = $score;}
			     }
   
			     
			  }
		}
		close (IN);
		
		if($max > 100){
			$max = 100; ## just for compass score
		}
		
		my $n_bins = 100;
		my $bin_size = ($max-$min)/$n_bins;
		#my $bin_size = ($max-$min)/($n_bins*20); #for fine-grained cutoffs only
		
		#print OUT "$infile\n";
		print OUT "Cutoff\tFP\tTP\tsensitivity\tprecision\tbaits_seen\tpreys_seen\tnetwork_size\tTPR\tFPR\tbait_recovery\tbait_weighted_precision\n";
		
		my @net_indices = ();
		my @net_all_indices = ();
		my $end_index_net = 0;
		my $end_index_all = 0;
		
		my @net_array_sorted = ();
		my @net_array_sorted_keys = sort{$net{$b}<=>$net{$a}}keys%net;
		foreach my $int (@net_array_sorted_keys){
			push (@net_array_sorted, $net{$int});
			#print T "$int\t$net{$int}\n";
		}
		my @net_all_array_sorted = ();
		my @net_all_array_sorted_keys = sort{$net_all{$b}<=>$net_all{$a}}keys%net_all;
		foreach my $int (@net_all_array_sorted_keys){
		    push (@net_all_array_sorted, $net_all{$int});
		    #print T "$int\t$net_all{$int}\n";
		}
		
		my $net_size = scalar(@net_array_sorted);
		my $net_all_size = scalar(@net_all_array_sorted);
		
		
		my $TTP = 0; # total test positive
		my $TTN = 0; # total test negative
		my $TP = 0; # true positive
		my $FP = 0; # false positive
		my $TN = 0; # true negative
		my $FN = 0; # false negative
		
		my %baits_seen = ();
		my %preys_seen = ();
		my %pair_seen = (); #pairs above a cutoff
		
		
		my $net_start = 0;
		my $net_end = 0;
		my $net_all_start = 0;
		my $net_all_end = 0;
		for(my $i=$max; $i>$min; $i-=$bin_size){
			if($net_array_sorted[$net_start] >= $i){
				 $net_start = $net_end;
			    
			    while ($net_end < $net_size-1 and $net_array_sorted[$net_end] >= $i){
			        $net_end++;
			    }
			    $net_end--; #last one does not satisfy the condition
			
			    $net_all_start = $net_all_end;
			   
			    while ($net_all_end < $net_all_size-1 and $net_all_array_sorted[$net_all_end] >= $i){
			        $net_all_end++;
			    }
			    $net_all_end--;
			   my $net_score = $net_array_sorted[$net_start];
			   my $net_end_score = $net_array_sorted[$net_end];
			   #print "$i\t\t$net_all_start\t$net_score\t\t$net_all_end\t$net_end_score\n";
			}
			
		   
		   for(my $k = $net_start; $k < $net_end; $k++){
		   	my $net_pair = $net_array_sorted_keys[$k];
		   	unless(defined $net_pair){die "i is $i, k is $k,  net_pair is not defined\n";}
		   	my $v = $net{$net_pair};
		   	$TTP++;
		      if (exists $intersection{$net_pair}){
		        $TP++;
		      }else{
		      	$FP++;
		      	if($FP <50){
		      	#	print "$net_pair\t$v\n";
		      	}
		      }
		   }
		   for(my $j = $net_all_start; $j <= $net_all_end; $j++){  
		      my $net_all_pair = $net_all_array_sorted_keys[$j];
		      unless(defined $net_all_pair){die "i is $i, j is $j,  net_all_pair is not defined\n";}
		      $pair_seen{$net_all_pair} = 1;
		      my @orfs = split(/\t/, $net_all_pair);
		      my $orf1 = $orfs[0];
		      my $orf2 = $orfs[1];
		          
		      if(exists $baits{$orf1}){
		        $baits_seen{$orf1} = 1;
		      }
		      if(exists $baits{$orf2}){
		        #$baits_seen{$orf2} = 1;
		      }
		          
		      if(exists $preys{$orf1}){
		        $preys_seen{$orf1} = 1;
		      }
		      if(exists $preys{$orf2}){
		         $preys_seen{$orf2} = 1;
		      }
		   }
		   
		   my @b_seen = keys(%baits_seen);
		   my @p_seen = keys(%preys_seen);
		   my @pr_seen = keys(%pair_seen);
		    
		    my $nb_seen = scalar(@b_seen);
		    my $np_seen = scalar(@p_seen);
		    my $seen_size = scalar(@pr_seen);
		
		     my $precision = 0;
		     if ($TTP > 0){
		          $precision = $TP/$TTP;
		     }
		     #my $sensitivity = $TP/$total_positive;
		     my $sensitivity = $TP/$size_positive_ref;
		     my $TPR = $TP/$size_positive_ref;
		     my $FPR = $FP/$size_negative_ref;
		     my $bait_recovery = $nb_seen/$number_baits;
		     my $bait_weighted_precision = sqrt($precision*$bait_recovery);
		     print OUT "$i\t$FP\t$TP\t$sensitivity\t$precision\t$nb_seen\t$np_seen\t$seen_size\t$TPR\t$FPR\t$bait_recovery\t$bait_weighted_precision\n";	     
			
			 while(my ($b, $v) = each %baits){
	            unless (exists $baits_seen{$b}){
	                if($i le 11){
	              #  	print "bait not included $b\n";
	                }
            	}
          	}
		}
		
		my $total_bait_not_seen = 0;
		while(my ($b, $v) = each %baits){
        unless (exists $baits_seen{$b}){
            $total_bait_not_seen++;
        }
     }
		
		close(OUT);
		print "min\t$min\nmax\t$max\ntotal baits not seen\t$total_bait_not_seen\nfinished $infile\n";
		print "total positive $total_positive\ntotal negative $total_negative\n";
	}
	close(R);
	close(NR);
}
close(B);
close(P);



print "DONE";
