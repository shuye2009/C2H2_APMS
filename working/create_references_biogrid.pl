#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Jan 27, 2012
#Last modified on: Feb 21, 2018

my $multiplier = 10; # size of positive references is multiplied by this number to get the size of negative references
# Only $highc_positive_reference is used for generating gold standard positives. 
# $all_positive_reference are excluded from gold standard negatives.
my $highc_positive_reference_file = "C:/RSYNC/ID_mapping/BIOGRID-MV-Physical-4.4.224.tab3_human.txt"; #Aug 2023
my $all_positive_reference_file = "C:/RSYNC/ID_mapping/BIOGRID-ORGANISM-Homo_sapiens-4.4.224.tab3_human.txt"; 


my $baitlist_file = $ARGV[0]; #"inter_Nov_10_2016_EM_bait.txt";
my $preylist_file = $ARGV[1]; #"inter_Nov_10_2016_EM_prey.txt";
my $dataset_name = $ARGV[2]; # the name of the data set to be analysed
my $bait_col = $ARGV[3]; # 0-based column number of bait name
my $prey_col = $ARGV[4]; # 0-based column number of prey name
my $id_col = $ARGV[5]; # if multiple ids are used like "Name|Uniprot|Ensembl", $id_col give the position of the id to be used after split the ids.
my $pos_small = $dataset_name."_positive_reference_PPI_small.tab";  #for spoke model PPI
my $pos_large = $dataset_name."_positive_reference_PPI_large.tab";	#for matrix model PPI
my $neg_small = $dataset_name."_negative_reference_PPI_small.tab";
my $neg_large = $dataset_name."_negative_reference_PPI_large.tab";


open (B,  "<$baitlist_file") or die "can not open file";
open (P,  "<$preylist_file") or die "can not open file";
open (HC, "<$highc_positive_reference_file") or die "can not open reference file";
open (ALL, "<$all_positive_reference_file") or die "can not open reference file";

# the following code is unneccessary, since contaminants are removed from prey list
if(0){
	# exclude any interactions involving contaminants from the positive reference, but keep them in negative reference
	open(C, "C:/RSYNC/Mass_pec_score/CM_TF_GFP/contaminants.txt") or die "can not open contaminants file";
	#open(C, "C:/RSYNC/Mass_pec_score/CM_data/Feb_22_2013_EM/top_nonspecific_binders_VAP.txt") or die "can not open contaminants file";
	
	my %contaminants = ();
	
	#<C>;
	while(<C>){
	    chomp($_);
	    my $temp = $_;
	    $temp =~ s/\n//;
	    $temp =~ s/\r//;
	
	    my @tary = split(/\t/,$temp);
	    my $ctm = $tary[0];
	    #my @c_ids = split(/\|/, $ctm);
	    #my $name = $c_ids[1];
	    #$contaminants{$name} = 1;
	    $contaminants{$ctm} = 1;
	}
}


open (PS, ">$pos_small");
open (PL, ">$pos_large");

my %all_pos = ();

my %baits = (); #gi only
my %preys = ();
my @bait_ary = (); #all ids
my @prey_ary = ();

while (<P>){
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    
    my $prey = $tary[$prey_col];
    my @p_ids = split(/\|/, $prey);
    my $name = $p_ids[$id_col];
    unless(exists $preys{$name}){
    	$preys{$name} = $name;
        push(@prey_ary, $name);
    }
    
}

while (<B>){
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    
    my $bait = $tary[$bait_col];
    my @b_ids = split(/\|/, $bait);
    my $name = $b_ids[$id_col];
    unless(exists $baits{$name}){
    	$baits{$name} = $name;
        push(@bait_ary, $name);
    }
    
}

my %PL = ();
my %PS = ();
my $num_hc_in_data = 0;
while (<HC>){
    
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    my $bait = $tary[0];
    my $prey = $tary[1];
    
    my @b_ids = split(/\|/, $bait);
    my $a = $bait;
    
    my @p_ids = split(/\|/, $prey);
    my $b = $prey;
    
    
    #only works if preys subsume baits
   if(exists $preys{$a} and exists $preys{$b}){
   	 my $p1 = $a;
      my $p2 = $b;
      
      if($p1 =~ m/^KRT.*/){
           $p1 = "KRT";
	    }
	    if($p2 =~ m/^KRT.*/){
	        $p2 = "KRT";
	    }
      #if(!exists $contaminants{$p1} and !exists $contaminants{$p2}){
		    if(exists $baits{$a} and exists $preys{$b}){
		    	unless(exists($PS{"$a\t$b"}) or exists($PS{"$b\t$a"})){
		    		print PS "$baits{$a}\t$preys{$b}\n";
		    		$PS{"$a\t$b"} = 1;
		    	}
		    	unless(exists($PL{"$a\t$b"}) or exists($PL{"$b\t$a"})){
		    	  print PL "$baits{$a}\t$preys{$b}\n";
		    	  $PL{"$a\t$b"} = 1;
		    	  $num_hc_in_data++;
		    	}
		    }elsif(exists $baits{$b} and exists $preys{$a}){
		      unless(exists($PS{"$a\t$b"}) or exists($PS{"$b\t$a"})){
                print PS "$preys{$a}\t$baits{$b}\n";
                $PS{"$a\t$b"} = 1;
              }
              unless(exists($PL{"$a\t$b"}) or exists($PL{"$b\t$a"})){
                 print PL "$preys{$a}\t$baits{$b}\n";
                 $PL{"$a\t$b"} = 1;
                 $num_hc_in_data++;
              }
		    	
		    }elsif(exists $baits{$b} and exists $baits{$a}){
		      unless(exists($PS{"$a\t$b"}) or exists($PS{"$b\t$a"})){
                 print PS "$baits{$a}\t$baits{$b}\n";
                 $PS{"$a\t$b"} = 1;
              }
              unless(exists($PL{"$a\t$b"}) or exists($PL{"$b\t$a"})){
                 print PL "$baits{$a}\t$baits{$b}\n";
                 $PL{"$a\t$b"} = 1;
                 $num_hc_in_data++;
              }
		        
		    }elsif(exists $preys{$b} and exists $preys{$a}){       
		        unless(exists($PL{"$a\t$b"}) or exists($PL{"$b\t$a"})){
                 print PL "$preys{$a}\t$preys{$b}\n";
                 $PL{"$a\t$b"} = 1;
                 $num_hc_in_data++;
              }
		    }
      #}
   }   
}

#<ALL>;

my $size_all_pos = 0;
my $size_bait_prey_pos = 0;
while (<ALL>){
    
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    $temp =~ s/^\d*\t\d*\t//; #remove first 2 columns, for GS only
      
      my @tary = split(/\t/,$temp);
    
    my $a = $tary[0];
    my $b = $tary[1];
    
    if(exists $baits{$a} and exists $preys{$b}){
        $size_bait_prey_pos++;
    }elsif(exists $baits{$b} and exists $preys{$a}){
        $size_bait_prey_pos++;
    }
    
    if((exists $baits{$a} or exists $preys{$a}) and (exists $baits{$b} or exists $preys{$b})){
        $all_pos{"$a\t$b"} = 1;
        $all_pos{"$b\t$a"} = 1;
        
        $size_all_pos++;
    }
	    
}

close B;
close P;
close HC;
close PS;
close PL;

### start negative reference



open (PS, "<$pos_small") or die "can not open reference file";
open (PL, "<$pos_large") or die "can not open reference file";


open (NS, ">$neg_small");
open (NL, ">$neg_large");


my %positive_reference_PPI_small = ();
my $size_positive_ref_small = 0;


while (<PS>){
    
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    my $bait = $tary[0];
    my $prey = $tary[1];
    
    
    $positive_reference_PPI_small{"$bait\t$prey"} = 1;
    $positive_reference_PPI_small{"$prey\t$bait"} = 1;
    
    $size_positive_ref_small++;
}

my %positive_reference_PPI_large = ();
my $size_positive_ref_large = 0;


while (<PL>){
    
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    my $bait = $tary[0];
    my $prey = $tary[1];
    
    
    $positive_reference_PPI_large{"$bait\t$prey"} = 1;
    $positive_reference_PPI_large{"$prey\t$bait"} = 1;
    
    $size_positive_ref_large++;
}


my $n_bait = scalar(@bait_ary);
my $n_prey = scalar(@prey_ary);

## large reference
my $all_pairs_large = $n_prey*($n_prey-1)/2;
$size_all_pos *= 1.5; # assuming that 1.5 times true interactions to be identified as the amount of know interactions
#$multiplier = ($all_pairs_large-$size_all_pos)/$size_all_pos; ### this is not used 
my $size_negative_large = $multiplier*$size_positive_ref_large; # times of positive PPI
my $fraction_large = $size_negative_large/($all_pairs_large-$size_all_pos); # the fraction of negative PPI in all prey-prey pairs
my $neg_count_large = 0;
for (my $i=0; $i<$n_prey-1; $i++){
    my $p1 = $prey_ary[$i];
    
   
    for (my $j=$i+1; $j<$n_prey; $j++){
        my $p2 = $prey_ary[$j];
        
        unless(exists $all_pos{"$p1\t$p2"}){
            my $r = rand();
            if($r < $fraction_large and $neg_count_large < $size_negative_large){
                $neg_count_large++;
                print NL "$p1\t$p2\n";
            }
        }
    }
}

# small reference

my $all_pairs_small = $n_prey*$n_bait;

my $size_negative_small = $multiplier*$size_positive_ref_small; #  times of positive PPI

my $fraction_small = $size_negative_small/($all_pairs_small-$size_bait_prey_pos); # the fraction of negative PPI in all prey-prey pairs
my $neg_count_small = 0;
for (my $i=0; $i<$n_bait; $i++){
    my $p1 = $bait_ary[$i];
    
    for (my $j=0; $j<$n_prey; $j++){
        my $p2 = $prey_ary[$j];
        
        unless(exists $all_pos{"$p1\t$p2"}){
            my $r = rand();
            if($r < $fraction_small  and $neg_count_small < $size_negative_small){
                $neg_count_small++;
                print NS "$p1\t$p2\n";
            }
        }
    }
}

close(PS);
close(PL);
close(NS);
close(NL);

print "multiplier\t$multiplier\n";
print "size hc positive in data\t$num_hc_in_data\n";
print "size all positive in data\t$size_all_pos\n";
print "size all possible pairs\t$all_pairs_large\n";
print "fraction small ref\t$fraction_small\n";
print "size positve samll\t$size_positive_ref_small\n";
print "size negative small\t$neg_count_small\n";
print "fraction large ref\t$fraction_large\n";
print "size positive large\t$size_positive_ref_large\n";
print "size negative large\t$neg_count_large\nfinish\n";