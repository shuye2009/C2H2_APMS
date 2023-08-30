#!/usr/bin/perl -w

# Author: Shuye Pu
# Created on: Aug 02, 2016
# Last modified on: Aug 27, 2018, changed 345 to 376, as 31 new baits are added
# Aug 28, 2018. Working directory: C:\RSYNC\AP_MS_C2H2\PPI_SAINT_2018\data180826

my $protein_length_file = $ARGV[0]; #"C:/RSYNC/AP_MS_C2H2/human_ids_length.txt";
my $spectral_count_file1 = $ARGV[1]; #"C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data/1_All_C2H2_ZNF_Purifications_raw_data.txt"; #"20190118_AP_MS_data_for_SAINT_analysis.csv";#"160726_C2H2_ZNF_Purifications_fixed_plus_180826_New_purification.txt";
my $spectral_count_file2 = $ARGV[2]; #"C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data/1_160726_TF_Negative_Control_Purifications.txt";
my $spectral_count_file3 = $ARGV[3]; #"C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data/2_GFP_Purifications_raw_data.txt";

my $analysis = $ARGV[4]; #"C2H2_ZFP";

my $remove_list = $ARGV[5]; #"C:/RSYNC/AP_MS_C2H2/PPI_SAINT_2018/data180826/160816ListofproteinstoremovfromSAINTnetwork.txt";
my $update_list = $ARGV[6]; #"C:/RSYNC/AP_MS_C2H2/human_name_update.txt";

open (L, "<$protein_length_file") or die "can not open length file";

open (S1,  "<$spectral_count_file1") or die "can not open spectral count file";
open (S2,  "<$spectral_count_file2") or die "can not open spectral count file";
open (S3,  "<$spectral_count_file3") or die "can not open spectral count file";

my $GFP = 1; #for GFP only, 
my $TF_neg = 1; #for TF_Negative_Control only
my $bait_file = "ZNF_baitfile_saint_".$analysis.".tab"; #default $GFP = 1, $TF_neg = 1
my $prey_file = "ZNF_preyfile_saint_".$analysis.".tab";
my $interaction_file = "ZNF_interactionfile_saint_".$analysis.".tab"; 
my $bait_coverage = "ZNF_bait_coverage_".$analysis.".tab";

if($GFP == 1 and $TF_neg == 0){
	$bait_file = "ZNF_baitfile_saint_".$analysis."_GFP.tab";
	$prey_file = "ZNF_preyfile_saint_".$analysis."_GFP.tab";
	$interaction_file = "ZNF_interactionfile_saint_".$analysis."_GFP.tab"; 
	$bait_coverage = "ZNF_bait_coverage_".$analysis."_GFP.tab";
}

if($GFP == 0 and $TF_neg == 1){
	$bait_file = "ZNF_baitfile_saint_".$analysis."_TF_neg.tab";
	$prey_file = "ZNF_preyfile_saint_".$analysis."_TF_neg.tab";
	$interaction_file = "ZNF_interactionfile_saint_".$analysis."_TF_neg.tab";
	$bait_coverage = "ZNF_bait_coverage_".$analysis."_TF_neg.tab";
}

if($GFP == 0 and $TF_neg == 0){
	$bait_file = "ZNF_interactionfile_".$analysis."_ihgscore0p5_bait_list.txt";  ### Hgscore file name is ZNF_interactionfile_".$analysis."_ihgscore0p5.txt
	$prey_file = "ZNF_interactionfile_".$analysis."_ihgscore0p5_prey_list.txt";  ### Will be used by Assessing_PPI_with_ROC_matrix_model.pl
	$interaction_file = "ZNF_interactionfile_".$analysis."_for_HGSCORE.tab";
	$bait_coverage = "ZNF_bait_coverage_".$analysis."_HGSCORE.tab";
}

open (I, ">$interaction_file");

open (B, ">$bait_file");
open (P, ">$prey_file");
open (BC, ">$bait_coverage");

open (R, "<$remove_list");
open (U, "<$update_list");


my %remove = ();
my %update = ();

<R>;
while (<R>){
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    my $gene = $tary[0];
    
    $remove{$gene} = 1;
}

<U>;
while (<U>){
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    
    my $old = $tary[0];
    my $new = $tary[1];
    
    $update{$old} = $new;
        
    
}
my %protein_length = (); #use uniprot as key
my %name_length = (); #use gene name as key

my %baits = (); #purification + name
my %bait_name = (); # name only
my %preys = ();
my %no_length = ();
<L>;
while (<L>){
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

    my @tary = split(/\t/,$temp);
    my $gene = $tary[0];
    my @ids = split(/\|/, $gene);
    my $gi = $ids[0];
    my $name = $ids[1];
    if(exists $update{$name}){
    	$name = $update{$name};
    	print "name update $name in length file\n";
    }
    
    
    my $uniprot = $ids[2];
    my $length = $tary[1];
  	if($uniprot){
  		$protein_length{$uniprot} = $length;
  	}
  	if($name){
  		$name_length{$name} = $length;
  		
  	}
    
}
#srand(time()^($$ + ($$ << 15)));
#my $r = rand();
if($GFP + $TF_neg == 0){
	print I "Purification\tBait\tPrey\tSpectral counts\n";
}
<S1>;
while (<S1>){
	 chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

   # $r = rand(); #for random sampling of 10% of the data
    #if($r > 0.1){
    #	next;
    #}
    my @tary = split(/\t/,$temp);
    my $batch = $tary[0];
    $batch =~ s/\s//g;
    my $bait = $tary[1];
    $bait =~ s/\s//g;
    $bait =~ s/\([N|C]\-EGFP\)//g;
    if($bait =~ m/^(.*)_[1|2]$/){
    	$bait = $1;
    }
    my $prey = $tary[2];
    $prey =~ s/\s//g;
    my $prey_uniprot = $tary[3];
    my $tsc = $tary[4];
    
    if(exists $remove{$prey}){
    	print "prey $prey removed\n";
    	next;
    }
    
    if(exists $update{$bait}){
    	$bait = $update{$bait};
    	print "name update $bait bait\n";
    }
    
    if(exists $update{$prey}){
    	$prey = $update{$prey};
    	print "name update $prey prey\n";
    }
    
    my $sequence_length = $protein_length{$prey_uniprot};
    my $name_length = $name_length{$prey};
    if($prey ne "rm" and $prey ne ""){
	    if($sequence_length){
	        $preys{"$prey"} = $sequence_length;
	    }elsif($name_length){
	    	$preys{"$prey"} = $name_length;
	    }else{
	        $preys{"$prey"} = 300; ## default length 300
	           $no_length{$prey} = 1;
	    }	
	    
    }
    #unless (exists $selected_baits{$bait}){next;}
    
    $baits{"$batch\t$bait"} = "T";   	
    $bait_name{$bait} = 1;
    if($bait and exists($preys{$prey})){
    	if(($bait ne $prey) or ($GFP + $TF_neg == 0)){
    		print I "$batch\t$bait\t$prey\t$tsc\n";
    	}
    	   
    }

     
     if($bait eq $prey){
     	print BC  "$batch\t$bait\t$prey\t$tsc\n";
     }
}

if($TF_neg){

	<S2>;
	while (<S2>){
		 chomp($_);
	    my $temp = $_;
	    $temp =~ s/\n//;
	    $temp =~ s/\r//;
	
	   # $r = rand(); #for random sampling of 10% of the data
	    #if($r > 0.1){
	    #	next;
	    #}
	    my @tary = split(/\t/,$temp);
	    my $batch = $tary[0];
	    $batch =~ s/\s//g;
	    my $bait = $tary[1];
	    $bait =~ s/\s//g;
	    $bait =~ s/\([N|C]\-EGFP\)//g;
	    
	    if($bait =~ m/^(.+?)_.+$/){
	    	$bait = $1;
	    }
	    my $prey = $tary[2];
	    $prey =~ s/\s//g;
	    my $prey_uniprot = $tary[3];
	    my $tsc = $tary[4];
	    
	    if(exists $remove{$prey}){
	    	print "prey $prey removed\n";
	    	next;
	    }
	    
	    if(exists $update{$bait}){
	    	$bait = $update{$bait};
	    	print "name update $bait\n";
	    }
	    
	    if(exists $update{$prey}){
	    	$prey = $update{$prey};
	    	print "name update $prey\n";
	    }
	    
	    my $sequence_length = $protein_length{$prey_uniprot};
	    my $name_length = $name_length{$prey};
	    if($prey ne "rm" and $prey ne ""){
		    if($sequence_length){
		        $preys{"$prey"} = $sequence_length;
		    }elsif($name_length){
		    	$preys{"$prey"} = $name_length;
		    }else{
		        $preys{"$prey"} = 300; ## default length 300
		           $no_length{$prey} = 1;
		    }	
		    
	    }
	    
	    
	    $baits{"$batch\t$bait"} = "C";   	
	    
	    if($bait and exists($preys{$prey})){
	    	if($bait ne $prey){
	    		print I "$batch\t$bait\t$prey\t$tsc\n";
	    	}
	    	   
	    }
	
	     
	     if($bait eq $prey){
	     	print BC  "$batch\t$bait\t$prey\t$tsc\n";
	     }
	}
}

if($GFP){
	<S3>;
	while (<S3>){
		 chomp($_);
	    my $temp = $_;
	    $temp =~ s/\n//;
	    $temp =~ s/\r//;
	
	   # $r = rand(); #for random sampling of 10% of the data
	    #if($r > 0.1){
	    #	next;
	    #}
	    my @tary = split(/\t/,$temp);
	    my $batch = $tary[0];
	    $batch =~ s/\s//g;
	    my $bait = $tary[1];
	    $bait =~ s/\s//g;
	    $bait =~ s/\([N|C]\-EGFP\)//g;
	    if($bait =~ m/^(.+?)_.+$/){
	    	$bait = $1;
	    }
	    my $prey = $tary[2];
	    $prey =~ s/\s//g;
	    my $prey_uniprot = $tary[3];
	    my $tsc = $tary[4];
	    
	    ####
	    if(exists $remove{$prey}){
	    	next;
	    	print "prey $prey removed\n";
	    }
	    if(exists $update{$bait}){
	    	$bait = $update{$bait};
	    	print "name update $bait\n";
	    }
	    if(exists $update{$prey}){
	    	$prey = $update{$prey};
	    	print "name update $prey\n";
	    }
	    ####
	    
	    my $sequence_length = $protein_length{$prey_uniprot};
	    my $name_length = $name_length{$prey};
	    if($prey ne "rm" and $prey ne ""){
		    if($sequence_length){
		        $preys{"$prey"} = $sequence_length;
		    }elsif($name_length){
		    	$preys{"$prey"} = $name_length;
		    }else{
		        $preys{"$prey"} = 300; ## default length 300
		           $no_length{$prey} = 1;
		    }	
		    
	    }
	    
	    
	    $baits{"$batch\t$bait"} = "C";   	
	    
	    if($bait and exists($preys{$prey})){
	    	if($bait ne $prey){
	    		print I "$batch\t$bait\t$prey\t$tsc\n";
	    	}
	    	   
	    }
	
	     
	     if($bait eq $prey){
	     	print BC  "$batch\t$bait\t$prey\t$tsc\n";
	     }
	}
}

if($GFP + $TF_neg > 0){
	while (my ($bait, $t) = each %baits){	
		print B "$bait\t$t\n";
	}
	while (my ($prey, $l) = each %preys){
    	print P "$prey\t$l\t$prey\n";
	}
	
}

if($GFP + $TF_neg == 0){
	while (my ($bait, $t) = each %bait_name){
		print B "$bait\n";
	}
	while (my ($prey, $l) = each %preys){
    	print P "$prey\n";
	}
	
}




my @nl = keys(%no_length);
my $nnl = scalar(@nl);
foreach my $n(@nl){
	print "$n has no length\n";
}

print "$nnl preys have no length\n\nFinish\n";


#### these outputs are used for the following commnands ####
# SAINTexpress-spc ZNF_interactionfile_saint_376.tab ZNF_preyfile_saint_376.tab ZNF_baitfile_saint_376.tab
# SAINTexpress-spc ZNF_interactionfile_saint_376_GFP.tab ZNF_preyfile_saint_376_GFP.tab ZNF_baitfile_saint_376_GFP.tab
# SAINTexpress-spc ZNF_interactionfile_saint_376_TF_neg.tab ZNF_preyfile_saint_376_TF_neg.tab ZNF_baitfile_saint_376_TF_neg.tab
#
# Aug 28, 2018. Working directory: C:\RSYNC\AP_MS_C2H2\PPI_SAINT_2018\data180826
# java -jar -Xmx6144m ..\..\hgscore2\iTapScore.jar -ppw 0.5 -s HG -l ..\..\human_name_length.txt ZNF_interactionfile_376_for_HGSCORE.tab ZNF_interactionfile_376_ihgscore0p5.txt

