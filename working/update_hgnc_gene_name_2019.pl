#!/usr/bin/perl -w

#Author: Shuye Pu
#Created on: Aug 16, 2016
#Last modified on: Sept 01, 2023

# hgnc_gene_map_2019.txt is produced by querying HGNC website on Feb 21, 2019

my $input_file = $ARGV[0]; #"C:/RSYNC/ID_mapping/hgnc_gene_map_2019.txt";
my $output_file = $ARGV[1]; #"C:/RSYNC/ID_mapping/hgnc_geneName_updates_2019.txt";

open(I, "<$input_file");
open(O, ">$output_file");

my %old_new = (); #name only
print O "Old\tNew\n";

my $colliding = 0;
my $unchanged = 0;

my $h = <I>;
while (<I>){
    chomp($_);
    my $temp = $_;
    $temp =~ s/\n//;
    $temp =~ s/\r//;

	$temp =~ s/\@//g;
	$temp =~ s/\#/_/g;
	
    my @tary = split(/\t/,$temp);
    my $size = scalar(@tary);
    
    my $new = $tary[0];
    #$new = uc($new);
    #$new =~ s/\s//g;
    if($size == 1){
    	print O "$new\t$new\n";
    	$unchanged++;
    	$old_new{$new} = $new;
    }else{
    	$old_new{$new} = $new;
    	print O "$new\t$new\n";
    	my @olds = split(/, /, $tary[1]);
    	foreach $old(@olds){
    		#$old = uc($old);
    		#$old =~ s/\s//g;
    		if(exists $old_new{$old}){
    			$colliding++;
    			my $existing = $old_new{$old};
    			print "$old\t$existing\t$new\n";
    		}
    		$old_new{$old} = $new;
    		print O "$old\t$new\n";
    	}
    	if(defined $tary[2]){
	    	my @synonyms = split(/, /, $tary[2]);
	    	
	    	foreach $synonym(@synonyms){
	    		$synonym = uc($synonym);
	    		$synonym =~ s/\s//g;
	    		if(exists $old_new{$synonym}){
	    			$colliding++;
	    			my $existing = $old_new{$synonym};
	    			print "$synonym\t$existing\t$new\n";
	    		}
	    		$old_new{$synonym} = $new;
	    		print O "$synonym\t$new\n";
	    	}
    	}
    }
    
}


print "number of colliding names $colliding\n";
print "number of unchanged names $unchanged\n";

