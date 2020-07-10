#! /usr/bin/perl

use warnings;
use strict;
use List::MoreUtils qw(uniq);

# a script to extract a gene list from big files
# modified to put 0s when no data found

unless (scalar @ARGV == 1) {
	die "\n clean_GEO_gmt.pl <Harmonizome GMT file> .\n";

}


# database file
my $database = shift @ARGV;

unless (open(DATA, $database)) {
	die "cannot open $database file. \n";
}

# read in list of genes in array
my %database = ();

while (my $line = <DATA>) {
	
	$line =~ s/\n//g; # remove newline so that gene name is recognized
	$line =~ s/\r//g; # remove newline so that gene name is recognized

	my @tmp = get_line_data ($line);
	
	my $setname = shift @tmp; # removes 1st element, will be used to parse TF name
	my $setinfo = shift @tmp; # removes 2nd element

	my @tfparse = split(/_/,$setname);
	#print uc($tfparse[0])."\n";
	
	my $TF = uc($tfparse[0]);
	
	push (@{$database{$TF}}, @tmp);
}

close DATA;

#print keys %database;
open(GMT, '>', "2020-05-20_GEO_Harmonizome_Target_summary.gmt");

foreach my $tf (keys %database) {

	print GMT $tf."\tGEO_TF_targets\t".join("\t", uniq @{$database{$tf}} )."\n";

}

close GMT;

exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}