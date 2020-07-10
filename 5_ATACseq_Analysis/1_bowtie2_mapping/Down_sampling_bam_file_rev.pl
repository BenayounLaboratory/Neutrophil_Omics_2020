#! /usr/bin/perl

use warnings;
use strict;
use IO::CaptureOutput qw(capture_exec);

# a script to downsample to desired size
unless (scalar @ARGV >= 2) {
	die "\nNot enough command line arguments.\n".
	"Usage : Down_sampling_bam_file.pl <target read number> <bam file list>.\n";

}

my $target = shift @ARGV; #  mapped reads

# create an array that contains the list of files to be treated
my @bam_files = @ARGV;


foreach my $elem (@bam_files) {

	# matched pattern containing the mark
	# $1 is set to the pattern in the parentheses
	$elem =~ m/(.+)\.bam/;
	my $filename = $1;
	print "Now processing $filename...\n";
	
	#sort bam file before call to samtools
	my $new_filename = "$filename"."_DS_$target.bam";
	
	my $MetricsFile = "$filename"."_metric_file.txt";
	
	# get number of mapped reads in bam file
	my $stdout = capture_exec("samtools view -c $elem");
    chomp $stdout;
    $stdout =~ s/\s//g;
    my $mapped =  $stdout;
	print "Mapped reads in $elem: $mapped\n";
	
	# calculate downsize factor for picard tool P option
	my $ds_proba = $target/$mapped;

	print "Downsizing randomly $elem by $ds_proba\n";
	
	srand(time|$$);
	my $seed = int(rand(1000));
	
	my $shell_line = "java -Xmx16g -jar /Users/benayoun/Softwares/picard-2.20.3/picard.jar DownsampleSam I=$elem O=$new_filename STRATEGY=HighAccuracy P=$ds_proba RANDOM_SEED=$seed M=$MetricsFile";
	
	print "sampling bam file $elem ...\n";
	system($shell_line);
	
}




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
