use strict;
use warnings;

# Go through VCF file and change only PL (genotype likelihood) fields by replacing "." with "255"

my $replaceWith = 255;	# This could be changed to another integer less than 255 

while (<>) {
	if (/^#/) { print and next }
	chomp;
	my @bits = split "\t";
	
	# Find out which field is PL field (separation by ":")
	my $format = $bits[8];
	my @formatFields = split ":", $format;
	my $plPos = "NA";
	for (my $ii = 0; $ii < scalar @formatFields; ++$ii) {
		$formatFields[$ii] eq "PL" and $plPos = $ii;
		}
	$plPos eq "NA" and die "Format does not contain PL field: $format\n";
	
	# Get PL data for each sample
	for (my $ii = 9; $ii < scalar @bits; ++$ii) {
		$bits[$ii] eq "." and next;	# Ignore those samples that have insufficient coverage and are therefore entirely "."
		my @sampleFields = split ":", $bits[$ii];
		my $samplePl = $sampleFields[$plPos];
		$samplePl =~ s/\./$replaceWith/g;
		$sampleFields[$plPos] = $samplePl;
		$bits[$ii] = join ":", @sampleFields;
		}
	
	print join "\t", @bits;
	print "\n";
	
	}
	
