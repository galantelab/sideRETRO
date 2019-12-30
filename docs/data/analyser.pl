#!/usr/bin/env perl
################################################################################
#
#         FILE: analyser.pl
#
#        USAGE: perl analyser.pl
#
#  DESCRIPTION: Filter relevant insertions from a VCF file.
#
################################################################################

use strict;
use warnings;

# Some variables.
my $count = 0;
my $count_val = 0;
my $count_ins = 0;
my $count_homo_zyg = 0;
my $line_sz = 0;
my $shr = 0;

# Test for input arguments from command line.
if ( @ARGV ) {
	my $vcf = $ARGV[0];

	# Test input existence.
	die "ERROR: Non existent file \"$vcf\"" unless -f $vcf;

	open my $fh, '<', $vcf or die "Can't open file $vcf";
	while ( <$fh> ) {
		$count++;

		# Skip header lines.
		next if /^#/;
		chomp;
		my @line = split;
		$count_val++;

		# Skip Unknown chromosomes.
		next if (length($line[0]) < 4)
			or (length($line[0]) > 5)
			or ($line[6] ne "PASS");

		$count_ins++;
		$line_sz = @line - 9;
		my $ct = 0;

		# Test haplotype.
		foreach ( 9 .. $#line ) {
			$ct++ if (split(":", $line[$_]))[0] =~ /1/;
			$count_homo_zyg++ if $line[$_] =~ /^1\/1:/;
		}
		$shr++ if $ct > 1;
	}
	close $fh;

	# Print results.
	my $str = <<"	EOS";
	GENERAL STATISTICS: $vcf
	================================================================================
	Total number of lines in file:                  $count
	Total non header lines:                         $count_val
	Total valid insertions:                         $count_ins
	Number of shared insertions:                    $shr
	Ratio of homozygous insertions in all samples:  $count_homo_zyg / ($count_ins*$line_sz)
	================================================================================
	EOS

	print $str;
}
else {
	# Print an usage for help.
	die "ERROR: No input arguments";
}
