#include "config.h"

#include <string.h>
#include <assert.h>
#include "utils.h"
#include "chr.h"

#define CHR_BUFSIZ 32

ChrStd *
chr_std_new (void)
{
	ChrStd *cs = hash_new (NULL, NULL);

	/*
	* Standardize human chromosome
	* acording to GENCODE
	*/

	// Autosomal chromosomes with no "chr"
	hash_insert (cs, "1", "chr1");
	hash_insert (cs, "2", "chr2");
	hash_insert (cs, "3", "chr3");
	hash_insert (cs, "4", "chr4");
	hash_insert (cs, "5", "chr5");
	hash_insert (cs, "6", "chr6");
	hash_insert (cs, "7", "chr7");
	hash_insert (cs, "8", "chr8");
	hash_insert (cs, "9", "chr9");
	hash_insert (cs, "10", "chr10");
	hash_insert (cs, "11", "chr11");
	hash_insert (cs, "12", "chr12");
	hash_insert (cs, "13", "chr13");
	hash_insert (cs, "14", "chr14");
	hash_insert (cs, "15", "chr15");
	hash_insert (cs, "16", "chr16");
	hash_insert (cs, "17", "chr17");
	hash_insert (cs, "18", "chr18");
	hash_insert (cs, "19", "chr19");
	hash_insert (cs, "20", "chr20");
	hash_insert (cs, "21", "chr21");
	hash_insert (cs, "22", "chr22");

	// Autosomal chromosomes with "chr"
	hash_insert (cs, "chr1", "chr1");
	hash_insert (cs, "chr2", "chr2");
	hash_insert (cs, "chr3", "chr3");
	hash_insert (cs, "chr4", "chr4");
	hash_insert (cs, "chr5", "chr5");
	hash_insert (cs, "chr6", "chr6");
	hash_insert (cs, "chr7", "chr7");
	hash_insert (cs, "chr8", "chr8");
	hash_insert (cs, "chr9", "chr9");
	hash_insert (cs, "chr10", "chr10");
	hash_insert (cs, "chr11", "chr11");
	hash_insert (cs, "chr12", "chr12");
	hash_insert (cs, "chr13", "chr13");
	hash_insert (cs, "chr14", "chr14");
	hash_insert (cs, "chr15", "chr15");
	hash_insert (cs, "chr16", "chr16");
	hash_insert (cs, "chr17", "chr17");
	hash_insert (cs, "chr18", "chr18");
	hash_insert (cs, "chr19", "chr19");
	hash_insert (cs, "chr20", "chr20");
	hash_insert (cs, "chr21", "chr21");
	hash_insert (cs, "chr22", "chr22");

	// Sexual chromosomes with no "chr"
	hash_insert (cs, "y", "chrY");
	hash_insert (cs, "x", "chrX");
	hash_insert (cs, "m", "chrM");
	hash_insert (cs, "mt", "chrM");

	// Sexual chromosomes with "chr"
	hash_insert (cs, "chry", "chrY");
	hash_insert (cs, "chrx", "chrX");
	hash_insert (cs, "chrm", "chrM");
	hash_insert (cs, "chrmt", "chrM");

	return cs;
}

void
chr_std_free (ChrStd *cs)
{
	hash_free (cs);
}

const char *
chr_std_lookup (ChrStd *cs, const char *chr)
{
	assert (cs != NULL && chr != NULL);

	char chr_copy[CHR_BUFSIZ];
	const char *chr_std = NULL;

	strncpy (chr_copy, chr, CHR_BUFSIZ - 1);
	chr_copy[CHR_BUFSIZ - 1] = '\0';

	chr_std = hash_lookup (cs, lower (chr_copy));

	if (chr_std == NULL)
		chr_std = chr;

	return chr_std;
}
