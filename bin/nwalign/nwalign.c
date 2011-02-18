/*
 * Program:
 *	nwalign
 *
 * Author:
 *	David Arenillas (dave@cmmt.ubc.ca)
 *
 * Copyright:
 *	Wasserman Lab
 *	Centre for Molecular Medicine and Therapeutics 
 *	University of British Columbia
 *
 * Purpose:
 *	Find best global alignment between two nucleotide sequences using
 *	the Needleman-Wunsch algorithm as descibed by Gotoh [1982] and with
 *	the additional modification of a gap extension penalty.
 *
 * Synopsis:
 *	nwalign [sequence1 sequence2] | [infile]  match_score mismatch_score
 *		gap_penalty gap_ext_penalty outfile
 *
 *	where:
 *		sequence1,
 *		sequence2 	= the two nucleotide sequences to be aligned
 *		infile		= name of file containing two nucleotide
 *				  sequences to be aligned in either FastA
 *				  or text (space separated) format
 *		match_score	= score applied additively when corresponding
 *				  nucleotides within the two sequences match
 *		mismatch_score	= score applied additively when corresponding
 *				  nucleotides within the two sequences do not
 *				  match
 *		gap_penalty	= score applied subtractively when a gap must 
 *				  be opened to improve alignment between the two
 *				  sequences
 *		gap_ext_penalty	= score applied subtractively when an existing
 *				  gap must be extended to improve alignment
 *				  between the two sequences
 *		outfile		= name of output file to which alignment is
 *				  written
 *
 * NOTES:
 * 	The nucleotide comparison scores either a match or mismatch.
 * 	All non-ACTG characters (ie 'N') are scored as an automatic mismatch.
 */
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "nwalign.h"
#include "matrix.h"
#include "alignment.h"

static const int	DEBUG = 0;

/*
 * An abitrary maximum sequence length, to prevent reading of erroneously
 * long sequences. This is independent of the maximum matrix size which
 * limits the combined length of sequences which can actually be aligned.
 */
static const int	MaxSeqLen = 1000000;


static int		seqFileRead(const char *fname, char **seq1,
				char **seq2);
static int		blankLine(const char *line);
static void		chomp(char *str);
static void		usage();

int			lastErrCode = 0;

Alignment *		NWAlign(const char *seq1, const char *seq2,
				int matchScore, int mismatchScore,
				int gapPenalty, int gapExtPenalty);


int
main(int argc, char **argv)
{
	char		*seq1 = NULL;
	char		*seq2 = NULL;
	char		*ifname, *ofname;
	int		matchScore, mismatchScore, gapPenalty, gapExtPenalty;
	FILE		*ofp;
	Alignment	*align = NULL;

	/*
	 * Simplistic argument reading/checking. More sophisticated argument
	 * checking using getopt is sacrificed in order to minimize program
	 * overhead in order to maximize execution speed.
	 */
	if (argc == 7) {
		/*
		 * sequences contained in a file
		 */
		ifname		= strdup(argv[1]);
		matchScore	= atoi(argv[2]);
		mismatchScore	= atoi(argv[3]);
		gapPenalty	= atoi(argv[4]);
		gapExtPenalty	= atoi(argv[5]);
		ofname		= strdup(argv[6]);
		/*
		 * Read sequences from file
		 */
		if (!seqFileRead(ifname, &seq1, &seq2)) {
			fputs("nwalign: error reading input sequences\n",
				stderr);
			exit(lastErrCode);
		}
	} else if (argc == 8) {
		/*
		 * sequences passed on the command line
		 */
		seq1 = strdup(argv[1]);
		seq2 = strdup(argv[2]);
		matchScore	= atoi(argv[3]);
		mismatchScore	= atoi(argv[4]);
		gapPenalty	= atoi(argv[5]);
		gapExtPenalty	= atoi(argv[6]);
		ofname		= strdup(argv[7]);
	} else {
		usage();
		exit(E_BADARGS);
	}

	if (!seq1) {
		fprintf(stderr, "nwalign: error obtaining sequence 1\n");
		exit(E_OTHER);
	}

	if (!seq2) {
		fputs("nwalign: error obtaining sequence 2\n", stderr);
		free(seq1);
		exit(E_OTHER);
	}

	if (DEBUG) {
		printf("\nseq1 = %s\n", seq1);
		printf("\nseq2 = %s\n", seq2);
	}

	/*
	 * perform Needleman-Wunsch alignment
	 */
	align = NWAlign(seq1, seq2, matchScore, mismatchScore, gapPenalty,
			gapExtPenalty);
	free(seq1);
	free(seq2);

	/*
	 * assume alignment is meaningless if less than some theoretical
	 * mimimum score
	 */
	if (align) {
		if (DEBUG)
			printf("\nAlignment score = %d\n", align->score);

		/*
		 * XXX
		 * Should alignment be considered a failure if it does not
		 * score above some minimum threshold value?
		 */
		/*
		if (align->score >= MinAlignmentScore) {
			AlignmentPrint(align, ofp);
		} else {
			fputs("\nNo significant alignment found\n", stderr);
			lastErrCode = E_ALIGNQUALITY;
		}
		*/

		if (!(ofp = fopen(ofname, "w"))) {
			perror("nwalign: error opening output file\n");
			exit(E_FILEIO);
		}
		AlignmentPrint(align, ofp);
		fclose(ofp);
		AlignmentDestroy(align);
	} else {
		fputs("nwalign: alignment failed\n", stderr);
		if (!lastErrCode)
		    lastErrCode = E_ALIGNQUALITY;
	}

	return lastErrCode;
}

/*
 * Performs Needleman-Wunsch alignment on the sequences seq1 and seq2. Based
 * on the provided matchScore, mismatchScore, gapPenalty and gapExtPenalty.
 * Returns the best alignment.
 */
Alignment *
NWAlign(const char *seq1, const char *seq2, int matchScore, int mismatchScore,
	int gapPenalty, int gapExtPenalty)
{
	int		seq1Len, seq2Len;
	Matrix		*scores = NULL;
	Alignment	*align = NULL;

	/*
	printf("\nnwalign called with:\n");
	printf("\tsequence1      = %s\n", seq1);
	printf("\tsequence2      = %s\n", seq2);
	printf("\tmatch score    = %d\n", matchScore);
	printf("\tmismatch score = %d\n", mismatchScore);
	printf("\tgap penalty    = %d\n\n", gapPenalty);
	*/

	if (!seq1 || !seq1[0]) {
	    	lastErrCode = E_BADARGS;
		return NULL;
	}

	if (!seq2 || !seq2[0]) {
	    	lastErrCode = E_BADARGS;
		return NULL;
	}

	seq1Len = strlen(seq1);
	seq2Len = strlen(seq2);

	/*
	 * Create the empty scoring matrix
	 */
	scores = MatrixCreate(seq1Len + 1, seq2Len + 1);
	if (!scores) {
	    	if (!lastErrCode)
		    	lastErrCode = E_OTHER;
		return NULL;
	}

	/*
	 * Fill in the scoring matrix according to the Needleman-Wunsch
	 * algorithm.
	 */
	if (!MatrixScore(scores, seq1, seq2, matchScore, mismatchScore,
			gapPenalty, gapExtPenalty))
	{
	    	if (!lastErrCode)
		    	lastErrCode = E_OTHER;
		return NULL;
	}

	/*
	 * Build alignment by tracing back the optimal path through the
	 * scoring matrix.
	 */
	align = MatrixBuildAlignment(scores, seq1, seq2);

	MatrixDestroy(scores);

	return align;
}

/*
 * Read two nucleotide sequences from a file.
 * The file should be either FastA format or have the two sequences separated
 * by at least a blank line.
 */
static int
seqFileRead(const char *fname, char **pseq1, char **pseq2)
{
	static const int	bufSize = 128;
	static const int	sizeInc = 32768;
	char			*seq1 = NULL;
	char			*seq2 = NULL;
	int			seqNum;
	int			seq1Len, seq2Len, bufLen;
	int			seq1Size, seq2Size;
	int			done;
	int			retVal = 1;
	int			lastLineBlank;
	char			buf[bufSize];
	FILE			*fp;

	seq1Len = seq2Len = seq1Size = seq2Size = 0;

	if (!fname) {
	    	lastErrCode = E_BADARGS;
		return 0;
	}

	if (!(fp = fopen(fname, "r"))) {
	        perror("seqFileRead: error opening sequence file\n");
	    	lastErrCode = E_FILEIO;
		return 0;
	}

	seq1 = (char *) malloc(sizeInc);
	if (!seq1) {
	        perror("seqFileRead: seq1 malloc error\n");
		lastErrCode = E_MEM;
		retVal = 0;
		goto cleanup;
	}
	seq1Size = sizeInc;

	seq2 = (char *) malloc(sizeInc);
	if (!seq2) {
	        perror("seqFileRead: seq2 malloc error\n");
		lastErrCode = E_MEM;
		retVal = 0;
		goto cleanup;
	}
	seq2Size = sizeInc;

	lastLineBlank = 1;
	seqNum = 0;
	done = 0;
	while (!done && fgets(buf, bufSize, fp)) {
		/*
		 * get rid of any newlines transferred by fgets
		 */
		chomp(buf);
		bufLen = strlen(buf);

		if (buf[0] == '>') {
			/*
			 * FastA format; indicates start of a new sequence
			 */
			seqNum++;
			lastLineBlank = 0;
			if (seqNum > 2)
				/*
				 * ignore extra sequences
				 */
				done = 1;
		} else if (blankLine(buf)) {
			lastLineBlank = 1;
		} else {
			/*
			 * non-blank line following a blank line indicates start
			 * of a new sequence (non-FastA format)
			 */
			if (lastLineBlank)
				seqNum++;
			lastLineBlank = 0;
			if (seqNum > 2) {
				/*
				 * ignore extra sequences
				 */
				done = 1;
			} else if (seqNum == 1) {
				if (seq1Len + bufLen > MaxSeqLen) {
					/*
					 * If sequence is too long indicate
					 * an error
					 */
					fprintf(stderr, "seqFileRead: sequence"
						" 1 exceeds maximum allowed"
						" length of %d\n", MaxSeqLen);
					lastErrCode = E_SEQLEN;
					done = 1;
					retVal = 0;
				}
				if (seq1Len + bufLen > seq1Size) {
					seq1 = realloc(seq1,
							seq1Size + sizeInc);
					if (!seq1) {
						perror("seqFileRead: seq1"
							" realloc error\n");
						lastErrCode = E_MEM;
						done = 1;
						retVal = 0;
					}
					seq1Size += sizeInc;
				}
				if (retVal) {
					/*
					 * add to sequence one
					 */
					strncpy(&seq1[seq1Len], buf, bufLen);
					seq1Len += bufLen;
				}
			} else if (seqNum == 2) {
				if (seq2Len + bufLen > MaxSeqLen) {
					/*
					 * If sequence is too long indicate
					 * an error
					 */
					fprintf(stderr, "seqFileRead: Sequence"
						" 2 exceeds maximum allowed"
						" length of %d\n", MaxSeqLen);
					lastErrCode = E_SEQLEN;
					done = 1;
					retVal = 0;
				}
				if (seq2Len + bufLen > seq2Size) {
					/*
					 * add to sequence two
					 */
					seq2 = realloc(seq2,
							seq2Size + sizeInc);
					if (!seq2) {
						perror("seqFileRead: seq2"
							" realloc error\n");
						lastErrCode = E_MEM;
						done = 1;
						retVal = 0;
					}
					seq2Size += sizeInc;
				}
				if (retVal) {
					strncpy(&seq2[seq2Len], buf, bufLen);
					seq2Len += bufLen;
				}
			}
		}
	}

	/*
	 * make sure we have read at least one character for each sequence
	 */
	if (retVal) {
		if (seq1Len < 1) {
			fprintf(stderr, "seqFileRead: sequence 1 short\n");
			lastErrCode = E_BADSEQ;
			retVal = 0;
			goto cleanup;
		}
		if (seq2Len < 1) {
			fprintf(stderr, "seqFileRead: sequence 2 short\n");
			lastErrCode = E_BADSEQ;
			retVal = 0;
			goto cleanup;
		}
	}
	*pseq1 = seq1;
	*pseq2 = seq2;

cleanup:
	fclose(fp);

	if (!retVal) {
		if (seq1)
			free(seq1);
		if (seq2)
			free(seq2);
	}

	return retVal;
}

/*
 * Same as the perl function of the same name.
 * If the last character of a string is a newline, replaces it with the NULL
 * character.
 */
static void
chomp(char *str)
{
	int	len;

	if (str) {
		len = strlen(str);
		if (len && str[len - 1] == '\n')
			str[len - 1] = '\0';
	}
}

/*
 * If line is 'blank' returns true (1), otherwise false (0).
 * Definition of blank line is simplistic. If line is NULL or begins with a
 * whitespace character (as determined by isblank()) it is considered blank.
 */
static int
blankLine(const char *line)
{
	if (!line)
		return 1;
	if (*line == '\0' || isspace((int) *line)) 
		return 1;

	return 0;
}

/*
 * Print a program usage message
 */
static void
usage()
{
	printf("Usage: nwalign [sequence1 sequence2] | [infile]  match_score"
		" mismatch_score gap_penalty gap_extension_penalty outfile\n");
}
