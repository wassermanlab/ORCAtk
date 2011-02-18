/*
 * Program:
 *	align_cons
 *
 * Date:
 * 	June 23rd, 2003
 *
 * Author:
 * 	David Arenillas
 *	Wasserman Lab
 *	Centre for Molecular Medicine and Therapeutics 
 *	University of British Columbia
 *	e-mail: dave@cmmt.ubc.ca
 *
 * Purpose:
 * 	Calculate the conservation profile/conserved regions of a pairwise
 * 	alignment.
 *
 * Synopsis:
 *	align_cons [-i in_file] [-w window_size] [-n increment] [-t threshold]
 *		[-s stringency] [-m method] [-g gff_file] [-r report_type]
 *		[-f format] [-o out_file] [-c] [-D]
 *
 *	where:
 *		in_file		= Name of file in multi-FastA format containing
 *				  two aligned sequences.
 *				  (default = stdin)
 *		window_size	= Size of the comparison window in nucleotides
 *				  on the query sequence (first sequence)
 *		increment	= Number of nucleotides by which to increment
 *				  the starting position of the window
 *		threshold	= Min. score of conserved regions to report.
 *		stringency	= Stringency used to dynamically compute the
 *				  threshold. If theshold is also specified and
 *				  dynamically computed threshold is less than
 *				  that specified, the specified threshold
 *				  overrides the computed one.
 *		gff_file	= GFF formatted file containing features to
 *				  filter out of the conserved regions.
 *		report_type	= Character indicating type of report to
 *				  output. A 'p' outputs position and score
 *				  suitable for making conservation plots.
 *				  A 'c' outputs the conserved region positions
 *				  and scores. Format of this conserved regions
 *				  output is:
 *				  seq_start seq_end length align_start
 *				  align_end score
 *				  (default 'p')
 *		format		= For conservation plot reports (-r p),
 *				  character indicating whether the positions
 *				  output are the start 's', end 'e' or center
 *				  'c' positions of the sliding window.
 *		out_file	= Optional output file name
 *				  (default = stdout)
 *		-c		= Indicates that second sequence was reverse
 *				  complemented during alignment process.
 *				  In such cases this flag is required so that
 *				  the conserved region coordinates are
 *				  correctly computed for the second sequence.
 *		-D		= Turn on debugging output.
 *		-h		= Print usage message.
 *		-?		= Print usage message.
 *
 * Algorithm:
 * 	Slide a window of size W across the alignment in increments of N.
 * 	For each window position compute conservation score as percentage ID.
 * 	Two different methods for computing percentage ID are employed,
 * 	standard and overall, controlled by the -m switch.
 *
 * 	Standard percentage ID for a given window of size W is computed as
 * 	number of nucleotides in second sequence which match the corresponding
 * 	W nucleotides in the first sequence, divided by the size of the base
 * 	sequence window in nucleotides, W.
 * 	  i.e., for a W of 16 nucleotides the sequences:
 * 	    CT-AGT---CGA--TGT--GC--ATG
 * 	    CT-AG--TGCT----AG-TGCT--T-
 * 	  would score 8/16 = 0.5
 *
 * 	Overall percentage ID for a given window of size W is computed as
 * 	number of nucleotides in both sequences which match, divided by the
 * 	size of the alignment window.
 * 	  i.e., for a W of 16 the sequences:
 * 	    CT-AGT---CGA--TGT--GC--ATG
 * 	    CT-AG--TGCT----AG-TGCT--T-
 * 	  would score 8/26 = 0.31
 *
 * 	NOTE: Window size is based on actual number of nucleotides in the base
 * 	sequence (sequence 1) That is, gaps in sequence 1 are ignored and not
 * 	considered part of the window.
 *
 * 	Optionally compute the conserved regions. A Conserved region of length
 * 	L and percentage identity score S where S is greater than or equal to
 * 	the given threshold, is computed by combining each overlapping window
 * 	within that region with score greater than or equal to the threshold.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include <tgmath.h>

typedef enum {
    TRIM_LEFT	= 0x01,
    TRIM_RIGHT	= 0x02,
    TRIM_BOTH	= TRIM_LEFT | TRIM_RIGHT
} WhichTrim;

typedef struct {
    size_t	size;		/* allocated size of alignment */
    int		length;		/* length of alignment */
    int		seq1Len;	/* (ungapped) length of sequence 1 */
    int		seq2Len;	/* (ungapped) length of sequence 2 */
    char	seq1ID[128];	/* sequence 1 identifier */
    char	seq2ID[128];	/* sequence 2 identifier */
    char	*seq1;		/* aligned sequence 1 */
    char	*seq2;		/* alignmed sequence 2 */
    int		*matches;	/* array containing 1 for each nucleotide
    				   position in seq1 which matches the aligned
				   nucleotide in seq2, 0 otherwize */
    int		*seq1ToAlnIdx;	/* seq1 to align coord (0-based) transform
    				   array */
    int		*seq2ToAlnIdx;	/* seq2 to align coord (0-based) transform
    				   array */
} Alignment;

typedef struct {
    int		alignStart;	/* Start of sim. win. within alignment */
    int		alignEnd;	/* End of sim. win. within alignment */
    int		queryStart;	/* Start of sim. win. on query (1st) sequence */
    int		queryEnd;	/* End of sim. win. on query (1st) sequence */
    float	score;		/* Similarity score for this window */
} SimilarityWindow;

typedef struct {
    int			winSize;
    int			winInc;
    int			nWins;
    SimilarityWindow	*wins;
} Similarities;

typedef struct {
    int		alignStart;	/* Start of conserved region within alignment */
    int		alignEnd;	/* End of conserved region within alignment */
    int		queryStart;	/* Start of conserved region on base (1st)
				 * sequence */
    int		queryEnd;	/* End of conserved region on base (1st)
				 * sequence */
    int		hitStart;	/* Start of conserved region on comparison (2nd)
				 * sequence */
    int		hitEnd;		/* End of conserved region on comparison (2nd)
				 * sequence */
    float	score;		/* Conservation score (% identity)*/
} ConservedRegion;

typedef struct {
    int			nRegs;
    int			maxRegs;
    ConservedRegion	*regs;
} Regions;

typedef struct {
    int		seqStart;	/* Start of feature on sequence */
    int		seqEnd;		/* End of feature on sequence */
} Feature;

typedef struct {
    int		nFeats;
    int		maxFeats;
    Feature	*feats;
} Features;

/*
 * Default argument values
 */
static const int	DFLT_WIN_SIZE	= 100;	 /* window size */
static const int	DFLT_WIN_INC	= 1;	 /* window slide incremement */
static const int	DFLT_MIN_REG_LEN = 20;	/* min conserved region
						    length to report */
static const int	REG_TAIL_LEN = 20;	/* size of region edges to check
						   for trimming */
static const float	DFLT_STRINGENCY	= 0.10; /* stringency for computing
						    dynamic threshold */
static const char	DFLT_REPORT_TYPE = 'c';	 /* report format
						    'c' = conserved regions
						    	  report */
static const char	DFLT_FORMAT	= 's';	 /* output format
						    's' = start position */
static const char	DFLT_PCT_ID_METHOD = 's'; /* method used for computing
						    percent identity */
static const char	GAP_CHAR	= '-';
static const char	MASK_CHAR	= 'N';

static int		DEBUG = 1;
static int		VERBOSE = 0;

inline static int	nucMatch(char nuc1, char nuc2);
static void		chomp(char *str);
static int		blankLine(const char *line);
static void		usage();

static Alignment *	readAlignment(const char *fname);
static void		buildAlignmentIndexes(Alignment *align);
static void		buildAlignmentMatchArray(Alignment *align);
static Features *	readFeatures(const char *fname);
static void		writeSimilarity(const Similarities *sim,
				const char *file, char format);
static void		writeSimilarityWithAligns(const Similarities *sim,
				const Alignment *align, const char *file);
static void		writeFeature(const Feature *feat, const char *msg);
static void		writeConservedRegion(const ConservedRegion *cr,
				const char * msg);
static void		writeConservedAlignment(const ConservedRegion *cr,
				const Alignment *aln);
static void		writeConservedRegionsReport(const Regions *reg,
				const Alignment *align, int winSize,
				int winInc, float stringency,
				float min_threshold, float dynamic_threshold,
				float threshold, int minRegLen,
				const char *file);
static void		writeConservedRegions(const Regions *reg);
static Similarities *	computeSimilarity(const Alignment *align, int winSize,
				int winInc, char pctIDMethod);
/*
static Similarities *	computeOverallSimilarity(const Alignment *align,
				int winSize, int winInc);
*/
/*
static int		alignToSeqPos(const char *seq, int alnPos,
				char match_type);
*/
static int		seqToAlignPos(const Alignment *align, int seqNum,
				int seqPos);
static float		computeDynamicThreshold(const Similarities *sim,
				const Features *feat, float stringency);
static int		similarityOverlapsFeatures(const SimilarityWindow *win,
				const Features *feat);
static int		scoreCompare(const void *score1, const void *score2);
static Regions *	computeConservedRegions(const Alignment *align,
				const Similarities *sim, float threshold,
				int winSize, const Features *feat,
				int complemented, char pctIDMethod, int truncate);
static Regions *	computeConservedRegionsOld(const Alignment *align,
				const Similarities *sim, float threshold,
				int winSize, const Features *feat,
				int complemented, char pctIDMethod);
static int		truncateRegions(Regions *reg, const Alignment *align,
				char pctIDMethod, float threshold);
static int		extendRegions(Regions *reg, const Alignment *align,
				char pctIDMethod, float threshold);
static int		truncateRegionComplex(ConservedRegion *reg,
				const Alignment *align, char pctIDMethod,
				float threshold);
static int		truncateRegionSimple(ConservedRegion *reg,
				const Alignment *align, char pctIDMethod,
				float threshold);
static int		extendRegion(ConservedRegion *reg,
				const Alignment *align, char pctIDMethod,
				float threshold);
static int		mergeRegions(Regions *reg, const Alignment * align,
				char pctIDMethod, int winSize, float threshold);
static void		computeConservedRegionsHitCoords(const Alignment *align,
				Regions *reg, int complemented);
static Regions *	filterRegions(Regions *reg, const Features *feat,
				Alignment *align, float threshold,
				char pctIDMethod);
static float		scoreRegion(ConservedRegion *reg,
				const Alignment *align, char pctIDMethod);
static float		percentageIdentity(const Alignment *align,
				int start, int end, char pctIDMethod);
static int		trimRegion(ConservedRegion *reg, const Alignment *align,
				WhichTrim whichTrim);
static Alignment *	AlignmentAlloc(size_t size);
static Alignment *	AlignmentIncreaseSize(Alignment *align, size_t size);
static void		AlignmentFree(Alignment *align);
static Similarities *	SimilaritiesAlloc(int nSimWins);
static void		SimilaritiesFree(Similarities *sim);
static Regions *	RegionsAlloc(int nRegs);
static Regions *	RegionsIncreaseSize(Regions *reg, int nRegs);
static void		RegionsFree(Regions *reg);
static Features *	FeaturesAlloc(int nFeats);
static Features *	FeaturesIncreaseSize(Features *feat, int nFeats);
static void		FeaturesFree(Features *feat);
static Regions *	RegionsDeleteRegionNum(Regions *reg, int regNum);
Regions *		RegionsSplitRegionNum(Regions *reg, int regNum,
				int start, int end, Alignment *align,
				float threshold, char pctIDMethod);
static void		maskFeatures(Alignment *align, Features *feat);

int
main(int argc, char **argv)
{
    extern char		*optarg;

    char		format = '\0'; 
    char		reportType = '\0';
    char		pctIDMethod = '\0';
    int			winSize = 0;
    int			winInc = 0;
    int			minRegLen = 0;
    int			complemented = 0;
    int                 truncate = 1;
    char		*ifname = NULL;
    char		*ofname = NULL;
    char		*gffname = NULL;
    float		threshold = 0;
    float		min_threshold = 0;
    float		stringency = 0;
    float		dynamic_threshold = 0;
    char		*endptr;
    Similarities	*sim = NULL;
    Regions		*reg = NULL;
    Features		*feat = NULL;
    Alignment		*align = NULL;

    /*
     * Get program arguments
     */
    int opt;
    while ((opt = getopt(argc, argv, "?cDWhf:g:i:l:m:n:o:r:s:t:w:")) != EOF) {
        switch (opt) {
        case '?':
        case 'h':
            usage();
            exit(0);
            break;
	case 'c':
	    complemented = 1;
	    break;
        case 'D':
            DEBUG = 1;
            break;
        case 'f':
            /*
             * Ouput is of the form:
             *     position  similarity
             * Where format determines position as:
             *     's' - start position (default)
             *     'e' - end position
             *     'c' - center position
             */
            format = optarg[0];
            break;
        case 'g':
	    gffname = strdup(optarg);
            break;
        case 'i':
            ifname = strdup(optarg);
            break;
        case 'l':
	    minRegLen = atoi(optarg);
            break;
        case 'm':
            pctIDMethod = optarg[0];
            break;
        case 'n':
            winInc = atoi(optarg);
            break;
        case 'o':
            ofname = strdup(optarg);
            break;
	case 'r':
	    reportType = optarg[0];
	    break;
	case 's':
	    stringency = (float) strtod(optarg, &endptr);
	    if (endptr == optarg) {
		fputs("\nInvalid stringency specified\n", stderr);
		usage();
		exit(-1);
	    }
	    if (stringency <= 0 || stringency > 1) {
		fputs("\nStringency should be > 0 and <= 1.0\n", stderr);
		usage();
		exit(-1);
	    }
	    break;
        case 't':
            min_threshold = (float) strtod(optarg, &endptr);
	    if (endptr == optarg) {
		fputs("\nInvalid minimum threshold specified\n", stderr);
		usage();
		exit(-1);
	    }
	    if (min_threshold <= 0 || min_threshold > 1) {
		fputs("\nThreshold should be > 0 and <= 1\n", stderr);
		usage();
		exit(-1);
	    }
            break;
        case 'w':
            winSize = atoi(optarg);
            break;
	case 'W':
	    truncate = 0;
	    break;
        default:
            break;
        }
    }

    if (winSize <= 0)
        winSize = DFLT_WIN_SIZE;

    if (winInc <= 0)
        winInc = DFLT_WIN_INC;

    if (!reportType)
	reportType = DFLT_REPORT_TYPE;

    if (!pctIDMethod)
	pctIDMethod = DFLT_PCT_ID_METHOD;

    if (pctIDMethod != 's' && pctIDMethod != 'o') {
	fprintf(stderr, "Unknown percent identity calculation method %c\n",
		pctIDMethod);
    }

    switch (reportType) {
    case 'p':
    case 'P':
	if (!format)
	    format = DFLT_FORMAT;
	break;
    case 'c':
    case 'C':
	/*
	 * If neither min_threshold or stringency is specified, computation is
	 * based on default stringency.
	 */
	if (!stringency && !min_threshold)
	    stringency = DFLT_STRINGENCY;

	if (!minRegLen)
	    minRegLen = DFLT_MIN_REG_LEN;

	if (gffname) {
	    feat = readFeatures(gffname);
	    if (!feat) {
		fprintf(stderr, "\nError reading features from %s\n", gffname);
		exit(-1);
	    }
	}
	break;
    default:
	fprintf(stderr, "Unknown report type %c\n", reportType);
	usage();
	exit(-1);
	break;
    }

    align = readAlignment(ifname);
    if (!align) {
        fprintf(stderr, "\nError reading alignment from %s\n", ifname);
        exit(-1);
    }

    /* build mapping of sequence to alignment positions */
    buildAlignmentIndexes(align);

    if (reportType == 'p' || reportType == 'P') {
	/* compute matching nucleotides in alignment */
	buildAlignmentMatchArray(align);
	sim = computeSimilarity(align, winSize, winInc, pctIDMethod);
	if (!sim) {
	    fprintf(stderr, "\nError computing similarity\n");
	    AlignmentFree(align);
	    FeaturesFree(feat);
	    exit(-1);
	}

	if (DEBUG && VERBOSE)
	    writeSimilarityWithAligns(sim, align, ofname);
	else
	    writeSimilarity(sim, ofname, format);
    } else if (reportType == 'c' || reportType == 'C') {
	/*
	 * Mask out features (exons) before computing similarity
	 * DJA 2006/02/28
	 */
	if (feat) {
	    maskFeatures(align, feat);
	}

	/* compute matching nucleotides in alignment */
	buildAlignmentMatchArray(align);

	sim = computeSimilarity(align, winSize, winInc, pctIDMethod);
	if (!sim) {
	    fprintf(stderr, "\nError computing similarity\n");
	    AlignmentFree(align);
	    FeaturesFree(feat);
	    exit(-1);
	}

	/*
	 * If stringency is specified, compute the threshold dynamically based
	 * on this stringency (specified or default).
	 * Note: if the threshold is also specified and the dynamically
	 * computed threshold falls below the specified threshold, use the
	 * specified threshold.
	 */
	if (stringency) {
	    dynamic_threshold = computeDynamicThreshold(sim, feat, stringency);

	    if (dynamic_threshold > min_threshold)
		threshold = dynamic_threshold;
	    else
		threshold = min_threshold;
	} else {
	    threshold = min_threshold;
	}


	if (!threshold) {
	    /*
	     * No longer consider a threshold of 0 to be an error. Just
	     * give a warning.
	     *
	    AlignmentFree(align);
	    SimilaritiesFree(sim);
	    FeaturesFree(feat);
	    */
	    fputs("\nWarning: Specified stringency resulted in a"
		    " dynamic threshold of 0.\nYou may want to specify a"
		    " lower stringency.\n\n", stderr);
	    /*exit(-1);*/
	}

	reg = computeConservedRegions(align, sim, threshold, winSize, feat,
					complemented, pctIDMethod, truncate);
	if (reg) {
	    /*
	     * Re-activated DJA 2005/06/07
	     */
	    if (feat) {
		reg = filterRegions(reg, feat, align, threshold, pctIDMethod);
	    }

	    computeConservedRegionsHitCoords(align, reg, complemented);
	}

	if (reg) {
	    if (DEBUG) {
		puts("\nFinal conserved regions:\n");
		writeConservedRegions(reg);
	    }

	    writeConservedRegionsReport(reg, align, winSize, winInc, stringency,
				    min_threshold, dynamic_threshold, threshold,
				    minRegLen, ofname);
	} else {
	    printf("\nNo conserved regions found\n");
	}
    }

    if (feat)
	FeaturesFree(feat);

    if (reg)
	RegionsFree(reg);

    if (align)
	AlignmentFree(align);

    if (sim)
	SimilaritiesFree(sim);

    return 0;
}

/*
 * Read sequence features from a GFF formatted file.
 */
static Features *
readFeatures(const char *fname)
{
    static const int	bufSize = 128;
    static const int	featSizeInc = 50;

    int		nFeats = 0;
    char	buf[bufSize];
    FILE	*fp = NULL;
    Features	*feat;

    if (!fname)
	return NULL;

    fp = fopen(fname, "r");

    if (!fp) {
	fprintf(stderr, "Error readFeatures: could not open GFF file %s\n",
		fname);
        return NULL;
    }

    feat = FeaturesAlloc(featSizeInc);

    while (fgets(buf, bufSize, fp)) {
	if (nFeats >= feat->maxFeats) {
	    feat = FeaturesIncreaseSize(feat, featSizeInc);
	}
	if (sscanf(buf, "%*s\t%*s\t%*s\t%d\t%d", &feat->feats[nFeats].seqStart,
		&feat->feats[nFeats].seqEnd) == 2)
	{
	    nFeats++;
	}

    }

    if (!nFeats)
	fputs("Warning readFeatures: no features read\n", stderr);

    feat->nFeats = nFeats;

    return feat;
}

/*
 * Read two aligned nucleotide sequences from a file and return the alignment
 * struct.
 * The file should be either FastA format or have the two sequences separated
 * by at least a blank line.
 */
static Alignment *
readAlignment(const char *fname)
{
    static const int	bufSize = 128;
    static const int	alignSizeInc = 200000;

    char	buf[bufSize];
    int		seqNum;
    int		seq1Len, seq2Len, bufLen;
    int		done;
    int		lastLineBlank;
    int		ok;
    int		idChars;
    char	*bufPtr;
    FILE	*fp;
    Alignment	*align = NULL;

    if (fname)
        fp = fopen(fname, "r");
    else
        fp = stdin;

    if (!fp)
        return NULL;

    align = AlignmentAlloc(alignSizeInc);

    seq1Len = seq2Len = 0;
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
	    if (seqNum == 1) {
		bufPtr = &buf[1];
		while (*bufPtr == ' ')
		    bufPtr++;
		idChars = 0;
		while (idChars < 127 && (isalnum(*bufPtr) || *bufPtr == '_')) {
		    align->seq1ID[idChars] = *bufPtr;
		    bufPtr++;
		    idChars++;
		}
		align->seq1ID[idChars] = '\0';
	    } else if (seqNum == 2) {
		bufPtr = &buf[1];
		while (*bufPtr == ' ')
		    bufPtr++;
		idChars = 0;
		while (idChars < 127 && (isalnum(*bufPtr) || *bufPtr == '_')) {
		    align->seq2ID[idChars] = *bufPtr;
		    bufPtr++;
		    idChars++;
		}
		align->seq2ID[idChars] = '\0';
	    } else {
                /*
                 * ignore extra sequences
                 */
                done = 1;
	    }
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
                 * Ignore extra sequences
                 */
                done = 1;
            } else if (seqNum == 1) {
                if (seq1Len + bufLen >= align->size)
                    /*
                     * If sequence is too long increase
                     * size.
                     */
		    align = AlignmentIncreaseSize(align, alignSizeInc);
                /*
                 * add to sequence one
                 */
                strncpy(&(align->seq1[seq1Len]), buf,
                        bufLen);
                seq1Len += bufLen;
                align->seq1[seq1Len] = '\0';
            } else if (seqNum == 2) {
                if (seq2Len + bufLen >= align->size)
                    /*
                     * If sequence is too long increase
                     * size.
                     */
		    align = AlignmentIncreaseSize(align, alignSizeInc);
                /*
                 * add to sequence two
                 */
                strncpy(&(align->seq2[seq2Len]), buf,
                        bufLen);
                seq2Len += bufLen;
                align->seq2[seq2Len] = '\0';
            }
        }
    }

    fclose(fp);

    ok = 0;
    if (align) {
        if (align->seq1 && align->seq2) {
            if (strlen(align->seq1) == strlen(align->seq2)) {
                    align->length = strlen(align->seq1);
                    ok = 1;
            }
        }
    }

    if (!ok) {
        if (align) {
            AlignmentFree(align);
            align = NULL;
        }
    }

    return align;
}

/*
 * Compute the sequence to alignment coordinate conversion arrays and the
 * ungapped sequence lengths.
 */
static void
buildAlignmentIndexes(Alignment *align)
{
    int		seq1Idx, seq2Idx, alnIdx;
    int		alnLen;

    if (!align)
	return;

    alnLen = align->length;
    if (!alnLen)
	return;

    align->seq1ToAlnIdx = (int *) calloc(alnLen, sizeof(int));
    memset(align->seq1ToAlnIdx, 0, alnLen * sizeof(int));

    align->seq2ToAlnIdx = (int *) calloc(alnLen, sizeof(int));
    memset(align->seq2ToAlnIdx, 0, alnLen * sizeof(int));

    seq1Idx = 0;
    seq2Idx = 0;
    for (alnIdx = 0; alnIdx < alnLen; alnIdx++) {
	if (align->seq1[alnIdx] != GAP_CHAR) {
	    align->seq1ToAlnIdx[seq1Idx] = alnIdx;
	    seq1Idx++;
	}

	if (align->seq2[alnIdx] != GAP_CHAR) {
	    align->seq2ToAlnIdx[seq2Idx] = alnIdx;
	    seq2Idx++;
	}
    }

    align->seq1Len = seq1Idx;
    align->seq2Len = seq2Idx;
}

/*
 * Compute the match array
 */
static void
buildAlignmentMatchArray(Alignment *align)
{
    int		seq1Len, alnLen, seq1Idx, alnIdx;

    if (!align)
	return;

    seq1Len = align->seq1Len;
    if (!seq1Len)
	return;

    alnLen = align->length;
    if (!alnLen)
	return;

    align->matches = (int *) calloc(seq1Len, sizeof(int));
    memset(align->matches, 0, seq1Len * sizeof(int));

    seq1Idx = 0;
    for (alnIdx = 0; alnIdx < alnLen; alnIdx++) {
	if (align->seq1[alnIdx] != GAP_CHAR) {
	    if (nucMatch(align->seq1[alnIdx], align->seq2[alnIdx])) {
	    	align->matches[seq1Idx] = 1;
	    }
	    seq1Idx++;
	}
    }
}

static void
maskFeatures(Alignment *align, Features *feat)
{
    int		i, j;
    int		featAlnStart, featAlnEnd;

    for (i = 0; i < feat->nFeats; i++) {
	featAlnStart = seqToAlignPos(align, 1, feat->feats[i].seqStart);
	featAlnEnd = seqToAlignPos(align, 1, feat->feats[i].seqEnd);
	if (DEBUG) {
	    printf("Masking feature %d - %d; %d - %d\n",
			feat->feats[i].seqStart, feat->feats[i].seqEnd,
	    		featAlnStart, featAlnEnd);
	}

	for (j = featAlnStart; j <= featAlnEnd; j++) {
	    if (align->seq1[j - 1] != GAP_CHAR)
		align->seq1[j - 1] = MASK_CHAR;
	}
    }
}

static Similarities *
computeSimilarity(const Alignment *align, int winSize, int winInc,
		    char pctIDMethod)
{
    int			winSeq1Start, winSeq1End;
    int			seq1Len, alignLen;
    int			winIdx;
    int			nSimWins;
    Similarities	*sim = NULL;

    if (!align)
        return NULL;

    if (winSize <= 0)
        return NULL;

    if (winInc <= 0)
        return NULL;

    seq1Len = align->seq1Len;
    if (seq1Len <= 0)
        return NULL;

    alignLen = align->length;
    if (alignLen <= 0)
        return NULL;

    if (winSize > seq1Len) {
	fprintf(stderr, "\nspecified window size (%d) is greater than"
			" length of sequence 1 (%d)\n", winSize, seq1Len);
	return NULL;
    }

    nSimWins = seq1Len - winSize + 1;
    sim = SimilaritiesAlloc(nSimWins);
    sim->winSize = winSize;
    sim->winInc = winInc;
    sim->nWins = 0;

    winIdx = 0;
    winSeq1Start = 1;
    winSeq1End = winSeq1Start + winSize - 1;
    while (winSeq1End <= seq1Len) {
	sim->wins[winIdx].queryStart = winSeq1Start;
	sim->wins[winIdx].queryEnd   = winSeq1End;
	sim->wins[winIdx].alignStart = seqToAlignPos(align, 1, winSeq1Start);
	sim->wins[winIdx].alignEnd   = seqToAlignPos(align, 1, winSeq1End);

	sim->wins[winIdx].score = percentageIdentity(align, winSeq1Start,
						    winSeq1End, pctIDMethod);

	sim->nWins++;

	winSeq1Start++;
	winSeq1End++;
	winIdx++;
    }

    return sim;
}

static float
computeDynamicThreshold(const Similarities *sim, const Features *feat,
			float stringency)
{
    int		i;
    int		nScores;
    int		thresholdIdx;
    float	*scores = NULL;
    float	threshold;

    scores = (float *) malloc(sim->nWins * sizeof(float));
    memset(scores, 0, sim->nWins * sizeof(float));

    nScores = 0;
    for (i = 0; i < sim->nWins; i++) {
    	/*
	 * We now mask features so this check is no longer needed.
	 * DJA 2006/02/28
	 *
	 * if (!similarityOverlapsFeatures(&(sim->wins[i]), feat)) {
	 */
	    scores[i] = sim->wins[i].score;
	    nScores++;
	/*
	 * } else {
	 *     scores[i] = 0;
	 * }
	 */
    }

    qsort(scores, sim->nWins, sizeof(float), scoreCompare);

    thresholdIdx = sim->nWins - (int) ceilf(nScores * stringency);
    if (thresholdIdx >= sim->nWins)
	threshold = 1.0;
    else
	threshold = scores[thresholdIdx];
    if (DEBUG) {
	fprintf(stderr, "Number of similarity windows %d\n", sim->nWins);
	fprintf(stderr, "Number of (non-exon) windows %d\n", nScores);
	fprintf(stderr, "Stringency %0.3f\n", stringency);
	fprintf(stderr, "Threshold index %d\n", thresholdIdx);
	fprintf(stderr, "Threshold %0.3f\n", threshold);
    }

    free(scores);

    return threshold;
}

static int
similarityOverlapsFeatures(const SimilarityWindow *win, const Features *feat)
{
    int		i;

    if (!win || !feat)
	return 0;

    for (i = 0; i < feat->nFeats; i++) {
	if (win->queryEnd >= feat->feats[i].seqStart
		&& win->queryStart <= feat->feats[i].seqEnd)
	    return 1;
    }

    return 0;
}

static int
scoreCompare(const void *score1, const void *score2)
{
    if (*((float *) score1) > *((float *) score2))
    	return 1;
    if (*((float *) score1) < *((float *) score2))
    	return -1;

    return 0;
}

static int
featureCompare(const void *feat1, const void *feat2)
{
    const Feature	*f1, *f2;

    f1 = (const Feature *) feat1;
    f2 = (const Feature *) feat2;

    if (f1->seqStart > f2->seqStart)
    	return 1;
    if (f1->seqStart < f2->seqStart)
    	return -1;

    return 0;
}

/*
 * Combine similarity windows to make regions of conservation
 * New version of computeConservedRegions. DJA 2006/03/02.
 */
static Regions *
computeConservedRegions(const Alignment *align, const Similarities *sim,
	float threshold, int winSize, const Features *feat, int complemented,
	char pctIDMethod, int truncate)
{
    static const int	regSizeInc = 100;

    int			nWins;
    int			curWinIdx;
    int			curRegIdx;
    int			merged;
    float		pctID;
    Regions		*reg = NULL;
    SimilarityWindow	*curWin = NULL;
    ConservedRegion	*curReg = NULL;

    if (!align)
        return NULL;

    if (!sim)
        return NULL;

    reg = RegionsAlloc(regSizeInc);

    nWins = sim->nWins;

    curWinIdx = 0;
    curRegIdx = -1;
    while (curWinIdx < nWins) {
	curWin = &sim->wins[curWinIdx];
	if (curWin->score >= threshold) {

	    /* start a new region */
	    curRegIdx++;
	    if (curRegIdx >= reg->maxRegs)
		reg = RegionsIncreaseSize(reg, regSizeInc);
	    curReg = &reg->regs[curRegIdx];
	    
	    curReg->queryStart	= curWin->queryStart;
	    curReg->queryEnd	= curWin->queryEnd;
	    curReg->alignStart	= curWin->alignStart;
	    curReg->alignEnd	= curWin->alignEnd;
	    curReg->score	= curWin->score;
	    reg->nRegs++;
	    
	    if (DEBUG)
		writeConservedRegion(curReg, "Initial conserved region");

	    /* extend region */
	    while(curWinIdx < nWins - 1 &&
		  ((pctIDMethod == 's' && sim->wins[curWinIdx+1].queryStart <= curReg->queryEnd + 1) ||
		   (pctIDMethod == 'o' && sim->wins[curWinIdx+1].alignStart <= curReg->alignEnd + 1))) {
		curWinIdx++;
		curWin = &sim->wins[curWinIdx];
		if (curWin->score >= threshold) {
		    /*
		      pctID = percentageIdentity(align, curReg->queryStart,
		      curWin->queryEnd, pctIDMethod);
		      if (pctID >= threshold) {
		    */
		    /* add this window to the current region */
		    curReg->queryEnd	= curWin->queryEnd;
		    curReg->alignEnd	= curWin->alignEnd;
		    /*
		      }
		    */
		}
	    }

	    /* end region */
	    scoreRegion(curReg, align, pctIDMethod);

	    if (DEBUG)
		writeConservedRegion(curReg, "Finished conserved region");

	    if (truncate && truncateRegionSimple(curReg, align, pctIDMethod,
						    threshold))
	    {
		if (DEBUG)
		    writeConservedRegion(curReg, "Truncated conserved region");

		/*
		 * We may have truncated region back beyond where a
		 * previous window scored above threshold. So backtrack
		 * windows to just after current region starts.
		 */
		while (curWin->queryStart > curReg->queryEnd + 2) {
		    curWinIdx--;
		    curWin = &sim->wins[curWinIdx];
		}
	    }
	}

	curWinIdx++;
    }
	

    if (DEBUG) {
	puts("\nInitial conserved regions:\n");
	writeConservedRegions(reg);
    }

    /*
    merged = 1;
    while (merged) {
    */
	/*
	 * Extend any windows which are above threshold to maximum span such
	 * that they are still above threshold
	 */
	/*
	extendRegions(reg, align, pctIDMethod, threshold);
	*/
	/*
	 * Truncate any windows which are below threshold such that they score
	 * above threshold
	 */
	/*
	truncateRegions(reg, align, pctIDMethod, threshold);
	*/
	/*
	 * Merge any windows which overlap or which are proximal and whose
	 * combined score is still above threshold
	 */
	/*
	merged = mergeRegions(reg, align, pctIDMethod, winSize, threshold);
	*/
    /*
    }
    */
    for (curRegIdx = 0; curRegIdx < reg->nRegs; curRegIdx++) {
    	writeConservedAlignment(&reg->regs[curRegIdx], align);
    }

    /*
     * Trim all regions so that each edge is a nucleotide match and re-score
     */
    for (curRegIdx = 0; curRegIdx < reg->nRegs; curRegIdx++) {
    	if (trimRegion(&reg->regs[curRegIdx], align, TRIM_BOTH))
	    scoreRegion(&reg->regs[curRegIdx], align, pctIDMethod);
    }

    if (DEBUG) {
	puts("\nExtended/truncated/merged conserved regions:\n");
	writeConservedRegions(reg);
    }

    if (reg) {
        if (reg->nRegs <= 0) {
	    RegionsFree(reg);
	    reg = NULL;
	}
    }

    return reg;
}

/*
 * Combine similarity windows to make regions of conservation
 * New version of computeConservedRegions. DJA 2005/06/06.
 */
static Regions *
computeConservedRegionsOld(const Alignment *align, const Similarities *sim,
	float threshold, int winSize, const Features *feat, int complemented,
	char pctIDMethod)
{
    static const int	regSizeInc = 100;

    int			nWins;
    int			nFeats;
    int			curWinIdx, lastMergedWinIdx, nextWinIdx;
    int			curRegIdx;
    int			curFeatIdx;
    int			queryMid;
    Regions		*reg = NULL;
    SimilarityWindow	*curWin = NULL;
    SimilarityWindow	*nextWin = NULL;
    ConservedRegion	*curReg = NULL;
    Feature		*curFeat = NULL;

    if (!align)
        return NULL;

    if (!sim)
        return NULL;

    /*
     * In order to improve exon filtering, we used to assign all windows which
     * overlap an exon a score of zero. Instead, we now mask exons before
     * computing similarity scores. We also improved the filterRegions routine
     * such that if an exonic part of a conserved region is removed, if this
     * causes the region score to fall below the threshold, the region is
     * truncated on the opposite side until the score of the region once again
     * falls above threshold.
     * DJA 2006/02/21
     */
    /*
    if (feat) {
	nFeats = feat->nFeats;
	curFeatIdx = 0;
	while (curFeatIdx < nFeats) {
	    curFeat = &feat->feats[curFeatIdx];
	    curWinIdx = 0;
	    nWins = sim->nWins;
	    while (curWinIdx < nWins) {
		curWin = &sim->wins[curWinIdx];
		queryMid = curWin->queryStart
			    + (int) nearbyint((double) (curWin->queryEnd
						    - curWin->queryStart) / 2);
		if (queryMid >= curFeat->seqStart
			&& queryMid <= curFeat->seqEnd)
		{
		    curWin->score = 0;
		}
		
		curWinIdx++;
	    }
	    curFeatIdx++;
	}
    }
    */

    reg = RegionsAlloc(regSizeInc);

    nWins = sim->nWins;

    curWinIdx = 0;
    curRegIdx = 0;
    while (curWinIdx < nWins) {
	/* Get next similarity window which scores above threshold. */
	if (curWinIdx < nWins)
	    curWin = &sim->wins[curWinIdx];
	while (curWinIdx < nWins && curWin->score < threshold) {
	    curWinIdx++;
	    if (curWinIdx < nWins)
		curWin = &sim->wins[curWinIdx];
	}

	if (curWinIdx < nWins) {
	    /* create tentative conserved region from this window */
	    curReg = &reg->regs[curRegIdx];
	    reg->nRegs++;
	    curReg->alignStart	= curWin->alignStart;
	    curReg->alignEnd	= curWin->alignEnd;
	    curReg->queryStart	= curWin->queryStart;
	    curReg->queryEnd	= curWin->queryEnd;
	    curReg->score	= curWin->score;

	    if (DEBUG) {
		writeConservedRegion(curReg, "New tentative conserved region");
	    }

	    /*
	     * As long as we have windows which overlap the current conserved
	     * region and they score above threshold, merge them with the
	     * current current region
	     */
	    nextWinIdx = curWinIdx + 1;
	    if (nextWinIdx < nWins)
		nextWin = &sim->wins[nextWinIdx];
	    while (nextWinIdx < nWins
		    && nextWin->queryStart <= curReg->queryEnd + 1)
	    {
		if (nextWin->score >= threshold) {
		    curReg->alignEnd = nextWin->alignEnd;
		    curReg->queryEnd = nextWin->queryEnd;
		    lastMergedWinIdx = nextWinIdx;
		}
		nextWinIdx++;
		if (nextWinIdx < nWins)
		    nextWin = &sim->wins[nextWinIdx];
	    }

	    scoreRegion(curReg, align, pctIDMethod);
	    if (DEBUG) {
		writeConservedRegion(curReg, "Merged conserved region");
	    }

	    /*
	     * If the overall merged region has fallen below threshold,
	     * cut if back on the right until it scores above threshold.
	     */
	    while (curReg->score < threshold) {
		lastMergedWinIdx--;
		nextWin = &sim->wins[lastMergedWinIdx];
		curReg->alignEnd = nextWin->alignEnd;
		curReg->queryEnd = nextWin->queryEnd;
		scoreRegion(curReg, align, pctIDMethod);
	    }

	    /*
	     * Trim region on both sides such that it starts and ends on a
	     * match and re-score.
	     */
	    if (trimRegion(curReg, align, TRIM_BOTH))
		scoreRegion(curReg, align, pctIDMethod);

	    if (DEBUG) {
		writeConservedRegion(curReg, "Trimmed conserved region");
	    }

	    /*
	     * Start processing next window after conserved region.
	     */
	    while (nextWinIdx < nWins
		    && nextWin->queryStart <= curReg->queryEnd)
	    {
		nextWinIdx++;
		if (nextWinIdx < nWins)
		    nextWin = &sim->wins[nextWinIdx];
	    }
	    curWinIdx = nextWinIdx;

	    curRegIdx++;
	    if (curRegIdx >= reg->maxRegs)
		reg = RegionsIncreaseSize(reg, regSizeInc);
	    curReg = &reg->regs[curRegIdx];
	}
    }

    if (reg) {
        if (reg->nRegs <= 0) {
	    RegionsFree(reg);
	    reg = NULL;
	}
    }

    return reg;
}

static int
truncateRegions(Regions *reg, const Alignment *align, char pctIDMethod,
		float threshold)
{
    int			curRegIdx, truncated; 
    ConservedRegion	*curReg = NULL;

    truncated = 0;

    curRegIdx = 0;
    while (curRegIdx < reg->nRegs - 1) {
    	curReg = &reg->regs[curRegIdx]; 	
	if (curReg->score < threshold) {
	    if (truncateRegionSimple(curReg, align, pctIDMethod, threshold))
	    	truncated = 1;
	}
	curRegIdx++;
    }

    return truncated;
}

static int
extendRegions(Regions *reg, const Alignment *align, char pctIDMethod,
		float threshold)
{
    int			curRegIdx, extended; 
    ConservedRegion	*curReg = NULL;

    extended = 0;

    curRegIdx = 0;
    while (curRegIdx < reg->nRegs - 1) {
    	curReg = &reg->regs[curRegIdx]; 	
	if (curReg->score > threshold) {
	    if (extendRegion(curReg, align, pctIDMethod, threshold))
	    	extended = 1;
	}
	curRegIdx++;
    }

    return extended;
}

static int
truncateRegionComplex(ConservedRegion *reg, const Alignment *align,
			char pctIDMethod, float threshold)
{
    char		lastTrunc;
    int			truncated, midLeft, midRight; 
    int			leftTailEnd, rightTailStart;
    float		leftPctID, rightPctID;
    float		mid;

    truncated = 0;

    if (trimRegion(reg, align, TRIM_BOTH))
    	scoreRegion(reg, align, pctIDMethod);

    while (reg->score < threshold && reg->queryStart < reg->queryEnd) {
	mid =  reg->queryStart + (float) (reg->queryEnd - reg->queryStart) / 2;
	midLeft = (int) floorf(mid);
	midRight = (int) ceilf(mid);
	if (midLeft == midRight) {
	    midLeft--;
	    midRight++;
	}

	leftTailEnd = reg->queryStart + REG_TAIL_LEN - 1;
	rightTailStart = reg->queryEnd - REG_TAIL_LEN + 1;
	if (leftTailEnd > midLeft)
	    leftTailEnd = midLeft;
	if (rightTailStart < midRight)
	    rightTailStart = midRight;

	leftPctID = 0;
	rightPctID = 0;
	leftPctID = percentageIdentity(align, reg->queryStart, leftTailEnd,
					pctIDMethod);

	rightPctID = percentageIdentity(align, rightTailStart, reg->queryEnd,
					pctIDMethod);

	lastTrunc = '\0';
	if (leftPctID > rightPctID) {
	    /* truncate on right */
	    reg->queryEnd--;
	    reg->alignEnd = seqToAlignPos(align, 1, reg->queryEnd);
	    trimRegion(reg, align, TRIM_RIGHT);
	    scoreRegion(reg, align, pctIDMethod);
	    truncated = 1;
	    lastTrunc = 'r';
	} else if (rightPctID > leftPctID) {
	    /* truncate on left */
	    reg->queryStart++;
	    reg->alignStart = seqToAlignPos(align, 1, reg->queryStart);
	    trimRegion(reg, align, TRIM_LEFT);
	    scoreRegion(reg, align, pctIDMethod);
	    truncated = 1;
	    lastTrunc = 'l';
	} else {
	    /*
	     * Left and right percentage identities are the same
	     * so truncate on opposite side from last time
	     */
	    if (lastTrunc == 'r') {
		reg->queryStart++;
		reg->alignStart = seqToAlignPos(align, 1, reg->queryStart);
		trimRegion(reg, align, TRIM_LEFT);
		scoreRegion(reg, align, pctIDMethod);
		truncated = 1;
		lastTrunc = 'l';
	    } else {
		reg->queryEnd--;
		reg->alignEnd = seqToAlignPos(align, 1, reg->queryEnd);
		trimRegion(reg, align, TRIM_RIGHT);
		scoreRegion(reg, align, pctIDMethod);
		truncated = 1;
		lastTrunc = 'r';
	    }
	}
    }

    return truncated;
}

static int
truncateRegionSimple(ConservedRegion *reg, const Alignment *align,
			char pctIDMethod, float threshold)
{
    int			truncated = 0;

    if (trimRegion(reg, align, TRIM_BOTH))
    	scoreRegion(reg, align, pctIDMethod);

    while (reg->score < threshold) {
	/* truncate on right */
	reg->queryEnd--;
	reg->alignEnd = seqToAlignPos(align, 1, reg->queryEnd);
	trimRegion(reg, align, TRIM_RIGHT);
	scoreRegion(reg, align, pctIDMethod);
	truncated = 1;
    }

    return truncated;
}

static int
extendRegion(ConservedRegion *reg, const Alignment *align, char pctIDMethod,
		float threshold)
{
    int			extended, this_extended; 
    float		leftPctID, rightPctID;

    extended = 0;

    this_extended = 1;
    while (this_extended) {
    	this_extended = 0;
	leftPctID = 0;
	rightPctID = 0;
	if (reg->queryStart > 1) {
	    leftPctID = percentageIdentity(align, reg->queryStart - 1,
					    reg->queryEnd, pctIDMethod);
	}

	if (reg->queryEnd < align->seq1Len) {
	    rightPctID = percentageIdentity(align, reg->queryStart,
					    reg->queryEnd + 1, pctIDMethod);
	}

	if (leftPctID > rightPctID) {
	    if (leftPctID >= threshold) {
		reg->queryStart--;
		reg->alignStart = seqToAlignPos(align, 1, reg->queryStart);
		reg->score = leftPctID;
		this_extended = 1;
		extended = 1;
	    }
	} else if (rightPctID > leftPctID) {
	    if (rightPctID >= threshold) {
		reg->queryEnd++;
		reg->alignEnd = seqToAlignPos(align, 1, reg->queryEnd);
		reg->score = rightPctID;
		this_extended = 1;
		extended = 1;
	    }
	} else {
	    /*
	     * Left and right percentage identities are the same (perhaps
	     * both 0). Bias toward extending on left.
	     */
	    if (leftPctID >= threshold) {
		reg->queryStart--;
		reg->alignStart = seqToAlignPos(align, 1, reg->queryStart);
		reg->score = leftPctID;
		this_extended = 1;
		extended = 1;
	    } else if (rightPctID >= threshold) {
		reg->queryEnd++;
		reg->alignEnd = seqToAlignPos(align, 1, reg->queryEnd);
		reg->score = rightPctID;
		this_extended = 1;
		extended = 1;
	    }
	}
    }

    return extended;
}

static int
mergeRegions(Regions *reg, const Alignment *align, char pctIDMethod,
		int winSize, float threshold)
{
    int			i, curRegIdx, merged; 
    float		pctID;
    ConservedRegion	*curReg = NULL, *nextReg = NULL;

    merged = 0;

    curRegIdx = 0;
    while (curRegIdx < reg->nRegs - 1) {
    	curReg = &reg->regs[curRegIdx]; 	
    	nextReg = &reg->regs[curRegIdx + 1]; 	

	if (nextReg->queryStart <= curReg->queryEnd + 1) {
	    /*
	     * Overlapping windows
	     */
	    if (DEBUG) {
	    	writeConservedRegion(curReg, "Merging overlapping region");
	    	writeConservedRegion(nextReg, "with region");
	    }

	    if (nextReg->queryEnd > curReg->queryEnd) {
		curReg->queryEnd = nextReg->queryEnd;
		curReg->alignEnd = nextReg->alignEnd;
		scoreRegion(curReg, align, pctIDMethod);
	    }

	    memset(nextReg, 0, sizeof(ConservedRegion));
	    for (i = curRegIdx + 2; i < reg->nRegs; i++) {
	    	memcpy(&reg->regs[i - 1], &reg->regs[i],
			sizeof(ConservedRegion));
	    }
	    reg->nRegs--;
	    merged = 1;
	/*
	} else if (nextReg->queryStart < curReg->queryEnd + winSize) {
	*/
	} else {
	    /*
	     * Proximal windows
	     */
	    pctID = percentageIdentity(align, curReg->queryStart,
	    				nextReg->queryEnd, pctIDMethod);
	    if (pctID >= threshold) {
	    	/*
		 * Merge if combined score is above threshold
		 */
		if (DEBUG) {
		    writeConservedRegion(curReg, "Merging region");
		    writeConservedRegion(nextReg, "with region");
		}

		curReg->queryEnd = nextReg->queryEnd;
		curReg->alignEnd = nextReg->alignEnd;
		curReg->score = pctID;

		memset(nextReg, 0, sizeof(ConservedRegion));
		for (i = curRegIdx + 2; i < reg->nRegs; i++) {
		    memcpy(&reg->regs[i - 1], &reg->regs[i],
			    sizeof(ConservedRegion));
		}
		reg->nRegs--;
		merged = 1;
	    } else {
		curRegIdx++;
	    }
	/*
	} else {
	    curRegIdx++;
	*/
	}
    }

    return merged;
}

static void
computeConservedRegionsHitCoords(const Alignment *align, Regions *reg,
	int complemented)
{
    int			regIdx; 
    int			alnIdx; 
    int			seq1Nucs, seq2Nucs, numSeq1Nucs, numSeq2Nucs;
    int			*mapPos;
    char		*alnSeq1, *alnSeq2;
    ConservedRegion	*cr;

    alnSeq1 = align->seq1;
    alnSeq2 = align->seq2;

    numSeq1Nucs = 0;
    numSeq2Nucs = 0;
    alnIdx = 0;
    while (alnIdx < align->length) {
	if (alnSeq1[alnIdx] != GAP_CHAR)
	    numSeq1Nucs++;
	if (alnSeq2[alnIdx] != GAP_CHAR)
	    numSeq2Nucs++;
	alnIdx++;
    }

    /*
     * Create array mapping sequence 1 positions to aligned sequence
     * 2 positions
     */
    mapPos = calloc(numSeq1Nucs + 1, sizeof(int));
    memset(mapPos, 0, (numSeq1Nucs + 1) * sizeof(int));

    seq1Nucs = 0;
    seq2Nucs = 0;
    alnIdx = 0;
    while (alnIdx < align->length) {
	if (alnSeq2[alnIdx] != GAP_CHAR) {
	    seq2Nucs++;
	}
	if (alnSeq1[alnIdx] != GAP_CHAR) {
	    seq1Nucs++;
	    if (alnSeq2[alnIdx] != GAP_CHAR) {
		if (complemented)
		    mapPos[seq1Nucs] = numSeq2Nucs - seq2Nucs + 1;
		else
		    mapPos[seq1Nucs] = seq2Nucs;
	    }
	}
	alnIdx++;
    }

    regIdx = 0;
    while (regIdx < reg->nRegs) {
	cr = &(reg->regs[regIdx]);
	if (complemented) {
	    cr->hitStart = mapPos[cr->queryEnd];
	    cr->hitEnd = mapPos[cr->queryStart];
	} else {
	    cr->hitStart = mapPos[cr->queryStart];
	    cr->hitEnd = mapPos[cr->queryEnd];
	}
	regIdx++;
    }
}

/*
 * Compute the percentage identity of sequence 1 between start and
 * end (specified in 1-based coords)
 */
static float
percentageIdentity(const Alignment *align, int start, int end, char pctIDMethod)
{
    int		temp, i, matchScore;
    int		alnStart, alnEnd;
    float	pctID = 0.0;

    if (!align)
	return 0;

    if (!align->matches)
	return 0;

    if (end < start) {
	temp = end;
	end = start;
	start = temp;
    }

    if (start < 1)
    	start = 1;
    if (end > align->seq1Len)
    	end = align->seq1Len;

    matchScore = 0;
    for (i = start - 1; i <= end - 1; i++)
    	if (align->matches[i])
	    matchScore++;

    if (pctIDMethod == 's') {
	pctID = (float) matchScore / (end - start + 1);
    } else if (pctIDMethod == 'o') {
	alnStart = seqToAlignPos(align, 1, start);
	alnEnd = seqToAlignPos(align, 1, end);
	pctID = (float) matchScore / (alnEnd - alnStart + 1);
    }

    return pctID;
}

static int
trimRegion(ConservedRegion *reg, const Alignment *align, WhichTrim whichTrim)
{
    int		alignPos, seqPos, trimmed;
    const char	*seq1, *seq2;

    if (!reg || !align)
	return 0;

    seq1 = align->seq1;
    seq2 = align->seq2;

    trimmed = 0;

    if (whichTrim & TRIM_LEFT) {
	seqPos = reg->queryStart;
	alignPos = reg->alignStart;
	while (!nucMatch(seq1[alignPos - 1], seq2[alignPos - 1])) {
	    seqPos++;
	    alignPos = seqToAlignPos(align, 1, seqPos);
	    trimmed = 1;
	}
	reg->queryStart = seqPos;
	reg->alignStart = alignPos;
    }

    if (whichTrim & TRIM_RIGHT) {
	seqPos = reg->queryEnd;
	alignPos = reg->alignEnd;
	while (!nucMatch(seq1[alignPos - 1], seq2[alignPos - 1])) {
	    seqPos--;
	    alignPos = seqToAlignPos(align, 1, seqPos);
	    trimmed = 1;
	}
	reg->queryEnd = seqPos;
	reg->alignEnd = alignPos;
    }

    return trimmed;
}

static float
scoreRegion(ConservedRegion *reg, const Alignment *align, char pctIDMethod)
{
    float	pctIdentity = 0.0;

    pctIdentity = percentageIdentity(align, reg->queryStart, reg->queryEnd,
					pctIDMethod);

    reg->score = pctIdentity;

    /*
    if (DEBUG) {
    	writeConservedRegion(reg, "Scored conserved region");
    }
    */

    return pctIdentity;
}

/*
static int
alignToSeqPos(const char *seq, int alnPos, char match_type)
{
    int			seqPos = -1;
    int			alnIdx = 0;

    if (alnPos < 0) {
	return -1;
    }

    if (!seq || strlen(seq) <= alnPos) {
	return -1;
    }

    if (!match_type) {
	match_type = '=';
    }

    if (match_type == '=' || match_type == '<') {
	while (alnIdx <= alnPos) {
	    if (seq[alnIdx++] != GAP_CHAR)
		seqPos++;
	}
    } else if (match_type == '>') {
	while (alnIdx < strlen(seq)) {
	    if (seq[alnIdx++] != GAP_CHAR) {
		seqPos++;
		if (alnIdx > alnPos)
		    break;
	    }
	}
	if (alnIdx >= strlen(seq))
	    seqPos = -1;
    }

    return seqPos;
}
*/

/*
 * Convert a sequence coordinate to an aligment coordinate. The sequence
 * coordinate (seqPos) and the returned aligment coordinate are in the
 * 1-based coordinate system.
 */
static int
seqToAlignPos(const Alignment *align, int seqNum, int seqPos)
{
    int alnPos = 0;

    if (!align)
	return 0;

    if (seqNum == 1) {
	if (seqPos >= 1 && seqPos <= align->seq1Len)
	    if (align->seq1ToAlnIdx)
		alnPos = align->seq1ToAlnIdx[seqPos - 1] + 1;

    } else if (seqNum == 2) {
	if (seqPos >= 1 || seqPos <= align->seq2Len)
	    if (align->seq2ToAlnIdx)
		alnPos = align->seq2ToAlnIdx[seqPos - 1] + 1;
    }

    return alnPos;
}


static Regions *
filterRegions(Regions *reg, const Features *feat, Alignment *align,
		float threshold, char pctIDMethod)
{
    int			curRegNum, curFeatNum;

    ConservedRegion	*curReg = NULL;
    Feature		*curFeat = NULL;

    if (!reg || !feat || !align)
	return NULL;

    /*
     * Sort features.
     * The regions will already be sorted.
     */
    qsort(feat->feats, feat->nFeats, sizeof(Feature), featureCompare);

    curRegNum = 0;
    curFeatNum = 0;
    while (curRegNum < reg->nRegs && curFeatNum < feat->nFeats) {
	curReg = &reg->regs[curRegNum];
	curFeat = &feat->feats[curFeatNum];

	if (curFeat->seqEnd < curReg->queryStart) {
	    curFeatNum++;
	} else if (curFeat->seqStart > curReg->queryEnd) {
	    curRegNum++;
	} else if (curFeat->seqStart <= curReg->queryStart) {
	    if (curFeat->seqEnd >= curReg->queryStart) {
		if (DEBUG) {
		    writeFeature(curFeat, "filterRegions: processing feature");
		}

		if (curFeat->seqEnd >= curReg->queryEnd) {
		    /* delete region altogether */
		    if (DEBUG) {
			writeConservedRegion(
			    curReg, "Deleting conserved region within exon");
		    }
		    reg = RegionsDeleteRegionNum(reg, curRegNum);
		} else {
		    /* truncate region on left at exon end */
		    curReg->queryStart = curFeat->seqEnd + 1;
		    curReg->alignStart = seqToAlignPos(
					    align, 1, curReg->queryStart);

		    /*
		     * After exon is removed from region it can result in
		     * score falling below threshold. Truncate the region
		     * on the opposite side until score increases back
		     * above threshold. DJA 2006/02/21
		     */
		    scoreRegion(curReg, align, pctIDMethod);
		    if (DEBUG) {
			writeConservedRegion(
			    curReg, "CR truncated at left exon");
		    }
		    while (curReg->score < threshold
			    && curReg->queryEnd >= curReg->queryStart)
		    {
			curReg->queryEnd--;
			curReg->alignEnd = seqToAlignPos(
						align, 1, curReg->queryEnd);
			scoreRegion(curReg, align, pctIDMethod);
		    }
		    if (trimRegion(curReg, align, TRIM_RIGHT))
			scoreRegion(curReg, align, pctIDMethod);
		    if (DEBUG) {
			writeConservedRegion(
			    curReg, "Trimmed CR truncated at left exon");
		    }

		    curFeatNum++;
		}
	    }
	} else {
	    if (curFeat->seqEnd >= curReg->queryEnd) {
		if (DEBUG) {
		    writeFeature(curFeat, "filterRegions: processing feature");
		}

	        /* truncate region on right at exon start */
		curReg->queryEnd = curFeat->seqStart - 1;
		curReg->alignEnd = seqToAlignPos(align, 1, curReg->queryEnd);

		/*
		 * After exon is removed from region it can result in
		 * score falling below threshold. Truncate the region
		 * on the opposite side until score increases back
		 * above threshold. DJA 2006/02/21
		 */
		scoreRegion(curReg, align, pctIDMethod);
		if (DEBUG) {
		    writeConservedRegion(
			curReg, "CR truncated at right exon");
		}
		while (curReg->score < threshold
			&& curReg->queryStart <= curReg->queryEnd)
		{
		    curReg->queryStart++;
		    curReg->alignStart = seqToAlignPos(
		    				align, 1, curReg->queryStart);
		    scoreRegion(curReg, align, pctIDMethod);
		}
		if (trimRegion(curReg, align, TRIM_LEFT))
		    scoreRegion(curReg, align, pctIDMethod);
		if (DEBUG) {
		    writeConservedRegion(
			curReg, "Trimmed CR truncated at right exon");
		}

		curRegNum++;
	    } else {
		/* split region in two */
		if (DEBUG) {
		    writeFeature(curFeat, "filterRegions: processing feature");
		}

		if (DEBUG) {
		    writeConservedRegion(curReg, "CR split at feature");
		}

		reg = RegionsSplitRegionNum(
				    reg, curRegNum,
				    curFeat->seqStart, curFeat->seqEnd,
				    align, threshold, pctIDMethod);
		curRegNum++;
		curFeatNum++;
	    }
	}
    }

    return reg;
}

static Regions *
RegionsDeleteRegionNum(Regions *reg, int regNum)
{
    int		i;

    if (!reg)
	return NULL;

    if (regNum < 0 || regNum >= reg->nRegs)
	return reg;

    for (i = regNum + 1; i < reg->nRegs; i++) {
	reg->regs[i - 1] = reg->regs[i];
    }
    reg->nRegs--;

    return reg;
}

Regions *
RegionsSplitRegionNum(Regions *reg, int regNum, int start, int end,
			Alignment *align, float threshold,
			char pctIDMethod)
{
    static const int	regSizeInc = 50;
    int			i;
    ConservedRegion	*reg1, *reg2;

    if (!reg)
	return NULL;

    if (regNum < 0 || regNum >= reg->nRegs)
	return reg;

    if (reg->nRegs >= reg->maxRegs)
	reg = RegionsIncreaseSize(reg, regSizeInc);

    /* Make room for right half of split region */
    for (i = reg->nRegs - 1; i > regNum; i--) {
	reg->regs[i + 1] = reg->regs[i];
    }

    /* make a copy of region to split into new slot */
    reg1 = &reg->regs[regNum];
    memcpy(&reg->regs[regNum + 1], reg1, sizeof(ConservedRegion));
    reg2 = &reg->regs[regNum + 1];

    reg1->queryEnd = start - 1;
    reg1->alignEnd = seqToAlignPos(align, 1, reg1->queryEnd);

    reg2->queryStart = end + 1;
    reg2->alignStart = seqToAlignPos(align, 1, reg2->queryStart);

    /*
     * After exon is removed from region it can result in
     * score falling below threshold. Truncate the region
     * on the opposite side until score increases back
     * above threshold. DJA 2006/02/21
     */
    scoreRegion(reg1, align, pctIDMethod);
    while (reg1->score < threshold && reg1->queryStart < reg1->queryEnd) {
	reg1->queryStart++;
	reg1->alignStart = seqToAlignPos(align, 1, reg1->queryStart);
	scoreRegion(reg1, align, pctIDMethod);
    }

    if (trimRegion(reg1, align, TRIM_LEFT))
	scoreRegion(reg1, align, pctIDMethod);
    if (DEBUG) {
	writeConservedRegion(reg1, "Trimmed left split CR");
    }

    scoreRegion(reg2, align, pctIDMethod);
    while (reg2->score < threshold && reg2->queryEnd > reg2->queryStart) {
	reg2->queryEnd--;
	reg2->alignEnd = seqToAlignPos(align, 1, reg2->queryEnd);
	scoreRegion(reg2, align, pctIDMethod);
    }

    if (trimRegion(reg2, align, TRIM_RIGHT))
	scoreRegion(reg2, align, pctIDMethod);
    if (DEBUG) {
	writeConservedRegion(reg2, "Trimmed right split CR");
    }

    reg->nRegs++;

    return reg;
}

/*
 * Simple nucleotide match check.
 * NOTES:
 * Returns either true (1) or false (0).
 * Only matches DNA nucleotides (not proteins).
 * Gaps and masked nucleotides are considered mis-matches.
 */
inline static int
nucMatch(char nuc1, char nuc2)
{
    char	uNuc1 = toupper(nuc1);
    char	uNuc2 = toupper(nuc2);

    if (uNuc1 == uNuc2 && (uNuc1 == 'A' || uNuc1 == 'C' || uNuc1 == 'G'
                || uNuc1 == 'T'))
        return 1;

    return 0;
}

static void
writeSimilarity(const Similarities *sim, const char *file, char format)
{
    int		winIdx;
    int		mid;
    FILE	*fp;

    if (!sim)
        return;

    if (!file)
        fp = stdout;
    else
        fp = fopen(file, "w");

    if (!fp) {
        fprintf(stderr,
	    "\nError writeSimilarity: could not open file %s for writing\n",
	    file);
        return;
    }

    winIdx = 0;
    for (winIdx = 0; winIdx < sim->nWins; winIdx++) {
        switch (format) {
        case 's':
        case 'S':
            fprintf(fp, "%d\t%7.4f\n",
                    sim->wins[winIdx].queryStart,
                    sim->wins[winIdx].score);
            break;
        case 'e':
        case 'E':
            fprintf(fp, "%d\t%7.4f\n",
                    sim->wins[winIdx].queryEnd,
                    sim->wins[winIdx].score);
            break;
        case 'c':
        case 'C':
            mid = (int) ((sim->wins[winIdx].queryStart
			+ sim->wins[winIdx].queryEnd) / 2);
            fprintf(fp, "%d\t%7.4f\n", mid, sim->wins[winIdx].score);
            break;
        default:
            fprintf(fp, "%d\t%7.4f\n",
		sim->wins[winIdx].queryStart, sim->wins[winIdx].score);
            break;
        }
    }

    fclose(fp);
}

static void
writeConservedRegion(const ConservedRegion *cr, const char *msg)
{
    printf("%s:\t%d - %d (%d - %d); %6.4f\n",
		msg,
    		cr->queryStart, cr->queryEnd,
    		cr->alignStart, cr->alignEnd,
		cr->score);
}

static void
writeConservedAlignment(const ConservedRegion *cr, const Alignment *aln)
{
    char	*seq1, *seq2;
    char	subSeq1[51], subSeq2[51], cigarLine[51];
    int		i, j;
    int		matches = 0, winSize = 0;

    printf("\n%d - %d (%d - %d); %6.4f\n",
    		cr->queryStart, cr->queryEnd,
    		cr->alignStart, cr->alignEnd,
		cr->score);

    seq1 = aln->seq1;
    seq2 = aln->seq2;

    memset(subSeq1, 0, 51);
    memset(subSeq2, 0, 51);
    memset(cigarLine, 0, 51);

    j = 0;
    for (i = cr->alignStart - 1; i <= cr->alignEnd - 1; i++) {
    	if (seq1[i] != '-') {
	    winSize++;
	}
	subSeq1[j] = seq1[i];
	subSeq2[j] = seq2[i];
	if (nucMatch(seq1[i], seq2[i])) {
	    cigarLine[j] = '|';
	    matches++;
	} else {
	    cigarLine[j] = ' ';
	}
    	j++;
	if (j == 50 || i == cr->alignEnd - 1) {
	    printf("\n%s\n%s\n%s\n", subSeq1, cigarLine, subSeq2);
	    memset(subSeq1, 0, 51);
	    memset(subSeq2, 0, 51);
	    memset(cigarLine, 0, 51);
	    j = 0;
	}
    }
    printf("\n%d  %d  %6.4f\n\n", winSize, matches, (double) matches / winSize);
}

static void
writeFeature(const Feature *feat, const char *msg)
{
    printf("%s:\n%d - %d\n", msg, feat->seqStart, feat->seqEnd);
}

static void
writeConservedRegionsReport(const Regions *reg, const Alignment *align,
	    int winSize, int winInc, float stringency, float min_threshold,
	    float dynamic_threshold, float threshold, int minLen,
	    const char *file)
{
    int			regIdx;
    FILE		*fp;
    ConservedRegion	*cr;

    if (!reg)
        return;

    if (!file)
        fp = stdout;
    else
        fp = fopen(file, "w");

    if (minLen < 1)
	minLen = 1;

    if (!fp) {
        fprintf(stderr, "\nError writeConservedRegionsReport: could not open"
			" file %s for writing\n", file);
        return;
    }

    fprintf(fp, "Sequence 1          : %s\n", align->seq1ID);
    fprintf(fp, "Sequence 2          : %s\n", align->seq2ID);
    fprintf(fp, "Window Size         : %d bp\n", winSize);
    fprintf(fp, "Window Increment    : %d bp\n", winInc);
    if (stringency) {
	fprintf(fp, "Top X%% Identities   : %.1f%%\n", stringency * 100);
	fprintf(fp, "Computed %% Identity : %.1f%%\n", dynamic_threshold * 100);
    }
    if (min_threshold)
	fprintf(fp, "Minimum %% Identity  : %.1f%%\n", min_threshold * 100);
    fprintf(fp, "Effective %% Identity: %.1f%%\n\n", threshold * 100);
    regIdx = 0;
    for (regIdx = 0; regIdx < reg->nRegs; regIdx++) {
	cr = &(reg->regs[regIdx]);
	if (cr->queryEnd - cr->queryStart + 1 >= minLen
		/* && cr->score >= threshold */)
	{
	    fprintf(fp, "%7d %7d %7d   %7d %7d %7d    %7d %7d    %6.4f\n",
		cr->queryStart,
		cr->queryEnd,
		cr->queryEnd - cr->queryStart + 1,
		cr->hitStart,
		cr->hitEnd,
		cr->hitEnd - cr->hitStart + 1,
		cr->alignStart,
		cr->alignEnd,
		cr->score);
	}
    }

    fclose(fp);
}

static void
writeConservedRegions(const Regions *reg)
{
    int			regIdx;
    ConservedRegion	*cr;

    if (!reg)
        return;

    regIdx = 0;
    for (regIdx = 0; regIdx < reg->nRegs; regIdx++) {
	cr = &(reg->regs[regIdx]);
	printf("%7d %7d %7d (%7d %7d) %8.6f\n",
	    cr->queryStart,
	    cr->queryEnd,
	    cr->queryEnd - cr->queryStart + 1,
	    cr->alignStart,
	    cr->alignEnd,
	    cr->score);
    }
}

static void
writeSimilarityWithAligns(const Similarities *sim, const Alignment *align,
        const char *file)
{
    int		winIdx;
    int		start, end, length;
    FILE	*fp;

    if (!sim)
        return;

    if (!file)
        fp = stdout;
    else
        fp = fopen(file, "w");

    if (!fp) {
        fprintf(stderr, "\nError writeSimilarityWithAligns: could not open"
			" file %s for writing\n", file);
        return;
    }

    winIdx = 0;
    for (winIdx = 0; winIdx < sim->nWins; winIdx++) {
        start = sim->wins[winIdx].alignStart - 1, 
        end = sim->wins[winIdx].alignEnd - 1,
        length = end - start + 1;
        fprintf(fp, "%.*s\n%.*s\n",
	    length, &align->seq1[start],
	    length, &align->seq2[start]);
        fprintf(fp, "%8d %8d %8d %8d %7.4f\n",
	    sim->wins[winIdx].alignStart, 
	    sim->wins[winIdx].alignEnd,
	    sim->wins[winIdx].queryStart, 
	    sim->wins[winIdx].queryEnd,
	    sim->wins[winIdx].score);
    }

    fclose(fp);
}

Alignment *
AlignmentAlloc(size_t size)
{
    Alignment	*align = NULL;

    if (size < 1)
	return NULL;

    align = (Alignment *) malloc(sizeof(Alignment));
    if (!align) {
	fputs("AlignmentAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(align, 0, sizeof(Alignment));

    align->seq1 = (char *) malloc(size);
    if (!align->seq1) {
	fputs("AlignmentAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    align->seq2 = (char *) malloc(size);
    if (!align->seq2) {
	fputs("AlignmentAlloc: memory allocation error\n", stderr);
	exit(-1);
    }

    align->size = size;

    return align;
}

Alignment *
AlignmentIncreaseSize(Alignment *align, size_t sizeInc)
{
    if (!align)
	return NULL;

    if (sizeInc < 1)
	return align;

    align->seq1 = realloc(align->seq1, align->size + sizeInc);
    if (!align->seq1) {
	fputs("AlignmentIncreaseSize: memory allocation error\n", stderr);
	exit(-1);
    }
    align->seq2 = realloc(align->seq2, align->size + sizeInc);
    if (!align->seq2) {
	fputs("AlignmentIncreaseSize: memory allocation error\n", stderr);
	exit(-1);
    }

    align->size += sizeInc;

    return align;
}

static void
AlignmentFree(Alignment *align)
{
    if (align) {
        if (align->seq1)
            free(align->seq1);
        if (align->seq2)
            free(align->seq2);
        free(align);
    }
}

Similarities *
SimilaritiesAlloc(int nSimWins)
{
    Similarities	*sim = NULL;

    if (nSimWins < 1)
	return NULL;

    sim = (Similarities *) malloc(sizeof(Similarities));
    if (!sim) {
	fputs("SimilaritiesAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(sim, 0, sizeof(Similarities));

    sim->wins = (SimilarityWindow *) malloc(nSimWins *
	    		sizeof(SimilarityWindow));
    if (!sim->wins) {
	fputs("SimilaritiesAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(sim->wins, 0, nSimWins * sizeof(SimilarityWindow));

    return sim;
}

static void
SimilaritiesFree(Similarities *sim)
{
    if (sim) {
        if (sim->wins)
            free(sim->wins);
        free(sim);
    }
}

Regions *
RegionsAlloc(int nRegs)
{
    Regions	*reg = NULL;

    if (nRegs < 1)
	return NULL;

    reg = (Regions *) malloc(sizeof(Regions));
    if (!reg) {
	fputs("RegionsAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(reg, 0, sizeof(Regions));

    reg->regs = (ConservedRegion *) malloc(nRegs * sizeof(ConservedRegion));
    if (!reg->regs) {
	fputs("RegionsAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(reg->regs, 0, nRegs * sizeof(ConservedRegion));

    reg->maxRegs = nRegs;

    return reg;
}

Regions *
RegionsIncreaseSize(Regions *reg, int nRegs)
{
    if (!reg)
	return NULL;

    if (nRegs < 1)
	return reg;

    reg->regs = (ConservedRegion *) realloc(reg->regs, (reg->maxRegs + nRegs)
					* sizeof(ConservedRegion));
    if (!reg->regs) {
	fputs("RegionsIncreaseSize: memory allocation error\n", stderr);
	exit(-1);
    }
    reg->maxRegs += nRegs;

    return reg;
}

static void
RegionsFree(Regions *reg)
{
    if (reg) {
        if (reg->regs)
            free(reg->regs);
        free(reg);
    }
}

Features *
FeaturesAlloc(int nFeats)
{
    Features	*feat = NULL;

    if (nFeats < 1)
	return NULL;

    feat = (Features *) malloc(sizeof(Features));
    if (!feat) {
	fputs("FeaturesAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(feat, 0, sizeof(Features));

    feat->feats = (Feature *) malloc(nFeats * sizeof(Feature));
    if (!feat->feats) {
	fputs("FeaturesAlloc: memory allocation error\n", stderr);
	exit(-1);
    }
    memset(feat->feats, 0, nFeats * sizeof(Feature));

    feat->maxFeats = nFeats;

    return feat;
}

Features *
FeaturesIncreaseSize(Features *feat, int nFeats)
{
    if (!feat)
	return NULL;

    if (nFeats < 1)
	return feat;

    feat->feats = (Feature *) realloc(feat->feats, (feat->maxFeats + nFeats)
					* sizeof(Feature));
    if (!feat->feats) {
	fputs("FeaturesIncreaseSize: memory allocation error\n", stderr);
	exit(-1);
    }
    feat->maxFeats += nFeats;

    return feat;
}

static void
FeaturesFree(Features *feat)
{
    if (feat) {
        if (feat->feats)
            free(feat->feats);
        free(feat);
    }
}

/*
 * Same as the perl function of the same name.
 * If the last character of a string is a newline, replaces it with the NULL
 * character.
 */
static void
chomp(char *str)
{
    int		len;

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
    printf("\nUsage:\n"
        "align_cons [-i in_file] [-w window_size] [-n increment]\n"
        "\t\t[-t threshold] [-s stringency] [-g gff_file] [-r report_type]\n"
	"\t\t[-l length] [-f format] [-o out_file] [-D]\n");
    printf(
        "where:\n"
        "\t-i in_file     = Input alignment file in multi-FastA format.\n"
        "\t                 (default = stdin)\n"
        "\t-w window_size = Size of %% identity scoring window.\n"
	"\t                 (default %d)\n"
        "\t-n increment   = Amount to slide scoring window.\n"
	"\t                 (default = %d)\n"
        "\t-t theshold    = Min. %% identity score of conserved regions to\n"
	"\t                 report.\n"
        "\t-s stringency  = Stringency for dynamically computing the\n"
	"\t                 threshold. Stringency is interpreted to mean the\n"
	"\t                 top X%% of scoring windows. Threshold is computed\n"
	"\t                 such that X%% of the windows score above that\n"
	"\t                 threshold. If the -t parameter if also specified,\n"
	"\t                 and the dynamically computed threshold falls\n"
	"\t                 below the specified threshold, the specified\n"
	"\t                 threshold, the specified threshold overrides the\n"
	"\t                 computed threshold.\n"
        "\t                 (default = %4.2f)\n"
        "\t-m method      = Method used to compute the percent identity\n"
	"\t                 's' = standard\n"
	"\t                 'o' = overall\n"
	"\t                 Standard - computes the percent identity based\n"
	"\t                 based on the number of nucleotides in the second\n"
	"\t                 sequence which match nucleotides in the first.\n"
	"\t                 (i.e. gaps in the first sequence are ignored).\n"
	"\t                 Overall - computes the percent identity based\n"
	"\t                 on the number of nucleotides in the two sequences\n"
	"\t                 which match within a fixed alignment window.\n"
	"\t                 (i.e. gaps in either sequence are counted as\n"
	"\t		    mis-matches).\n"
        "\t                 (default = '%c')\n"
        "\t-g gff_file    = GFF file containing features to filter out of\n"
	"\t                 conserved regions\n"
	"\t-r report_type = A 'p' or 'c' indicating type of report to\n"
	"\t                 output. A 'p' outputs position and score\n"
	"\t                 suitable for making conservation plots.\n"
	"\t                 A 'c' outputs conserved regions. Format of this\n"
	"\t                 file is:\n"
	"\t                 seq_start seq_end length align_start align_end\n"
	"\t		    score\n"
	"\t                 (default = '%c')\n"
        "\t-l length      = Minimum length of conserved regions to report\n"
        "\t                 (default = %d)\n"
	"\t-W             = Report conserved regions resulting from merging\n"
	"\t                 of significant windows. The percent-identity\n"
	"\t                 threshold is applied to the windows, but not to the\n"
	"\t                 final conserved regions.\n"
	"\t                 Note: options -W and -g  may be incompatible.\n"
        "\t-f format      = For conservation plot reports (-r p), Character\n"
	"\t                 indicating whether the position output is start\n"
	"\t                 's', end 'e' or center 'c'\n"
        "\t                 (default = '%c')\n"
        "\t-o out_file    = Output file\n"
        "\t                 (default = stdout)\n"
        "\t-c             = Indicate sequence 2 was reverse complemented\n"
	"\t                 during alignment.\n"
	"\t                 In such cases this flag is required so that\n"
     	"\t                 the conserved region coordinates are\n"
	"\t                 correctly computed for the second sequence.\n"
        "\t-D             = Turn on debugging output\n"
        "\t-h             = Print this usage message\n"
        "\t-?             = Print this usage message\n",
        DFLT_WIN_SIZE, DFLT_WIN_INC, DFLT_STRINGENCY, DFLT_PCT_ID_METHOD,
	DFLT_REPORT_TYPE, DFLT_MIN_REG_LEN, DFLT_FORMAT);
}
