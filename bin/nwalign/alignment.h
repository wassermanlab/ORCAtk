/*
 * Alignment.h
 *
 * Defines the Alignment data type and associated functions.
 *
 */

#ifndef	_alignment_h_INCLUDED_
#define	_alignment_h_INCLUDED_

#include <stdio.h>

static const char	HomolChar = '*';

/*
 * Minimum alignment score.
 * Scores below this indicate no significant alignment.
 */
static const int	MinAlignmentScore = -9999;

typedef struct Alignment {
	size_t		size;		/* allocated size of seq1/seq2 */
	int		length;		/* actual length of seq1/seq2 */
	int		score;		/* the alignment score */
	char		*seq1;		/* nucleotide sequence 1 */
	char		*seq2;		/* nucleotide sequence 2 */
	char		*homology;	/* homology sequence */
} Alignment;

Alignment *	AlignmentCreate(size_t size);
void		AlignmentDestroy(Alignment *align);
int		AlignmentAdd(Alignment *align, char seq1char, char seq2char,
			char homolChar);
void		AlignmentPrint(const Alignment *align, FILE *fp);

#endif	/* _alignment_h_INCLUDED_ */
