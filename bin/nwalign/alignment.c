#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "nwalign.h"
#include "alignment.h"

extern int	lastErrCode;

/*
 * Allocate and initialize memory for an Alignment structure
 */
Alignment *
AlignmentCreate(size_t size)
{
	Alignment	*align = NULL;

	if (size < 1)
		return NULL;

	align = (Alignment *) malloc(sizeof(Alignment));
	if (!align) {
		fprintf(stderr, "Could not allocate memory for alignment\n");
		lastErrCode = E_MEM;
		return NULL;
	}

	memset(align, 0, sizeof(Alignment));

	align->seq1 = (char *) malloc(size);
	if (!align->seq1) {
		fprintf(stderr, "Could not allocate %d bytes of memory for"
			" alignment sequence 1\n", size);
		AlignmentDestroy(align);
		lastErrCode = E_MEM;
		return NULL;
	}
	memset(align->seq1, 0, size);

	align->seq2 = (char *) malloc(size);
	if (!align->seq2) {
		fprintf(stderr, "Could not allocate %d bytes of memory for"
			" alignment sequence 2\n", size);
		AlignmentDestroy(align);
		lastErrCode = E_MEM;
		return NULL;
	}
	memset(align->seq2, 0, size);

	align->homology = (char *) malloc(size);
	if (!align->homology) {
		fprintf(stderr, "Could not allocate %d bytes of memory for"
			" alignment homology sequence\n", size);
		AlignmentDestroy(align);
		lastErrCode = E_MEM;
		return NULL;
	}
	memset(align->homology, 0, size);

	align->size = size;
	align->length = 0;

	return align;
}

/*
 * Free memory allocated to an Alignment structure
 */
void
AlignmentDestroy(Alignment *align)
{
	if (!align)
		return;

	if (align->seq1)
		free(align->seq1);
	if (align->seq2)
		free(align->seq2);
	if (align->homology)
		free(align->homology);

	free(align);
}

/*
 * Add a sequence pair to the Alignment
 */
int
AlignmentAdd(Alignment *align, char seq1char, char seq2char, char homolChar)
{
	int	length;

	length = align->length;

	if (length < 0 || length >= (int) align->size)
		return 0;

	align->seq1[length] = seq1char;
	align->seq2[length] = seq2char;
	align->homology[length] = homolChar;

	++(align->length);

	return 1;
}

/*
 * Print an aligned pair of sequences and their corresponding alignment score
 */
void
AlignmentPrint(const Alignment *align, FILE *fp)
{
	static const int	charsPerLine = 60;
	register int		i, j;
	int			length;
	int			endPos;

	length = align->length;

	if (!align)
		return;

	if (!align->seq1 || !align->seq2)
		return;

	/*
	 * alignment is built in reverse so print from back to front
	 */
	for (j = length - 1; j >= 0; j -= charsPerLine) {
		endPos = j - charsPerLine + 1;
		if (endPos < 0)
			endPos = 0;

		for (i = j; i >= endPos; i--) {
			putc(align->seq1[i], fp);
		}
		putc('\n', fp);

		for (i = j; i >= endPos; i--) {
			putc(align->seq2[i], fp);
		}
		putc('\n', fp);
		
		for (i = j; i >= endPos; i--) {
			putc(align->homology[i], fp);
		}
		putc('\n', fp);
		putc('\n', fp);
	}

	/* fprintf(fp, "score = %d\n", align->score); */
}
