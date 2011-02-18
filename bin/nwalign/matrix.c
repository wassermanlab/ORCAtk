#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "nwalign.h"
#include "matrix.h"

#define MAX2(a, b)	((a) > (b) ? (a) : (b))
#define MAX3(a, b, c)	MAX2(MAX2((a), (b)), (c))

extern int	lastErrCode;

/*
 * Bit masks for the flags indicating that the current cell's maximum score
 * was derived from the upper, left, and/or upper-left cell's score. For
 * algorithm 1 we have to keep track of three different matrices and which
 * previous matrix score each of those was derived from, i.e. 9 possible
 * traceback pointers.
 */
#if ALGORITHM==1
static const cell_t VfromV = 0x1;
static const cell_t VfromG = 0x2;
static const cell_t VfromH = 0x4;
static const cell_t GfromV = 0x8;
static const cell_t GfromG = 0x10;
static const cell_t GfromH = 0x20;
static const cell_t HfromV = 0x40;
static const cell_t HfromG = 0x80;
static const cell_t HfromH = 0x100;
#elif ALGORITHM==2
static const cell_t UMask = 0x1;
static const cell_t LMask = 0x2;
static const cell_t ULMask = 0x4;
#endif

inline static int	nucMatch(register char nuc1, register char nuc2);
/*
 * for debugging only
 * static void		printMatrixRow(int *scoreV, int *scoreG, int *scoreH,
 *				    cell_t *cell, int row, int cols);
 */

/*
 * Allocate and initialize memory for a Matrix.
 */
Matrix *
MatrixCreate(int n, int m)
{
    size_t		matrixSize;
    Matrix		*matrix = NULL;
    
    if (n < 1 || m < 1) {
	lastErrCode = E_BADARGS;
	return NULL;
    }

    /*
     * Check to make sure we don't exceed some maximum size limits
     */
    matrixSize = n * m * sizeof(cell_t);
    if (matrixSize > MaxMatrixSize) {
	fprintf(stderr, "Required matrix size of %u bytes exceeds"
		" maximum allowed %u bytes\n", matrixSize, MaxMatrixSize);
	lastErrCode = E_SEQLEN;
	return NULL;
    }

    matrix = (Matrix *) malloc(sizeof(Matrix));
    if (!matrix) {
	fprintf(stderr, "Could not allocate memory for scoring matrix\n");
	lastErrCode = E_MEM;
	return NULL;
    }
    memset(matrix, 0, sizeof(Matrix));

    matrix->cells = (cell_t *) malloc(matrixSize);
    if (!matrix->cells) {
	fprintf(stderr, "Could not allocate %u bytes of memory for"
		" scoring matrix\n", matrixSize);
	MatrixDestroy(matrix);
	lastErrCode = E_MEM;
	return NULL;
    }
    memset(matrix->cells, 0, matrixSize);

    matrix->cols = n;
    matrix->rows = m;

    return matrix;
}

/*
 * Free memory allocated to a matrix
 */
void
MatrixDestroy(Matrix *matrix)
{
    if (!matrix)
	return;

    if (matrix->cells)
	free(matrix->cells);

    free(matrix);
}

#if ALGORITHM==1
/*
 * Fill in matrix scores as specified by the Needleman-Wunsch algorithm.
 * The matrix is scored from top-left to bottom-right. For NW with affine
 * gaps there seem to be two different recursive algorithms in use.
 *
 * This routine implements algorithm 1.
 *
 * Algorithm 1 (e.g. Durbin book):
 *
 *		      | V(i-1, j-1) + s(xi, yj)
 *	V(i, j) = MAX | G(i-1, j-1) + s(xi, yj)
 *		      | H(i-1, j-1) + s(xi, yj)
 *
 *		      | V(i-1, j) - d
 *	G(i, j) = MAX | G(i-1, j) - e
 *		      | H(i-1, j) - d (*)
 *
 *		      | V(i, j-1) - d
 *	H(i, j) = MAX | G(i, j-1) - d (*)
 *		      | H(i, j-1) - e
 *
 * where:
 *	d		= gap open penalty
 *	e		= gap extension penalty
 *	s(xi, yj) 	= match score, if seq1(i) = seq2(j),
 *		  	  mismatch score, if seq1(i) != seq2(i)
 *	
 * NOTES:
 *	(*) The Durbin book eliminates these terms on the assumption that
 *	an insertion is not immediately followed by a deletion.
 *	
 * It is necessary to keep track of traceback pointers for each of the
 * matrices V, G and H. It is possible to jump between matrices during
 * traceback if for example we are currently tracing back in the V matrix
 * and V(i,J) was derived from G(i-1,j-1) then we start tracing back from
 * the G matrix at (i-1, j-1).
 * See http://www.dina.dk/~sestoft/bsa/graphalign.html for a graphical
 * example of this.
 *
 */
int
MatrixScore(Matrix *matrix, const char *seq1, const char *seq2,
		int match, int mismatch, int gapOpen, int gapExt)
{
    register int	i, j;		/* try to max. speed */
    int			idx;
    int			rows, cols;
    int			s;
    int			F;
    int			scoreV, scoreG, scoreH;
    int			*curV, *prevV, *curG, *prevG, *curH, *prevH, *temp;
    cell_t		*cell;
    WhichMatrix		which;

    if (!matrix)
	return 0;

    if (!matrix->cells)
	return 0;

    if (!seq1 || !seq2)
	return 0;

    cols = matrix->cols;
    cell = matrix->cells;
    rows = matrix->rows;

    curV = (int *) calloc(cols, sizeof(int));
    if (!curV) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix V current row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(curV, 0, cols * sizeof(int));

    prevV = (int *) calloc(cols, sizeof(int));
    if (!prevV) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix V previous row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(prevV, 0, cols * sizeof(int));

    curG = (int *) calloc(cols, sizeof(int));
    if (!curG) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix G current row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(curG, 0, cols * sizeof(int));

    prevG = (int *) calloc(cols, sizeof(int));
    if (!prevG) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix G previous row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(prevG, 0, cols * sizeof(int));

    curH = (int *) calloc(cols, sizeof(int));
    if (!curH) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix H current row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(curH, 0, cols * sizeof(int));

    prevH = (int *) calloc(cols, sizeof(int));
    if (!prevH) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix H previous row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(prevH, 0, cols * sizeof(int));

    /*
     * Initialize upper left. Was using MINUS_INFINITY instead of BIG_NEGATIVE
     * but the problem is when you subtract from it, the number rolls over and
     * produces a huge positive number which is not the desired affect!
     */
    curV[0] = 0;
    curG[0] = BIG_NEGATIVE;
    curH[0] = BIG_NEGATIVE;
    cell[0] = 0;
    F = 0,
    which = V;

    /* first row initialization */
    for (i = 1; i < cols; i++) {
	curV[i] = BIG_NEGATIVE;
	curG[i] = -(gapOpen + gapExt * (i - 1));
	curH[i] = BIG_NEGATIVE;
	cell[i] = GfromG;
	F = curG[i];
	which = G;
    }
    /*printMatrixRow(curV, curG, curH, cell, 0, cols);*/

    for (j = 1; j < rows; j++) {
	temp = prevV;
	prevV = curV;
	curV = temp;
	memset(curV, 0, cols);

	temp = prevG;
	prevG = curG;
	curG = temp;
	memset(curG, 0, cols);

	temp = prevH;
	prevH = curH;
	curH = temp;
	memset(curH, 0, cols);

	/* first column initialization */
	curV[0] = BIG_NEGATIVE;
	curG[0] = BIG_NEGATIVE;
	curH[0] = -(gapOpen + gapExt * (j - 1));
	cell[j * cols] = HfromH;
	F = curH[0];
	which = H;

	for (i = 1; i < cols; i++) {
	    s = nucMatch(seq1[i - 1], seq2[j - 1]) ? match : mismatch;

	    scoreV = MAX3(prevV[i-1], prevG[i-1], prevH[i-1]) + s;
	    scoreG = MAX3(curV[i-1] - gapOpen, curG[i-1] - gapExt,
			    curH[i-1] - gapOpen);
	    scoreH = MAX3(prevV[i] - gapOpen, prevG[i] - gapOpen,
			    prevH[i] - gapExt);

	    /*
	     * Set overall score of alignment and which matrix resulted in
	     * the highest overall alignment score up to this point
	     */
	    F = MAX3(scoreV, scoreG, scoreH); 
	    if (F == scoreV)
		which = V;
	    else if (F == scoreG)
		which = G;
	    else if (F == scoreH)
		which = H;

	    /* set current scores for the V, G and H matrix */
	    curV[i] = scoreV;
	    curG[i] = scoreG;
	    curH[i] = scoreH;

	    /* index of current cell of traceback matrix */
	    idx = j * cols + i;

	    /* set traceback pointers for each of the matrices */
	    if (scoreV == prevV[i-1] + s)
	        cell[idx] |= VfromV;
	    if (scoreV == prevG[i-1] + s)
	        cell[idx] |= VfromG;
	    if (scoreV == prevH[i-1] + s)
	        cell[idx] |= VfromH;

	    if (scoreG == curV[i-1] - gapOpen)
	        cell[idx] |= GfromV;
	    if (scoreG == curG[i-1] - gapExt)
	        cell[idx] |= GfromG;
	    if (scoreG == curH[i-1] - gapOpen)
	        cell[idx] |= GfromH;

	    if (scoreH == prevV[i] - gapOpen)
	        cell[idx] |= HfromV;
	    if (scoreH == prevG[i] - gapOpen)
	        cell[idx] |= HfromG;
	    if (scoreH == prevH[i] - gapExt)
	        cell[idx] |= HfromH;
	}

	/*printMatrixRow(curV, curG, curH, cell, j, cols);*/
    }

    /*
     * Set overall alignment score and which matrix resulted in the highest
     * overall score.
     */
    matrix->score = F;
    matrix->which = which;

    free(prevV);
    free(curV);
    free(prevG);
    free(curG);
    free(prevH);
    free(curH);
    
    return 1;
}

#elif ALGORITHM==2
/*
 * Fill in matrix scores as specified by the Needleman-Wunsch algorithm.
 * The matrix is scored from top-left to bottom-right. For NW with affine
 * gaps there seem to be two different recursive algorithms in use.
 *
 * This routine implements algorithm 2.
 *
 * Algorithm 2 (e.g. Gusfield book):
 *
 *		      | F(i, j)
 *	V(i, j) = MAX | G(i, j)
 *		      | H(i, j)
 *
 *	F(i, j) = V(i-1, j-1) + s(xi, yj) 
 *
 *	G(i, j) = MAX | V(i-1, j) - d
 *		      | G(i-1, j) - e
 *
 *	H(i, j) = MAX | V(i, j-1) - d
 *		      | H(i, j-1) - e;
 *
 * where:
 *	d		= gap open penalty
 *	e		= gap extension penalty
 *	s(xi, yj) 	= match score, if seq1(i) = seq2(j),
 *		  	  mismatch score, if seq1(i) != seq2(i)
 *
 * The reference materials state that V(i,j) is the optimal alignment. Although
 * it is never explicitly stated, this would seem to imply that it is only
 * necessary to keep track of traceback pointers for V. However, when I tested
 * this algorithm it did not work correctly.
 */
int
MatrixScore(Matrix *matrix, const char *seq1, const char *seq2,
		int match, int mismatch, int gapOpen, int gapExt)
{
    register int	i, j;		/* try to max. speed */
    int			idx;
    int			rows, cols;
    int			s;
    int			F;
    int			scoreV, scoreG, scoreH;
    int			*curV, *prevV, *curG, *prevG, *curH, *prevH, *temp;
    cell_t		*cell;

    if (!matrix)
	return 0;

    if (!matrix->cells)
	return 0;

    if (!seq1 || !seq2)
	return 0;

    cols = matrix->cols;
    cell = matrix->cells;
    rows = matrix->rows;

    curV = (int *) calloc(cols, sizeof(int));
    if (!curV) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix V current row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(curV, 0, cols * sizeof(int));

    prevV = (int *) calloc(cols, sizeof(int));
    if (!prevV) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix V previous row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(prevV, 0, cols * sizeof(int));

    curG = (int *) calloc(cols, sizeof(int));
    if (!curG) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix G current row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(curG, 0, cols * sizeof(int));

    prevG = (int *) calloc(cols, sizeof(int));
    if (!prevG) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix G previous row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(prevG, 0, cols * sizeof(int));

    curH = (int *) calloc(cols, sizeof(int));
    if (!curH) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix H current row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(curH, 0, cols * sizeof(int));

    prevH = (int *) calloc(cols, sizeof(int));
    if (!prevH) {
	fprintf(stderr, "Could not allocate %d bytes of memory for"
		" matrix H previous row\n", cols * sizeof(int));
	lastErrCode = E_MEM;
	return 0;
    }
    memset(prevH, 0, cols * sizeof(int));

    /*
     * Initialize upper left. Was using MINUS_INFINITY instead of BIG_NEGATIVE
     * but the problem is when you subtract from it, the number rolls over and
     * produces a huge positive number which is not the desired affect!
     */
    curV[0] = 0;
    curG[0] = BIG_NEGATIVE;
    curH[0] = BIG_NEGATIVE;
    cell[0] = 0;
    scoreV = curV[0];

    /* first row initialization */
    for (i = 1; i < cols; i++) {
	curV[i] = curG[i] = -(gapOpen + gapExt * (i - 1));
	curH[i] = BIG_NEGATIVE;
	cell[i] = LMask;
	scoreV = curV[i];
    }

    for (j = 1; j < rows; j++) {
	temp = prevV;
	prevV = curV;
	curV = temp;
	memset(curV, 0, cols);

	temp = prevG;
	prevG = curG;
	curG = temp;
	memset(curG, 0, cols);

	temp = prevH;
	prevH = curH;
	curH = temp;
	memset(curH, 0, cols);

	/* first column initialization */
	curV[0] = curH[0] = -(gapOpen + gapExt * (j - 1));
	curG[0] = BIG_NEGATIVE;
	cell[j * cols] = UMask;
	scoreV = curV[0];

	for (i = 1; i < cols; i++) {
	    s = nucMatch(seq1[i - 1], seq2[j - 1]) ? match : mismatch;

	    F = prevV[i-1] + s;
	    scoreG = MAX2(curV[i-1] - gapOpen, curG[i-1] - gapExt);
	    scoreH = MAX2(prevV[i] - gapOpen, prevH[i] - gapExt);

	    /* set overall score of alignment and which matrix resulted in the
	       highest overall alignment score */
	    scoreV = MAX3(F, scoreG, scoreH); 

	    /* set current scores for the V, G and H matrix */
	    curV[i] = scoreV;
	    curG[i] = scoreG;
	    curH[i] = scoreH;

	    /* index of current cell of traceback matrix */
	    idx = j * cols + i;

	    /*
	     * set pointer(s) back to the cell(s) from which the highest
	     * score was derived
	     */
	    if (scoreV == F)
	    	cell[idx] = ULMask;
	    if (scoreV == scoreG)
		cell[idx] |= LMask;
	    if (scoreV == scoreH)
		cell[idx] |= UMask;
	}
    }

    matrix->score = scoreV;

    free(prevV);
    free(curV);
    free(prevG);
    free(curG);
    free(prevH);
    free(curH);
    
    return 1;
}
#endif

/*
 * Build up an alignment by traversing the scoring matrix
 */
Alignment *
MatrixBuildAlignment(const Matrix *scores, const char *seq1, const char *seq2)
{
    int			rows, cols;
    size_t		alignSize;
    Alignment		*align;

    rows = scores->rows;
    cols = scores->cols;

    /*
     * Maximum possible size for total non-alignment.
     * What is a more realistic max. poss. size???
     */
    alignSize = rows + cols;
    align = AlignmentCreate(alignSize);
    if (!align)
	return NULL;

    MatrixTraverse(scores, seq1, seq2, align);

    align->score = MatrixGetScore(scores);

    return align;
}

#if ALGORITHM==1
/*
 * Iteratively traverse scoring matrix building up alignment structure.
 *
 * Note: this routine only traces back a single path, not all possible
 * paths.
 */
void
MatrixTraverse(const Matrix *scores, const char *seq1, const char *seq2,
		Alignment *align)
{
    register int		row, col, idx;
    register cell_t		*cell;
    char			homolChar;
    int				rows, cols;
    WhichMatrix			which;

    rows = scores->rows;
    cols = scores->cols;
    cell = scores->cells;
    which = scores->which;

    /*
     * Start with the lower-right most cell which by definition,
     * corresponds to the highest possible alignment score
     */
    idx = rows * cols - 1;
    row = rows - 1;
    col = cols - 1;
    while (idx) {
	switch (which) {
	case V:
	    if (cell[idx] & VfromV)
		which = V;
	    else if (cell[idx] & VfromG)
		which = G;
	    else if (cell[idx] & VfromH)
		which = H;

	    homolChar = nucMatch(seq1[col-1], seq2[row-1]) ? HomolChar : ' ';
	    AlignmentAdd(align, seq1[col - 1], seq2[row - 1], homolChar);
	    row--;
	    col--;
	    idx -= (cols + 1);
	    break;
	case G:
	    if (cell[idx] & GfromV)
		which = V;
	    else if (cell[idx] & GfromG)
		which = G;
	    else if (cell[idx] & GfromH)
		which = H;

	    AlignmentAdd(align, seq1[col - 1], '-', ' ');
	    col--;
	    idx--;
	    break;
	case H:
	    if (cell[idx] & HfromV)
		which = V;
	    else if (cell[idx] & HfromG)
		which = G;
	    else if (cell[idx] & HfromH)
		which = H;

	    AlignmentAdd(align, '-', seq2[row - 1], ' ');
	    row--;
	    idx -= cols;
	    break;
	default:
	    fprintf(stderr, "Error determining which matrix in traceback\n");
	    exit(-1);
	    break;
	}
    }
}

#elif ALGORITHM==2
/*
 * Iteratively traverse scoring matrix building up alignment structure.
 *
 * Note: this routine only traces back a single path, not all possible
 * paths.
 */
void
MatrixTraverse(const Matrix *scores, const char *seq1, const char *seq2,
		Alignment *align)
{
    register int		row, col, idx;
    register cell_t		*cell;
    char			homolChar;
    int				rows, cols;

    rows = scores->rows;
    cols = scores->cols;
    cell = scores->cells;

    /*
     * Start with the lower-right most cell which by definition,
     * corresponds to the highest possible alignment score
     */
    idx = rows * cols - 1;
    row = rows - 1;
    col = cols - 1;

    while (idx) {
	if (cell[idx] & ULMask) {
	    homolChar = nucMatch(seq1[col-1], seq2[row-1])
			    ? HomolChar : ' ';
	    AlignmentAdd(align, seq1[col - 1], seq2[row - 1], homolChar);
	    idx -= (cols + 1);
	    row--;
	    col--;
	} else if (cell[idx] & UMask) {
	    AlignmentAdd(align, '-', seq2[row - 1], ' ');
	    idx -= cols;
	    row--;
	} else if (cell[idx] & LMask) {
	    AlignmentAdd(align, seq1[col - 1], '-', ' ');
	    idx--;
	    col--;
	}
    }
}
#endif

int
MatrixGetScore(const Matrix *matrix)
{
    if (!matrix)
	return MinMatrixScore;

    return matrix->score;
}

/*
 * Simple nucleotide match check.
 *
 * ASSUMPTION: an N or X is considered to be a masking character and is
 *	automatically counted as a mismatch. (Unfortunately N is also the
 *	single character representation for the amino acid Asparagine, so
 *	this will not work so well for amino acid sequences)
 *
 * NOTE: Returns either true (1) or false (0). Does not compute a weighted
 * score based on degree of match or mismatch (i.e. for amino acids).
 */
inline static int
nucMatch(register char nuc1, register char nuc2)
{
    register char	uNuc1 = toupper(nuc1);
    register char	uNuc2 = toupper(nuc2);

    if (uNuc1 == uNuc2 && uNuc1 != 'N' && uNuc2 != 'N' && uNuc1 != 'X'
	    && uNuc2 != 'X')
	return 1;

    return 0;
}

/*
static void
printMatrixRow(int *scoreV, int *scoreG, int *scoreH, cell_t *cell,
		int row, int cols)
{
    char		fromStr[4];
    int			i, idx;

    for (i = 0; i < cols; i++) {
        idx = row * cols + i;

	if (scoreH[i] <= BIG_NEGATIVE) {
	    printf("        %4s", "-Inf");
	} else {
	    printf("        %4d", scoreH[i]);
	}

	strcpy(fromStr, "   ");
	if (cell[idx] & HfromV)
	    fromStr[0] = 'V';
	if (cell[idx] & HfromG)
	    fromStr[1] = 'G';
	if (cell[idx] & HfromH)
	    fromStr[2] = 'H';
	printf("%3s|", fromStr);
    }
    printf("\n");

    for (i = 0; i < cols; i++) {
        idx = row * cols + i;

	if (scoreG[i] <= BIG_NEGATIVE) {
	    printf("%4s", "-Inf");
	} else {
	    printf("%4d", scoreG[i]);
	}

	strcpy(fromStr, "   ");
	if (cell[idx] & GfromV)
	    fromStr[0] = 'V';
	if (cell[idx] & GfromG)
	    fromStr[1] = 'G';
	if (cell[idx] & GfromH)
	    fromStr[2] = 'H';
	printf("%3s", fromStr);

	if (scoreV[i] <= BIG_NEGATIVE) {
	    printf(" %4s", "-Inf");
	} else {
	    printf(" %4d", scoreV[i]);
	}

	strcpy(fromStr, "   ");
	if (cell[idx] & VfromV)
	    fromStr[0] = 'V';
	if (cell[idx] & VfromG)
	    fromStr[1] = 'G';
	if (cell[idx] & VfromH)
	    fromStr[2] = 'H';
	printf("%3s|", fromStr);
    }
    printf("\n");
    for (i = 0; i < cols; i++) {
	printf("----------------");
    }
    printf("\n");
}
*/
