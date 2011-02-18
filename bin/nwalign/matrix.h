/*
 * matrix.h
 *
 * Defines the scoring matrix data type and associated functions
 *
 */
#ifndef	_matrix_h_INCLUDED_
#define	_matrix_h_INCLUDED_

#include "alignment.h"

/* most negative number representable in four bytes (machine independent?) */
#define MINUS_INFINITY	0x80000000
#define BIG_NEGATIVE	-999999

#define ALGORITHM 1	/* must be either 1 or 2 */

#if ALGORITHM==1
#define cell_t		unsigned short
#elif ALGORITHM==2
#define cell_t		unsigned char
#else
#error Algorithm must be 1 or 2!
#endif

/*
 * maximum allowable size in bytes of the matrix
 */
static const size_t			MaxMatrixSize = 256 * 1024 * 1024;

#if ALGORITHM==1
typedef enum {V, G, H}	WhichMatrix;
#endif

typedef struct Matrix {
	int		rows;	/* number of rows */
	int		cols;	/* number of columns */
	int		score;	/* overall matrix (alignment) score */
	cell_t		*cells;	/* cells of the matrix, bit field containing
				 * bit flags to previous cells from which
				 * current cell's score was derived */
#if ALGORITHM==1
	WhichMatrix	which;	/* defines which matrix resulted in the
				 * largest overall alignment score */
#endif
} Matrix;

/* some minimum possible score */
static const int	MinMatrixScore = MINUS_INFINITY;

Matrix *	MatrixCreate(int rows, int cols);
void		MatrixDestroy(Matrix *);
int		MatrixScore(Matrix *matrix, const char *seq1, const char *seq2,
			int matchScore, int mismatchScore, int gapPenalty,
			int gapExtPenalty);
Alignment *	MatrixBuildAlignment(const Matrix *scores, const char *seq1,
			const char *seq2);
void		MatrixTraverse(const Matrix *scores, const char *seq1,
			const char *seq2, Alignment *align); 
int		MatrixGetScore(const Matrix *matrix);
void		MatrixPrint(const Matrix *matrix, const char *seq1,
			const char *seq2);

#endif	/* _matrix_h_INCLUDED_ */
