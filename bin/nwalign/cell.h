/*
 * cell.h
 *
 * Defines the scoring matrix cell data type and associated functions
 *
 */
#ifndef	_cell_h_INCLUDED_
#define	_cell_h_INCLUDED_

typedef struct Cell {
	int		score;		/* cell score */
	/*
	 * The pointers below are set to indicate from which of the previous
	 * cell(s) the current cell's score was derived
	 */
	struct Cell	*cellUL;	/* upper-left cell */
	struct Cell	*cellU;		/* upper cell */
	struct Cell	*cellL;		/* left cell */
} Cell;

Cell *		CellCreate(void);
void		CellDestroy(Cell *);

#endif	/* _cell_h_INCLUDED_ */
