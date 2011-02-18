#include <stdlib.h>
#include <string.h>

#include "cell.h"

/*
 * Create and initialize a new cell structure
 */
Cell *
CellCreate()
{
	Cell 		*cell;

	cell = (Cell *) malloc(sizeof(Cell));

	if (!cell)
		return NULL;

	memset(cell, 0, sizeof(Cell));

	return cell;
}

/*
 * Free memory allocated to a cell structure 
 */
void
CellDestroy(Cell *cell)
{
	if (cell)
		free(cell);
}
