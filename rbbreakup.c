#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define COMBINE 1
#define BREAKUP 0
#define NBOUNDARIES 4
/* create a macro to calculate the index into the 1D array of a point with coordinates i, j */
#define Idx1D(row, col, rowlen) ((row) * rowlen + (col))
#define MAX(x, y) ((x)>(y) ? (x) : (y))

/* structure to hold arguments passed to program */

struct opt_struct {
	int operation;   /* COMBINE or BREAKUP.  Currently only BREAKUP is supported */
	char *inpfile;   /* input grid file.  defaults to sample.txt*/
	char *resfile;   /* base name for broken up files.  defaults to input file */
	int ibreak;      /* number of columnwise chunks */
	int jbreak;      /* number of rowwise chunks */
};

/* initialize argument values.  This routine should initialize everything in opt_struct  */
static void init_opts  (struct opt_struct *opts_p)  {
	static char default_filename[] = "sample.txt";

	opts_p->operation = BREAKUP;
	opts_p->inpfile = default_filename;
	opts_p->resfile = NULL;
	opts_p->ibreak = 1;
	opts_p->jbreak = 1;

	return;
}


static void write_usage(char *progname)
{
	char const usage[]=
		"Usage:  %s -ncores n [-inpfile filename] [-resfile filename]\n";
	char const opts_help[] =
		"\t-help           print this message\n"
		"\t-ichunk n       number of columnwise chunks to break\n"
		"\t-jchunk n       number of rowwise chunks to break\n"
		"\t-in filename     name of input file  to break or combine. if combine, files name filename.rank will be combined\n"
		"\t-out filename   name of results file to write if break, files will be written to filename.rank\n"
		;  /* don't lose this */

	fprintf( stdout, usage, progname );
	fputs( opts_help, stdout );

	return;
}


/* parse command line arguments. */
static int parse_args (int argc, char **argv, struct opt_struct *opts_p)
{
	register int i;

	init_opts(opts_p);

	for ( i=1; i< argc; ++i ) {

		if ( !strcmp("-help", argv[i] ) || !strcmp( "/?", argv[i] ) ) {      	/* print help msg */
			write_usage(argv[0]);
			exit (EXIT_SUCCESS);
		}

		else if ( !strcmp( "-in", argv[i] ) && (++i < argc))  {
			opts_p->inpfile = argv[i];
		}

		else if ( !strcmp( "-out", argv[i] ) && (++i < argc))  {
			opts_p->resfile = argv[i];
		}

		else if ( !strcmp( "-combine", argv[i] ))  {
			opts_p->operation = COMBINE;
		}

		else if ( !strcmp( "-ichunk", argv[i] ) && (++i < argc))  {
			opts_p->ibreak = strtol(argv[i], NULL, 10);
		}

		else if ( !strcmp( "-jchunk", argv[i] ) && (++i < argc))  {
			opts_p->jbreak = strtol(argv[i], NULL, 10);
		}

		else  {
			/* argument is unrecognized.  Barf. */
			fprintf(stdout, "Unrecognized argument \"%s\"\n", argv[i]);
			write_usage(argv[0]);
			return -1;
		}
	}

	/* if an output file name was not specified, use the input file name */
	if (opts_p->resfile == NULL)
		opts_p->resfile = opts_p->inpfile;

	return 0;
}


/* ReadGridHdr - read the header from a grid file with error checking */
static int ReadGridHdr(FILE *fp, int *rows_p, int *cols_p, int *niter_p, double *eps_p)
{
	if ( (fscanf(fp, "%d", rows_p) < 1) ||
		 (fscanf(fp, "%d", cols_p) < 1) ||
		 (fscanf(fp, "%lf", eps_p) < 1) ||
		 (fscanf(fp, "%d", niter_p) < 1) )
	{
		fprintf(stderr, "Error reading header of input file\n");
		return -1;
	}

	if (*rows_p <= 0 || *cols_p <= 0)  {

		fprintf(stderr, "Invalid row or column size specified in input file.  Must be greater than 0.\n");
		return -1;
	}

	if (*niter_p <= 0 )  {

		fprintf(stderr, "Invalid max iteration specified in input file.  Must be greater than 0.\n");
		return -1;
	}

	if (*eps_p < 0.0 )  {

		fprintf(stderr, "Invalid epsilon specified in input file.  Must be 0.0 or greater.\n");
		return -1;
	}

	return 0;
}


/* ReadGridData - Read gridpoint values from grid file */
static int ReadGridData (FILE *fp, int rows, int cols, double *g)
{
	register int r, c, n;

	for (n = 0, r = 0; r < rows; ++r)  {
		for (c = 0; c < cols; ++c)  {
			if (fscanf(fp, "%lf", g + n++) < 1)  {
				fprintf(stdout, "error reading grid data at (%d, %d)\n", r, c);
				return -1;
			}
		}
	}
	return 0;
}


/* LoadGridFile - read grid file into memory.  returns pointer to allocated array, and extents, etc. via arguments */
int LoadGridFile (char *filename, int *rows_p, int *cols_p, int *niter_p, double *eps_p, double **g_p)
{
	FILE *fp;
	double *g;

	fp = fopen(filename, "r");
	if (fp == NULL)  {
		fprintf(stdout, "Unable to open input file %s: %s\n", filename, strerror(errno));
		return -1;
	}

	if (ReadGridHdr(fp, rows_p, cols_p, niter_p, eps_p))
		return -1;

	/* allocate memory for grid.  The grid is stored in a 1D array, row by row */
	g = malloc(*cols_p * *rows_p * sizeof *g);
	if (g == NULL)  {
		fprintf (stdout, "Unable to allocate memory for %d X %d elements\n", *rows_p, *cols_p);
		(void) fclose(fp);
		return -1;
	}

	if (ReadGridData(fp, *rows_p, *cols_p, g))  {
		(void) free(g);
		(void) fclose(fp);
		return -1;
	}

	*g_p = g;  /* return pointer to allocated grid array to caller */

	(void) fclose(fp);
	return 0;
}


/* open_file - simple file open with error check */
static FILE *open_file( char *filename, const char *status)
{
	FILE *fp;

	fp = fopen(filename, status);
	if (fp == NULL)  {
		fprintf(stdout, "Unable to open file %s: %s", filename, strerror(errno));
	}

	return fp;
}


/* OpenBreakFile - open chunk file for given base name and rank */
FILE *OpenBreakFile(char *filebase, int rank)  {
	char *filepath;
	FILE *fp;

	/* allocate memory for temporary filename */
	filepath = malloc(strlen(filebase) + 12);
	if (filepath == NULL)  {
		fprintf (stderr, "Unable to allocate memory for temporary file name\n");
		return NULL;
	}

	/* open a file with the name filebase.N, where n is the rank */
    sprintf(filepath, "%s.%d", filebase, rank);

  	fp = open_file(filepath, "w");

  	(void) printf ("writing data for rank %d to %s\n", rank, filepath);
  	(void) fflush(stdout);

    free (filepath);
    return fp;
}


/* GetChunkExtents - calculate local extents of a chunk based on global extents and rank */
static void GetChunkExtents(int rank, int global_nrows, int global_ncols, int nrowchunks, int ncolchunks,
                            int *local_startrow_p, int *local_endrow_p, int *local_startcol_p, int *local_endcol_p)
{

	int mychunkrow, mychunkcol; /* row and column index of chunk, starting with 0 */

	mychunkrow = rank / ncolchunks;    /* reminder: integer division */
	mychunkcol = rank % ncolchunks;

	*local_startrow_p = mychunkrow * global_nrows / nrowchunks;
	*local_endrow_p = (mychunkrow + 1) * global_nrows / nrowchunks -1;

	*local_startcol_p = mychunkcol * global_ncols / ncolchunks;
	*local_endcol_p = (mychunkcol + 1) * global_ncols / ncolchunks -1;

	/* Unless we're the top row of chunks, add a ghost line to the north */
	if (mychunkrow > 0)  {
		*local_startrow_p -= 1;
	}

	/* Unless we're the easternmost row of chunks, add a ghost line to the east */
	if (mychunkcol < (ncolchunks-1))
		*local_endcol_p += 1;

	/* Unless we're the bottom row of chunks, add a ghost line to the south */
	if (mychunkrow < (nrowchunks-1))
		*local_endrow_p +=1;

	/* Unless we're the westernmost row of chunks, add a ghost line to the west */
	if (mychunkcol > 0)
		*local_startcol_p -= 1;

#ifdef DEBUG
	(void) printf("chunk for rank %d is %d-%d, %d-%d\n",
	              rank, *local_startcol_p, *local_endcol_p, *local_startrow_p, *local_endrow_p);
	fflush(stdout);
#endif

	return;
}


/* PrintSubgridToFile
   Print a subgrid specified by rank, IChunk, JChunk to a file.  The pointer can be stdout.
*/
static int PrintSubgridToFile(int rank, int global_ncols, int global_nrows, int IChunk, int JChunk,
	                          double eps, int niter, double *g, FILE *fp)
{
	register int r, c;
	int local_startrow, local_endrow, local_startcol, local_endcol;

	GetChunkExtents(rank, global_nrows, global_ncols, JChunk, IChunk,
                    &local_startrow, &local_endrow, &local_startcol, &local_endcol);

	if (((local_endrow-local_startrow) < 3) || ((local_endcol-local_startcol) < 3))  {
		fprintf(stderr, "Invalid breakup specified for this grid - resulting chunks must be at least 3 cells in all dimensions\n");
		fflush(stderr);
		return -1;
	}

	if (fprintf (fp, "%d\n%d\n%g\n%d\n", local_endcol-local_startcol+1, local_endrow-local_startrow+1, eps, niter) <= 0)  {
		fprintf(stderr, "Error writing header to chunk file for rank %d\n", rank);
		return -1;
	}

	for (r = local_startrow; r <= local_endrow; ++r)  {
		for (c = local_startcol; c <= local_endcol; ++c)  {

			if (fprintf (fp, "%7g ", g[Idx1D(r, c, global_ncols)] ) <= 0)  {
				fprintf(stderr, "Error writing data to chunk file for rank %d\n", rank);
				return -1;
			}
		}
		fputs ("\n", fp);  /* write newline at end of each row */
	}

	return 0;
}


/* BreakGridFile - for given grid file break into IChunk X JChunk chunks using original name as a base for output */
int BreakGridFile (char *infile, char *filebase, int IChunk, int JChunk)
{
	double *gg; /* global grid */
	int global_nrows, global_ncols, niter;
	double eps;
	register int rank;
	FILE *fp;

	if (LoadGridFile (infile, &global_ncols, &global_nrows, &niter, &eps, &gg))
		return -1;

	for (rank = 0; rank < (IChunk * JChunk); ++rank)  {

		if ((fp = OpenBreakFile(filebase, rank)) == NULL)
			break;

		if (PrintSubgridToFile(rank, global_ncols, global_nrows, IChunk, JChunk,
	                           eps, niter, gg, fp))
			break;

		(void) fclose(fp);
	}

	free(gg);

	return rank >= (IChunk * JChunk);  /* if we didn't do all the loop iterations, set bad status */
}


int main(int argc, char *argv[])
{
	static struct opt_struct opts;

	if (parse_args (argc, argv, &opts))
		return EXIT_FAILURE;

	(void) printf ("Red/Black Breakup - processing file %s\n", opts.inpfile);

	if (BreakGridFile (opts.inpfile, opts.resfile, opts.ibreak, opts.jbreak))
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
