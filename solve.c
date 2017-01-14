#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "cyclicreduction.h"

int main(int argc, char **argv)
{
  char *inputfile, *outfile;
  int i;
  int myid, nprocs;
  int n, ntot, logntot, rowsperproc;	// n is padded with redundant rows if n!=2^i, i>=0
  int lastrows, nofillprocs;		
  // 1 when 1st and last rows on proc are 1st and last lines in a tridiagonal matrix and will not be swapped out
  int hasfirst = 0, haslast = 0;		
  int prevrowproc, nextrowproc;		// processor storing row before/after 1st/last row on myid
  int offset = 1;			// doubles after each iteration; denotes the difference between the ranks of the processors in each pair
  // +1 after each iteration; denotes the difference between consecutive rows in each system on a processor (after all inter-proc row exchanges are done) 
  int rowoffset = 1;			
  int sendodds;			// procs are divided (on each iteration) into pairs which trade rows, with one sending even rows and the other odds
  int *myrownums;
  double lognprocs;
  RowInfo *myrowinfo;
  Xsol *xsol, *x_all;
  MPI_Datatype MPI_Rowinfo;
  MPI_Datatype MPI_Xsol;

  if (argc != 4) { fprintf(stderr, "usage: %s infile outfile nrows\n", argv[0]); exit(1); }
  inputfile = argv[1];
  outfile = argv[2];
  n = atoi(argv[3]);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  lognprocs = log2(nprocs);
  logntot = ceil(log2(n));
  ntot = exp2(logntot);
  rowsperproc = ntot / nprocs;
  if (lognprocs != ceil(lognprocs))	
  {
    fprintf(stderr, "nprocs must be 2^i, for some integer i in {0, 1, ..., 10}\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (lognprocs >= logntot)	
  {
    fprintf(stderr, "nprocs must be at most half of n\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  myrowinfo = (RowInfo *) malloc((rowsperproc + 2) * sizeof(RowInfo));
  myrownums = (int *) malloc((rowsperproc + 2) * sizeof(int));
  xsol = (Xsol *) malloc(rowsperproc * sizeof(RowInfo));
  rows_init(myrowinfo, myrownums, rowsperproc + 2, myid * rowsperproc);
  MPI_Rowinfo = create_rowinfotype();
  MPI_Xsol = create_xsoltype();

  if (log2(n) != logntot)
  {
    int tmpread;
    int tmpub = rowsperproc + 2;
    if (myid == nprocs - 1) tmpub = rowsperproc + 1;
    nofillprocs = n / rowsperproc;
    lastrows = n % rowsperproc;
    if (myid < nofillprocs - 1)
    {
      if (!readtridiag(inputfile, myrowinfo, rowsperproc, rowsperproc, 12, myid, nprocs))
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else if (myid == nofillprocs - 1)
    {
      if (!lastrows) tmpread = rowsperproc - 1;
      else tmpread = rowsperproc;
      if (!readtridiag(inputfile, myrowinfo, rowsperproc, tmpread, 12, myid, nprocs))
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else if (myid == nofillprocs)
    {
      tmpread = lastrows - 1;
      if (myid == nprocs - 1) tmpread = lastrows;
      if (!readtridiag(inputfile, myrowinfo, rowsperproc, tmpread, 12, myid, nprocs))
        MPI_Abort(MPI_COMM_WORLD, 1);
      myrowinfo[lastrows+1].a = 0.0;
      for (i = lastrows + 1; i < tmpub; ++i)
        myrowinfo[i].y = 0.0;
    }
    else
      for (i = 0; i < tmpub; ++i)
        myrowinfo[i].y = 0.0;
  }
  else
  {
    if (!readtridiag(inputfile, myrowinfo, rowsperproc, rowsperproc, 12, myid, nprocs))
      MPI_Abort(MPI_COMM_WORLD, 1);
  }

  double starttime = MPI_Wtime();
  // 1 and n are special cases  
  if (myid == 0) hasfirst = 1;
  if (myid == nprocs - 1) haslast = 1;
  if (myid != 0) prevrowproc = myid - 1;
  if (myid != nprocs - 1) nextrowproc = myid + 1;

  if (nprocs > 1)
  {
    for (i = 0; i < lognprocs; ++i)
    {
      if (update_rowinfo(myrowinfo, rowsperproc, hasfirst, haslast) < 0)
        MPI_Abort(MPI_COMM_WORLD, 1);
      sendodds = (myid / offset) % 2;	// the processor sending odd rows (receiving evens) is always the larger in each (exchanging) pair
      pack(myrowinfo, rowsperproc, 1, rowsperproc, sendodds);
      xchg_oddevens(myrowinfo, rowsperproc/2, 1, sendodds, myid, offset, MPI_Rowinfo);
      offset *= 2;
      update_rownums(myrownums, 1, rowsperproc, offset, sendodds);
      if (myrownums[1] == myid + 1) hasfirst = 1;
      if (myrownums[rowsperproc] == ntot - (nprocs - 1 - myid)) haslast = 1;
      xchg_bounds(myrowinfo, rowsperproc, myid, nprocs, offset, hasfirst, haslast, MPI_Rowinfo);
    }
    for (; i < logntot - 1; ++i, ++rowoffset)
      if (update_rowinfo_serial(myrowinfo, rowsperproc, rowoffset) < 0)
        MPI_Abort(MPI_COMM_WORLD, 1);
  }
  else
    for (i = 0; i < logntot - 1; ++i, ++rowoffset)
      if (update_rowinfo_serial(myrowinfo, rowsperproc, rowoffset) < 0)
        MPI_Abort(MPI_COMM_WORLD, 1);

  if (solve_2by2(&myrowinfo[1], xsol, rowsperproc) < 0)
    MPI_Abort(MPI_COMM_WORLD, 1);

  if (!myid) { printf("%f,", MPI_Wtime() - starttime); fflush(stdout); }
  if (!myid) x_all = (Xsol *) malloc(ntot * sizeof(Xsol));
  MPI_Gather(xsol, rowsperproc, MPI_Xsol, x_all, rowsperproc, MPI_Xsol, 0, MPI_COMM_WORLD);
  if (!myid) if (!write_x(x_all, outfile, ntot)) MPI_Abort(MPI_COMM_WORLD, 1);

  MPI_Finalize();
  return 0;
}


