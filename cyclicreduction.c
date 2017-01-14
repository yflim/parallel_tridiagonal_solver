#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "cyclicreduction.h"

#define MASTERID 0

void rows_init(RowInfo *rinfo, int *rnums, int nrl, int firstrow)
{
  int i, j;
  for (i = 0, j = firstrow; i < nrl; ++i, ++j)
  {
    rnums[i] = rinfo[i].row = j;
    rinfo[i].a = 1.0; rinfo[i].b = 5.0; rinfo[i].c = 1.0; 
    rinfo[i].y = 0.0;
  }
}

void copy_rowinfo(RowInfo *src, RowInfo *dest)
{
  dest->row = src->row;
  dest->a = src->a; dest->b = src->b; dest->c = src->c; 
  dest->y = src->y;
}

void copy_rowinfo_array(RowInfo *src, RowInfo *dest, int nrows)
{
  int i;
  for (i = 0; i < nrows; ++i)
  {
    dest[i].row = src[i].row;
    dest[i].a = src[i].a; dest[i].b = src[i].b; dest[i].c = src[i].c; 
    dest[i].y = src[i].y;
  }
}

// print a specified range of entries:2D
void print_rowinfo(RowInfo *rinfo, int nrl, int myid)
{
  int i;
  for (i = 0; i < nrl; ++i)
    printf("proc %d row %d: a b c y %f %f %f %f\n", myid, rinfo[i].row, rinfo[i].a, rinfo[i].b, rinfo[i].c, rinfo[i].y);
}

void print_rownums(int *rnums, int nrl, int myid)
{
  int i;
  printf("proc %d: ", myid);
  for (i = 0; i < nrl; ++i) printf("%d ", rnums[i]);
  printf("\n");
}

MPI_Datatype create_rowinfotype()
{
  int blocklengths[2]={4,1};
  MPI_Datatype types[2]={MPI_DOUBLE, MPI_INT};
  MPI_Aint disp[2];
  MPI_Datatype RowInfoType;
  MPI_Aint int_ex, double_ex;
  MPI_Type_extent(MPI_DOUBLE, &double_ex);
  MPI_Type_extent(MPI_INT, &int_ex);
  disp[0] = (MPI_Aint) 0; disp[1] = 4 * double_ex;
  MPI_Type_struct(2, blocklengths, disp, types, &RowInfoType);
  MPI_Type_commit(&RowInfoType);
  return RowInfoType;
}

MPI_Datatype create_xsoltype()
{
  int blocklengths[2]={1,1};
  MPI_Datatype types[2]={MPI_DOUBLE, MPI_INT};
  MPI_Aint disp[2];
  MPI_Datatype XsolType;
  MPI_Aint int_ex, double_ex;
  MPI_Type_extent(MPI_DOUBLE, &double_ex);
  MPI_Type_extent(MPI_INT, &int_ex);
  disp[0] = (MPI_Aint) 0; disp[1] = double_ex;
  MPI_Type_struct(2, blocklengths, disp, types, &XsolType);
  MPI_Type_commit(&XsolType);
  return XsolType;
}

RowInfo *readtridiag(char *initfile, RowInfo *rinfo, int nrl, int mynrows, int toklen, int myid, int nprocs)
{
  int i, j;
  int lb = 0, ub = mynrows + 2;
  FILE *file;
  int colchars = 4 * toklen + 3;
  char *s = (char *) malloc(colchars * sizeof(char));
  char **num = (char **) calloc(4, (toklen + 1) * sizeof(char));
  char *delim = " ";

  if (!(file = fopen(initfile, "r")))
  {
    fprintf(stderr, "failed to open %s for read.\n", initfile);
    return NULL;
  }

  for (i = 0; i < nrl * myid - 1; ++i)
    if (!fgets(s, colchars, file))
    {
      fprintf(stderr, "error reading row %d in %s on proc %d. nrows given may be incorrect.\n", i, initfile, myid);
      return NULL;
    }

  if (myid == 0) lb = 1;
  if (myid == nprocs - 1) ub = mynrows + 1;
  for (i = lb; i < ub; ++i)
  {
    if (!fgets(s, colchars, file))
    {
      fprintf(stderr, "error reading row %d in %s on proc %d. nrows given may be incorrect.\n", i, initfile, myid);
      return NULL;
    }
    if (!(num[0] = strtok(s, delim)))
    {
      fprintf(stderr, "error parsing %s with strtok\n", initfile);
      return NULL;
    }
    for (j = 1; j < 4; ++j)
      if (!(num[j] = strtok(NULL, delim)))
      {
        fprintf(stderr, "error parsing %s with strtok\n", initfile);
        return NULL;
      }
    rinfo[i].a = atof(num[0]);	rinfo[i].b = atof(num[1]);	
    rinfo[i].c = atof(num[2]);	rinfo[i].y = atof(num[3]);	
  }

  fclose(file);
  return rinfo;
}

void update_rows1and2(RowInfo *r1, RowInfo *r2, RowInfo *dest)
{
  double gamma = (-1) * r1->c / r2->b;
  dest->row = r1->row;
  dest->b = r1->b + gamma * r2->a;
  dest->b = r1->b + gamma * r2->a;
  dest->c = gamma * r2->c;
  dest->y = r1->y + gamma * r2->y;
}

void update_rows(RowInfo *r1, RowInfo *r2, RowInfo *r3, RowInfo *dest)
{
  double alpha = (-1) * r2->a / r1->b;
  double gamma = (-1) * r2->c / r3->b;
  dest->row = r2->row;
  dest->a = alpha * r1->a;
  dest->b = alpha * r1->c + r2->b + gamma * r3->a;
  dest->c = gamma * r3->c;
  dest->y = alpha * r1->y + r2->y + gamma * r3->y;
}

void update_last2rows(RowInfo *r1, RowInfo *r2, RowInfo *dest)
{
  double alpha = (-1) * r2->a / r1->b;
  dest->row = r2->row;
  dest->a = alpha * r1->a;
  dest->b = r2->b + alpha * r1->c;
  dest->y = r2->y + alpha * r1->y;
}

int update_rowinfo(RowInfo *myrows, int nrows, int hasfirst, int haslast)
{
  int i, j;
  int lb = 1, ub = nrows;
  RowInfo newrows[nrows];

  if (hasfirst)
  {
    lb = 2;
    if (!myrows[1].b)
    {
      fprintf(stderr, "b in row 1 is 0: x2 is already determined. please reformulate system without x2.\n");
      return -1;
    }
    if (!myrows[2].b)
    {
      fprintf(stderr, "zero on diagonal (row 2): row information update for cyclic reduction fails.\n");
      return -1;
    }
    update_rows1and2(&myrows[1], &myrows[2], &newrows[0]);
  }
  if (haslast) ub = nrows - 1;
  for (i = lb, j = lb - 1; i <= ub; ++i, ++j)
  {
    if (!myrows[i+1].b)
    {
      fprintf(stderr, "zero on diagonal: row information update for cyclic reduction fails.\n");
      if (i == nrows - 1) fprintf(stderr, "b in last row is 0: x%d is already determined. please reformulate system without x%d.\n", nrows - 1, nrows - 1);
      return -1;
    }
    update_rows(&myrows[i-1], &myrows[i], &myrows[i+1], &newrows[j]);
  }
  if (haslast)
    update_last2rows(&myrows[nrows-1], &myrows[nrows], &newrows[nrows-1]);

  copy_rowinfo_array(newrows, &myrows[1], nrows);
  return nrows;
}

int update_rowinfo_serial(RowInfo *myrows, int nrows, int rowoffset)
{
  int i;
  RowInfo newrows[nrows];

  for (i = 1; i <= nrows; ++i)
  {
    if (i + rowoffset <= nrows)
    {
      if (!myrows[i+rowoffset].b)
      {
        fprintf(stderr, "zero on diagonal (row %d): row information update for cyclic reduction fails.\n", i + rowoffset);
        return -1;
      }
      if (i - rowoffset < 1)
      {
        if (!myrows[i].b)
        {
          fprintf(stderr, "zero on diagonal (row %d): row information update for cyclic reduction fails.\n", i);
          return -1;
        }
        update_rows1and2(&myrows[i], &myrows[i+rowoffset], &newrows[i-1]);
      }
      else update_rows(&myrows[i-rowoffset], &myrows[i], &myrows[i+rowoffset], &newrows[i-1]);
    }
    else 
      update_last2rows(&myrows[i-rowoffset], &myrows[i], &newrows[i-1]);
  }

  copy_rowinfo_array(newrows, &myrows[1], nrows);
  return nrows;
}

// assumes nrows is even
void pack(RowInfo *rinfo, int nrows, int first, int last, int sendodds)
{
  int i, j, k;
  int half = nrows / 2;
  RowInfo tmprow[half];
  if (sendodds)
  {
    for (i = j = last - 1, k = half - 1; j >= first; --i, j-=2, --k)
    {
      copy_rowinfo(&rinfo[j], &tmprow[k]);
      copy_rowinfo(&rinfo[j-1], &rinfo[i]);
    }
    copy_rowinfo_array(tmprow, &rinfo[first], half);
    //for (i = last / 2, k = half - 1; k >= 0; --i, --k) copy_rowinfo(&tmprow[k], &rinfo[i]);
  }
  else
  {
    for (i = j = first + 1, k = 0; j <= last; ++i, j+=2, ++k)
    {
      copy_rowinfo(&rinfo[j], &tmprow[k]);
      copy_rowinfo(&rinfo[j+1], &rinfo[i]);
    }
    copy_rowinfo_array(tmprow, &rinfo[last/2+1], half);
    //for (i = last / 2 + 1, k = 0; k < half; ++i, ++k) copy_rowinfo(&tmprow[k], &rinfo[i]);
  }
}

void xchg_oddevens(RowInfo *rinfo, int xchgrows, int firstidx, int sendodds, int myid, int offset, MPI_Datatype MPI_Rowinfo)
{
  MPI_Request recvreq, sendreq;
  MPI_Status recvstat, sendstat;
  int recvflag, sendflag;
  recvflag = sendflag = 0;
  int peer; 
  RowInfo recvbuf[xchgrows];

  if (sendodds)
  {
    // proc sending odd rows (receiving evens) always the larger in each (exchanging) pair
    peer = myid - offset;	
    /* odd rows (e.g. 1, 3, ...) are collected in the first half of the processor's array of RowInfo structs, 
    then sent to its peer where it forms the second half of the peer's RowInfo array. Symmetrically, the peer 
    collects even rows in the second half of its array, then sends them to this processor, where they are placed 
    in the first half, where the collected odd rows used to be.							*/
    MPI_Isend(&rinfo[firstidx], xchgrows, MPI_Rowinfo, peer, 0, MPI_COMM_WORLD, &sendreq);
    MPI_Irecv(recvbuf, xchgrows, MPI_Rowinfo, peer, 0, MPI_COMM_WORLD, &recvreq);
  }
  else
  {
    peer = myid + offset;	
    // what this processor does (sending even rows) is symmetric with what its peer (sending odd rows) does
    MPI_Isend(&rinfo[firstidx+xchgrows], xchgrows, MPI_Rowinfo, peer, 0, MPI_COMM_WORLD, &sendreq);
    MPI_Irecv(recvbuf, xchgrows, MPI_Rowinfo, peer, 0, MPI_COMM_WORLD, &recvreq);
  }

  while (!recvflag || !sendflag)
  {
    MPI_Test(&recvreq, &recvflag, &recvstat);
    MPI_Test(&sendreq, &sendflag, &sendstat);
  }
  if (sendodds)
    copy_rowinfo_array(recvbuf, &rinfo[firstidx], xchgrows);
  else
    copy_rowinfo_array(recvbuf, &rinfo[firstidx+xchgrows], xchgrows);
}

void update_rownums(int *rnums, int firstidx, int totrows, int offset, int sendodds)
{
  int i, j;
  if (sendodds)
  {
    i = totrows - 2; j = rnums[i]; rnums[i+1] = j + offset;
    for (; i >= 0; --i, j -= offset) rnums[i] = j;
  }
  else
  {
    i = 1; j = rnums[i]; rnums[0] = j - offset;
    for (; i < totrows; ++i, j += offset) rnums[i] = j;
  }
}

void xchg_bounds(RowInfo *rinfo, int nrows, int myid, int nprocs, int offset, int hasfirst, int haslast, MPI_Datatype MPI_Rowinfo)
{
  MPI_Request prev_rreq, prev_sreq, next_rreq, next_sreq;
  MPI_Status prev_rstat, prev_sstat, next_rstat, next_sstat;
  int prev_rflag, prev_sflag, next_rflag, next_sflag;
  int prevhaslast = 1, nexthasfirst = 1;
  prev_rflag = next_rflag = 0;
  int prevpeer = myid - offset, nextpeer = myid + offset; 

  // may not be necessary but I couldn't figure out rules I was sure of,
  // so the metadata exchange is done for safety
  if (prevpeer >= 0)
  {
    MPI_Isend(&hasfirst, 1, MPI_INT, prevpeer, 0, MPI_COMM_WORLD, &prev_sreq);
    MPI_Irecv(&prevhaslast, 1, MPI_INT, prevpeer, 0, MPI_COMM_WORLD, &prev_rreq);
  }
  if (nextpeer < nprocs)
  {
    MPI_Isend(&haslast, 1, MPI_INT, nextpeer, 0, MPI_COMM_WORLD, &next_sreq);
    MPI_Irecv(&nexthasfirst, 1, MPI_INT, nextpeer, 0, MPI_COMM_WORLD, &next_rreq);
  }
  if (prevpeer >= 0)
  {
    while (!prev_rflag || !prev_sflag)
    {
      MPI_Test(&prev_rreq, &prev_rflag, &prev_rstat);
      MPI_Test(&prev_sreq, &prev_sflag, &prev_sstat);
    }
    if (!prevhaslast) MPI_Isend(&rinfo[1], 1, MPI_Rowinfo, prevpeer, 0, MPI_COMM_WORLD, &prev_sreq);
    if (!hasfirst) MPI_Irecv(&rinfo[0], 1, MPI_Rowinfo, prevpeer, 0, MPI_COMM_WORLD, &prev_rreq);
  }
  if (nextpeer < nprocs)
  {
    while (!next_rflag || !next_sflag)
    {
      MPI_Test(&next_rreq, &next_rflag, &next_rstat);
      MPI_Test(&next_sreq, &next_sflag, &next_sstat);
    }
    if (!nexthasfirst) MPI_Isend(&rinfo[nrows], 1, MPI_Rowinfo, nextpeer, 0, MPI_COMM_WORLD, &next_sreq);
    if (!haslast) MPI_Irecv(&rinfo[nrows+1], 1, MPI_Rowinfo, nextpeer, 0, MPI_COMM_WORLD, &next_rreq);
  }

  prev_rflag = next_rflag = 0;
  if (!hasfirst && !haslast && prevpeer >= 0 && nextpeer < nprocs) 
    while (!prev_rflag || !next_rflag)
    {
      MPI_Test(&prev_rreq, &prev_rflag, &prev_rstat);
      MPI_Test(&next_rreq, &next_rflag, &next_rstat);
    }
  else if (!hasfirst && prevpeer >= 0)
    while (!prev_rflag)
      MPI_Test(&prev_rreq, &prev_rflag, &prev_rstat);
  else if (!haslast && nextpeer < nprocs) 
    while (!next_rflag)
      MPI_Test(&next_rreq, &next_rflag, &next_rstat);
}

int solve_2by2(RowInfo *myrows, Xsol *xs, int nrows)
{
  int i;
  int half = nrows / 2;
  int second;
  double det;
  for (i = 0; i < half; ++i)
  {
    second = i + half;
    xs[second].row = myrows[second].row;
    xs[i].row = myrows[i].row;
    if (!(det = myrows[i].b * myrows[second].b - myrows[second].a * myrows[i].c))
    {
      fprintf(stderr, "singular 2-by-2 matrix - cyclic reduction fails.\n");
      return -1;
    }
    xs[second].x = (myrows[i].b * myrows[second].y - myrows[second].a * myrows[i].y) / det;
    xs[i].x = (myrows[second].b * myrows[i].y - myrows[i].c * myrows[second].y) / det;
    //printf("%d %f; %d %f\n", xs[i].row, xs[i].x, xs[second].row, xs[second].x);
  }
  return nrows;
}

Xsol *write_x(Xsol *xs, char *outfile, int nrows)
{
  FILE *f;
  int i;
  if (!(f = fopen(outfile, "w")))
  {
    fprintf(stderr, "error opening file %s\n", outfile);
    return 0;
  }
  for (i = 0; i < nrows; ++i) fprintf(f, "%d %f\n", xs[i].row, xs[i].x);
  fclose(f);
  return xs;
}

void print_x(Xsol *xs, int nrows)
{
  int i;
  for (i = 0; i < nrows; ++i) printf("%d %f\n", xs[i].row, xs[i].x);
}

