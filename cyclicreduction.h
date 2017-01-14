#include "mpi.h"

typedef struct rowinfo
{
  //int ax, bx, cx;
  double a, b, c;
  double y;
  int row;
} RowInfo;

typedef struct xsol
{
  double x;
  int row;
} Xsol;

void rows_init(RowInfo*, int*, int, int);
void copy_rowinfo(RowInfo*, RowInfo*);
void copy_rowinfo_array(RowInfo*, RowInfo*, int);
void print_rowinfo(RowInfo*, int, int);
void print_rownums(int*, int, int);
MPI_Datatype create_rowinfotype();
MPI_Datatype create_xsoltype();
RowInfo *readtridiag(char*, RowInfo*, int, int, int, int, int);
void update_rows1and2(RowInfo*, RowInfo*, RowInfo*);
void update_rows(RowInfo*, RowInfo*, RowInfo*, RowInfo*);
void update_last2rows(RowInfo*, RowInfo*, RowInfo*);
int update_rowinfo(RowInfo*, int, int, int);
int update_rowinfo_serial(RowInfo*, int, int);
void pack(RowInfo*, int, int, int, int);
void xchg_oddevens(RowInfo*, int, int, int, int, int, MPI_Datatype);
void update_rownums(int*, int, int, int, int);
void xchg_bounds(RowInfo*, int, int, int, int, int, int, MPI_Datatype);
int solve_2by2(RowInfo*, Xsol*, int);
Xsol *write_x(Xsol*, char*, int);
void print_x(Xsol*, int);

