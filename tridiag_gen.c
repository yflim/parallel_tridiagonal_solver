#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <math.h>

// generates a, b, c of a strictly diagonally dominant tridiagonal matrix 
// with identical values on each diagonal, as well as y (solution)
int main(int argc, char **argv)
{
  FILE *file;
  int nrows;
  int i;
  double aval, bval, cval;

  if (argc != 6) { printf("usage: inputgen writefile nrows aval bval cval\n"); exit(1); }
  nrows = atoi(argv[2]);
  aval = atof(argv[3]);
  bval = atof(argv[4]);
  cval = atof(argv[5]);
  if (bval <= fabs(aval) + fabs(cval)) 
  {
    printf("bval should be greater than the absolute values of aval and bval summed\n");
    exit(1);
  }

  file = fopen(argv[1], "w");
  for (i = 0; i < nrows; ++i)
  {
    fprintf(file, "%.2f ", aval);
    fprintf(file, "%.2f ", bval);
    fprintf(file, "%.2f ", cval);
    fprintf(file, "%.2f ", (double) i);
    fprintf(file, "\n");
  }
  fclose(file);

  return 0;
}

