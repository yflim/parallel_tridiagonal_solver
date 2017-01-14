Note: Code written in late 2010 or early 2011!

(Example) Usage (on machines/clusters with MPI installations only):
$ make all
$ inputgen nums 32 1 3 1 
$ solve num outfile 32 

cyclicreduction.c defines the functions used in this program, and solve.c contains main.

This assignment is based on a version of the parallel cylic reduction algorithm for solving tridiagonal matrices. It is described under the headings "cyclic reduction" and "parallel cyclic reduction" on pages 2 and (the beginning of) 3 of http://www.cse.illinois.edu/courses/cs554/notes/09_tridiagonal_8up.pdf

I'm not sure it's the method of parallel cyclic reduction usually used for solving tridiagonal matrices, which may be the one described on pages 481 to 489 of http://books.google.com/books?id=EnqgXI1AXAQC&pg=PA481&lpg=PA481&dq=karniadakis+kirby+cyclic+reduction&source=bl&ots=XFXoqidVpo&sig=4u1qNInF1RijIHKAvHb0xVSEeJ8&hl=en&ei=MoiGTceXD-y00QH3lvjECA&sa=X&oi=book_result&ct=result&resnum=1&ved=0CBUQ6AEwAA#v=onepage&q&f=false.
However, it seems to have more parallelism than the one described by Karniadakis and Kirby (in which odd rows are recursively eliminated, so more processors drop out at each level, before being used again for back substitution), as well as having the advantage of not calling for back substitution.

However, it wasn't clear how to correctly compute the solution with multiple triples (of rows/equations) assigned to each processor. I found a way to do so for strictly diagonally dominant tridiagonal matrices*, which I shall try to describe roughly:

1. if n (number of rows) isn't 2^i for some i = 0 or a natural number, pad out the given system with redundant equations having x = 0 (and forming a strictly diagonally dominant submatrix of A'). Call the resulting tridiagonal matrix A', which is a square matrix of size ntot.

2. nprocs must also be some power of 2; in particular, 0 <= log(nprocs) <= log(ntot) - 1

3. So each processor has ntot/nprocs equations at every step. Initially, for example, for ntot = 16 and nprocs = 4:
Step 1: 1,2,3,4; 5,6,7,8; 9,10,11,12; 13,14,15,16;

4. The equations are recursively decoupled (as described in slides) at each level into 2 systems each consisting of either all odd equations or all even equations. For example:
Step 2: 1,3,5,7; 2,4,6,8; 9,11,13,15; 10,12,14,16;
Step 3: 1,5,9,13; 2,6,10,14; 3,7,11,15; 4,8,12,16;
Step 3: 1,9,5,13; 2,10,6,14; 3,11,7,15; 4,12,8,16;
This is done by dividing the processors, at each level of recursion, into pairs which trade equations (rows), with one sending its even rows to make up the first half of its peer's rows at the next level, and receiving its peer's odd rows to form the second half of its rows at the next level. For example,
the pairs at each level are:
Step 1: 1,2; 3,4;
Step 2: 1,3; 2,4;
5. After that, each processor "has all it needs", and the process of decoupling
the system into more and more smaller systems happens within each processor.
6. The base case occurs when n = 2, at which point the xs in each system are computed directly.

*b > |a| + |c| for all rows. It seems to guarantee that zeros on the diagonal will never appear during computation.
