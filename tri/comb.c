/*
This code is associated with "Computational Geometry in C" (Second Edition),
Chapter 1.  It generates a comb polygon suitable for input to
tri.c

Input: n
Output: n vertex coordinates for a comb polygon, in ccw order

Compile: gcc -o cone cone.c

Written by Joseph O'Rourke, with contributions by Min Xu.
Last modified: October 1998
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/

#include <stdio.h>
#define	X	0
#define	Y	1

#define DIM     2               /* Dimension of points */
typedef int     tPointi[DIM];   /* type integer point */
#define PMAX    10000           /* Max # of pts in polygon */
   
typedef tPointi tPolygoni[PMAX];/* type integer polygon */

main()
{
  tPolygoni P;
  int i, j, n;

  printf("n="); scanf("%d",&n);
  for ( i = 0; i < n/2; i++ ) {
    P[2*i][X] = n-2 - (2*i);
    P[2*i][Y] = 0;
    P[2*i+1][X] = n-2 - (2*i + 1);
    P[2*i+1][Y] = 10;
  }
  
  P[n-1][X] = (n-2) / 2;
  P[n-1][Y] = -2;

  printf("%d\n", n);
  for ( i = 0; i< n; i++)
    printf("%d %d\n", P[i][X], P[i][Y]);
}
