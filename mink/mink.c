/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 8.  It is not written to be comprehensible without the
explanation in that book.

Compile: gcc -o mink mink.c (or simply: make)

Written by Joseph O'Rourke.
Last modified: December 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/

#include   <stdio.h>
#include   <math.h>
#include   <stdlib.h>

#define EXIT_FAILURE   1
#define X 0
#define Y 1
typedef enum { FALSE, TRUE }   bool;

#define DIM 2               /* Dimension of points */
typedef int tPointi[DIM];   /* Type integer point */


/*----------Point(s) Structure-------------*/
typedef struct tPointStructure tsPoint;
typedef tsPoint *tPoint;
struct tPointStructure {
   int     vnum;
   tPointi v;
   bool    primary;
};

/* Global variables */
#define PMAX    1000               /* Max # of points */
typedef tsPoint tPointArray[PMAX];
static tPointArray P;

int m;  /* Total number of points in both polygons */
int n;  /* Number of points in primary polygon */
int s;  /* Number of points in secondary polygon */



/*----------Function Prototypes-------------*/
int     AreaSign( tPointi a, tPointi b, tPointi c );
int     ReadPoints( tPointi p0 );
void    PrintPoints( void );
void    SubVec( tPointi a, tPointi b, tPointi c );
void    AddVec( tPointi a, tPointi b, tPointi c );
void    Vectorize( void );
int     Compare( const void *tp1, const void *tp2 );
void    OutputPolygons();
void    Convolve( int j0, tPointi p0 );
/*------------------------------------------*/


main()
{
   tPointi p0 = {0,0};
   int j0;               /* index of start point */

   j0 = ReadPoints( p0 );
   PrintPoints();
   OutputPolygons();
   Vectorize();
   PrintPoints();
   qsort(
      &P[0],             /* pointer to 1st elem */
      m,                 /* number of elems */
      sizeof( tsPoint ), /* size of each elem */
      Compare            /* -1,0,+1 compare function */
   );
   printf("%%After qsort of vectors:");
   PrintPoints();
   Convolve( j0, p0 );
}

int     AreaSign( tPointi a, tPointi b, tPointi c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}

/*---------------------------------------------------------------------
Reads in the coordinates of the points from stdin, filling in P,
and setting the three global counters m=n+s.
Returns a start point p0.
---------------------------------------------------------------------*/
int    ReadPoints( tPointi p0 )
{
   int i;
   int xmin, ymin, xmax, ymax;     /* Primary min & max */
   int sxmin, symin, sxmax, symax; /* Secondary min & max */
   int mp, ms; /* i index of max (u-r) primary and secondary points */

   m = 0;
   scanf("%d", &n );
   if ( n > PMAX )
      printf("Error in ReadPoints:  > %d\n", PMAX), exit(EXIT_FAILURE);
   for( i = 0; i < n; i++ ) {
      scanf("%d %d",&P[i].v[X],&P[i].v[Y]);
      P[i].vnum = i;
      P[i].primary = TRUE;
      m++;
   }
   /*printf("%%ReadPoints: %d primary points read\n", n);*/

   scanf("%d", &s );
   if ( n+s > PMAX )
      printf("Error in ReadPoints:  > %d\n", PMAX), exit(EXIT_FAILURE);
   for( i = 0; i < s; i++ ) {
      scanf("%d %d",&P[n+i].v[X],&P[n+i].v[Y]);
      /* Reflect secondary polygon */
      P[n+i].v[X] = -P[n+i].v[X];
      P[n+i].v[Y] = -P[n+i].v[Y];
      P[n+i].vnum = i;
      P[n+i].primary = FALSE;
      m++;
   }
   /*printf("%%ReadPoints: %d secondary points read\n", s);*/

   /* Compute Bounding Box and output Postscript header. */
   xmin = xmax = P[0].v[X];
   ymin = ymax = P[0].v[Y];
   mp = 0;
   for (i = 1; i < n; i++) {
      if      ( P[i].v[X] > xmax ) xmax = P[i].v[X];
      else if ( P[i].v[X] < xmin ) xmin = P[i].v[X];
      if      ( P[i].v[Y] > ymax ) {ymax = P[i].v[Y]; mp = i;}
      else if ( P[i].v[Y] == ymax && (P[i].v[X] > P[mp].v[X]) ) mp = i;
      else if ( P[i].v[Y] < ymin ) ymin = P[i].v[Y];
   }
   printf("%%Index of upper rightmost primary, i=mp = %d\n", mp);
   printf("%%highest primary: %d %d\n", P[mp].v[X], P[mp].v[Y]);

   sxmin = sxmax = P[n].v[X];
   symin = symax = P[n].v[Y];
   ms = n;
   /* ms update lines corrected by Rishikesh.Parthasarathi@inrialpes.fr */
   for (i = 1; i < s; i++) {
      if      ( P[n+i].v[X] > sxmax ) sxmax = P[n+i].v[X];
      else if ( P[n+i].v[X] < sxmin ) sxmin = P[n+i].v[X];
      if      ( P[n+i].v[Y] > symax ) {symax = P[n+i].v[Y]; ms+= i;}
      else if ( P[n+i].v[Y] == symax && (P[n+i].v[X] > P[ms].v[X]) ) ms+= i;
      else if ( P[n+i].v[Y] < symin ) symin = P[n+i].v[Y];
   }
   printf("%%Index of upper rightmost secondary, i=ms = %d\n", ms);
   printf("%%highest secondary: %d %d\n", P[ms].v[X], P[ms].v[Y]);

   /* Compute the start point: upper rightmost of both. */
   printf("%%p0=(%d,%d)\n", p0[X], p0[Y]);
   AddVec( p0, P[mp].v, p0 );
   printf("%%p0=(%d,%d)\n", p0[X], p0[Y]);
   AddVec( p0, P[ms].v, p0 );
   printf("%%p0=(%d,%d)\n", p0[X], p0[Y]);

   /* PostScript header */
   printf("%%!PS\n");
   printf("%%%%Creator: mink.c (Joseph O'Rourke)\n");
   printf("%%%%BoundingBox: %d %d %d %d\n", 
      ((xmin < sxmin) ? xmin : sxmin),
      ((ymin < symin) ? ymin : symin),
      ((xmax > sxmax) ? xmax : sxmax),
      ((xmax > sxmax) ? xmax : sxmax)
      );
   printf("%%%%EndComments\n");
   printf(".00 .00 setlinewidth\n");
   printf("%d %d translate\n", 
      -((xmin < sxmin) ? xmin : sxmin) + 72,
      -((ymin < symin) ? ymin : symin) + 72 );
   /* The +72 shifts the figure one inch from the lower left corner */

   printf("%%j0 = %d\n", mp);
   return mp; /* j0: starting index. */
}

void   PrintPoints( void )
{
   int i;

   printf("%%Points:\n");
   for( i = 0; i < m; i++ )
      printf("%%i=%d: primary=%d | vnum=%3d, x=%4d, y=%4d\n", 
             i, P[i].primary, P[i].vnum, P[i].v[X], P[i].v[Y]);
}

void   OutputPolygons( void )
{
   int i;

   printf("\n%%Primary Polygon:\n");
   printf("newpath\n");
   printf("%d\t%d\tmoveto\n", P[0].v[X], P[0].v[Y]);
   for( i = 1; i <= n; i++ )
      printf("%d\t%d\tlineto\n", P[i%n].v[X], P[i%n].v[Y]);
   printf("closepath\n0.8 setgray fill stroke\n0 setgray");

   printf("\n%%Secondary Polygon (reflected):\n");
   printf("newpath\n");
   printf("%d\t%d\tmoveto\n", P[n].v[X], P[n].v[Y]);
   for( i = 1; i <= s; i++ )
      printf("%d\t%d\tlineto\n", P[n+(i%s)].v[X], P[n+(i%s)].v[Y]);
   printf("closepath stroke\n");
 }

/*---------------------------------------------------------------------
a - b ==> c.
---------------------------------------------------------------------*/
void    SubVec( tPointi a, tPointi b, tPointi c )
{
   int i;

   /*printf("a=(%d,%d); b=(%d,%d)\n",a[X],a[Y],b[X],b[Y]);*/
   for( i = 0; i < DIM; i++ )
      c[i] = a[i] - b[i];
}

/*---------------------------------------------------------------------
a + b ==> c.
---------------------------------------------------------------------*/
void    AddVec( tPointi a, tPointi b, tPointi c )
{
   int i;

   /*printf("a=(%d,%d); b=(%d,%d)\n",a[X],a[Y],b[X],b[Y]);*/
   for( i = 0; i < DIM; i++ )
      c[i] = a[i] + b[i];
}

void   Vectorize( void )
{
   int i;
   tPointi last;  /* Holds the last vector difference. */

   printf("%%Vectorize\n");
   SubVec( P[0].v, P[n-1].v, last );
   for( i = 0; i < n-1; i++ ) {
      /*printf("i/i+1=%d,%d\n", i, i+1 );*/
      SubVec( P[i+1].v, P[i].v, P[i].v );
   }
   P[n-1].v[X] = last[X];
   P[n-1].v[Y] = last[Y];

   SubVec( P[n].v, P[n+s-1].v, last );
   for( i = 0; i < s-1; i++ ) {
      /*printf("i/i+1=%d,%d\n", n+i, n+i+1 );*/
      SubVec( P[n+i+1].v, P[n+i].v, P[n+i].v );
   }
   P[n+s-1].v[X] = last[X];
   P[n+s-1].v[Y] = last[Y];

}

int   Compare( const void *tpi, const void *tpj )
{
   int a;             /* AreaSign result */
   int x, y;          /* projections in 1st quadrant */
   tPoint pi, pj;     /* Recasted points */
   tPointi Origin = {0,0};
   pi = (tPoint)tpi;
   pj = (tPoint)tpj;

   /* A vector in the open   upper halfplane is after
      a vector in the closed lower halfplane. */
   if      ( ( pi->v[Y] > 0 ) && (pj->v[Y] <= 0 ) )
      return  1;
   else if ( ( pi->v[Y] <= 0 ) && (pj->v[Y] > 0 ) )
      return -1;

   /* A vector on the x-axis and one in the lower halfplane
      are handled by the Left computation below. */

   /* Both vectors on the x-axis requires special handling. */
   else if ( ( pi->v[Y] == 0 ) && (pj->v[Y] == 0 ) ) { 
      if      ( ( pi->v[X] < 0 ) && ( pj->v[X] > 0 ) )
         return -1;
      if      ( ( pi->v[X] > 0 ) && ( pj->v[X] < 0 ) )
         return  1;
      else if ( abs(pi->v[X]) < abs(pj->v[X]) )
         return  -1;
      else if ( abs(pi->v[X]) > abs(pj->v[X]) )
         return   1;
      else
         return  0;
   }

   /* Otherwise, both in open upper halfplane, 
      or both in closed lower halfplane, but not both on x-axis. */
   else { 
   
      a = AreaSign( Origin, pi->v, pj->v );
      if      (a > 0)
         return -1;
      else if (a < 0)
         return 1;
      else { /* Begin collinear */
         x =  abs( pi->v[X] ) - abs( pj->v[X] );
         y =  abs( pi->v[Y] ) - abs( pj->v[Y] );
   
         if      ( (x < 0) || (y < 0) )
            return -1;
         else if ( (x > 0) || (y > 0) )
            return 1;
         else /* points are coincident */
            return 0;
      } /* End collinear */
   }
}


void    Convolve( int j0, tPointi p )
{
   int i;  /* Index into sorted edge vectors P */
   int j;  /* Primary polygon indices */

   printf("%%Convolve: Start array i = %d, primary j0= %d\n", i, j0);
   printf("1 1 setlinewidth\n");
   printf("newpath\n");
   printf("%d\t%d\tmoveto\n", p[X], p[Y]);

   i = 0;  /* Start at angle -pi, rightward vector. */
   j = j0; /* Start searching for j0. */
   do {

      /* Advance around secondary edges until next j reached. */
      while ( !(P[i].primary && P[i].vnum == j) ) {
	 if ( !P[i].primary ) {
            AddVec( p, P[i].v, p );
            printf("%d\t%d\tlineto\n", p[X], p[Y]);
	 }
         i = (i+1)%m;
	 printf("%%X: i incremented to %d\n", i);
      }

      /* Advance one primary edge. */
      printf("%%X: j=%d found at i=%d\n", j, i);
      AddVec( p, P[i].v, p );
      printf("%d\t%d\tlineto\n", p[X], p[Y]);
      j = (j+1)%n;
      printf("%%X: j incremented to %d\n", j);

   } while ( j != j0);

   /* Finally, complete circuit on secondary/robot polygon. */
   while (i != 0) {
      if ( !P[i].primary ) {
         AddVec( p, P[i].v, p );
         printf("%d\t%d\tlineto\n", p[X], p[Y]);
      }
      i = (i+1)%m;
      printf("%%X: i incremented to %d in final circuit\n", i);
   }

   printf("closepath stroke\n");
   printf("showpage\n%%%%EOF\n");
}
