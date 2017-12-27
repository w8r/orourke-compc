/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 8.  It is not written to be comprehensible without the
explanation in that book.

Prints out one arm configuration to reach given target.
Assumes number of links >= 3.
Input:
   nlinks           Number of links
   L1 L2 ... Ln     Link lengths
   x0 y0            target0
   x1 x2            target1
   ...

Written by Joseph O'Rourke.
Last modified: December 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/

#include <stdio.h>
#include <math.h>
#define EXIT_SUCCESS 0
#define	X 0
#define	Y 1
#define	NLINKS	20
typedef	enum { FALSE, TRUE }	bool;

#define DIM     2               /* Dimension of points */
typedef int     tPointi[DIM];   /* type integer point */
typedef double  tPointd[DIM];   /* type double point */

/*----------Function Prototypes-------------*/
void	ReadTarget( tPointi target );
int	ReadLinks( void );
void	PrintLinks( void );
void	LineTo_i( tPointi p );
void	MoveTo_i( tPointi p );
void	LineTo_d( tPointd p );
bool	Solve3( int L1, int L2, int L3, tPointi target );
int	TwoCircles( tPointi c1, int r1, tPointi c2, int r2, tPointd p );
int	TwoCircles0a( int r1, tPointi c2, int r2, tPointd p );
int	TwoCircles0b( int r1, tPointi c2, int r2, tPointd p );
void	TwoCircles00( int r1, double a2, int r2, tPointd p );
void	SubVec( tPointi a, tPointi b, tPointi c );
bool	EqPointi( tPointi a, tPointi b );
double	Length2( tPointi v );
bool	Solve2( int L1, int L2, tPointi target, tPointd J );
bool	Nested( int L1, int L2, int L3, tPointi target );
bool	Solven( int nlinks );
double	RadDeg( double x );
void	PostscriptOpen( void );
void	PostscriptClose( void );
/*------------------------------------------*/

/* Global variables. */
int	linklen[NLINKS];	/* link lengths */
int	nlinks;			/* number of links */
tPointi	target;			/* target point */

main()
{
   tPointi origin = {0,0};

   nlinks = ReadLinks();
   PostscriptOpen();
   while (TRUE) { /* loop broken by EOF in ReadTarget */
      ReadTarget( target );
      printf("newpath\n");
      MoveTo_i( origin );
      if ( !Solven( nlinks ) ) 
         printf("Solven: no solutions!\n");
      LineTo_i( target );
      printf("closepath stroke\n");
   }
}


bool    Solven( int nlinks )
{
   int i;
   int m;              /* index of median link */
   int L1, L2, L3;     /* length of links between kinks */
   int totlength;      /* total length of all links */
   int halflength;     /* floor of half of total */

   /* Compute total & half length. */
   totlength = 0;
   for ( i = 0; i < nlinks; i++ ) 
           totlength += linklen[i];
   halflength = totlength / 2;

   /* Find median link. */
   L1 = 0;
   for ( m = 0; m < nlinks; m++ ) {
      if ( (L1 + linklen[m]) > halflength)
         break;
      L1 += linklen[m];
   }

   L2 = linklen[m];
   L3 = totlength - L1 - L2;
   if ( Solve3( L1, L2, L3, target ) )
        return  TRUE;
   else return  FALSE;
}

bool    Solve3( int L1, int L2, int L3, tPointi target )
{
   tPointd Jk;      /* coords of kinked joint returned by Solve2 */
   tPointi J1;      /* Joint1 on x-axis */
   tPointi Ttarget; /* translated target */

   printf("%%==>Solve3: links = %d,%d,%d\n", L1, L2, L3);
   if ( Solve2( L1 + L2, L3, target, Jk ) ) {
      printf("%%Solve3: link1=%d, link2=%d, joint=\n", L1 + L2, L3);
      LineTo_d( Jk );
      return TRUE;
   }
   else if ( Solve2( L1, L2 + L3, target, Jk ) ) {
      printf("%%Solve3: link1=%d, link2=%d, joint=\n", L1, L2 + L3);
      LineTo_d( Jk );
      return TRUE;
   }
   else {   /* pin J0 to 0. */
      /* Shift so J1 is origin. */
      J1[X] = L1;   J1[Y] = 0;
      SubVec( target, J1, Ttarget );
      if ( Solve2( L2, L3, Ttarget, Jk ) ) {
         /* Shift solution back to origin. */
         Jk[X] += L1;
         printf("%%Solve3: link1=%d, link2=%d, link3=%d, joints=\n", L1, L2, L3);
         LineTo_i( J1 );
         LineTo_d( Jk );
         return TRUE;
      }
      else 
         return FALSE;
   }
}
bool    Solve2( int L1, int L2, tPointi target, tPointd J )
{
   tPointi c1 = {0,0};  /* center of circle 1 */
   int nsoln;           /* # of solns: 0,1,2,3(infinite) */

   nsoln = TwoCircles( c1, L1, target, L2, J );
   printf("%%<==Solve2: L1=%d, L2=%d; nsoln=%d, target=(%d,%d),J=(%g,%g)\n", 
      L1, L2, nsoln,target[X],target[Y],J[X],J[Y]); 
   return   nsoln != 0;
}


/*---------------------------------------------------------------------
TwoCircles finds an intersection point between two circles.
General routine: no assumptions. Returns # of intersections; point in p.
---------------------------------------------------------------------*/
int     TwoCircles( tPointi c1, int r1, tPointi c2, int r2, tPointd p)
{
   tPointi c;
   tPointd q;
   int nsoln = -1;

   /* Translate so that c1={0,0}. */
   SubVec( c2, c1, c );
   nsoln = TwoCircles0a( r1, c, r2, q );
   /* Translate back. */
   p[X] = q[X] + c1[X];
   p[Y] = q[Y] + c1[Y];
   return nsoln;
}

/*---------------------------------------------------------------------
TwoCircles0a assumes that the first circle is centered on the origin.
Returns # of intersections: 0, 1, 2, 3 (inf); point in p.
---------------------------------------------------------------------*/
int     TwoCircles0a( int r1, tPointi c2, int r2, tPointd p )
{
   double dc2;              /* dist to center 2 squared */
   double rplus2, rminus2;  /* (r1 +/- r2)^2 */
   double f;                /* fraction along c2 for nsoln=1 */

   /* Handle special cases. */
   dc2 = Length2( c2 );
   rplus2  = (r1 + r2) * (r1 + r2);
   rminus2 = (r1 - r2) * (r1 - r2);

   /* No solution if c2 out of reach + or -. */
   if ( ( dc2 > rplus2 ) || ( dc2 < rminus2 ) )
      return   0;

   /* One solution if c2 just reached. */
   /* Then solution is r1-of-the-way (f) to c2. */
   if ( dc2 == rplus2 ) {
      f = r1 / (double)(r1 + r2);
      p[X] = f * c2[X];   p[Y] = f * c2[Y];
      return 1;
   }
   if ( dc2 == rminus2 ) {
      if ( rminus2 == 0 ) {   /* Circles coincide. */
         p[X] = r1;    p[Y] = 0;
         return 3;
      }
      f = r1 / (double)(r1 - r2);
      p[X] = f * c2[X];   p[Y] = f * c2[Y];
      return 1;
   }

   /* Two intersections. */
   return TwoCircles0b( r1, c2, r2, p );
}

/*---------------------------------------------------------------------
TwoCircles0b also assumes that the 1st circle is origin-centered.
---------------------------------------------------------------------*/
int     TwoCircles0b( int r1, tPointi c2, int r2, tPointd p )
{
   double a2;          /* center of 2nd circle when rotated to x-axis */
   tPointd q;          /* one solution when c2 on x-axis */
   double cost, sint;  /* sine and cosine of angle of c2 */

   /* Rotate c2 to a2 on x-axis. */
   a2 = sqrt( Length2( c2 ) );
   cost = c2[X] / a2;
   sint = c2[Y] / a2;

   TwoCircles00( r1, a2, r2, q );

   /* Rotate back */
   p[X] =  cost * q[X] + -sint * q[Y];
   p[Y] =  sint * q[X] +  cost * q[Y];
   
   return 2;
}

/*---------------------------------------------------------------------
TwoCircles00 assumes circle centers are (0,0) and (a2,0).
---------------------------------------------------------------------*/
void     TwoCircles00( int r1, double a2, int r2, tPointd p )
{
   double r1sq, r2sq;
   r1sq = r1*r1;
   r2sq = r2*r2;

   /* Return only positive-y soln in p. */
   p[X] = ( a2 + ( r1sq - r2sq ) / a2 ) / 2;
   p[Y] = sqrt( r1sq - p[X]*p[X] );
   printf("%%TwoCircles00: p=(%g,%g)\n", p[X], p[Y]);
}
void   ReadTarget( tPointi target )
{
   int i, readerror;

   for( i = 0; i < DIM; i++) {
      readerror = scanf("%d", &target[i] );
      if ( readerror == EOF ) {
    PostscriptClose();
    exit( EXIT_SUCCESS );
      }
   }
   printf("%%Target point = (%d,%d)\n", target[X], target[Y]);
}

void   LineTo_i( tPointi p )
{
   printf("%5d    %5d    lineto\n", p[X], p[Y] );
}
void   MoveTo_i( tPointi p )
{
   printf("%5d    %5d    moveto\n", p[X], p[Y] );
}

void   LineTo_d( tPointd p )
{
   printf("%8.2lf %8.2lf lineto\n", p[X], p[Y] );
}

int   ReadLinks( void )
{
   int i, length;

   scanf("%d", &nlinks);
   for( i=0; i < nlinks; i++) {
      scanf("%d", &length);
      linklen[i] = length;
   }
   /*PrintLinks();*/
   return nlinks;
}
      
void   PrintLinks()
{
   int   i;

   printf("%%");
   for ( i=0; i < nlinks; i++)
      printf("%d:%d\t", i, linklen[i]);
   putchar('\n');
}
bool   EqPointi( tPointi a, tPointi b )
{
   int   i;

   for ( i=0; i < DIM; i++ )
      if ( a[i] != b[i] )
    return FALSE;
   return TRUE;
}

double   Length2( tPointi v )
{
   int   i;
   double   ss;

   ss = 0;
   for ( i=0; i < DIM; i++ )
      ss += v[i] * (double)v[i];
   return ss;
}



double   RadDeg( double x )
{
   return 180 * x / M_PI;
}
void   SubVec( tPointi a, tPointi b, tPointi c )
{
   int   i;

   for ( i=0; i < DIM; i++ ) {
      c[i] = a[i] - b[i];
   }
}
void	PostscriptOpen( void )
{
   int xmin, xmax, ymin, ymax;
   int i, totlength=0;

   for ( i = 0; i < nlinks; i++ )
      totlength += linklen[i];
   xmin = ymin = -totlength;
   xmax = ymax = totlength;

   printf("%%!PS\n");
   printf("%%%%Creator: arm.c (Joseph O'Rourke)\n");
   printf("%%%%BoundingBox: %d %d %d %d\n",
      xmin, ymin, xmax, ymax);
   printf("%%%%EndComments\n");
   printf(".00 .00 setlinewidth\n");
   printf("%d %d translate\n", -xmin+100, -ymin+100 );
   /* The +100 shifts the figure from the lower left corner. */
}
void	PostscriptClose( void )
{
   printf("showpage\n%%%%EOF\n");
}
