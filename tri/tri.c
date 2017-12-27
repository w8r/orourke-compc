/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 1.  It is not written to be comprehensible without the 
explanation in that book.

Input: 2n integer coordinates for vertices of a simple polygon,
       in counterclockwise order.  NB: the code will not work
       for points in clockwise order!
Output: the diagonals of a triangulation, in PostScript.

Compile: gcc -o tri tri.c (or simply: make)

Written by Joseph O'Rourke, with contributions by Min Xu.
Last modified: October 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#define	EXIT_FAILURE 1
#define	X 0
#define	Y 1
typedef	enum { FALSE, TRUE } bool;

#define	DIM 2               /* Dimension of points */
typedef	int tPointi[DIM];   /* Type integer point */

typedef	struct tVertexStructure tsVertex;	/* Used only in NEW(). */
typedef	tsVertex *tVertex;
struct	tVertexStructure {
   int		vnum;		/* Index */
   tPointi	v;		/* Coordinates */
   bool 	ear;		/* TRUE iff an ear */
   tVertex 	next,prev;
};

/* Global variable definitions */
tVertex	vertices  = NULL;	/* "Head" of circular list. */
int	nvertices = 0;		/* Total number of polygon vertices. */

#include "macros.h"	/* NEW, ADD */

/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/
int	AreaPoly2( void );
void	Triangulate( void );
void	EarInit( void );
bool	Diagonal( tVertex a, tVertex b );
bool	Diagonalie( tVertex a, tVertex b );
bool	InCone( tVertex a, tVertex b );

int     Area2( tPointi a, tPointi b, tPointi c );
int     AreaSign( tPointi a, tPointi b, tPointi c );
bool	Xor( bool x, bool y );
bool	Left( tPointi a, tPointi b, tPointi c );
bool	LeftOn( tPointi a, tPointi b, tPointi c );
bool	Collinear( tPointi a, tPointi b, tPointi c );
bool	Between( tPointi a, tPointi b, tPointi c );
bool	Intersect( tPointi a, tPointi b, tPointi c, tPointi d );
bool	IntersectProp( tPointi a, tPointi b, tPointi c, tPointi d );

tVertex	MakeNullVertex( void );
void	ReadVertices( void );
void	PrintVertices( void );
void	PrintDiagonal( tVertex a, tVertex b );
void	PrintPoly( void );
/*-------------------------------------------------------------------*/

main()
{
   ReadVertices();
   PrintVertices();
   printf("%%Area of polygon = %g\n", 0.5 * AreaPoly2() );
   Triangulate();
   printf("showpage\n%%%%EOF\n");
}

/*---------------------------------------------------------------------
Returns twice the signed area of the triangle determined by a,b,c.
The area is positive if a,b,c are oriented ccw, negative if cw,
and zero if the points are collinear.
---------------------------------------------------------------------*/
int     Area2( tPointi a, tPointi b, tPointi c )
{
   return
          (b[X] - a[X]) * (c[Y] - a[Y]) -
          (c[X] - a[X]) * (b[Y] - a[Y]);
}


/*---------------------------------------------------------------------
Exclusive or: TRUE iff exactly one argument is true.
---------------------------------------------------------------------*/
bool	Xor( bool x, bool y )
{
   /* The arguments are negated to ensure that they are 0/1 values. */
   /* (Idea due to Michael Baldwin.) */
   return   !x ^ !y;
}

/*---------------------------------------------------------------------
Returns true iff ab properly intersects cd: they share
a point interior to both segments.  The properness of the
intersection is ensured by using strict leftness.
---------------------------------------------------------------------*/
bool	IntersectProp( tPointi a, tPointi b, tPointi c, tPointi d )
{
   /* Eliminate improper cases. */
   if (
      Collinear(a,b,c) ||
      Collinear(a,b,d) ||
      Collinear(c,d,a) ||
      Collinear(c,d,b)
      )
      return FALSE;

   return
         Xor( Left(a,b,c), Left(a,b,d) )
      && Xor( Left(c,d,a), Left(c,d,b) );
}

/*---------------------------------------------------------------------
Returns true iff c is strictly to the left of the directed
line through a to b.
---------------------------------------------------------------------*/
bool	Left( tPointi a, tPointi b, tPointi c )
{ 
   return  AreaSign( a, b, c ) > 0;
}

bool	LeftOn( tPointi a, tPointi b, tPointi c )
{
   return  AreaSign( a, b, c ) >= 0;
}

bool	Collinear( tPointi a, tPointi b, tPointi c )
{
   return  AreaSign( a, b, c ) == 0;
}

/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
First checks that c is collinear with a and b.
---------------------------------------------------------------------*/
bool	Between( tPointi a, tPointi b, tPointi c )
{
   tPointi	ba, ca;

   if ( ! Collinear( a, b, c ) )
      return  FALSE;

   /* If ab not vertical, check betweenness on x; else on y. */
   if ( a[X] != b[X] ) 
      return ((a[X] <= c[X]) && (c[X] <= b[X])) ||
             ((a[X] >= c[X]) && (c[X] >= b[X]));
   else
      return ((a[Y] <= c[Y]) && (c[Y] <= b[Y])) ||
             ((a[Y] >= c[Y]) && (c[Y] >= b[Y]));
}

/*---------------------------------------------------------------------
Returns TRUE iff segments ab and cd intersect, properly or improperly.
---------------------------------------------------------------------*/
bool	Intersect( tPointi a, tPointi b, tPointi c, tPointi d )
{
   if      ( IntersectProp( a, b, c, d ) )
      return  TRUE;
   else if (   Between( a, b, c )
            || Between( a, b, d )
            || Between( c, d, a )
            || Between( c, d, b )
           )
      return  TRUE;
   else    return  FALSE;
}

/*---------------------------------------------------------------------
Returns TRUE iff (a,b) is a proper internal *or* external
diagonal of P, *ignoring edges incident to a and b*.
---------------------------------------------------------------------*/
bool   Diagonalie( tVertex a, tVertex b )
{
   tVertex c, c1;

   /* For each edge (c,c1) of P */
   c = vertices;
   do {
      c1 = c->next;
      /* Skip edges incident to a or b */
      if (    ( c != a ) && ( c1 != a )
           && ( c != b ) && ( c1 != b )
           && Intersect( a->v, b->v, c->v, c1->v )
         )
         return FALSE;
      c = c->next;
   } while ( c != vertices );
   return TRUE;
}

/*---------------------------------------------------------------------
This function initializes the data structures, and calls
Triangulate2 to clip off the ears one by one.
---------------------------------------------------------------------*/
void   EarInit( void )
{
   tVertex v0, v1, v2;   /* three consecutive vertices */

   /* Initialize v1->ear for all vertices. */
   v1 = vertices;
   printf("newpath\n");
   do {
      v2 = v1->next;
      v0 = v1->prev;
      v1->ear = Diagonal( v0, v2 );
      v1 = v1->next;
   } while ( v1 != vertices );
   printf("closepath stroke\n\n");

}


/*---------------------------------------------------------------------
Prints out n-3 diagonals (as pairs of integer indices)
which form a triangulation of P.
---------------------------------------------------------------------*/

void   Triangulate( void )
{
   tVertex v0, v1, v2, v3, v4;	/* five consecutive vertices */
   int   n = nvertices;		/* number of vertices; shrinks to 3. */
   bool earfound;		/* for debugging and error detection only. */

   EarInit();
   printf("\nnewpath\n");
   /* Each step of outer loop removes one ear. */
   while ( n > 3 ) {     
      /* Inner loop searches for an ear. */
      v2 = vertices;
      earfound = FALSE;
      do {
         if (v2->ear) {
            earfound = TRUE;
            /* Ear found. Fill variables. */
            v3 = v2->next; v4 = v3->next;
            v1 = v2->prev; v0 = v1->prev;

            /* (v1,v3) is a diagonal */
            PrintDiagonal( v1, v3 );
            
            /* Update earity of diagonal endpoints */
            v1->ear = Diagonal( v0, v3 );
            v3->ear = Diagonal( v1, v4 );
            
            /* Cut off the ear v2 */
            v1->next = v3;
            v3->prev = v1;
            vertices = v3;	/* In case the head was v2. */
            n--;
            break;   /* out of inner loop; resume outer loop */
         } /* end if ear found */
         v2 = v2->next;
      } while ( v2 != vertices );

      if ( !earfound ) {
         printf("%%Error in Triangulate:  No ear found.\n");
         PrintPoly();
         printf("showpage\n%%%%EOF\n");
         exit(EXIT_FAILURE);
      }
   } /* end outer while loop */
   printf("closepath stroke\n\n");
}


/*---------------------------------------------------------------------
Returns TRUE iff the diagonal (a,b) is strictly internal to the 
polygon in the neighborhood of the a endpoint.  
---------------------------------------------------------------------*/
bool   InCone( tVertex a, tVertex b )
{
   tVertex a0,a1;	/* a0,a,a1 are consecutive vertices. */

   a1 = a->next;
   a0 = a->prev;

   /* If a is a convex vertex ... */
   if( LeftOn( a->v, a1->v, a0->v ) )
       return    Left( a->v, b->v, a0->v )
              && Left( b->v, a->v, a1->v );

   /* Else a is reflex: */
       return !(    LeftOn( a->v, b->v, a1->v )
                 && LeftOn( b->v, a->v, a0->v ) );
}

/*---------------------------------------------------------------------
Returns TRUE iff (a,b) is a proper internal diagonal.
---------------------------------------------------------------------*/
bool	Diagonal( tVertex a, tVertex b )
{
   return InCone( a, b ) && InCone( b, a ) && Diagonalie( a, b );
}


/*---------------------------------------------------------------------
ReadVertices: Reads in the vertices, and links them into a circular
list with MakeNullVertex.  There is no need for the # of vertices to be
the first line: the function looks for EOF instead.
---------------------------------------------------------------------*/
void   ReadVertices( void )
{
   tVertex	v;
   int		x, y;
   int		vnum = 0;

   while ( scanf ("%d %d", &x, &y ) != EOF )  {
      v = MakeNullVertex();
      v->v[X] = x;
      v->v[Y] = y;
      v->vnum = vnum++;
   }
   nvertices = vnum;
   if (nvertices < 3) 
      printf("%%Error in ReadVertices: nvertices=%d<3\n", nvertices),
      exit(EXIT_FAILURE);
}
/*---------------------------------------------------------------------
MakeNullVertex: Makes a vertex.
---------------------------------------------------------------------*/
tVertex   MakeNullVertex( void )
{
   tVertex  v;
   
   NEW( v, tsVertex );
   ADD( vertices, v );
   return v;
}

/*---------------------------------------------------------------------
For debugging; not currently invoked.
---------------------------------------------------------------------*/
void   PrintPoly( void )
{
   tVertex  v;
   printf("%%Polygon circular list:\n");
   v = vertices;
   do {                                 
      printf( "%% vnum=%5d:\tear=%d\n", v->vnum, v->ear );
      v = v->next;
   } while ( v != vertices );
}

/*---------------------------------------------------------------------
Print: Prints out the vertices.  Uses the vnum indices 
corresponding to the order in which the vertices were input.
Output is in PostScript format.
---------------------------------------------------------------------*/
void   PrintVertices( void )
{
   /* Pointers to vertices, edges, faces. */
   tVertex  v;
   int xmin, ymin, xmax, ymax;

   /* Compute bounding box for Encapsulated PostScript. */
   v = vertices;
   xmin = xmax = v->v[X];
   ymin = ymax = v->v[Y];
   do {
      if      ( v->v[X] > xmax ) xmax = v->v[X];
      else if ( v->v[X] < xmin ) xmin = v->v[X];
      if      ( v->v[Y] > ymax ) ymax = v->v[Y];
      else if ( v->v[Y] < ymin ) ymin = v->v[Y];
      v = v->next;
   } while ( v != vertices );
   
   /* PostScript header */
   printf("%%!PS\n");
   printf("%%%%Creator: tri.c (Joseph O'Rourke)\n");
   printf("%%%%BoundingBox: %d %d %d %d\n", xmin, ymin, xmax, ymax);
   printf("%%%%EndComments\n");
   printf(".00 .00 setlinewidth\n");
   printf("%d %d translate\n", -xmin+72, -ymin+72 );
   /* The +72 shifts the figure one inch from the lower left corner */

   /* Output vertex info as a PostScript comment. */
   printf("\n%% number of vertices = %d\n", nvertices);
   v = vertices;
   do {                                 
      printf( "%% vnum=%5d:\tx=%5d\ty=%5d\n", v->vnum, v->v[X], v->v[Y] );
      v = v->next;
   } while ( v != vertices );

   /* Draw the polygon. */
   printf("\n%%Polygon:\n");
   printf("newpath\n");
   v = vertices;
   printf("%d\t%d\tmoveto\n", v->v[X], v->v[Y] );
   v = v->next;
   do {                                 
      printf("%d\t%d\tlineto\n", v->v[X], v->v[Y] );
      v = v->next;
   } while ( v != vertices );
   printf("closepath stroke\n");

}

void	PrintDiagonal( tVertex a, tVertex b )
{
   printf("%%Diagonal: (%d,%d)\n", a->vnum, b->vnum );
   printf("%d\t%d\tmoveto\n", a->v[X], a->v[Y] );
   printf("%d\t%d\tlineto\n", b->v[X], b->v[Y] );
}

int	AreaPoly2( void )
{
   int     sum = 0;
   tVertex p, a;

   p = vertices;   /* Fixed. */
   a = p->next;    /* Moving. */
   do {
      sum += Area2( p->v, a->v, a->next->v );
      a = a->next;
   } while ( a->next != vertices );
   return sum;

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
