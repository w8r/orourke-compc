/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the
explanation in that book.

Compile:    gcc -o inhedron inhedron.c -lm (or simply: make)
Run (e.g.): inhedron < i.8

Written by Joseph O'Rourke, with contributions by Min Xu.
Last modified: April 1998
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/
#include <stdio.h>
#include <math.h>
#define EXIT_FAILURE 1
#define X 0
#define Y 1
#define Z 2
#define MAX_INT   2147483647 
typedef enum { FALSE, TRUE } bool;

#define DIM 3                  /* Dimension of points */
typedef int    tPointi[DIM];   /* Type integer point */
typedef double tPointd[DIM];   /* Type double point */
#define PMAX 10000             /* Max # of pts */
tPointi Vertices[PMAX];        /* All the points */
tPointi Faces[PMAX];           /* Each triangle face is 3 indices */
int check = 0;
tPointi Box[PMAX][2];          /* Box around each face */

/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/
char 	InPolyhedron( int F, tPointi q, tPointi bmin, tPointi bmax, int radius );
char    SegPlaneInt( tPointi Triangle, tPointi q, tPointi r, tPointd p, int *m )
;
int     PlaneCoeff( tPointi T, tPointd N, double *D );
void    Assigndi( tPointd p, tPointi a );
int     ReadVertices( void );
int     ReadFaces( void );
void    NormalVec( tPointi q, tPointi b, tPointi c, tPointd N );
double  Dot( tPointi q, tPointd d );
void    SubVec( tPointi q, tPointi b, tPointi c );
char    InTri3D( tPointi T, int m, tPointi p );
char    InTri2D( tPointi Tp[3], tPointi pp );
int     AreaSign( tPointi q, tPointi b, tPointi c );
char    SegTriInt( tPointi Triangle, tPointi q, tPointi r, tPointd p );
char    InPlane( tPointi Triangle, int m, tPointi q, tPointi r, tPointd p);
int     VolumeSign( tPointi a, tPointi b, tPointi c, tPointi d );
char    SegTriCross( tPointi Triangle, tPointi q, tPointi r );
int  	ComputeBox( int F, tPointi bmin, tPointi bmax );
void 	RandomRay( tPointi ray, int radius );
void 	AddVec( tPointi q, tPointi ray );
int  	InBox( tPointi q, tPointi bmin, tPointi bmax );
char 	BoxTest ( int n, tPointi a, tPointi b );
void 	PrintPoint( tPointi q );
int	irint( double x);
/*-------------------------------------------------------------------*/

main()
{
  int n, F, i;
  tPointi q, bmin, bmax;
  int radius;

  srandom( (int) time( (long *) 0 ) ); 

  n = ReadVertices();
  F = ReadFaces();

  /* Initialize the bounding box */
  for ( i = 0; i < DIM; i++ )
    bmin[i] = bmax[i] = Vertices[0][i];
  radius = ComputeBox( n, bmin, bmax );
  printf("radius=%d\n", radius);

  while( scanf( "%d %d %d", &q[X], &q[Y], &q[Z] ) != EOF ) {
    printf( "\n----------->q = %d %d %d\n", 
       q[X], q[Y], q[Z] );
    printf( "In = %c\n", 
       InPolyhedron( F, q, bmin, bmax, radius ) );
  }
}

/*
  This function returns a char:
    'V': the query point a coincides with a Vertex of polyhedron P.
    'E': the query point a is in the relative interior of an Edge of polyhedron P.
    'F': the query point a is in the relative interior of a Face of polyhedron P.
    'i': the query point a is strictly interior to polyhedron P.
    'o': the query point a is strictly exterior to( or outside of) polyhedron P.
*/
char InPolyhedron( int F, tPointi q, tPointi bmin, tPointi bmax, int radius )
{
   tPointi r;  /* Ray endpoint. */
   tPointd p;  /* Intersection point; not used. */
   int f, k = 0, crossings = 0;
   char code = '?';
 
   /* If query point is outside bounding box, finished. */
   if ( !InBox( q, bmin, bmax ) )
      return 'o';
  
   LOOP:
   while( k++ < F ) {
      crossings = 0;
  
      RandomRay( r, radius ); 
      AddVec( q, r ); 
      printf("Ray endpoint: (%d,%d,%d)\n", r[0],r[1],r[2] );
  
      for ( f = 0; f < F; f++ ) {  /* Begin check each face */
         if ( BoxTest( f, q, r ) == '0' ) {
              code = '0';
              printf("BoxTest = 0!\n");
         }
         else code = SegTriInt( Faces[f], q, r, p );
         printf( "Face = %d: BoxTest/SegTriInt returns %c\n\n", f, code );

         /* If ray is degenerate, then goto outer while to generate another. */
         if ( code == 'p' || code == 'v' || code == 'e' ) {
            printf("Degenerate ray\n");
            goto LOOP;
         }
   
         /* If ray hits face at interior point, increment crossings. */
         else if ( code == 'f' ) {
            crossings++;
            printf( "crossings = %d\n", crossings );
         }

         /* If query endpoint q sits on a V/E/F, return that code. */
         else if ( code == 'V' || code == 'E' || code == 'F' )
            return( code );

         /* If ray misses triangle, do nothing. */
         else if ( code == '0' )
            ;

         else 
            fprintf( stderr, "Error, exit(EXIT_FAILURE)\n" ), exit(1);

      } /* End check each face */

      /* No degeneracies encountered: ray is generic, so finished. */
      break;

   } /* End while loop */
 
   printf( "Crossings = %d\n", crossings );
   /* q strictly interior to polyhedron iff an odd number of crossings. */
   if( ( crossings % 2 ) == 1 )
      return   'i';
   else return 'o';
}

int ComputeBox( int F, tPointi bmin, tPointi bmax )
{
  int i, j, k;
  double radius;
  
  for( i = 0; i < F; i++ )
    for( j = 0; j < DIM; j++ ) {
      if( Vertices[i][j] < bmin[j] )
	bmin[j] = Vertices[i][j];
      if( Vertices[i][j] > bmax[j] ) 
	bmax[j] = Vertices[i][j];
    }
  
  radius = sqrt( pow( (double)(bmax[X] - bmin[X]), 2.0 ) +
                 pow( (double)(bmax[Y] - bmin[Y]), 2.0 ) +
                 pow( (double)(bmax[Z] - bmin[Z]), 2.0 ) );
  printf("radius = %lf\n", radius);

  return irint( radius +1 ) + 1;
}

/* Return a random ray endpoint */
void RandomRay( tPointi ray, int radius )
{
  double x, y, z, w, t;

  /* Generate a random point on a sphere of radius 1. */
  /* the sphere is sliced at z, and a random point at angle t
     generated on the circle of intersection. */
  z = 2.0 * (double) random() / MAX_INT - 1.0;
  t = 2.0 * M_PI * (double) random() / MAX_INT;
  w = sqrt( 1 - z*z );
  x = w * cos( t );
  y = w * sin( t );
  
  ray[X] = irint ( radius * x );
  ray[Y] = irint ( radius * y );
  ray[Z] = irint ( radius * z );
  
  /*printf( "RandomRay returns %6d %6d %6d\n", ray[X], ray[Y], ray[Z] );*/
}

void AddVec( tPointi q, tPointi ray )
{
  int i;
  
  for( i = 0; i < DIM; i++ )
    ray[i] = q[i] + ray[i];
}

int InBox( tPointi q, tPointi bmin, tPointi bmax )
{
  int i;

  if( ( bmin[X] <= q[X] ) && ( q[X] <= bmax[X] ) &&
      ( bmin[Y] <= q[Y] ) && ( q[Y] <= bmax[Y] ) &&
      ( bmin[Z] <= q[Z] ) && ( q[Z] <= bmax[Z] ) )
    return TRUE;
  return FALSE;
}
    

/*---------------------------------------------------------------------
    'p': The segment lies wholly within the plane.
    'q': The q endpoint is on the plane (but not 'p').
    'r': The r endpoint is on the plane (but not 'p').
    '0': The segment lies strictly to one side or the other of the plane.
    '1': The segement intersects the plane, and 'p' does not hold.
---------------------------------------------------------------------*/
char	SegPlaneInt( tPointi T, tPointi q, tPointi r, tPointd p, int *m)
{
    tPointd N; double D;
    tPointi rq;
    double num, denom, t;
    int i;

    *m = PlaneCoeff( T, N, &D );
    /*printf("m=%d; plane=(%lf,%lf,%lf,%lf)\n", m, N[X],N[Y],N[Z],D);*/
    num = D - Dot( q, N );
    SubVec( r, q, rq );
    denom = Dot( rq, N );
    /*printf("SegPlaneInt: num=%lf, denom=%lf\n", num, denom );*/

    if ( denom == 0.0 ) {  /* Segment is parallel to plane. */
       if ( num == 0.0 )   /* q is on plane. */
           return 'p';
       else
           return '0';
    }
    else
       t = num / denom;
    /*printf("SegPlaneInt: t=%lf \n", t );*/

    for( i = 0; i < DIM; i++ )
       p[i] = q[i] + t * ( r[i] - q[i] );

    if ( (0.0 < t) && (t < 1.0) )
         return '1';
    else if ( num == 0.0 )   /* t == 0 */
         return 'q';
    else if ( num == denom ) /* t == 1 */
         return 'r';
    else return '0';
}
/*---------------------------------------------------------------------
Computes N & D and returns index m of largest component.
---------------------------------------------------------------------*/
int	PlaneCoeff( tPointi T, tPointd N, double *D )
{
    int i;
    double t;              /* Temp storage */
    double biggest = 0.0;  /* Largest component of normal vector. */
    int m = 0;             /* Index of largest component. */

    NormalVec( Vertices[T[0]], Vertices[T[1]], Vertices[T[2]], N );
    /*printf("PlaneCoeff: N=(%lf,%lf,%lf)\n", N[X],N[Y],N[Z]);*/
    *D = Dot( Vertices[T[0]], N );

    /* Find the largest component of N. */
    for ( i = 0; i < DIM; i++ ) {
      t = fabs( N[i] );
      if ( t > biggest ) {
        biggest = t;
        m = i;
      }
    }
    return m;
}
/*---------------------------------------------------------------------
Reads in the number and coordinates of the vertices of a polyhedron
from stdin, and returns n, the number of vertices.
---------------------------------------------------------------------*/
int 	ReadVertices( void )
{
   int   i, n;

   do {
     scanf( "%d", &n );
     if ( n <= PMAX )
       break;
     printf("Error in read_vertex:  too many points; max is %d\n", PMAX);
   }
   while ( 1 );

   printf( "Polyhedron Vertices:\n" );
   printf( "  i:   x   y   z\n");
   for ( i = 0; i < n; i++ ) {
     scanf( "%d %d %d", &Vertices[i][X], &Vertices[i][Y], &Vertices[i][Z] );
     printf( "%3d:%4d%4d%4d\n", i, Vertices[i][X], Vertices[i][Y], Vertices[i][Z
] );
   }
   printf("n = %3d vertices read\n",n);
   putchar('\n');

   return n;
}

/*---------------------------------------------------------------------
a - b ==> c.
---------------------------------------------------------------------*/
void    SubVec( tPointi a, tPointi b, tPointi c )
{
   int i;

   for( i = 0; i < DIM; i++ )
      c[i] = a[i] - b[i];
}

/*---------------------------------------------------------------------
Returns the dot product of the two input vectors.
---------------------------------------------------------------------*/
double	Dot( tPointi a, tPointd b )
{
    int i;
    double sum = 0.0;

    for( i = 0; i < DIM; i++ )
       sum += a[i] * b[i];

    return  sum;
}

/*---------------------------------------------------------------------
Compute the cross product of (b-a)x(c-a) and place into N.
---------------------------------------------------------------------*/
void	NormalVec( tPointi a, tPointi b, tPointi c, tPointd N )
{
    N[X] = ( c[Z] - a[Z] ) * ( b[Y] - a[Y] ) -
           ( b[Z] - a[Z] ) * ( c[Y] - a[Y] );
    N[Y] = ( b[Z] - a[Z] ) * ( c[X] - a[X] ) -
           ( b[X] - a[X] ) * ( c[Z] - a[Z] );
    N[Z] = ( b[X] - a[X] ) * ( c[Y] - a[Y] ) -
           ( b[Y] - a[Y] ) * ( c[X] - a[X] );
}

/* Reads in the number of faces of the polyhedron and their indices from stdin,
    and returns the number n. */
int ReadFaces( void )
{
  int	i,j,k, n;
  int   w; /* temp storage for coordinate. */
  
  do {
    scanf( "%d", &n );
    if ( n <= PMAX )
      break;
    printf("Error in read_vertex:  too many points; max is %d\n", PMAX);
  }
  while ( 1 );

  printf( "Faces:\n" );
  printf( "  i   i0  i1  i2\n");
  for ( i = 0; i < n; i++ ) {
    scanf( "%d %d %d", &Faces[i][0], &Faces[i][1], &Faces[i][2] );
    printf( "%3d:%3d%4d%4d\n", i, Faces[i][0], Faces[i][1], Faces[i][2] );
    /* Compute bounding box. */
    /* Initialize to first vertex. */
    for ( j=0; j < 3; j++ ) {
       Box[i][0][j] = Vertices[ Faces[i][0] ][j];
       Box[i][1][j] = Vertices[ Faces[i][0] ][j];
    }
    /* Check k=1,2 vertices of face. */
    for ( k=1; k < 3; k++ )
    for ( j=0; j < 3; j++ ) {
       w = Vertices[ Faces[i][k] ][j];
       if ( w < Box[i][0][j] ) Box[i][0][j] = w;
       if ( w > Box[i][1][j] ) Box[i][1][j] = w;
    }
    /* printf("Bounding box: (%d,%d,%d);(%d,%d,%d)\n",
       Box[i][0][0],
       Box[i][0][1],
       Box[i][0][2],
       Box[i][1][0],
       Box[i][1][1],
       Box[i][1][2] );
    */
  }
  printf("n = %3d faces read\n",n);
  putchar('\n');

  return n;
}

/* Assumption: p lies in the plane containing T.
    Returns a char:
     'V': the query point p coincides with a Vertex of triangle T.
     'E': the query point p is in the relative interior of an Edge of triangle T.
     'F': the query point p is in the relative interior of a Face of triangle T.
     '0': the query point p does not intersect (misses) triangle T.
*/

char 	InTri3D( tPointi T, int m, tPointi p )
{
   int i;           /* Index for X,Y,Z           */
   int j;           /* Index for X,Y             */
   int k;           /* Index for triangle vertex */
   tPointi pp;      /* projected p */
   tPointi Tp[3];   /* projected T: three new vertices */

   /* Project out coordinate m in both p and the triangular face */
   j = 0;
   for ( i = 0; i < DIM; i++ ) {
     if ( i != m ) {    /* skip largest coordinate */
       pp[j] = p[i];
       for ( k = 0; k < 3; k++ )
	Tp[k][j] = Vertices[T[k]][i];
       j++;
     }
   }
   return( InTri2D( Tp, pp ) );
}

char 	InTri2D( tPointi Tp[3], tPointi pp )
{
   int area0, area1, area2;

   /* compute three AreaSign() values for pp w.r.t. each edge of the face in 2D */
   area0 = AreaSign( pp, Tp[0], Tp[1] );
   area1 = AreaSign( pp, Tp[1], Tp[2] );
   area2 = AreaSign( pp, Tp[2], Tp[0] );
   printf("area0=%d  area1=%d  area2=%d\n",area0,area1,area2);

   if ( ( area0 == 0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
        ( area1 == 0 ) && ( area0 > 0 ) && ( area2 > 0 ) ||
        ( area2 == 0 ) && ( area0 > 0 ) && ( area1 > 0 ) ) 
     return 'E';

   if ( ( area0 == 0 ) && ( area1 < 0 ) && ( area2 < 0 ) ||
        ( area1 == 0 ) && ( area0 < 0 ) && ( area2 < 0 ) ||
        ( area2 == 0 ) && ( area0 < 0 ) && ( area1 < 0 ) )
     return 'E';                 
   
   if ( ( area0 >  0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
        ( area0 <  0 ) && ( area1 < 0 ) && ( area2 < 0 ) )
     return 'F';

   if ( ( area0 == 0 ) && ( area1 == 0 ) && ( area2 == 0 ) )
     fprintf( stderr, "Error in InTriD\n" ), exit(EXIT_FAILURE);

   if ( ( area0 == 0 ) && ( area1 == 0 ) ||
        ( area0 == 0 ) && ( area2 == 0 ) ||
        ( area1 == 0 ) && ( area2 == 0 ) )
     return 'V';

   else  
     return '0';  
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

char    SegTriInt( tPointi T, tPointi q, tPointi r, tPointd p )
{
    int code = '?';
    int m = -1;

    code = SegPlaneInt( T, q, r, p, &m );
    printf("SegPlaneInt code=%c, m=%d; p=(%lf,%lf,%lf)\n", code,m,p[X],p[Y],p[Z]
);

    if      ( code == '0')
       return '0';
    else if ( code == 'q')
       return InTri3D( T, m, q );
    else if ( code == 'r')
       return InTri3D( T, m, r );
    else if ( code == 'p' )
       return InPlane( T, m, q, r, p );
    else if ( code == '1' )
       return SegTriCross( T, q, r );
    else /* Error */
       return code;
}

char	InPlane( tPointi T, int m, tPointi q, tPointi r, tPointd p)
{
    /* NOT IMPLEMENTED */
    return 'p';
}

/*---------------------------------------------------------------------
The signed volumes of three tetrahedra are computed, determined
by the segment qr, and each edge of the triangle.  
Returns a char:
   'v': the open segment includes a vertex of T.
   'e': the open segment includes a point in the relative interior of an edge
   of T.
   'f': the open segment includes a point in the relative interior of a face
   of T.
   '0': the open segment does not intersect triangle T.
---------------------------------------------------------------------*/

char SegTriCross( tPointi T, tPointi q, tPointi r )
{
   int vol0, vol1, vol2;
   
   vol0 = VolumeSign( q, Vertices[ T[0] ], Vertices[ T[1] ], r ); 
   vol1 = VolumeSign( q, Vertices[ T[1] ], Vertices[ T[2] ], r ); 
   vol2 = VolumeSign( q, Vertices[ T[2] ], Vertices[ T[0] ], r );
 
   printf( "SegTriCross:  vol0 = %d; vol1 = %d; vol2 = %d\n", 
      vol0, vol1, vol2 ); 
     
   /* Same sign: segment intersects interior of triangle. */
   if ( ( ( vol0 > 0 ) && ( vol1 > 0 ) && ( vol2 > 0 ) ) || 
        ( ( vol0 < 0 ) && ( vol1 < 0 ) && ( vol2 < 0 ) ) )
      return 'f';
   
   /* Opposite sign: no intersection between segment and triangle */
   if ( ( ( vol0 > 0 ) || ( vol1 > 0 ) || ( vol2 > 0 ) ) &&
        ( ( vol0 < 0 ) || ( vol1 < 0 ) || ( vol2 < 0 ) ) )
      return '0';

   else if ( ( vol0 == 0 ) && ( vol1 == 0 ) && ( vol2 == 0 ) )
     fprintf( stderr, "Error 1 in SegTriCross\n" ), exit(EXIT_FAILURE);
   
   /* Two zeros: segment intersects vertex. */
   else if ( ( ( vol0 == 0 ) && ( vol1 == 0 ) ) || 
             ( ( vol0 == 0 ) && ( vol2 == 0 ) ) || 
             ( ( vol1 == 0 ) && ( vol2 == 0 ) ) )
      return 'v';

   /* One zero: segment intersects edge. */
   else if ( ( vol0 == 0 ) || ( vol1 == 0 ) || ( vol2 == 0 ) )
      return 'e';
   
   else
     fprintf( stderr, "Error 2 in SegTriCross\n" ), exit(EXIT_FAILURE);
}

int 	VolumeSign( tPointi a, tPointi b, tPointi c, tPointi d )
{ 
   double vol;
   double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
   double bxdx, bydy, bzdz, cxdx, cydy, czdz;

   ax = a[X];
   ay = a[Y];
   az = a[Z];
   bx = b[X];
   by = b[Y];
   bz = b[Z];
   cx = c[X]; 
   cy = c[Y];
   cz = c[Z];
   dx = d[X];
   dy = d[Y];
   dz = d[Z];

   bxdx=bx-dx;
   bydy=by-dy;
   bzdz=bz-dz;
   cxdx=cx-dx;
   cydy=cy-dy;
   czdz=cz-dz;
   vol =   (az-dz) * (bxdx*cydy - bydy*cxdx)
         + (ay-dy) * (bzdz*cxdx - bxdx*czdz)
         + (ax-dx) * (bydy*czdz - bzdz*cydy);


   /* The volume should be an integer. */
   if      ( vol > 0.5 )   return  1;
   else if ( vol < -0.5 )  return -1;
   else                    return  0;
}

/*
  This function returns a char:
    '0': the segment [ab] does not intersect (completely misses) the 
         bounding box surrounding the n-th triangle T.  It lies
         strictly to one side of one of the six supporting planes.
    '?': status unknown: the segment may or may not intersect T.
*/
char BoxTest ( int n, tPointi a, tPointi b )
{
   int i; /* Coordinate index */
   int w;

   for ( i=0; i < DIM; i++ ) {
       w = Box[ n ][0][i]; /* min: lower left */
       if ( (a[i] < w) && (b[i] < w) ) return '0';
       w = Box[ n ][1][i]; /* max: upper right */
       if ( (a[i] > w) && (b[i] > w) ) return '0';
   }
   return '?';
}


/* irint not available in some libraries, so... */

int	irint( double x )
{
	return (int) rint( x );
}
