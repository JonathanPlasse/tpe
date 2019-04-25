/***************************************************************************
*	Mathlib.h: définition des fonctions de la Mathlib.
****************************************************************************/
double  **dmatrice( int  p , int q );
double *dvect( int  p );
void  Detruitdvect( double *v ) ;
void  Detruitdmatrice( double **m, int p );
void mxv( double **m , double *v1 , double *v2 , int p , int q );
void mxm( double  **m1, double  **m2 , double  **m3 , int p , int q ,int r );
int pinvGreville(double **A,int m, int n,double **B);
