#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathlib.h"

#include "tpe.h"

#define NB_FRAME 199

//const char *imagefront = "C:\\Users\\root\\Documents\\Develop\\TPE\\TPE-3D\\3D\\im3D-";	/* Prefixe pour les noms des images */
//const char *imageend = ".pgm";		/* Suffixe */

#ifdef LINUX
const char *imagepath = "../Seq"; /* chemin vers la sequence d'image */
const char *imagefront = "testf";	/* Prefixe pour les noms des images */
const char *imageend = ".pgm";		/* Suffixe */
#else
const char *imagepath = "../Code";
const char *imagefront = "testf";	/* Prefixe pour les noms des images */
const char *imageend = ".pgm";		/* Suffixe */

#endif

double z_des = 200.0/1000.0;


/****************************************************************************************
*
*	Update_aoi : mise  �  jour de la position de l'AOI
*
*****************************************************************************************/
int Update_aoi ( struct tpe_t *tpe, int obj)	{
	double det,alpha,a1,a2;

	cible_t* c;

	c = &tpe->cible[obj];

	// Calcul les caracteristique de l'ellipse
 	det = sqrt((c->mu02 - c->mu20)*(c->mu02 - c->mu20) + 4 * c->mu11 * c->mu11);
	alpha = atan2(c->mu02 - c->mu20 + det, 2. * c->mu11);
	a1 = sqrt(2*(c->mu02 + c->mu20 + det) / c->m00);
	a2 = sqrt(2*(c->mu02 + c->mu20 - det) / c->m00);

	// Calcul la largeur et la longueur de l'imagette
	// On considere que l'ellipse est un rectangle
	c->aoi_dx = 2 * a1 * fabs(cos(alpha)) + 2 * a2 * fabs(sin(alpha)) + 2 * TPE_MARGE;
	c->aoi_dy = 2 * a1 * fabs(sin(alpha)) + 2 * a2 * fabs(cos(alpha)) + 2 * TPE_MARGE;

	// Calcul de l'origine de l'imagette
	c->aoi_x0 = c->cx - c->aoi_dx / 2;
	c->aoi_y0 = c->cy - c->aoi_dy / 2;

	// Recale si l'origine de l'imagette n'est pas dans l'image
	if (c->aoi_x0 < 0) {
		c->aoi_dx += c->aoi_x0;
		c->aoi_x0 = 0;
	}
	if (c->aoi_y0 < 0) {
		c->aoi_dy += c->aoi_y0;
		c->aoi_y0 = 0;
	}

	// Recale si l'imagette sort de l'image
	if (c->aoi_x0 + c->aoi_dx > tpe->im.width)
	{
		c->aoi_dx = tpe->im.width - c->aoi_x0;
	}
	if (c->aoi_y0 + c->aoi_dy > tpe->im.height)
	{
		c->aoi_dy = tpe->im.height - c->aoi_y0;
	}

	return 0;
}

/****************************************************************************************
*
*	Barycentre : calcul des centres de gravite des marqueurs d'une cible
*
*****************************************************************************************/
int Moment (struct tpe_t *tpe, int obj  ){
	int i, j;

	struct cible_t* c;

	c = &tpe->cible[obj];

	// Initialise les moments
	c->m00 = 0;
	c->m10 = 0;
	c->m01 = 0;
	c->m20 = 0;
	c->m11 = 0;
	c->m02 = 0;

	// Calcul les moments
	for (i = c->aoi_x0 ; i < c->aoi_x0 + c->aoi_dx ; i++) {
		for (j = c->aoi_y0 ; j < c->aoi_y0 + c->aoi_dy ; j++) {
			if (tpe->im.coord[j][i] < TPE_SEUIL) {
				c->m00 += 1;
				c->m10 += i;
				c->m01 += j;
				c->m20 += i*i;
				c->m11 += i*j;
				c->m02 += j*j;
			}
		}
	}
	// Calcul les centroids
	c->cx = c->m10 / c->m00;
	c->cy = c->m01 / c->m00;

	// Calcul les moments centres d'ordre 2
	c->mu20 = c->m20 - c->m10*c->m10/c->m00;
	c->mu11 = c->m11 - c->m10*c->m01/c->m00;
	c->mu02 = c->m02 - c->m01*c->m01/c->m00;

	return 0 ;
}

/****************************************************************************************
*
*	Update_mesure : reconstruction des indices visuels
*
*****************************************************************************************/
int Update_mesure(struct tpe_t *tpe, int obj_des, int obj_cur) 	{
	double a1,a2,a3,a4;
	double u_cur,v_cur,cx_cur,cy_cur,u_ref,v_ref,cx_ref,cy_ref,alpha,det,a1carre,a2carre,z;
	double **M;
	double **M_inv;
	double a_vec[4];
	double u_vec[4];

	struct cible_t* c_des;
	struct cible_t* c_cur;
	double* info_image;

	c_des = &tpe->cible[obj_des];
	c_cur = &tpe->cible[obj_cur];
	info_image = tpe->info_image[1];

	// Calcul des caracteristiques de l'ellipse sur l'image desiree
	det = sqrt((c_des->mu20 - c_des->mu02)*(c_des->mu20 - c_des->mu02) + 4*c_des->mu11*c_des->mu11);
	alpha = atan2(c_des->mu02 - c_des->mu20 + det, 2. * c_des->mu11);
	a1carre=2*(c_des->mu02 + c_des->mu20+det)/c_des->m00;
	a2carre=2*(c_des->mu02 + c_des->mu20-det)/c_des->m00;

	// Definition des points caracteristiques de l'ellipse sur l'image desiree
	// Le centre de l'ellipse
	cy_ref= c_des->cy;
	cx_ref= c_des->cx;
	// Un sommet de l'ellipse
	u_ref = cx_ref+sqrt(a1carre)*cos(alpha);
	v_ref = cy_ref+sqrt(a1carre)*sin(alpha);

	// Calcul des caracteristiques de l'ellipse sur l'image courante
 	det = sqrt((c_cur->mu20 - c_cur->mu02)*(c_cur->mu20 - c_cur->mu02) + 4*c_cur->mu11*c_cur->mu11);
	alpha = atan2(c_cur->mu02 - c_cur->mu20 + det, 2. * c_cur->mu11);
	a1carre=2*(c_cur->mu02 + c_cur->mu20+det)/c_cur->m00;
	a2carre=2*(c_cur->mu02 + c_cur->mu20-det)/c_cur->m00;

	// Definition des points caracteristiques de l'ellipse sur l'image courante
	cx_cur=c_cur->cx;
	cy_cur=c_cur->cy;
	u_cur=cx_cur+sqrt(a1carre)*cos(alpha);
	v_cur=cy_cur+sqrt(a1carre)*sin(alpha);

	// Calcul des solutions du systeme
	M = dmatrice(4,4);
	M_inv = dmatrice(4,4);

	M[0][0] = cx_ref-U_0; M[0][1] = -(cy_ref-V_0); M[0][2] = ALPHA_U; M[0][3] = 0.0;
	M[1][0] = cy_ref-V_0; M[1][1] = cx_ref-U_0; M[1][2] = 0.0; M[1][3] = ALPHA_V;
	M[2][0] = u_ref-U_0; M[2][1] = -(v_ref-V_0); M[2][2] = ALPHA_U; M[2][3] = 0.0;
	M[3][0] = v_ref-V_0; M[3][1] = u_ref-U_0; M[3][2] = 0.0; M[3][3] = ALPHA_V;

	pinvGreville(M,4,4,M_inv);

	u_vec[0] = cx_cur - U_0;
	u_vec[1] = cy_cur - V_0;
	u_vec[2] = u_cur - U_0;
	u_vec[3] = v_cur - V_0;

	mxv(M_inv,u_vec,a_vec,4,4);

	Detruitdmatrice(M,4);
	Detruitdmatrice(M_inv,4);

	a1 = a_vec[0];
	a2 = a_vec[1];
	a3 = a_vec[2];
	a4 = a_vec[3];

	z = Z_EST / sqrt(a1*a1 + a2*a2);
	tpe->info_image[1][2] = z - Z_EST; //tz
	tpe->info_image[1][0] = z * a3; //tx
	tpe->info_image[1][1] = z * a4; //ty
	tpe->info_image[1][3] = atan2(a2, a1); //alpha

	return 0;
}
/****************************************************************************************
*
*	Commande : calcul de la commande du robot dans l'espace de la camera
*
*****************************************************************************************/
int Commande ( struct tpe_t tpe, double *control)	{
	double* info_i;
	info_i = tpe.info_image[1];

	// Mise a jour de la commande
	control[0] = info_i[0] * GAIN_T;
	control[1] = info_i[1] * GAIN_T;
	control[2] = info_i[2] * GAIN_T;
	control[3] = info_i[3] * GAIN_R;

	return 0;
}

/****************************************************************************************
*
*	Init_vision : initialisation de la structure image
*
*****************************************************************************************/
int Init_vision ( struct image_t *image, short w, short h )	{

	int i;

	image->width = w;
	image->height = h;

	/*initiamisation de l'image courante */
	if ( ! (image->buf = (unsigned char *)malloc(w*h) ) )	{
		printf("Unable to allocate memory for current image.\n");
		return 1;
		}

	if ( ! (image->coord = (unsigned char **)malloc(sizeof(unsigned char *) * h ) ) )	{
		printf("Unable to allocate memory for current image.\n");
		return 2;

		}
	else	{

		for (i=0; i< h; i++)
			image->coord[i] = &image->buf[i*w];

		}

	if ( ! (image->buf_save = (unsigned char *)malloc(w*h) ) )	{
		printf("Unable to allocate memory for current image.\n");
		return 1;
		}

	if ( ! (image->coord_save = (unsigned char **)malloc(sizeof(unsigned char *) * h ) ) )	{
		printf("Unable to allocate memory for current image.\n");
		return 2;

		}
	else	{

		for (i=0; i< h; i++)
			image->coord_save[i] = &image->buf_save[i*w];

		}

	return 0;

	}
/****************************************************************************************
*
*	Detruit_vision : Liberation de la memoire allouee pour la structure image
*
*****************************************************************************************/
int Detruit_vision ( struct image_t image )	{

	free(image.buf);
	free(image.coord);

	free(image.buf_save);
	free(image.coord_save);

	return 0;
	}

/****************************************************************************************
*
*	Init_vision : initialisation de la structure image
*
*****************************************************************************************/
int Init_cible ( struct cible_t *cible, int *pos, int *siz  )	{

	cible->aoi_x0 = pos[0];
	cible->aoi_y0 = pos[1];


	cible->aoi_dx = siz[0];
	cible->aoi_dy = siz[1];



	return 0;


	}

/****************************************************************************************
*
*	SaveFile : sauvegarde d'un image N&B
*
*****************************************************************************************/
int SaveFile( char* name, int width, int height, unsigned char **pt )	{

	FILE		*Stream;
	unsigned char	buf[width*height];
	int		i,j;

	if ( (Stream = fopen( name, "w" )) )	{

		/* Header PGM */

		fprintf( Stream, "P5\n%d %d\n255\n", width, height );

		/* Changement de l'ordre BGR en RGB */

		for ( i = 0; i < width  ; i ++ )	{
			for ( j=0; j <height; j++)	{

				buf[j*width + i] = pt[j][i];

				}
			}

		/* Sauvegarde */

		if ( fwrite( buf, 1, width * height , Stream ) != width * height  )	{
			fclose( Stream );
			return -1;
			}

		fclose( Stream );
		}
	else
		return -1;

	return 0;
	}

/****************************************************************************************
*
*	SavePartFile : sauvegarde d'un image N&B
*
*****************************************************************************************/
int SavePartFile( char* name, int x,int y, int dx, int dy, unsigned char **pt )	{

	FILE		*Stream;
	unsigned char	buf[dx*dy];
	int		i,j;

	if ( (Stream = fopen( name, "w" )) )	{

		/* Header PGM */

		fprintf( Stream, "P5\n%d %d\n255\n", dx, dy );

		/* Changement de l'ordre BGR en RGB */

		for ( i = x; i < x+dx  ; i ++ )	{
			for ( j = y; j <y+dy ; j++)	{

				buf[(j-y)*dx + (i-x)] = pt[j][i];

				}
			}

		/* Sauvegarde */

		if ( fwrite( buf, 1, dx * dy , Stream ) != dx * dy  )	{
			fclose( Stream );
			return -1;
			}

		fclose( Stream );
		}
	else
		return -1;

	return 0;
	}


/**********************************************************************
*		initialisation d'une image pgm noir&blanc	      *
*								      *
**********************************************************************/
int InitImage(int width, int height, unsigned char **pt, int value)	{

  int i,j;
  for (i = 0; i < width; i ++)
    for (j = 0; j < height; j ++)
      pt[j][i] = value;
  return 0;

}
/**********************************************************************
*		Lecture d'une image pgm noir&blanc		      *
*								      *
**********************************************************************/
int ReadPGM (const char *imageFileName, struct image_t *image)	{

	FILE          *fd;
	char          MagicValue[3];
	short int	      MaxGrey;
	long          nb_pix;

	if (( fd = fopen( imageFileName , "r" ) ) == NULL ) 	{
		printf("File %s does not exist\n", imageFileName);
		return(1);
		}

	fgets( MagicValue , 3 , fd );

	/*--- PGM BINARY (256 grey-levels) ---*/
	if ( strcmp( MagicValue , "P5" ) == 0 )	{
   		fscanf( fd , "%hd %hd\n%hd\n",&(image->width),&(image->height),&MaxGrey);
   		nb_pix = (image->width) * (image->height);
		Init_vision( image, image->width, image->height );
   		fread( image->buf , 1 , nb_pix , fd );
   		}

	fclose(fd);

	return(0);

	}
/**********************************************************************
*		Trace un point dans l'image pt 		      	      *
*								      *
**********************************************************************/
void DisplayPoint(int x, int y, int color, unsigned char **pt, int w, int h)	{

    if (x < 0)
      return;
    if (y < 0)
      return;
    if (x > w)
      return;
    if (y > h)
      return;

      pt[y][x] = color;



}

/**********************************************************************
*		Trace une croix a l'emplacement x,y dans l'image pt   *
*								      *
**********************************************************************/
void DisplayCross(int x, int y, int size, int color, unsigned char **pt, int w, int h)	{

    int i,j;
    if (x-size/2 < 0)
      return;
    if (y-size/2 < 0)
      return;
    if (x+size/2 > w)
      return;
    if (y+size/2 > h)
      return;

    for (i = x-size/2; i <= x+size/2; i++)
      	    pt[y][i] = color;


    for (j = y-size/2; j <= y+size/2; j++)
      	    pt[j][x] = color;




}
/**********************************************************************************
*		Trace un carre d'origine x,y et de taille sx,sy dans l'image pt   *
*								                  *
***********************************************************************************/

void DisplaySquare(int x, int y, int sx, int sy, int color, unsigned char **pt, int w, int h)	{

    int i,j;

    if (x < 0)
      return;
    if (y < 0)
      return;
    if (x+sx > w)
      return;
    if (y+sy > h)
      return;

    for (i = x; i <= x+sx; i++)	{
      	    pt[y][i] = color;
      	    pt[y+sy][i] = color;

    }


    for (j = y; j <= y+sy; j++)	{
      	    pt[j][x] = color;
      	    pt[j][x+sx] = color;
    }




}

/**********************************************************************
*		Programme principal		      *
*								      *
**********************************************************************/
int main( int argc, char *argv[] )	{

	//const char Name[255] = "im3D-1.pgm";
	//char NameSave[255] = "im3d-1_save.pgm";

  char Name[255] ;
  char NameSave[255];

  int i;
  double control[4];


	/* zone et taille de l'imagette �  traiter */
	// initialisation image

	int pos[2][2]={{40,290},{50,120}};

	int siz[2][2]={{650,170},{350,235}};


  int crossSize=10;
  double alpha, det, a2;



	struct tpe_t T;
	i=199;
	snprintf( Name, 255, "%s/In/%s%03d%s", 	imagepath,imagefront, i, imageend );
	snprintf( NameSave, 255, "%s/Out/%s%03d%s", 	imagepath,imagefront, i, imageend );

	ReadPGM(Name,&T.im);

	Init_cible ( &T.cible[0], pos[0], siz[0]  );
	Moment ( &T, 0);
	Update_aoi ( &T, 0);
	//det = sqrt((tpe->cible[0].mu20 - tpe->cible[0].mu02)*(tpe->cible[0].mu20 - tpe->cible[0].mu02) + 4*tpe->cible[0].mu11*tpe->cible[0].mu11);
	//alpha = atan( (tpe->cible[0].mu02 - tpe->cible[0].mu20 + det)/(2.0*tpe->cible[0].mu11) );

	// Update_mesure(&T, 0, 1);
	Detruit_vision( T.im );

	for (i=0; i < NB_FRAME; i++)	{

	// Attention fonction specifique LINUX : _snprintf semble �tre l'equivalent sous windows
    snprintf( Name, 255, "%s/In/%s%03d%s", 	imagepath,imagefront, i, imageend );
    snprintf( NameSave, 255, "%s/Out/%s%03d%s", 	imagepath,imagefront, i, imageend );

    ReadPGM(Name,&T.im);

		if ( i== 0){
		  Init_cible ( &T.cible[1], pos[1], siz[1]  );

		}


		Moment ( &T, 1);
		Update_aoi ( &T, 1);
		Update_mesure(&T,0,1);
		Commande ( T, control);

		#ifdef LINUX
		//printf(" moment image : %f %f %f %f %f %f \n",T.cible[1].m00,T.cible[1].m01,T.cible[1].m10,T.cible[1].m02,T.cible[1].m20,T.cible[1].m11);
		 printf("vitesse camera : %f %f %f %f \n",control[0],control[1],control[2],control[3]);
		#endif

		DisplayCross((int)T.cible[1].cx, (int) T.cible[1].cy, crossSize, 255 , T.im.coord,T.im.width,T.im.height);
		DisplaySquare((int)T.cible[1].aoi_x0, (int) T.cible[1].aoi_y0, T.cible[1].aoi_dx,T.cible[1].aoi_dy, 0 , T.im.coord,T.im.width,T.im.height);

		SaveFile(NameSave,T.im.width,T.im.height,T.im.coord);

		Detruit_vision( T.im );

	}


	return 0;

	}
