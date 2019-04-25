#define TPE_NB_OBJET		2  /* Nb d'objet �detecter (cible desiree : 0 cible courante : 1 */
#define TPE_NB_MARQUEUR		1  /* Nombre de marqueurs dans l'objet */ 
#define TPE_SEUIL		100
#define TPE_MARGE		20

/* Estimation de la profondeur de la cible dans sa position de reference */
#define Z_EST			0.9

/* Gain de la commande */
#define GAIN_T			0.1
#define GAIN_R			0.1

/* Param�res intrinseque de la camera */
#define ALPHA_U			1000.0
#define ALPHA_V			1000.0

#define U_0			380.0
#define V_0			285.0

#define TPE_PAS			0.07 /*pas de la cible*/


/* structure contenant la description de la tache */
struct cible_t	{
	
	int aoi_x0; /*origine de l'imagette haut gauche */
	int aoi_y0; /*origine de l'imagette haut gauche */
	
	int aoi_dx; /* largeur de l'imagette */
	int aoi_dy; /* hauteur de l'imagette haut gauche */
	
	int nb_pix; /*Nb de pixels detectes */ 
	
	double cx,cy; /*coordonnees du centre de gravite de la tache */
	
	int nb_marque; /* nombre de marqueurs initialis�s de l'objet */
	
	double m00;
	double m10,m01;
	double m20,m11,m02;
	
	double mu20,mu11,mu02;
	 
	};
	
/* structure contenant les infos de l'image courante �traiter */

struct image_t	{
	
	unsigned short width; /* largeur de l'image */
	unsigned short height; /* hauteur de l'image */
	unsigned char *buf; /* pointeur sur l'image */
	unsigned char **coord; /*pointeur sur les coordonn�s de l'image */
	unsigned char *buf_save; /* pointeur sur l'image */
	unsigned char **coord_save; /*pointeur sur les coordonn�s de l'image */
	
	};
	

struct tpe_t		{
	
	struct image_t im; /* image courante */
	struct cible_t cible[TPE_NB_OBJET]; /* information sur la cible de reference : 0 et le cible courante : 1*/
	double info_image[TPE_NB_OBJET][6]; /*informations visuelles reconstruite �partir des moments de l'objet */
	int nb_objet; /* nombre d'objet initialis�*/
	};

int Update_aoi ( struct tpe_t *tpe, int obj ); /* fonction de mise �jour de l'imagette autour de l'objet */
int Moment(struct tpe_t *tpe, int obj );/* calcul des moments de la cible */
int Update_mesure(struct tpe_t *tpe, int obj_des, int obj_cur); /*fonction de mise �jour des indices visuelles utilis� pour la commande*/
int Commande ( struct tpe_t tpe, double *control); /* fonction de calcul de la commande */


int Init_vision ( struct image_t *image, short w, short h ); 
int Init_cible ( struct cible_t *cible,  int *pos, int *siz  );
int rgb2gray ( unsigned char ** src, short w, short h, unsigned char **dst);
int Detruit_vision ( struct image_t image );

int SaveFile( char* name, int width, int height, unsigned char **pt );
int SavePartFile( char* name, int x,int y, int dx, int dy, unsigned char **pt );

int ReadPGM (const char *imageFileName, struct image_t *image);
