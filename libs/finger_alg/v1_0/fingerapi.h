/*
* fingerapi.h
*/

#ifndef FINGERAPI_H
#define FINGERAPI_H

#ifdef __cplusplus
extern "C" {
#endif

#define HIGH_LEVEL					2
#define MEDIUM_LEVEL				1
#define LOW_LEVEL					0

#define ERR_NO						 0	
#define ERR_OK						 1
#define	ERR_GENERAL_ERROR			-1

#define	ERR_MATCH_FAILED			-11
#define	ERR_CAN_NOT_ALLOC_MEMORY	-21
#define	ERR_VECT_FAILED				-31  
#define	ERR_INVALID_IMAGESIZE		-41

#define BOOL				int 
#define TRUE				1
#define FALSE				0
#define BYTE				unsigned char
#define UINT				unsigned int
#define DWORD				unsigned long
#define USHORT				unsigned short

#define MAX_FEATUREVECT_LEN			512

typedef struct
{
	unsigned char Data[MAX_FEATUREVECT_LEN];
} MINUTIAVECT, *LPFPFEATURE;

#define CURRENT_VERSION			0x41

#define MAX_HEADER_SIZE			4

#define MAX_IMG_WIDTH			640
#define MAX_IMG_HEIGHT			640

#define END_MINUTIA             	0
#define BIF_MINUTIA             	1

#define MAXMINUTIATAGNUM        	300
#define MAX_MINUTIA_NUMBER      	50
#define MIN_MINUTIA_NUMBER      	3

#define MAX_CORE_NUMBER         	4
#define MAX_SP_NUMBER           	4
#define CORE		            	1
#define DELTA			        0

#define BLOCK_SIZE              	16
#define MAX_BLOCK_ROW           	20
#define MAX_BLOCK_COL           	17
#define MAX_BLOCK_NUMBER        	(MAX_BLOCK_ROW*MAX_BLOCK_COL)

#define END_POINT			1
#define BIF_POINT			3

#define MAX_SECTION_NUMBER      	50
#define NULLVALUE               	1000
#define	NO_STOP			        0
#define	EXCEPTION_STOP			(-1)
#define	LOOP_STOP		        (-2)

#define MAX_MAINLINE_NUM		100
#define MAX_POINTNUM_OF_POLYGON 	50

#define MP_SIZE_V22			5

#define MAX_TYPE_NUM			11

#define	MAX_BAR_NUM        		600
#define	MAX_SEARCH_BAR_NUM 		500

#define	MATCH_TH_HIGH			100
#define	MATCH_TH_MEDIUM			50
#define	MATCH_TH_LOW			40

#define	ANGLE                   	10
#define	LENGTH                  	13

#define	MAX_SAMEDIFF_NUM        	10
#define MAX_BAR_LENGTH			300

#define MAX_NEIGH_NUM 			8
#define MAX_CHUNK_NUM 			15

#define MAX_IDENT_NUM 			4

#define MAX_NUM_NEIGHBOR		16	


#define IANGLE_60(dir)		( ( (dir) >= 60 ) ? ( 120 - (dir) ) : ( (dir) ) )
#define IANGLE_120(dir)		( ( (dir) >= 120 ) ? ( 240 - (dir) ) : ( (dir) ) )
#define ANGLE_0(dir)		( ( (dir) < 0 ) ? ( (dir) + 240 ) : ( (dir) ) )
#define ANGLE_120(dir)		( ( (dir) >= 120 ) ? ( (dir) - 120 ) : ( (dir) ) )
#define ANGLE_240(dir)		( ( (dir) >= 240 ) ? ( (dir) - 240 ) : ( (dir) ) )
#define ANGLE_0_240(dir)	( ( (dir) >= 240 ) ? ( (dir) - 240 ) : ( (dir) < 0 ? ( (dir) + 240 ) : (dir) ) )

#define ROUND(x)		( ( (x) > 0 ) ? ( x + 8192 ) : ( (x) ) )


typedef struct
{
	short x, y;
} POINTTAG, *LPPOINTTAG;

typedef struct
{
	short nNumber;
	POINTTAG point[MAX_SECTION_NUMBER];
	BYTE nValue[MAX_SECTION_NUMBER];
} SECTION, *LPSECTION;

typedef struct
{
	short nNumber;
	short nX[MAX_SP_NUMBER];
	short nY[MAX_SP_NUMBER];
	short nOrient[MAX_SP_NUMBER];
	short nType[MAX_SP_NUMBER];
} SINGULAR;

typedef struct
{
	short x, y, dir;
	BYTE score, kind;
} REALMINUTIA;

typedef struct {
	short nNumber;
	REALMINUTIA item[MAXMINUTIATAGNUM];
} REALPVECT, *LPREALPVECT;

typedef struct {
	BYTE nNumbers[4];
	short points_x[4][MAX_MAINLINE_NUM];
	short points_y[4][MAX_MAINLINE_NUM];
} MAINLINE;

typedef struct {
	short points_x[4];
	short points_y[4];
} DESMAINLINE;

typedef struct
{
	short x, y;
	BYTE dir, kind;
} COREITEMEX;

typedef struct {
	char nNumber;
	COREITEMEX item[MAX_CORE_NUMBER];
} COREVECTEX, *LPCOREVECTEX;


typedef struct
{
	BYTE radius, alpha, delta;
} NEIGHBOR;
typedef struct
{
	BYTE nCount;
	NEIGHBOR item[MAX_NUM_NEIGHBOR];
} NEIGHBORS;

typedef struct
{
	BYTE score, sid, fid;
} PAIRLIST;


typedef struct
{
	short x, y;
	BYTE dir, curv, score, kind;
} MINUTIAEX; 
 
typedef struct {
	char nNumber;
	BYTE quality;
	MINUTIAEX item[MAX_MINUTIA_NUMBER];
} MPVECTEX, *LPMPVECTEX;

typedef struct {
	BYTE nCol;
	BYTE nRow;
	BYTE Data[MAX_BLOCK_NUMBER];
} BLOCKVECT, *LPBLOCKVECT;

typedef struct {
	BYTE	nVersion;
	BYTE	nHeaderLength;
	short nImageX;
	short nImageY;
} HEADERVECT, *LPHEADERVECT;

typedef struct
{
	HEADERVECT	nHeader;
	BYTE nType;
	BYTE nRidgeDensity;
	BYTE nFrequency;
	DESMAINLINE	MainLine;
	BLOCKVECT Block;
	COREVECTEX Core;
	MPVECTEX Mp;
} FPVECTEX, *LPFPVECTEX;


typedef struct
{
	short nLen;
	short nSlope;
	short nDiff1;
	short nDiff2;
	short nID1;
	short nID2;
} BAR;

typedef struct
{
	short nNumber;
	BAR item[MAX_BAR_NUM];
} BARVECT;

typedef struct
{
	short nLen, nID1, nID2;
} BARTEMP;

typedef struct
{
	short score, sid, fid;
} PAIRBAR;

typedef struct
{
	short nNumber;
	short nRot;
	short nXOffset, nYOffset;
	short nXC, nYC;
	short nSearchID[MAX_MINUTIA_NUMBER];
	short nFileID[MAX_MINUTIA_NUMBER];
} PAIRVECT;


typedef struct
{
	short orderNum;
	short x[MAX_POINTNUM_OF_POLYGON];
	short y[MAX_POINTNUM_OF_POLYGON];
} POLYGON;

typedef struct
{
	short x, y, dir, val;
} TEMPCORE;


typedef struct
{
	BYTE cdir;
	int radius[MAX_NEIGH_NUM];
	BYTE angle[MAX_NEIGH_NUM];
	BYTE delta[MAX_NEIGH_NUM];
} DATACHUNK;

typedef struct
{
	char number;
	DATACHUNK chunk[MAX_CHUNK_NUM];
} FPTEMPLATE;

int finger_match(BYTE* pFeature1,BYTE* pFeature2,int securitylevel);

int finger_search(BYTE* pFeature,BYTE* pDBFeature,int nDBSize,int securitylevel);

int create_template(BYTE* pImage,int nWidth,int nHeight,BYTE* pFeature); 

#ifdef __cplusplus
}
#endif

#endif // FINGERAPI_H
