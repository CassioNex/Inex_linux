/*
* fingerapi.h
*/

#ifndef FINGERAPI_H
#define FINGERAPI_H

#ifdef __cplusplus
extern "C" {
#endif

#define BOOL				int 
#define TRUE				1
#define FALSE				0
#define BYTE				unsigned char
#define UINT				unsigned int
#define DWORD				unsigned long
#define USHORT				unsigned short

#define MAX_FEATUREVECT_LEN			512

#define HIGH_LEVEL					2
#define MEDIUM_LEVEL				1
#define LOW_LEVEL					0

#define ERR_NO						 0	
#define ERR_OK						 1
#define	ERR_GENERAL_ERROR			-1

#define	ERR_MATCH_FAILED			-11
#define	ERR_CAN_NOT_ALLOC_MEMORY	-21
#define	ERR_VECT_FAILED			    -31  
#define	ERR_INVALID_IMAGESIZE		-41

#define MAX_IMG_WIDTH			200 //640
#define MAX_IMG_HEIGHT			200 //640
#define MAX_IMG_SIZE			MAX_IMG_WIDTH*MAX_IMG_HEIGHT

#define MAX_HEADER_SIZE			8

#define VERSION_NUMBER			0x40

#define END_MINUTIA             0
#define BIF_MINUTIA             1

#define MAXMINUTIATAGNUM        300
#define MAX_MINUTIA_NUMBER      50 
#define MIN_MINUTIA_NUMBER      3

#define CORE		            1
#define DELTA			        0

#define MAX_CORE_NUMBER         4
#define MAX_SP_NUMBER           4 //30

#define BLOCK_SIZE              16
#define MAX_BLOCK_ROW           18
#define MAX_BLOCK_COL           16
#define MAX_BLOCK_NUMBER        (MAX_BLOCK_ROW*MAX_BLOCK_COL)

#define	MAX_BAR_NUM				600 //1000
#define	MAX_SEARCH_BAR_NUM		500

#define	MATCH_TH_HIGH			100
#define	MATCH_TH_MEDIUM			50
#define	MATCH_TH_LOW			40

#define	ANGLE                   10
#define	LENGTH                  13

#define MAX_MAINLINE_NUM		100
#define MAX_POINTNUM_OF_POLYGON 50

#define MP_SIZE					6

#define MAX_SHIFT_X_SIZE        (MAX_IMG_WIDTH*2)
#define MAX_SHIFT_Y_SIZE        (MAX_IMG_HEIGHT*2)
#define MAX_SHIFT_SIZE          (MAX_SHIFT_X_SIZE)

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

typedef struct
{
	short x, y, dir;
	BYTE score, kind;
} TAG_POINT;

typedef struct {
	short nNumber;
	TAG_POINT item[MAXMINUTIATAGNUM];
} REALPVECT, *LPREALPVECT;

typedef struct
{
	short nNumber;
	short nX[MAX_SP_NUMBER];
	short nY[MAX_SP_NUMBER];
	short nOrient[MAX_SP_NUMBER];
	short nType[MAX_SP_NUMBER];
} SINGULAR;

typedef struct {
	short points_x[4];
	short points_y[4];
} TYPELINE;

typedef struct
{
	short x, y;
	BYTE dir, kind;
} COREITEMEX;
typedef struct {
	BYTE nNumber;
	COREITEMEX item[MAX_CORE_NUMBER];
} COREVECTEX, *LPCOREVECTEX;

typedef struct
{
	short x, y;
	BYTE dir;
	BYTE curv;
	BYTE score;
	BYTE kind;
} MP_POINT;

typedef struct {
	BYTE nNumber;
	MP_POINT item[MAX_MINUTIA_NUMBER];
} MPVECTEX, *LPMPVECTEX;

typedef struct {
	BYTE nCol;
	BYTE nRow;
	BYTE Data[MAX_BLOCK_NUMBER];
} BLOCKVECT, *LPBLOCKVECT;

typedef struct {
	BYTE	nVersion;
	BYTE	nHeaderLength;
	short	nImageX;
	short	nImageY;
} HEADERVECT, *LPHEADERVECT;

typedef struct
{
	HEADERVECT	nHeader;
	BYTE		nType;
	BYTE		nRidgeDensity;
	BYTE		nFrequency;
	TYPELINE	MainLine;
	BLOCKVECT	Block;
	COREVECTEX	Core;
	MPVECTEX	Mp;
} FPVECTEX, *LPFPVECTEX;

typedef struct
{
	short nLen, nSlope;
	short nDiff1, nDiff2;
	short nID1, nID2;
} MATCH_TAG;

typedef struct
{
	short nNumber;
	MATCH_TAG item[MAX_BAR_NUM];
} BARVECT;

typedef struct
{
	short nLen,	nID1, nID2;
} BARTEMP;

typedef struct
{
	short score, sid, fid;
} PAIRBAR;

typedef struct
{
	short nNumber;
	short nRot, nXOffset, nYOffset;
	short nXC, nYC;
	int nSearchID[MAX_MINUTIA_NUMBER];
	int nFileID[MAX_MINUTIA_NUMBER];
} PAIRVECT;

typedef struct
{
	short orderNum;
	short x[MAX_POINTNUM_OF_POLYGON];
	short y[MAX_POINTNUM_OF_POLYGON];
} POLYGON;

typedef struct {
	BYTE	nNumbers[4];
	short	points_x[4][MAX_MAINLINE_NUM];
	short	points_y[4][MAX_MAINLINE_NUM];
} MAINLINE;

typedef struct
{
	short x, y, dir, val;
} TEMPCORE;

typedef struct
{
	BYTE Data[MAX_FEATUREVECT_LEN];
} FPFEATURE, *LPFPFEATURE;

int finger_match(BYTE* pFeature1,BYTE* pFeature2,int securitylevel);

int finger_search(BYTE* pFeature,BYTE* pDBFeature,int nDBSize,int securitylevel);

int create_template(BYTE* pImage,int nWidth,int nHeight,BYTE* pFeature);

#ifdef __cplusplus
}
#endif

#endif // FINGERAPI_H
 
