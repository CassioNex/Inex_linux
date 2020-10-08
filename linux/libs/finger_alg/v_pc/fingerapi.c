// fingerapi.c: Defines the functions for recognition.
//

#include <stdlib.h>
#include "math.h"
#include "memory.h"
#include "consttable.h"
#include "fingerapi.h"

FPVECTEX ViewFPEx[2];

int op_func_01(int ex,int ey,int sx,int sy)
{
	int nDiffX = abs(sx-ex), nDiffY = abs(sy-ey);
	int nAngle;

	while(nDiffX>=50 || nDiffY>=50) {
		nDiffX >>= 1; nDiffY >>= 1;
	}
	nAngle = _table_01[nDiffY*50+nDiffX];
	if(ex > sx) {
		if(ey < sy) nAngle = 240 - nAngle;
	} 
	else {
		if(ey > sy) nAngle = 120 - nAngle;
		else nAngle += 120;
	}
	if(nAngle < 0) nAngle += 240;
	if(nAngle >= 240) nAngle -= 240;

	return nAngle;
}

int op_func_02(int sqr_val)
{
	int tmp1, sqrt_val = 0x100, tmp2 = 0x100;
	
	if(sqr_val <= 0) return 0;
	if(sqr_val >= 0x40000) return 0x200;
	do {
		tmp1 = sqrt_val*sqrt_val;
		tmp2 >>= 1;
		if(sqr_val == tmp1) return sqrt_val;
		if(sqr_val > tmp1) sqrt_val += tmp2;
		else sqrt_val -= tmp2;
	} while(tmp2 > 1);
	tmp1 = sqrt_val*sqrt_val;
	if(tmp1 == sqr_val) return sqrt_val;
	if(tmp1 > sqr_val) {
		tmp2 = tmp1 - sqrt_val;
		if(tmp2 >= sqr_val) return sqrt_val-1;
	} else {
		tmp2 = tmp1 + sqrt_val;
		if(tmp2 < sqr_val) return sqrt_val+1;
	}
	return sqrt_val;
}

BOOL check_in_polygon(int point_x,int point_y,POLYGON *polygon,int limit)
{
	int i, i1, i2, num, x1, x2, x3, y1, y2, y3, lim2, nSign, mx0, mx1, my0, my1;

	if ( polygon->orderNum < 3 ) return FALSE;
	num = polygon->orderNum; nSign = 1;
	if ( limit < 0 ) nSign = -1;
	lim2 = limit*limit;
	x1 = polygon->x[0]; y1 = polygon->y[0];
	x3 = polygon->x[num-1]; y3 = polygon->y[num-1];

	for ( i=0; i<num; i++ ){
		if ( i+1 < num ){
			x2 = polygon->x[i+1]; y2 = polygon->y[i+1];
		}
		else{
			x2 = polygon->x[0]; y2 = polygon->y[0];
		}
		i1 = (y1-y2)*(point_x-x2) - (x1-x2)*(point_y-y2);
		if ( i1 == 0 ){
			if ( x1 < x2 ){	mx0 = x1; mx1 =	x2;	}
			else{ mx0 = x2; mx1 = x1; }
			if ( y1 < y2 ){	my0 = y1; my1 =	y2;	}
			else{ my0 = y2; my1 = y1; }
			if ( mx0 <= point_x && point_x <= mx1 && my0 <= point_y && point_y <= my1 ){
				if ( nSign > 0 ) return TRUE;
				return FALSE;
			}
		}
		i2 = (y1-y2)*(x3-x2) - (x1-x2)*(y3-y2);
		if ( (i1>0 && i2<0) || (i1<0 && i2>0) ){
			if ( limit == 0 || nSign < 0 ) return FALSE;
			if ( x1 == x2 ) i1 = (point_x-x2)*(point_x-x2);
			else{
				i2 = op_func_02((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)) * 100;
				if ( i2 == 0 ) return FALSE;
				i2 = (i1 * 100)/i2;
				i1 = i2*i2;
			}
			if ( i1 > lim2 ) return FALSE;
		}
		if ( nSign < 0 ){
			if ( x1 == x2 ) i1 = (point_x-x2)*(point_x-x2);
			else{
				i2 = op_func_02((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)) * 100;
				if ( i2 == 0 ) return FALSE;
				i2 = (i1 * 100)/i2;
				i1 = i2*i2;
			}
			if ( i1 < lim2 ) return FALSE;
		}
		x3 = x1; y3 = y1; x1 = x2; y1 = y2;
	}
	return TRUE;
}

BOOL get_polygon_points_sub(int *point_x,int *point_y,int pNum,POLYGON *polygon)
{
	int i,j,i1,i2,i3,j1,j2,k1,k2;
	int boundNum = 1;
	int endFlag = 0,repeatFlag;

	if(pNum < 3) return FALSE;

	i1 = 1000;
	for(i = 0;i < pNum;i++)
	{
		if(point_x[i] < i1) 
		{
			i1 = point_x[i];
			polygon->x[0] = i1;
			polygon->y[0] = point_y[i];
		}
		else if(point_x[i] == i1 && polygon->y[0] > point_y[i]) 
		{
			i1 = point_x[i];
			polygon->x[0] = i1;
			polygon->y[0] = point_y[i];
		}
	}
	while(1)
	{
		repeatFlag = 0;
		for(i = 0;i < pNum;i++)
		{
			if(point_x[i] == polygon->x[boundNum-1] && 
				point_y[i] == polygon->y[boundNum-1]) continue;
			if(boundNum > 1 && point_x[i] == polygon->x[boundNum-2] && 
				point_y[i] == polygon->y[boundNum-2]) continue;

			j1 = point_x[i];
			j2 = point_y[i];

			i1 = i2 = 0;
			for(j = 0;j < pNum;j++)
			{
				if(i == j) continue;
				if(point_x[j] == polygon->x[boundNum-1] && 
					point_y[j] == polygon->y[boundNum-1]) continue;
				i3 = (point_y[j] - polygon->y[boundNum-1])*(j1-polygon->x[boundNum-1]);
				i3 -= (point_x[j] - polygon->x[boundNum-1])*(j2-polygon->y[boundNum-1]);
				if(i3 < 0) i1++;
				if(i3 > 0) i2++;
				if(i3 == 0)
				{
					k1 = point_x[j] - j1;
					k2 = point_x[j] - polygon->x[boundNum-1];
					if(k1*k2 > 0 && abs(k1) < abs(k2))
					{
						i1++;	i2++;
					}
					k1 = point_y[j] - j2;
					k2 = point_y[j] - polygon->y[boundNum-1];
					if(k1*k2 > 0 && abs(k1) < abs(k2))
					{
						i1++;	i2++;
					}
					if(i1 == 0 || i2 == 0)
					{
						if(point_x[j] == polygon->x[0] && point_y[j] == polygon->y[0])
						{
							endFlag = 1;	break;
						}
					}
				}
				if(i1 > 0 && i2 > 0) break;
			}
			if(i1 > 0 && i2 > 0) continue;
			if(polygon->x[0] == j1 && polygon->y[0] == j2) endFlag = 1;
			if(endFlag == 1) break;
			polygon->x[boundNum] = j1;
			polygon->y[boundNum] = j2;
			boundNum++;
			repeatFlag = 1;
			break;
		}
		if(endFlag == 1 || repeatFlag == 0) break;
	}

	if(endFlag == 0) 
	{
		polygon->orderNum = 0;
		return FALSE;
	}
	polygon->orderNum = boundNum;
	if(boundNum < 3) return FALSE;
	return TRUE;
}

BOOL get_polygon_points(LPMPVECTEX pVect,POLYGON *polygon)
{
	int i,xx[100],yy[100],num;

	num = pVect->nNumber;
	for(i = 0;i < num;i++){
		xx[i] = pVect->item[i].x;
		yy[i] = pVect->item[i].y;
	}
	return(get_polygon_points_sub(xx,yy,num,polygon));
}

void get_smoothed_image(BYTE* img,int cxDIB,int cyDIB)
{
	int i, j, i1, nSum, idx = 0, idd = 0;
	int id = -cxDIB, id1, idd1, kk;
	int *pSum;
	BYTE *pTmp;

	pTmp = (BYTE*)malloc(sizeof(BYTE)*3*cxDIB);
	if ( pTmp == NULL ) return;
	pSum = (int*)calloc(cxDIB,sizeof(int));
	if ( pSum == NULL ){ free(pTmp); return; }
	
	for ( i=0,i1=-1; i<cyDIB+1; i++,idd+=cxDIB ){
		if ( i1 >= 2 ){ i1 = 0; id = 0; }
		else{ i1++; id += cxDIB; }
		if ( i >= 3 ){
			for ( j=0,id1=id; j<cxDIB; j++ ) pSum[j] -= pTmp[id1++];
		}
		if ( i < cyDIB ){
			idd1 = idd;
			for ( j=0,id1=id; j<cxDIB; j++ ){
				pTmp[id1] = img[idd1++]; pSum[j] += pTmp[id1++];
			}
		}
		if ( i < 1 ) continue;
		kk = idx; idx += cxDIB; nSum = 0;
		for ( j=0; j<cxDIB+1; j++ ){
			if ( j < cxDIB ) nSum += pSum[j];
			if ( j < 1 ) continue;
			if ( j >= 3 ) nSum -= pSum[j-3];
			img[kk++] = (BYTE)((nSum*114) >> 10);
		}
	}
	free(pSum); free(pTmp);
}

void get_smoothed_image2(BYTE* img,int cxDIB,int cyDIB,int nStep)
{
	int i, i1, j, nNum0 = 0, nNum, nSum, idx, idx1 = -1, idx2;
	int nWindow = 2*nStep+1, *pSum;
	BYTE *pTmp;

	pTmp = (BYTE*)malloc(sizeof(BYTE)*nWindow*cxDIB);
	if ( pTmp == NULL ) return;
	pSum = (int*)calloc(cxDIB,sizeof(int));
	if ( pSum == NULL ){ free(pTmp); return; }

	for ( i=0,i1=-1; i<cyDIB+nStep; i++ ){
		if ( i1 >= nWindow-1 ) i1 = 0;
		else i1++;
		if ( i >= nWindow ){
			for ( j=0; j<cxDIB; j++ ) pSum[j] -= pTmp[i1*cxDIB+j];
			nNum0--;
		}
		if ( i < cyDIB ){
			for ( j=0; j<cxDIB; j++ ){
				pTmp[i1*cxDIB+j] = img[i*cxDIB+j];
				pSum[j] += pTmp[i1*cxDIB+j];
			}
			nNum0++;
		}
		if ( i < nStep ) continue;
		idx = (i-nStep)*cxDIB;
		if ( idx1 >= nWindow-1 ) idx1 = 0;
		else idx1++;
		idx2 = idx1*cxDIB;
		nNum = nSum = 0;
		for ( j=0; j<cxDIB+nStep; j++ ){
			if ( j < cxDIB ){
				nNum += nNum0; nSum += pSum[j];
			}
			if ( j < nStep ) continue;
			if ( j >= nWindow ){
				nNum -= nNum0; nSum -= pSum[j-nWindow];
			}
			img[idx++] = (BYTE)(nSum/nNum);
		}
	}
	free(pSum); free(pTmp);
}

void get_sharpend_image(BYTE *Img,BYTE *RefImg,int cxDIB,int cyDIB,int nStep)
{
	int i, j, *pSum, nNum0 = 0, nNum, nSum;
	int nWindow = min(cyDIB, 2*nStep+1), idx;
	int nMean, nImgVal, nLower, nUpper, nDiff, nRefVal, nVal;
	BYTE **pTmp;

	pSum = (int*)calloc(cxDIB, sizeof(int));
	if ( pSum == NULL ) return;
	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if ( pTmp == NULL ){ free(pSum); return; }
	for ( i=0; i<nWindow; i++ ){
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if ( pTmp[i] == NULL ){
			free(pSum);
			for ( j=0; j<i; j++ ) free(pTmp[j]);
			free(pTmp); return;
		}
	}

	for(i=0; i<cyDIB+nStep; i++) {
		if(i >= nWindow) {
			for(j=0; j<cxDIB; j++) { 
				pSum[j] -= abs(pTmp[i-nWindow][j]-RefImg[(i-nWindow)*cxDIB+j]);
			}
			nNum0--;
		}
		if(i < cyDIB) {
			if(i >= nWindow) pTmp[i] = pTmp[i-nWindow];
			for(j=0; j<cxDIB; j++) {
				pTmp[i][j] = Img[i*cxDIB+j];
				pSum[j] += abs(pTmp[i][j]-RefImg[i*cxDIB+j]);
			}
			nNum0++;
		}
		if(i < nStep) continue;
		nSum = nNum = 0;
		for(j=0; j<cxDIB+nStep; j++) {
			if(j < cxDIB) {
				nNum += nNum0; nSum += pSum[j];
			}
			if(j < nStep) continue;
			if(j >= nWindow) {
				nNum -= nNum0; nSum -= pSum[j-nWindow];
			}
			nMean = nSum/nNum;
			idx = (i-nStep)*cxDIB + (j-nStep);
			nImgVal = Img[idx];
			nLower = ( nImgVal >= nMean ) ? nImgVal-nMean:0;
			nUpper = (nImgVal+nMean > 255) ? 255: nImgVal+nMean;
			nDiff = nUpper - nLower;
			nRefVal = RefImg[idx];

			if(nDiff == 0) nVal = nRefVal;
			else {
				if(nRefVal <= nLower) nVal = 0;
				else {
					if(nRefVal >= nUpper) nVal = 255;
					else nVal = ((nRefVal-nLower)*255)/nDiff;
				}
			}
 			Img[idx] = (BYTE)nVal;
		}
	}
	free(pSum);
	for ( i=0; i<nWindow; i++ ) free(pTmp[i]);
	free(pTmp);
}

void get_block_image(BYTE *Block1,BYTE *Block2,int nCol,int nRow,
				  BYTE *Img,int cxDIB,int cyDIB,int nMinGrayDiff)
{
	int i, j, m, n, mStart, nStart, idx;
	int nNum = 0, nNum100 = 0, nDiv100, nDiv;
	int t0, t1, t2, t3, var0, var1, var2, var3, nMax;
	int nSum[4], nDir[4], *pSum[4], *pTmp, sum1, sum2;
	int  nMinGrayDiffSum = 2*BLOCK_SIZE*BLOCK_SIZE*nMinGrayDiff;
	BYTE value;

	for ( i=0; i<4; i++ ){
		pSum[i] = (int*)calloc(4*nCol,sizeof(int));
		if ( pSum[i] == NULL ){
			for ( j=0; j<i; j++ ) free(pSum[i]);
		}
	}
	for ( i=0; i<nRow+1; i++ ){
		if ( i >= 3 ){
			for ( j=0; j<nCol; j++ ){
				pSum[3][4*j+0] -= pSum[0][4*j+0];
				pSum[3][4*j+1] -= pSum[0][4*j+1];
				pSum[3][4*j+2] -= pSum[0][4*j+2];
				pSum[3][4*j+3] -= pSum[0][4*j+3];
			}
			nNum--; nNum100 -= 100;
		}
		if ( i < nRow ){
			mStart = (i*BLOCK_SIZE == 0) ? 1 : i*BLOCK_SIZE;
			for ( j=0; j<nCol; j++ ){
				nSum[0] = nSum[1] = nSum[2] = nSum[3] = 0;
				nStart = (j*BLOCK_SIZE == 0) ? 1 : j*BLOCK_SIZE;
				for ( m=mStart; m<(i+1)*BLOCK_SIZE; m++ ){
					if ( m >= cyDIB-1 ) break;
					for ( n=nStart; n<(j+1)*BLOCK_SIZE; n++ ){
						if ( n >= cxDIB-1 ) break;
						idx = m*cxDIB+n; value = Img[idx];
						nSum[0] += (abs(value-Img[idx-1]) + abs(value-Img[idx+1]));
						nSum[1] += (abs(value-Img[idx+cxDIB+1]) + abs(value-Img[idx-cxDIB-1]));
						nSum[2] += (abs(value-Img[idx-cxDIB]) + abs(value-Img[idx+cxDIB]));
						nSum[3] += (abs(value-Img[idx+cxDIB-1]) + abs(value-Img[idx-cxDIB+1]));
					}
				}
				pSum[3][4*j+0] += nSum[0];
				pSum[3][4*j+1] += nSum[1];
				pSum[3][4*j+2] += nSum[2];
				pSum[3][4*j+3] += nSum[3];

				pSum[0][4*j+0] =  nSum[0];
				pSum[0][4*j+1] =  nSum[1];
				pSum[0][4*j+2] =  nSum[2];
				pSum[0][4*j+3] =  nSum[3];
			}
			nNum++; nNum100 += 100;
		}
		pTmp = pSum[0]; pSum[0] = pSum[1]; pSum[1] = pSum[2]; pSum[2] = pTmp;
		if ( i < 1 ) continue;
		nDiv100 = 0; nDiv = 0;
		nSum[0] = nSum[1] = nSum[2] = nSum[3] = 0;
		for ( j=0; j<nCol+1; j++ ){
			if ( j < nCol ){
				nSum[0] += pSum[3][4*j+0];
				nSum[1] += pSum[3][4*j+1];
				nSum[2] += pSum[3][4*j+2];
				nSum[3] += pSum[3][4*j+3];
				nDiv100 += nNum100; nDiv += nNum;
			}
			if ( j < 1 ) continue;
			if ( j >= 3 ){
				nSum[0] -= pSum[3][4*(j-3)+0];
				nSum[1] -= pSum[3][4*(j-3)+1];
				nSum[2] -= pSum[3][4*(j-3)+2];
				nSum[3] -= pSum[3][4*(j-3)+3];
				nDiv100 -= nNum100; nDiv -= nNum;
			}
			var0 = nDir[0] = nSum[0] / nDiv;
			var1 = nDir[1] = 71*nSum[1] / nDiv100;
			var2 = nDir[2] = nSum[2] / nDiv;
			var3 = nDir[3] = 71*nSum[3] / nDiv100;

			idx = (i-1)*nCol+j-1; value = Block2[idx] & 0x80;
			nMax = nDir[0];
			if ( nMax < nDir[1] ) nMax = nDir[1];
			if ( nMax < nDir[2] ) nMax = nDir[2];
			if ( nMax < nDir[3] ) nMax = nDir[3];

			Block2[idx] = 45;
			sum1 = nDir[2] + nDir[3];
			sum2 = nDir[2] + nDir[1];
			if ( sum1 < sum2 ) {
				sum2 = sum1;
				Block2[idx] = 75;
				var0 = nDir[1];	var1 = nDir[2];
				var2 = nDir[3]; var3 = nDir[0];
			}
			sum1 = nDir[0] + nDir[3];
			if ( sum1 < sum2 ){
				sum2 = sum1;
				Block2[idx] = 105;
				var0 = nDir[2];	var1 = nDir[3];
				var2 = nDir[0];	var3 = nDir[1];
			}
			sum1 = nDir[0] + nDir[1];
			if ( sum1 < sum2 ){
				Block2[idx] = 15;
				var0 = nDir[3]; var1 = nDir[0];
				var2 = nDir[1];	var3 = nDir[2];
			}
			t0 = (var0+var1+var2+var3) - 4*nMax;
			if ( t0 == 0 ){
				Block2[idx] = 0x7F;
				Block1[idx] = 255;
			}
			else{
				t1 = 15*(3*(var3-var0)-(var1-var2))/t0;
				Block2[idx] += t1;
				if ( Block2[idx] == 120 ) Block2[idx] = 0; 
				t0 = var2; t2 = var3;
				if ( t0 < var1 ) t2 = var0; else t0 = var1;
				t1 = abs(t1);
				t3 = ((t2-t0)*(15-t1))/30;
				if ( t3 <= t0 ) t0 -= t3; 
				else t0 = 0;
				t2 += t3;
				if ( t2 != 0 )
					Block1[idx] = (BYTE)(255*t0/t2);
				else 
					Block1[idx] = 255;
				if ( nMax <= nMinGrayDiffSum ) Block1[idx] = 255;
				Block2[idx] |= value;
			}
		}
	}
	for ( i=0; i<4; i++ ) free(pSum[i]);
}

void filter_block_image(BYTE *buffer1,BYTE *buffer2,int nCol,int nRow,BYTE value)
{
	int i, j, m, n, count;

	for ( i=0; i<nRow; i++ ){
		for ( j=0; j<nCol; j++ ){
			count = 0;
			for ( m=(i<1) ? 0:i-1; m<=i+1; m++ ){
				if ( m >= nRow ) break;
				for ( n=(j<1) ? 0:j-1; n<=j+1; n++ ){
					if ( n >= nCol ) break;
					if ( buffer1[m*nCol+n] >= value ) continue;
					count++;
				}
			}
			if ( count > 4 ) continue;
			buffer2[i*nCol+j] = 255;
		}
	}
}

int get_poincare_value(BYTE* pLoop,int nLoopSize)
{
	int i, diff1, diff2, buf;
	int sum=0, nPrevVal = pLoop[nLoopSize-1];

	for ( i=0; i<nLoopSize; i++ ){
		diff1 = abs(pLoop[i] - nPrevVal);
		if ( diff1 > 120 ) diff1 = 240 - diff1;
		diff2 = abs((pLoop[i]+120) - nPrevVal);
		if ( diff2 > 120 ) diff2 = 240 - diff2;
		if ( diff1 == diff2 ) return 0;
		buf = pLoop[i];
		if ( diff2 < diff1 ) buf += 120; 
		buf = buf - nPrevVal;
		if ( buf > 120 ) buf -= 240;
		else if ( buf < -120 ) buf += 240;
		sum += buf;
		nPrevVal = pLoop[i];
		if ( diff2 < diff1 ) nPrevVal += 120;
	}
	return (sum/120);
}

void get_sp_cand(SINGULAR *pSingular,BYTE *Block,int nCol,int nRow)
{
	int i, j, ret1, ret2, idx;
	BYTE loop[20];

	pSingular->nNumber = 0;
	for ( i=2; i<nRow-2; i++ ){
		for ( j=2; j<nCol-2; j++ ){
			idx = i*nCol + j;
			if ( Block[idx-nCol] == 0xFF ) continue;
			if ( Block[idx-nCol-1] == 0xFF ) continue;
			if ( Block[idx-nCol+1] == 0xFF ) continue;
			if ( Block[idx+1] == 0xFF ) continue;
			if ( Block[idx-1] == 0xFF ) continue;
			if ( Block[idx+nCol] == 0xFF ) continue;
			if ( Block[idx+nCol-1] == 0xFF ) continue;
			if ( Block[idx+nCol+1] == 0xFF ) continue;

			if ( Block[idx-2*nCol] == 0xFF ) continue;
			if ( Block[idx-2*nCol-1] == 0xFF ) continue;
			if ( Block[idx-2*nCol+1] == 0xFF ) continue;
			if ( Block[idx-nCol+2] == 0xFF ) continue;
			if ( Block[idx+2] == 0xFF ) continue;
			if ( Block[idx+nCol+2] == 0xFF ) continue;
			if ( Block[idx+2*nCol] == 0xFF ) continue;
			if ( Block[idx+2*nCol+1] == 0xFF ) continue;
			if ( Block[idx+2*nCol-1] == 0xFF ) continue;
			if ( Block[idx+nCol-2] == 0xFF ) continue;
			if ( Block[idx-2] == 0xFF ) continue;
			if ( Block[idx-nCol-2] == 0xFF ) continue;

			loop[0]  = Block[idx-nCol-1];
			loop[1]  = Block[idx-nCol];
			loop[2]  = Block[idx-nCol+1];
			loop[3]  = Block[idx+1];
			loop[4]  = Block[idx+nCol+1];
			loop[5]  = Block[idx+nCol];
			loop[6]  = Block[idx+nCol-1];
			loop[7]  = Block[idx-1];
			ret1 = get_poincare_value(&loop[0], 8);

			loop[8]  = Block[idx-2*nCol-1];
			loop[9]  = Block[idx-2*nCol];
			loop[10] = Block[idx-2*nCol+1];
			loop[11] = Block[idx-nCol+2];
			loop[12] = Block[idx+2];
			loop[13] = Block[idx+nCol+2];
			loop[14] = Block[idx+2*nCol+1];
			loop[15] = Block[idx+2*nCol];
			loop[16] = Block[idx+2*nCol-1];
			loop[17] = Block[idx+nCol-2];
			loop[18] = Block[idx-2];
			loop[19] = Block[idx-nCol-2];
			ret2 = get_poincare_value(&loop[8], 12);

			if ( ret1 != ret2 ){
				if ( ret1*ret2 <= 0 ) continue;
				ret1 = ret2;
			}
			if ( ret1 == 0 ) continue;
			pSingular->nX[pSingular->nNumber] = j;
			pSingular->nY[pSingular->nNumber] = i;
			pSingular->nOrient[pSingular->nNumber] = -1;
			pSingular->nType[pSingular->nNumber] = ret1;
			pSingular->nNumber++;
		}
	}
}

void get_orient_image(BYTE *Img,BYTE *OrntImg,int cxDIB,int cyDIB,
					int nStep,int nMinGrayDiff,BYTE* image_buffer4)
{
	int i, j, idx0, idx1, idx2, t1, t2, t3, nAvg, diff;
	int nSum[4], nDir[4], var0, var1, var2, var3, nMax, sum1, sum2;
	int *pSum, nWindow = 2*nStep + 4, idx;
	int nMinGrayDiffSum = 2*(2*nStep+1)*(2*nStep+1)*nMinGrayDiff;
	int nVar1,nVar2, nMean, nVar;
	BYTE **pTmp, tmp;

	pSum = (int*)calloc(4*cxDIB,sizeof(int));
	if ( pSum == NULL ) return;
	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if ( pTmp == NULL ){ free(pSum); return; }
	for ( i=0; i<nWindow; i++ ){
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if ( pTmp[i] == NULL ){
			free(pSum);
			for ( j=0; j<i; j++ ) free(pTmp[j]);
			free(pTmp); return;
		}
	}
	for ( i=0; i<cyDIB+nStep+1; i++ ){
		if ( i < cyDIB ){
			if ( i >= nWindow ) pTmp[i] = pTmp[i-nWindow];
			memcpy( pTmp[i],&Img[i*cxDIB],sizeof(BYTE)*cxDIB );
		}
		if ( i < 2 ) continue;
		if ( i < cyDIB ){
			idx0 = i - 2; idx1 = i - 1; idx2 = i;
			for ( j=1; j<cxDIB-1; j++ ){
				tmp = pTmp[idx1][j+0];
				pSum[4*j+0] += abs(tmp-pTmp[idx1][j-1]) + abs(tmp-pTmp[idx1][j+1]);
				pSum[4*j+1] += abs(tmp-pTmp[idx0][j-1])	+ abs(tmp-pTmp[idx2][j+1]);
				pSum[4*j+2] += abs(tmp-pTmp[idx0][j+0])	+ abs(tmp-pTmp[idx2][j+0]);
				pSum[4*j+3] += abs(tmp-pTmp[idx2][j-1])	+ abs(tmp-pTmp[idx0][j+1]);
			}
		}
		if ( i < nStep+1 ) continue;
		if ( i >= 2*nStep+3 ){
			idx0 = i - (2*nStep+3); idx1 = idx0 + 1; idx2 = idx0 + 2;
			for ( j=1; j<cxDIB-1; j++ ){
				tmp = pTmp[idx1][j+0];
				pSum[4*j+0] -= abs(tmp-pTmp[idx1][j-1]) + abs(tmp-pTmp[idx1][j+1]);
				pSum[4*j+1] -= abs(tmp-pTmp[idx0][j-1])	+ abs(tmp-pTmp[idx2][j+1]);
				pSum[4*j+2] -= abs(tmp-pTmp[idx0][j+0])	+ abs(tmp-pTmp[idx2][j+0]);
				pSum[4*j+3] -= abs(tmp-pTmp[idx2][j-1])	+ abs(tmp-pTmp[idx0][j+1]);
			}
		}
		nSum[0] = nSum[1] = nSum[2] = nSum[3] = 0;
		for ( j=0; j<cxDIB+nStep; j++ ){
			if ( j < cxDIB ){
				nSum[0] += pSum[4*j+0];
				nSum[1] += pSum[4*j+1];
				nSum[2] += pSum[4*j+2];
				nSum[3] += pSum[4*j+3];
			}
			if ( j < nStep ) continue;
			if ( j >= 2*nStep+1 ){
				nSum[0] -= pSum[4*(j-(2*nStep+1))+0];
				nSum[1] -= pSum[4*(j-(2*nStep+1))+1];
				nSum[2] -= pSum[4*(j-(2*nStep+1))+2];
				nSum[3] -= pSum[4*(j-(2*nStep+1))+3];
			}
			var0 = nDir[0] = nSum[0];
			var1 = nDir[1] = 71*nSum[1]/100;
			var2 = nDir[2] = nSum[2];
			var3 = nDir[3] = 71*nSum[3]/100;
			idx = (i-nStep-1)*cxDIB + j-nStep; 

			nMean = (var0+var1+var2+var3)/4;
			nVar = (abs(nMean-var0)+abs(nMean-var1)+abs(nMean-var2)+abs(nMean-var3))/4;
			nVar1 = nMean; nVar2 = nVar;
			nVar2 *= 160;
			if ( nVar1 < 1 ){
				if ( nVar2 < 1 ) nVar2 = 0;
				else nVar2 = 63;
			}
			else nVar2 = nVar2/nVar1;
			if ( nVar2 > 63 ) nVar2 = 63;
			image_buffer4[idx] = nVar2;

			tmp = OrntImg[idx] & 0x80;
			nMax = nDir[0];
			if ( nMax < nDir[1] ) nMax = nDir[1];
			if ( nMax < nDir[2] ) nMax = nDir[2];
			if ( nMax < nDir[3] ) nMax = nDir[3];

			OrntImg[idx] = 45;
			sum1 = var2 + var3;
			sum2 = var2 + var1;
			if ( sum1 < sum2 ){
				sum2 = sum1;
				OrntImg[idx] = 75;
				nDir[0] = var1;	nDir[1] = var2;
				nDir[2] = var3;	nDir[3] = var0;
			}
			sum1 = var0 + var3;
			if ( sum1 < sum2 ){
				sum2 = sum1;
				OrntImg[idx] = 105;
				nDir[0] = var2;	nDir[1] = var3;
				nDir[2] = var0;	nDir[3] = var1;
			}
			sum1 = var0 + var1;
			if ( sum1 < sum2 ){
				OrntImg[idx] = 15;
				nDir[0] = var3;	nDir[1] = var0;
				nDir[2] = var1;	nDir[3] = var2;
			}
			nAvg = (nDir[0]+nDir[1]+nDir[2]+nDir[3]) - nMax*4;
			if ( nAvg == 0 ){
				OrntImg[idx] = 0x7F;
				Img[idx] = 0xFF;
			}
			else{
				diff = 15*(3*(nDir[3]-nDir[0])-nDir[1]+nDir[2])/nAvg;
				OrntImg[idx] += diff;
				if ( OrntImg[idx] == 120 ) OrntImg[idx] = 0;
				if ( nDir[2] < nDir[1] ){
					t1 = nDir[2]; t2 = nDir[0];
				}
				else{
					t1 = nDir[1]; t2 = nDir[3];
				}
				t3 = ((t2-t1)*(15-abs(diff))) / 30;
				if(t3 <= t1) t1 -= t3; else t1 = 0;
				t2 += t3;
				if ( t2 != 0 ) Img[idx] = (BYTE)((255*t1)/t2);
				else Img[idx] = 0xFF;
			}
			if ( nMax <= nMinGrayDiffSum ) Img[idx] = 0xFF;
			OrntImg[idx] |= tmp;
		}
	}
	free(pSum);
	for ( i=0; i<nWindow; i++ ) free(pTmp[i]);
	free(pTmp);
}

BOOL check_in_singular(SINGULAR *pSingular,int x,int y,int nStep)
{
	int i, dx, dy, nTh = nStep*nStep;

	for ( i=0; i<pSingular->nNumber; i++ ){
		dx = pSingular->nX[i] - x; dy = pSingular->nY[i] - y;
		if ( dx*dx+dy*dy < nTh ) return TRUE;
	}
	return FALSE;
}

int get_difference(int x,int y,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, j, nSum = 0;
	BYTE data1, data2;

	if ( x < 2 || x >= cxDIB-2 || y < 2 || y >= cyDIB-2 ) return 0;

	data1 = OrntImg[y*cxDIB+x] & 0x7F;
	if ( data1 == 0x7F ) return 0;

	for ( i=y-2; i<=y+2; i++ ){
		for ( j=x-2; j<=x+2; j++ ){
			data2 = OrntImg[i*cxDIB+j] & 0x7F;
			if ( data2 == 0x7F ) continue;
			if ( data1 <= data2 ) data2 -= data1;
			else data2 = data1 - data2;
			if ( data2 > 60 ) data2 = 120 - data2;
			nSum += data2;
		}
	}
	return (nSum);
}

void get_singular_points(SINGULAR *pSingular,SINGULAR *pTmpSP,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	SINGULAR tmpSP;
	int j, k, m, nType, nMaxOriSum, nNumber = pTmpSP->nNumber;
	int nSPX, nSPY, sx, sy, nX, nY, nOriSum;
	BOOL flag;

	pSingular->nNumber = 0;
	if ( nNumber == 0 ) return;

	do{
		tmpSP.nNumber = 0;
		do {
			flag = FALSE;
			if ( pTmpSP->nNumber <= tmpSP.nNumber ) break;
			for ( j=0; j<pTmpSP->nNumber; j++ ){
				if ( pTmpSP->nType[j] == 0 ) continue;
				if ( tmpSP.nNumber != 0 ){
					if ( check_in_singular(&tmpSP,pTmpSP->nX[j],pTmpSP->nY[j],2) == FALSE ) continue;
				}
				flag = TRUE;
				tmpSP.nX[tmpSP.nNumber] = pTmpSP->nX[j];
				tmpSP.nY[tmpSP.nNumber] = pTmpSP->nY[j];
				tmpSP.nType[tmpSP.nNumber++] = pTmpSP->nType[j];
				pTmpSP->nType[j] = 0; nNumber--;
			}
		} while ( flag != FALSE );
		nType = tmpSP.nType[0];
		for ( j=0; j<tmpSP.nNumber-1; j++ ){
			if ( tmpSP.nType[j+1] == nType ) continue;
			nType = 0;
		}
		if ( nType == 0 ) continue;
		nMaxOriSum = -1;
		for ( j=0; j<tmpSP.nNumber; j++ ){
			sx = tmpSP.nX[j] * BLOCK_SIZE; sy = tmpSP.nY[j] * BLOCK_SIZE;
			for ( k=sy; k<sy+BLOCK_SIZE; k++ ){
				for(m=sx; m<sx+BLOCK_SIZE; m++) {
					nOriSum = get_difference(m,k,OrntImg,cxDIB,cyDIB);
					if ( nOriSum <= nMaxOriSum ) continue;
					nMaxOriSum = nOriSum; nSPX = m; nSPY = k;
				}
			}
		}
		nX = nSPX; nY = nSPY;
		for ( k=nSPY-BLOCK_SIZE; k<=nSPY+BLOCK_SIZE; k++ ){
			for ( m=nSPX-BLOCK_SIZE; m<=nSPX+BLOCK_SIZE; m++ ){
				nOriSum = get_difference(m,k,OrntImg,cxDIB,cyDIB);
				if ( nOriSum <= nMaxOriSum ) continue;
				nMaxOriSum = nOriSum; nX = m; nY = k;
			}
		}
		pSingular->nX[pSingular->nNumber] = nX;
		pSingular->nY[pSingular->nNumber] = nY;
		pSingular->nNumber++;
	} while(nNumber != 0);
}

void get_bkgrnd(BYTE *Img,BYTE *OrntImg,int cxDIB,int cyDIB,
					  SINGULAR *pSingular,int nStep,int nThValue)
{
	int i, j, nWindow = 2*nStep+1, nSum;
	int nMeanSquare = nWindow*nWindow/2, *pSum;

	pSum = (int*)calloc(cxDIB, sizeof(int));
	if ( pSum == NULL ) return;

	for ( i=0; i<cyDIB+nStep; i++ ){
		if ( i < cyDIB ){
			for ( j=0; j<cxDIB; j++ ){
				if ( Img[i*cxDIB+j] >= nThValue ) continue;
				pSum[j]++;
			}
		}
		if ( i < nStep ) continue;
		if ( i >= nWindow ){
			for ( j=0; j<cxDIB; j++ ){
				if ( Img[(i-nWindow)*cxDIB+j] >= nThValue ) continue;
				pSum[j]--;
			}
		}
		nSum = 0;
		for ( j=0; j<cxDIB+nStep; j++ ){
			if ( j < cxDIB ) nSum += pSum[j];
			if ( j < nStep ) continue;
			if ( j >= nWindow ) nSum -= pSum[j-nWindow];
			if ( nSum <= nMeanSquare ){
				if ( check_in_singular(pSingular,j-nStep,i-nStep,16) == FALSE ){
					OrntImg[(i-nStep)*cxDIB+j-nStep] |= 0x80; continue;
				}
			}
		}
	}
	free(pSum);
}

void image_proc_01(BYTE *OrntImg,BYTE *Img,int cxDIB,int cyDIB)
{
	int i, j, nWeight, nSum, nNum, nIdx, nStep = 6;
	int k, n, tx, ty, nWindow = 2*nStep+1, idx;
	BYTE **pTmp, nOrnt;
	BOOL bBound;

	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if ( pTmp == NULL ) return;
	for ( i=0; i<nWindow; i++ ){
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if ( pTmp[i] == NULL ){
			for ( j=0; j<i; j++ ) free(pTmp[j]);
			free(pTmp); return;
		}
	}
	for ( i=0; i<cyDIB+nStep; i++ ){
		if ( i < cyDIB ){
			if ( i >= nWindow ) pTmp[i] = pTmp[i-nWindow];
			memcpy( pTmp[i],&Img[i*cxDIB],sizeof(BYTE)*cxDIB );
		}
		if ( i < nStep ) continue;
		if ( i-nStep>=nStep && i<cyDIB ) bBound = FALSE;
		else bBound = TRUE;
		for ( j=nStep; j<cxDIB-nStep; j++ ){
			if ( bBound ){
				idx = (i-nStep)*cxDIB + j;
				nOrnt = OrntImg[idx] & 0x7F;
				if ( nOrnt == 0x7F ){
					nSum = pTmp[i-nStep][j]; nNum = 1;
					if ( j-1 >= 0 && j-1 < cxDIB ){
						nSum += pTmp[i-nStep][j-1]; nNum++;
					}
					if ( j+1 >= 0 && j+1 < cxDIB ){
						nSum += pTmp[i-nStep][j+1]; nNum++;
					}
					if ( i-nStep-1 >= 0 && i-nStep-1 < cyDIB ){
						if ( j-1 >= 0 && j-1 < cxDIB ){
							nSum += pTmp[i-nStep-1][j-1]; nNum++;
						}
						if ( j >= 0 && j < cxDIB ){
							nSum += pTmp[i-nStep-1][j]; nNum++;
						}
						if ( j+1 >= 0 && j+1 < cxDIB ){
							nSum += pTmp[i-nStep-1][j+1]; nNum++;
						}
					}
					if ( i-nStep+1 >= 0 && i-nStep+1 < cyDIB ){
						if ( j-1 >= 0 && j-1 < cxDIB ){
							nSum += pTmp[i-nStep+1][j-1]; nNum++;
						}
						if ( j >= 0 && j < cxDIB ){
							nSum += pTmp[i-nStep+1][j]; nNum++;
						}
						if ( j+1 >= 0 && j+1 < cxDIB ){
							nSum += pTmp[i-nStep+1][j+1]; nNum++;
						}
					}
					Img[idx] = (BYTE)(nSum/nNum);
				}
				else{
					nSum = pTmp[i-nStep][j] * _table1[nOrnt];
					nWeight = _table1[nOrnt];
					nIdx = nOrnt*18;
					for ( k=0; k<_table2[nOrnt]; k++ ){
						ty = i - nStep + _table3[nIdx+k];
						tx = j + _table4[nIdx+k];
						if ( ty>=0 && ty<cyDIB && tx>=0 && tx<cxDIB ){
							nSum += pTmp[ty][tx] * _table5[nIdx+k];
							nWeight += _table5[nIdx+k];
						}
						ty = i - nStep - _table3[nIdx+k];
						tx = j - _table4[nIdx+k];
						if ( ty>=0 && ty<cyDIB && tx>=0 && tx<cxDIB ){
							nSum += pTmp[ty][tx] * _table5[nIdx+k];
							nWeight += _table5[nIdx+k];
						}
					}
					Img[idx] = (BYTE)(nSum/nWeight);
				}
			}
			else{
				idx = (i-nStep)*cxDIB + j;
				nOrnt = OrntImg[idx] & 0x7F;
				if ( nOrnt == 0x7F ){
					nSum  = pTmp[i-nStep-1][j-1];
					nSum += pTmp[i-nStep-1][j  ];
					nSum += pTmp[i-nStep-1][j+1];
					nSum += pTmp[i-nStep  ][j-1];
					nSum += pTmp[i-nStep  ][j  ];
					nSum += pTmp[i-nStep  ][j+1];
					nSum += pTmp[i-nStep+1][j-1];
					nSum += pTmp[i-nStep+1][j  ];
					nSum += pTmp[i-nStep+1][j+1];
					Img[idx] = (BYTE)(nSum/9);
				}
				else{
					nSum = pTmp[i-nStep][j] * _table1[nOrnt];
					nIdx = nOrnt*18;
					for(k=0; k<_table2[nOrnt]; k++) {
						ty = i - nStep + _table3[nIdx+k];
						tx = j + _table4[nIdx+k];
						nSum += pTmp[ty][tx] * _table5[nIdx+k];
						ty = i - nStep - _table3[nIdx+k];
						tx = j - _table4[nIdx+k];
						nSum += pTmp[ty][tx] * _table5[nIdx+k];
					}
					Img[idx] = (BYTE)(nSum/112211);
				}
			}
		}
		for ( j=0; j<nStep; j++ ){
			idx = (i-nStep)*cxDIB + j;
			nOrnt = OrntImg[idx] & 0x7F;
			if ( nOrnt == 0x7F ){
				nSum = pTmp[i-nStep][j]; nNum = 1;
				if ( i-nStep >= 0 && i-nStep < cyDIB ){
					if ( j-1 >= 0 && j-1 < cxDIB ){
						nSum += pTmp[i-nStep][j-1]; nNum++;
					}
					if ( j+1 >= 0 && j+1 < cxDIB ){
						nSum += pTmp[i-nStep][j+1]; nNum++;
					}
				}
				if ( i-nStep-1 >= 0 && i-nStep-1 < cyDIB ){
					if ( j-1 >= 0 && j-1 < cxDIB ){
						nSum += pTmp[i-nStep-1][j-1]; nNum++;
					}
					if ( j >= 0 && j < cxDIB ){
						nSum += pTmp[i-nStep-1][j]; nNum++;
					}
					if ( j+1 >= 0 && j+1 < cxDIB ){
						nSum += pTmp[i-nStep-1][j+1]; nNum++;
					}
				}
				if ( i-nStep+1 >= 0 && i-nStep+1 < cyDIB ){
					if ( j-1 >= 0 && j-1 < cxDIB ){
						nSum += pTmp[i-nStep+1][j-1]; nNum++;
					}
					if ( j >= 0 && j < cxDIB ){
						nSum += pTmp[i-nStep+1][j]; nNum++;
					}
					if ( j+1 >= 0 && j+1 < cxDIB ){
						nSum += pTmp[i-nStep+1][j+1]; nNum++;
					}
				}
				Img[idx] = (BYTE)(nSum/nNum);
			}
			else{
				nSum = pTmp[i-nStep][j] * _table1[nOrnt];
				nWeight = _table1[nOrnt];
				nIdx = nOrnt*18;
				for ( k=0; k<_table2[nOrnt]; k++ ){
					ty = i - nStep + _table3[nIdx+k];
					tx = j + _table4[nIdx+k];
					if ( ty>=0 && ty<cyDIB && tx>=0 && tx<cxDIB ){
						nSum += pTmp[ty][tx] * _table5[nIdx+k];
						nWeight += _table5[nIdx+k];
					}
					ty = i - nStep - _table3[nIdx+k];
					tx = j - _table4[nIdx+k];
					if ( ty>=0 && ty<cyDIB && tx>=0 && tx<cxDIB ){
						nSum += pTmp[ty][tx] * _table5[nIdx+k];
						nWeight += _table5[nIdx+k];
					}
				}
				Img[idx] = (BYTE)(nSum/nWeight);
			}
			n = cxDIB - 1 - j; idx = (i-nStep)*cxDIB + n;
			nOrnt = OrntImg[idx] & 0x7F;
			if ( nOrnt == 0x7F ){
				nSum = pTmp[i-nStep][n]; nNum = 1;
				if ( i-nStep >= 0 && i-nStep < cyDIB ){
					if ( n-1 >= 0 && n-1 < cxDIB ){
						nSum += pTmp[i-nStep][n-1]; nNum++;
					}
					if ( n+1 >= 0 && n+1 < cxDIB ){
						nSum += pTmp[i-nStep][n+1]; nNum++;
					}
				}
				if ( i-nStep-1 >= 0 && i-nStep-1 < cyDIB ){
					if ( n-1 >= 0 && n-1 < cxDIB ){
						nSum += pTmp[i-nStep-1][n-1]; nNum++;
					}
					if ( n >= 0 && n < cxDIB ){
						nSum += pTmp[i-nStep-1][n]; nNum++;
					}
					if ( n+1 >= 0 && n+1 < cxDIB ){
						nSum += pTmp[i-nStep-1][n+1]; nNum++;
					}
				}
				if ( i-nStep+1 >= 0 && i-nStep+1 < cyDIB ){
					if ( n-1 >= 0 && n-1 < cxDIB ){
						nSum += pTmp[i-nStep+1][n-1]; nNum++;
					}
					if ( n >= 0 && n < cxDIB ){
						nSum += pTmp[i-nStep+1][n]; nNum++;
					}
					if ( n+1 >= 0 && n+1 < cxDIB ){
						nSum += pTmp[i-nStep+1][n+1]; nNum++;
					}
				}
				Img[idx] = (BYTE)(nSum/nNum);
			}
			else{
				nSum = pTmp[i-nStep][n] * _table1[nOrnt];
				nWeight = _table1[nOrnt];
				nIdx = nOrnt*18;
				for(k=0; k<_table2[nOrnt]; k++) {
					ty = i - nStep + _table3[nIdx+k];
					tx = n + _table4[nIdx+k];
					if ( ty>=0 && ty<cyDIB && tx>=0 && tx<cxDIB ){
						nSum += pTmp[ty][tx] * _table5[nIdx+k];
						nWeight += _table5[nIdx+k];
					}
					ty = i - nStep - _table3[nIdx+k];
					tx = n - _table4[nIdx+k];
					if ( ty>=0 && ty<cyDIB && tx>=0 && tx<cxDIB ){
						nSum += pTmp[ty][tx] * _table5[nIdx+k];
						nWeight += _table5[nIdx+k];
					}
				}
				Img[idx] = (BYTE)(nSum/nWeight);
			}
		}
	}
	for ( i=0; i<nWindow; i++ ) free(pTmp[i]);
	free(pTmp);
}

void get_binary_image2(BYTE *OrntImg,BYTE *BinImg,BYTE *RefImg,int cxDIB,int cyDIB,int nStep1,int nStep2)
{
	int i, j, nTh, nWindow1 = 2*nStep1+1, nWindow2 = 2*nStep2+1;
	int *pSum1, *pSum2, nSum, nNum, nNum1 = 0, nNum2 = 0;

	pSum1 = (int*)calloc(cxDIB, sizeof(int));
	if(pSum1 == NULL) return;
	pSum2 = (int*)calloc(cxDIB, sizeof(int));
	if(pSum2 == NULL) { free(pSum1); return; }

	for(i=0; i<cyDIB+nStep2; i++) {
		if(i < cyDIB) {
			for(j=0; j<cxDIB; j++) {
				pSum1[j] += RefImg[i*cxDIB+j];
				pSum2[j] += RefImg[i*cxDIB+j];
			}
			nNum1++; nNum2++;
		}
		if(i-nStep1>=0 && i-nStep1<cyDIB) {
			if(i >= nWindow1) {
				for(j=0; j<cxDIB; j++) {
					pSum1[j] -= RefImg[(i-nWindow1)*cxDIB+j];
				}
				nNum1--;
			}
			nSum =0; nNum = 0;
			for(j=0; j<cxDIB+nStep1; j++) {
				if(j < cxDIB) {
					nSum += pSum1[j]; nNum += nNum1;
				}
				if(j < nStep1) continue;
				if(j >= nWindow1) {
					nSum -= pSum1[j-nWindow1]; nNum -= nNum1;
				}
				if ( (OrntImg[(i-nStep1)*cxDIB+(j-nStep1)]&0x80) == 0 ){
					BinImg[(i-nStep1)*cxDIB+(j-nStep1)] = (BYTE)(nSum/nNum);
				}
			}
		}
		if(i < nStep2) continue;
		if(i >= nWindow2) {
			for(j=0; j<cxDIB; j++) {
				pSum2[j] -= RefImg[(i-nWindow2)*cxDIB+j];
			}
			nNum2--;
		}
		nSum =0; nNum = 0;
		for(j=0; j<cxDIB+nStep2; j++) {
			if(j < cxDIB) {
				nSum += pSum2[j]; nNum += nNum2;
			}
			if(j < nStep2) continue;
			if(j >= nWindow2) {
				nSum -= pSum2[j-nWindow2]; nNum -= nNum2;
			}
			if ( (OrntImg[(i-nStep2)*cxDIB+(j-nStep2)]&0x80) == 0 ){
				nTh = (nSum/nNum + BinImg[(i-nStep2)*cxDIB+(j-nStep2)]) / 2;
				if(RefImg[(i-nStep2)*cxDIB+(j-nStep2)] >= nTh )
					BinImg[(i-nStep2)*cxDIB+(j-nStep2)] = 0xFF;
				else 
					BinImg[(i-nStep2)*cxDIB+(j-nStep2)] = 0x00;
			}
		}
	}
	free(pSum1); free(pSum2);
}

void image_proc_05(BYTE *Img,BYTE *OrntImg,int cxDIB,int cyDIB,int nStep)
{
	int i, j, nTotalSum, nSum, nWindow = min(cyDIB, 2*nStep+4);
	int *pTotalSum, *pSum, idx0, idx1, idx2;
	BYTE **pTmp, nOrnt, tmp;

	pTotalSum = (int*)calloc(cxDIB, sizeof(int));
	if(pTotalSum == NULL) return;
	pSum = (int*)calloc(cxDIB, sizeof(int));
	if(pSum == NULL) { free(pTotalSum); return; }
	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if(pTmp == NULL) {
		free(pTotalSum); free(pSum); return;
	}
	for(i=0; i<nWindow; i++) {
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if(pTmp[i] == NULL) {
			free(pTotalSum); free(pSum);
			for(j=0; j<i; j++) free(pTmp[j]);
			free(pTmp); return;
		}
	}
	for(i=0; i<cyDIB+nStep+1; i++) {
		if(i < cyDIB) {
			if(i >= nWindow) pTmp[i] = pTmp[i-nWindow];
			memcpy( pTmp[i],&Img[i*cxDIB],sizeof(BYTE)*cxDIB );
		}
		if(i>=2 && i<cyDIB) {
			for(j=1; j<cxDIB-1; j++) {
				if((OrntImg[(i-1)*cxDIB+j] & 0x80) != 0 ) continue;
				pTotalSum[j]++;
				nSum  = pTmp[i-1][j-1];
				nSum += pTmp[i-1][j+0];
				nSum += pTmp[i-1][j+1];
				nSum += pTmp[i-2][j-1];
				nSum += pTmp[i-2][j+0];
				nSum += pTmp[i-2][j+1];
				nSum += pTmp[i+0][j-1];
				nSum += pTmp[i+0][j+0];
				nSum += pTmp[i+0][j+1];
				if(nSum<9 || nSum>2286) continue;
				pSum[j]++;
			}
		}
		if(i < nStep+1 ) continue;
		if(i >= 2*nStep+3) {
			idx0 = i - (2*nStep+3);
			idx1 = i - (2*nStep+2);
			idx2 = i - (2*nStep+1);
			for(j=1; j<cxDIB-1; j++) {
				if((OrntImg[idx1*cxDIB+j] & 0x80) != 0) continue;
				pTotalSum[j]--;
				nSum  = pTmp[idx1][j-1];
				nSum += pTmp[idx1][j+0];
				nSum += pTmp[idx1][j+1];
				nSum += pTmp[idx0][j-1];
				nSum += pTmp[idx0][j+0];
				nSum += pTmp[idx0][j+1];
				nSum += pTmp[idx2][j-1];
				nSum += pTmp[idx2][j+0];
				nSum += pTmp[idx2][j+1];
				if(nSum<9 || nSum>2286) continue;
				pSum[j]--;
			}
		}
		nTotalSum = 0; nSum = 0;
		for(j=0; j<cxDIB+nStep; j++) {
			if(j < cxDIB) {
				nTotalSum += pTotalSum[j];
				nSum += pSum[j];
			}
			if(j < nStep) continue;
			if(j >= 2*nStep+1) {
				nTotalSum -= pTotalSum[j-(2*nStep+1)];
				nSum -= pSum[j-(2*nStep+1)];
			}
			if ( (OrntImg[(i-nStep-1)*cxDIB+(j-nStep)]&0x80) == 0 ){
				if(nTotalSum == 0) tmp = 0;
				else {
					tmp = (255*nSum)/nTotalSum;
					nOrnt = OrntImg[(i-nStep-1)*cxDIB+(j-nStep)] & 0x7F;
					if(nOrnt < 120)
						tmp = (tmp * _table_05[nOrnt]) / 100;
				}
				Img[(i-nStep-1)*cxDIB+(j-nStep)] = (BYTE)tmp;
			}
		}
	}
	free(pTotalSum); free(pSum);
	for(i=0; i<nWindow; i++) free(pTmp[i]);
	free(pTmp);
}

void get_binary_image1(
							BYTE *OrntImg,BYTE *BinImg,BYTE *RefImg,
							int cxDIB, int cyDIB,
							int nMaxValue, int nMinValue,
							int nMinStep, int nMidStep, int nMaxStep
								  )
{
	int i, j, nMinNum = 0, nMidNum = 0, nMaxNum = 0, nSum, nNumber;
	int nMinWindow = 2*nMinStep+1, nMidWindow = 2*nMidStep+1, nMaxWindow = 2*nMaxStep+1;
	int *pMinSum, *pMidSum, *pMaxSum, nWindow = min(cyDIB, nMaxWindow+1);
	BYTE **pTmp, nTh;

	pMinSum = (int*)calloc(cxDIB, sizeof(int));
	if(pMinSum == NULL) return;
	pMidSum = (int*)calloc(cxDIB, sizeof(int));
	if(pMidSum == NULL) { free(pMinSum); return; }
	pMaxSum = (int*)calloc(cxDIB, sizeof(int));
	if(pMaxSum == NULL) { free(pMinSum); free(pMidSum); return; }
	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if(pTmp == NULL) { free(pMinSum); free(pMidSum); free(pMaxSum); return; }
	for(i=0; i<nWindow; i++) {
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if(pTmp[i] == NULL) {
			free(pMinSum); free(pMidSum); free(pMaxSum);
			for(j=0; j<i; j++) free(pTmp[j]);
			free(pTmp); return;
		}
	}

	for(i=0; i<cyDIB+nMaxStep; i++ ){
		if(i<cyDIB){
			if(i >= nWindow) pTmp[i] = pTmp[i-nWindow];
			for(j=0; j<cxDIB; j++) {
				pTmp[i][j]  = BinImg[i*cxDIB+j];
				pMinSum[j]  += BinImg[i*cxDIB+j];
				pMidSum[j]  += BinImg[i*cxDIB+j];
				pMaxSum[j]  += BinImg[i*cxDIB+j];
			}
			nMinNum++; nMidNum++; nMaxNum++;
		}
		if(i-nMinStep>=0 && i-nMinStep<cyDIB) {
			if(i >= nMinWindow) {
				for(j=0; j<cxDIB; j++) 
					pMinSum[j] -= pTmp[i-nMinWindow][j];
				nMinNum--;
			}
			nNumber = nSum = 0;
			for(j=0; j<cxDIB+nMinStep; j++) {
				if(j < cxDIB) {
					nSum += pMinSum[j]; nNumber += nMinNum;
				}
				if(j < nMinStep) continue;
				if(j >= nMinWindow) {
					nSum -= pMinSum[j-nMinWindow]; nNumber -= nMinNum;
				}
				if ( (OrntImg[(i-nMinStep)*cxDIB+(j-nMinStep)]&0x80) == 0 ){
					if(RefImg[(i-nMinStep)*cxDIB+(j-nMinStep)] < nMaxValue) continue;
					nTh = (BYTE)(nSum/nNumber);
					if(pTmp[i-nMinStep][j-nMinStep] >= nTh) 
						BinImg[(i-nMinStep)*cxDIB+(j-nMinStep)] = 0xFF;
					else 
						BinImg[(i-nMinStep)*cxDIB+(j-nMinStep)] = 0x00;
				}
			}
		}
		if(i-nMidStep>=0 && i-nMidStep<cyDIB) {
			if(i >= nMidWindow) {
				for(j=0; j<cxDIB; j++) 
					pMidSum[j] -= pTmp[i-nMidWindow][j];
				nMidNum--;
			}
			nNumber = nSum = 0;
			for(j=0; j<cxDIB+nMidStep; j++) {
				if(j < cxDIB) {
					nSum += pMidSum[j]; nNumber += nMidNum;
				}
				if(j < nMidStep) continue;
				if(j >= nMidWindow) {
					nSum -= pMidSum[j-nMidWindow]; nNumber -= nMidNum;
				}
				if ( (OrntImg[(i-nMidStep)*cxDIB+(j-nMidStep)]&0x80) == 0 ){
					if(RefImg[(i-nMidStep)*cxDIB+(j-nMidStep)] < nMinValue) continue;
					if(RefImg[(i-nMidStep)*cxDIB+(j-nMidStep)] >= nMaxValue) continue;
					nTh = (BYTE)( nSum/nNumber);
					if(pTmp[i-nMidStep][j-nMidStep] >= nTh)
						BinImg[(i-nMidStep)*cxDIB+(j-nMidStep)] = 0xFF;
					else
						BinImg[(i-nMidStep)*cxDIB+(j-nMidStep)] = 0x00;
				}
			}
		}
		if(i-nMaxStep>=0 && i-nMaxStep<cyDIB) {
			if(i >= nMaxWindow) {
				for(j=0; j<cxDIB; j++) 
					pMaxSum[j] -= pTmp[i-nMaxWindow][j];
				nMaxNum--;
			}
			nNumber = nSum = 0;
			for(j=0; j<cxDIB+nMaxStep; j++) {
				if(j < cxDIB) {
					nSum += pMaxSum[j]; nNumber += nMaxNum;
				}
				if(j < nMaxStep) continue;
				if(j >= nMaxWindow) {
					nSum -= pMaxSum[j-nMaxWindow]; nNumber -= nMaxNum;
				}
				if ( (OrntImg[(i-nMaxStep)*cxDIB+(j-nMaxStep)]&0x80) == 0 ){
					if(RefImg[(i-nMaxStep)*cxDIB+(j-nMaxStep)] >= nMinValue) continue;
					nTh = (BYTE)( nSum/nNumber);
					if(pTmp[i-nMaxStep][j-nMaxStep] >= nTh )
						BinImg[(i-nMaxStep)*cxDIB+(j-nMaxStep)] = 0xFF;
					else
						BinImg[(i-nMaxStep)*cxDIB+(j-nMaxStep)] = 0x00;
				}
			}
		}
	}
	free(pMinSum); free(pMidSum); free(pMaxSum);
	for(i=0; i<nWindow; i++) free(pTmp[i]);
	free(pTmp);
}

void image_proc_02(BYTE *Img,BYTE* OrntImg,int cxDIB,int cyDIB,int nStep,int nCount)
{
	int i, j, k, m, x, y, idx0, idx1, idx2, sum1, sum2, num;
	int nSum[4], nDir[4], nMax, var0, var1, var2, var3, *pSum;
	int nOrnt, nAvg, diff, nWindow = 2*nStep+4, idx;
	BYTE tmp, **pTmp;

	pSum = (int*)malloc(4*cxDIB*sizeof(int));
	if ( pSum == NULL ) return;
	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if ( pTmp == NULL ){ free(pSum); return; }
	for ( i=0; i<nWindow; i++ ){
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if ( pTmp[i] == NULL ){
			free(pSum);
			for ( j=0; j<i; j++ ) free(pTmp[j]);
			free(pTmp); return;
		}
	}
	for ( k=0; k<nCount; k++ ){
		memset( pSum, 0, sizeof(int)*4*cxDIB );
		for ( i=0; i<cyDIB+nStep+1; i++ ){
			if(i < cyDIB) {
				if(i >= nWindow) pTmp[i] = pTmp[i-nWindow];
				memcpy( pTmp[i], &Img[i*cxDIB], sizeof(BYTE)*cxDIB );
			}
			if(i>=2 && i<cyDIB) {
				for(j=1; j<cxDIB-1; j++) {
					tmp = pTmp[i-1][j+0];
					pSum[4*j+0] += (abs(tmp-pTmp[i-1][j-1]) + abs(tmp-pTmp[i-1][j+1]));
					pSum[4*j+1] += (abs(tmp-pTmp[i-2][j-1]) + abs(tmp-pTmp[i+0][j+1]));
					pSum[4*j+2] += (abs(tmp-pTmp[i-2][j+0]) + abs(tmp-pTmp[i+0][j+0]));
					pSum[4*j+3] += (abs(tmp-pTmp[i-2][j+1]) + abs(tmp-pTmp[i+0][j-1]));
				}
			}
			if(i < nStep+1) continue;
			if(i >= 2*nStep+3) {
				idx0 = i - (2*nStep+3);
				idx1 = i - (2*nStep+2);
				idx2 = i - (2*nStep+1);
				for(j=1; j<cxDIB-1; j++) {
					tmp = pTmp[idx1][j+0];
					pSum[4*j+0] -= (abs(tmp-pTmp[idx1][j-1]) + abs(tmp-pTmp[idx1][j+1]));
					pSum[4*j+1] -= (abs(tmp-pTmp[idx0][j-1]) + abs(tmp-pTmp[idx2][j+1]));
					pSum[4*j+2] -= (abs(tmp-pTmp[idx0][j+0]) + abs(tmp-pTmp[idx2][j+0]));
					pSum[4*j+3] -= (abs(tmp-pTmp[idx0][j+1]) + abs(tmp-pTmp[idx2][j-1]));
				}
			}
			nSum[0] = nSum[1] = nSum[2] = nSum[3] = 0;
			for(j=0; j<cxDIB+nStep; j++) {
				if(j < cxDIB) {
					nSum[0] += pSum[4*j+0]; nSum[1] += pSum[4*j+1];
					nSum[2] += pSum[4*j+2];	nSum[3] += pSum[4*j+3];
				}
				if(j < nStep) continue;
				idx0 = j - (2*nStep+1);
				if(idx0 >= 0) {
					nSum[0] -= pSum[4*idx0+0]; nSum[1] -= pSum[4*idx0+1];
					nSum[2] -= pSum[4*idx0+2]; nSum[3] -= pSum[4*idx0+3];
				}

				idx = (i-nStep-1)*cxDIB + j-nStep;
				if ( (OrntImg[idx]&0x80) != 0 ) continue;					

				var0 = nDir[0] = nSum[0];
				var1 = nDir[1] = 71*nSum[1]/100;
				var2 = nDir[2] = nSum[2];
				var3 = nDir[3] = 71*nSum[3]/100;

				nMax = nDir[0];
				if ( nMax < nDir[1] ) nMax = nDir[1];
				if ( nMax < nDir[2] ) nMax = nDir[2];
				if ( nMax < nDir[3] ) nMax = nDir[3];
				
				nOrnt = 45;
				sum1 = var2 + var3;
				sum2 = var2 + var1;
				if(sum1 < sum2) {
					sum2 = sum1;
					nOrnt = 75;
					nDir[0] = var1; nDir[1] = var2;
					nDir[2] = var3;	nDir[3] = var0;
				}
				sum1 = var0 + var3;
				if(sum1 < sum2) {
					sum2 = sum1;
					nOrnt = 105;
					nDir[0] = var2;	nDir[1] = var3;
					nDir[2] = var0;	nDir[3] = var1;
				}
				sum1 = var0 + var1;
				if(sum1 < sum2) {
					nOrnt = 15;
					nDir[0] = var3;	nDir[1] = var0;
					nDir[2] = var1;	nDir[3] = var2;
				}
				nAvg = (nDir[0]+nDir[1]+nDir[2]+nDir[3]) - nMax*4;
				if(nAvg == 0) nOrnt = 0x7F;
				else {
					diff = 15*(3*(nDir[3]-nDir[0])-nDir[1]+nDir[2])/nAvg;
					nOrnt += diff;
					if(nOrnt == 120) nOrnt = 0;
				}
				nOrnt &= 0x7F;
				if(i-nStep-2>=6 && j-nStep-1>=5 && i-nStep<cyDIB-6 && j-nStep+1<cxDIB-5) {
					if(nOrnt == 0x7F) {
						idx2  = j - nStep;
						idx1  = i - nStep - 1;
						sum1  = (pTmp[idx1][idx2-1] + pTmp[idx1][idx2+0] + pTmp[idx1][idx2+1]);
						idx1  = i - nStep - 2;
						sum1 += (pTmp[idx1][idx2-1] + pTmp[idx1][idx2+0] + pTmp[idx1][idx2+1]);
						idx1  = i - nStep;
						sum1 += (pTmp[idx1][idx2-1] + pTmp[idx1][idx2+0] + pTmp[idx1][idx2+1]);
						
						tmp = (BYTE)(sum1/9);
					}
					else{
						idx1 = i- nStep - 1; idx2 = j - nStep;
						sum1 = pTmp[idx1][idx2]*_table1[nOrnt];
						idx0 = 18*nOrnt;
						for(m=0; m<_table2[nOrnt]; m++) {
							y = idx1 + _table3[idx0+m];
							x = idx2 + _table4[idx0+m];
							sum1 += pTmp[y][x]*_table5[idx0+m];
						
							y = idx1 - _table3[idx0+m];
							x = idx2 - _table4[idx0+m];
							sum1 += pTmp[y][x]*_table5[idx0+m];
						}
						tmp = (BYTE)(sum1/112211);
					}
					Img[idx] = tmp;
				}
				else{
					if ( nOrnt == 0x7F ){
						idx1 = i- nStep - 1; idx2 = j - nStep;
						sum1 = pTmp[idx1][idx2]; num = 1;
						if(idx2-1>=0 && idx2-1<cxDIB) {
							sum1 += pTmp[idx1][idx2-1]; num++;
						}
						if(idx2+1>=0 && idx2+1<cxDIB) {
							sum1 += pTmp[idx1][idx2+1]; num++;
						}
						
						idx1 = i - nStep - 2;
						if(idx1>=0 && idx1<cyDIB) {
							sum1 += pTmp[idx1][idx2]; num++;
							if(idx2-1>=0 && idx2-1<cxDIB) {
								sum1 += pTmp[idx1][idx2-1]; num++;
							}
							if(idx2+1>=0 && idx2+1<cxDIB) {
								sum1 += pTmp[idx1][idx2+1]; num++;
							}
						}
						idx1 = i - nStep;
						if(idx1>=0 && idx1<cyDIB) {
							sum1 += pTmp[idx1][idx2]; num++;
							if(idx2-1>=0 && idx2-1<cxDIB) {
								sum1 += pTmp[idx1][idx2-1]; num++;
							}
							if(idx2+1>=0 && idx2+1<cxDIB) {
								sum1 += pTmp[idx1][idx2+1]; num++;
							}
						}
						tmp = (BYTE)(sum1/num);
					}
					else{
						idx1 = i - nStep - 1; idx2 = j - nStep;
						sum1 = pTmp[idx1][idx2]*_table1[nOrnt];
						num  = _table1[nOrnt];
						idx0 = 18*nOrnt;
						for(m=0; m<_table2[nOrnt]; m++) {
							y = idx1 + _table3[idx0+m];
							x = idx2 + _table4[idx0+m];
							if(y>=0 && y<cyDIB && x>=0 && x<cxDIB) {
								sum1 += pTmp[y][x] * _table5[idx0+m];
								num += _table5[idx0+m];
							}
							y = idx1 - _table3[idx0+m];
							x = idx2 - _table4[idx0+m];
							if(y>=0 && y<cyDIB && x>=0 && x<cxDIB) {
								sum1 += pTmp[y][x]*_table5[idx0+m];
								num += _table5[idx0+m];
							}
						}
						tmp = (BYTE)(sum1/num);
					}
					Img[idx] = tmp;
				}
			}
		}
	}
	free(pSum);
	for(i=0; i<nWindow; i++) free(pTmp[i]);
	free(pTmp);
}

void re_get_orient_image(BYTE *Img,BYTE *image_buffer2,BYTE *OrntImg,int cxDIB,int cyDIB,
					SINGULAR *pSingular,int nStep,int nMaxValue,
					int nMinValue, int nStartValue, int nMultiValue)
{
	int i, j, idx0, idx1, idx2, nSum[4], nDir[4];
	int *pSum, nWindow = 2*nStep+1, idx;
	int var0, var1, var2, var3, sum1, sum2;
	int buf, avg, nMax, t1, t2, t3, diff;
	BYTE tmp, nPrevVal;

	pSum = (int*)calloc(4*cxDIB, sizeof(int));
	if(pSum == NULL) return;

	for(i=0; i<cyDIB+nStep; i++) {
		if(i-1>=0 && i+1<cyDIB) {
			idx0 = i - 1;	idx1 = i;	idx2 = i + 1;
			for ( j=1; j<cxDIB-1; j++ ) {
				tmp = Img[idx1*cxDIB+j];
				pSum[j*4+0] += (abs(tmp-Img[idx1*cxDIB+j-1]) + abs(tmp-Img[idx1*cxDIB+j+1]));
				pSum[j*4+1] += (abs(tmp-Img[idx0*cxDIB+j-1]) + abs(tmp-Img[idx2*cxDIB+j+1]));
				pSum[j*4+2] += (abs(tmp-Img[idx0*cxDIB+j]) + abs(tmp-Img[idx2*cxDIB+j]));
				pSum[j*4+3] += (abs(tmp-Img[idx2*cxDIB+j-1]) + abs(tmp-Img[idx0*cxDIB+j+1]));
			}
		}
		if(i < nStep+1) continue;
		if(i >= nWindow+1) {
			idx0 = i-nWindow - 1; idx1 = i-nWindow; idx2 = i-nWindow + 1;
			for(j=1; j<cxDIB-1; j++) {
				tmp = Img[idx1*cxDIB+j];
				pSum[j*4+0] -= (abs(tmp-Img[idx1*cxDIB+j-1]) + abs(tmp-Img[idx1*cxDIB+j+1]));
				pSum[j*4+1] -= (abs(tmp-Img[idx0*cxDIB+j-1]) + abs(tmp-Img[idx2*cxDIB+j+1]));
				pSum[j*4+2] -= (abs(tmp-Img[idx0*cxDIB+j]) + abs(tmp-Img[idx2*cxDIB+j]));
				pSum[j*4+3] -= (abs(tmp-Img[idx2*cxDIB+j-1]) + abs(tmp-Img[idx0*cxDIB+j+1]));
			}
		}
		nSum[0] = nSum[1] = nSum[2] = nSum[3] = 0;
		for(j=0; j<cxDIB+nStep; j++) {
			if(j < cxDIB) {
				nSum[0] += pSum[4*j+0];
				nSum[1] += pSum[4*j+1];
				nSum[2] += pSum[4*j+2];
				nSum[3] += pSum[4*j+3];
			}
			if(j < nStep) continue;
			if(j >= nWindow) {
				nSum[0] -= pSum[4*(j-nWindow)+0];
				nSum[1] -= pSum[4*(j-nWindow)+1];
				nSum[2] -= pSum[4*(j-nWindow)+2];
				nSum[3] -= pSum[4*(j-nWindow)+3];
			}
			var0 = nDir[0] = nSum[0];
			var1 = nDir[1] = 71*nSum[1]/100;
			var2 = nDir[2] = nSum[2];
			var3 = nDir[3] = 71*nSum[3]/100;

			idx = (i-nStep)*cxDIB + j-nStep;
			nPrevVal = OrntImg[idx] & 0x80;
			nMax = nDir[0];
			if ( nMax < nDir[1] ) nMax = nDir[1];
			if ( nMax < nDir[2] ) nMax = nDir[2];
			if ( nMax < nDir[3] ) nMax = nDir[3];

			OrntImg[idx] = 45;
			sum1 = var2 + var3; sum2 = var2 + var1;
			if(sum1 < sum2) {
				sum2 = sum1;
				OrntImg[idx] = 75;
				nDir[0] = var1;	nDir[1] = var2;
				nDir[2] = var3;	nDir[3] = var0;
			}
			sum1 = var0 + var3;
			if(sum1 < sum2) {
				sum2 = sum1;
				OrntImg[idx] = 105;
				nDir[0] = var2;	nDir[1] = var3;
				nDir[2] = var0;	nDir[3] = var1;
			}
			sum1 = var0 + var1;
			if(sum1 < sum2) {
				OrntImg[idx] = 15;
				nDir[0] = var3;	nDir[1] = var0;
				nDir[2] = var1;	nDir[3] = var2;
			}
			avg = (nDir[0]+nDir[1]+nDir[2]+nDir[3]) - nMax*4;
			if(avg == 0) {
				OrntImg[idx] = 0x7F;
				buf = 0xFF;
			} else {
				diff = 15*(3*(nDir[3]-nDir[0])-nDir[1]+nDir[2])/avg;
				OrntImg[idx] += diff;
				if(OrntImg[idx] == 120) OrntImg[idx] = 0; 
				t1 = nDir[2]; t2 = nDir[3];
				if(t1 < nDir[1]) t2 = nDir[0]; 
				else t1 = nDir[1];
				t3 = ((t2-t1)*(15-abs(diff)))/30;
				if(t3 <= t1) t1 -= t3; else t1 = 0;
				t2 += t3;
				if(t2 == 0) buf = 0xFF;
				else buf = (int)(255*t1)/t2;
			}
			image_buffer2[idx] &= 0xF8;
			if(nPrevVal != 0) {
				OrntImg[idx] |= 0x80; continue;
			}
			if ( check_in_singular(pSingular,j-nStep,i-nStep,16) ) continue;
			if(buf > nMaxValue) buf = nMaxValue;
			if(buf < nMinValue) buf = 0;
			else buf -= nMinValue;
			buf *= nMultiValue;
			tmp = (BYTE)nStartValue;
			while(buf >= (nMaxValue-nMinValue)/2) {
				buf -= (nMaxValue-nMinValue); tmp++;
			}
			image_buffer2[idx] |= tmp;
		}
	}
	free(pSum);
}

void image_proc_03(BYTE *image_buffer2,BYTE *OrntImg,BYTE *Img,int cxDIB,int cyDIB)
{
	int nStep = 6, nWindow = 2*nStep+1, idx;
	int i, j, k, jj, nCount, idx0, idx1, idx2, nSum, nNum;
	BYTE orValue, tmp, **pTmp;
	BOOL bOut;
	
	pTmp = (BYTE**)malloc(sizeof(BYTE*)*cyDIB);
	if(pTmp == NULL) return;
	for(i=0; i<nWindow; i++) {
		pTmp[i] = (BYTE*)malloc(sizeof(BYTE)*cxDIB);
		if(pTmp[i] == NULL) {
			for(j=0; j<i; j++) free(pTmp[j]);
			free(pTmp); return;
		}
	}
	
	do {
		nCount = 0;
		for(i=0; i<cyDIB+nStep; i++) {
			if(i < cyDIB) {
				if(i >= nWindow) pTmp[i] = pTmp[i-nWindow];
				memcpy( pTmp[i],&Img[i*cxDIB],sizeof(BYTE)*cxDIB );
			}
			if(i < nStep) continue;
			if(i-nStep>=nStep && i<cyDIB) bOut = FALSE; else bOut = TRUE;
			for(j=nStep; j<cxDIB-nStep; j++) {
				idx = (i-nStep)*cxDIB + j;
				orValue = image_buffer2[idx] & 0x07;
				if(orValue == 0) continue;
				tmp = OrntImg[idx] & 0x7F;
				if(bOut == TRUE) {
					if(tmp == 0x7F) {
						idx1 = i - nStep;
						nSum  = pTmp[idx1][j+0];
						nSum += pTmp[idx1][j-1];
						nSum += pTmp[idx1][j+1]; nNum = 3;
						if(i-nStep-1>=0 && i-nStep-1<cyDIB) {
							idx1 = i - nStep - 1;
							nSum += pTmp[idx1][j-1];
							nSum += pTmp[idx1][j+0];
							nSum += pTmp[idx1][j+1]; nNum += 3;
						}
						if(i-nStep+1>=0 && i-nStep+1<cyDIB) {
							idx1 = i - nStep + 1;
							nSum += pTmp[idx1][j-1];
							nSum += pTmp[idx1][j+0];
							nSum += pTmp[idx1][j+1]; nNum += 3; 
						}
						Img[idx] = (BYTE)(nSum/nNum);
					}
					else{
						nSum = pTmp[i-nStep][j]*_table1[tmp];
						nNum = _table1[tmp];
						idx0 = tmp*18;
						for(k=0; k<_table2[tmp]; k++) {
							idx1 = i - nStep + _table3[idx0+k];
							idx2 = j + _table4[idx0+k];
							if(idx1>=0 && idx2>=0 && idx1<cyDIB && idx2<cxDIB) {
								nSum += pTmp[idx1][idx2]*_table5[idx0+k];
								nNum += _table5[idx0+k];
							}
							idx1 = i - nStep - _table3[idx0+k];
							idx2 = j - _table4[idx0+k];
							if(idx1>=0 && idx2>=0 && idx1<cyDIB && idx2<cxDIB) {
								nSum += pTmp[idx1][idx2]*_table5[idx0+k];
								nNum += _table5[idx0+k];
							}
						}
						Img[idx] = (BYTE)(nSum/nNum);
					}
				}
				else{
					if(tmp == 0x7F) {
						idx1 = i - nStep + 0;
						nSum = pTmp[idx1][j-1] + pTmp[idx1][j+0] + pTmp[idx1][j+1];
						idx1 = i - nStep - 1;
						nSum = pTmp[idx1][j-1] + pTmp[idx1][j+0] + pTmp[idx1][j+1];
						idx1 = i - nStep + 1;
						nSum = pTmp[idx1][j-1] + pTmp[idx1][j+0] + pTmp[idx1][j+1];
						Img[idx] = (BYTE)(nSum/9);
					}
					else{
						nSum = pTmp[i-nStep][j]*_table1[tmp];
						idx0 = tmp*18;
						for(k=0; k<_table2[tmp]; k++) {
							idx1 = i - nStep + _table3[idx0+k];
							idx2 = j + _table4[idx0+k];
							nSum += pTmp[idx1][idx2]*_table5[idx0+k];
							idx1 = i - nStep - _table3[idx0+k];
							idx2 = j - _table4[idx0+k];
							nSum += pTmp[idx1][idx2]*_table5[idx0+k];
						}
						Img[idx] = (BYTE)(nSum/112211);
					}
				}
				image_buffer2[idx] = (image_buffer2[idx] & 0xF8) | (orValue-1);
				nCount++;
			}
			for(j=0,jj=cxDIB-1; j<nStep; j++,jj--) {
				idx = (i-nStep)*cxDIB + j;
				orValue = image_buffer2[idx] & 0x07;
				if(orValue > 0) {
					tmp = OrntImg[idx] & 0x7F;
					if(tmp == 0x7F) {
						idx1 = i - nStep;
						nSum = pTmp[idx1][j]; nNum = 1;
						if(j-1>=0 && j-1<cxDIB) {
							nSum += pTmp[idx1][j-1]; nNum++;
						}
						if(j+1>=0 && j+1<cxDIB) {
							nSum += pTmp[idx1][j+1]; nNum++;
						}
						if(i-nStep-1>=0 && i-nStep-1<cyDIB) {
							idx1 = i - nStep - 1;
							nSum += pTmp[idx1][j]; nNum++;
							if(j-1>=0 && j-1<cxDIB) {
								nSum += pTmp[idx1][j-1]; nNum++;
							}
							if(j+1>=0 && j+1<cxDIB) {
								nSum += pTmp[idx1][j+1]; nNum++;
							}
						}
						if(i-nStep+1>=0 && i-nStep+1<cyDIB) {
							idx1 = i - nStep + 1;
							nSum += pTmp[idx1][j]; nNum++;
							if(j-1>=0 && j-1<cxDIB) {
								nSum += pTmp[idx1][j-1]; nNum++;
							}
							if(j+1>=0 && j+1<cxDIB) {
								nSum += pTmp[idx1][j+1]; nNum++;
							}
						}
						Img[idx] = (BYTE)(nSum/nNum);
					} else {
						idx1 = i - nStep;
						nSum = pTmp[idx1][j]*_table1[tmp];
						nNum = _table1[tmp];
						idx0 = tmp*18;
						for(k=0; k<_table2[tmp]; k++) {
							idx1 = i - nStep + _table3[idx0+k];
							idx2 = j + _table4[idx0+k];
							if(idx1>=0 && idx2>=0 && idx1<cyDIB && idx2<cxDIB) {
								nSum += pTmp[idx1][idx2]*_table5[idx0+k];
								nNum += _table5[idx0+k];
							}
							idx1 = i - nStep - _table3[idx0+k];
							idx2 = j - _table4[idx0+k];
							if(idx1>=0 && idx2>=0 && idx1<cyDIB && idx2<cxDIB) {
								nSum += pTmp[idx1][idx2]*_table5[idx0+k];
								nNum += _table5[idx0+k];
							}
						}
						Img[idx] = (BYTE)(nSum/nNum);
					}
					image_buffer2[idx] = (image_buffer2[idx] & 0xF8) | (orValue-1);
					nCount++;
				}
				idx = (i-nStep)*cxDIB + jj;
				orValue = image_buffer2[idx] & 0x07;
				if(orValue == 0) continue;
				tmp = OrntImg[idx] & 0x7F;
				if(tmp == 0x7F) {
					idx1 = i - nStep;
					nSum = pTmp[idx1][jj]; nNum = 1;
					if(jj-1>=0 && jj-1<cxDIB) {
						nSum += pTmp[idx1][jj-1]; nNum++;
					}
					if(jj+1>=0 && jj+1<cxDIB) {
						nSum += pTmp[idx1][jj+1]; nNum++;
					}
					if(i-nStep-1 >= 0) {
						idx1 = i - nStep - 1;
						nSum += pTmp[idx1][jj]; nNum++;
						if(jj-1>=0 && jj-1<cxDIB) {
							nSum += pTmp[idx1][jj-1]; nNum++;
						}
						if(jj+1>=0 && jj+1<cxDIB) {
							nSum += pTmp[idx1][jj+1]; nNum++;
						}
					}
					if(i-nStep+1>=0 && i-nStep+1<cyDIB) {
						idx1 = i - nStep + 1;
						nSum += pTmp[i-nStep+1][jj]; nNum++;
						if(jj-1>=0 && jj-1<cxDIB) {
							nSum += pTmp[idx1][jj-1]; nNum++;
						}
						if(jj+1>=0 && jj+1<cxDIB) {
							nSum += pTmp[idx1][jj+1]; nNum++;
						}
					}
					Img[idx] = (BYTE)(nSum/nNum);
				} else {
					nSum = pTmp[i-nStep][jj]*_table1[tmp];
					nNum = _table1[tmp];
					idx0 = tmp*18;
					for(k=0; k<_table2[tmp]; k++) {
						idx1 = i - nStep + _table3[idx0+k];
						idx2 = jj + _table4[idx0+k];
						if(idx1>=0 && idx1<cyDIB && idx2>=0 && idx2<cxDIB) {
							nSum += pTmp[idx1][idx2]*_table5[idx0+k];
							nNum += _table5[idx0+k];
						}
						idx1 = i - nStep - _table3[idx0+k];
						idx2 = jj - _table4[idx0+k];
						if(idx1>=0 && idx1<cyDIB && idx2>=0 && idx2<cxDIB) {
							nSum += pTmp[idx1][idx2]*_table5[idx0+k];
							nNum += _table5[idx0+k];
						}
					}
					Img[idx] = (BYTE)(nSum/nNum);
				}
				image_buffer2[idx] = (image_buffer2[idx] & 0xF8) | (orValue-1);
				nCount++;
			}
		}
	} while(nCount != 0);
	for(i=0; i<nWindow; i++) free(pTmp[i]);
	free(pTmp);
}

void get_block_data(FPVECTEX* pFPEx,int nCol,int nRow,BYTE* image_buffer3,int cxDIB,int cyDIB)
{
	int i, j, x, y;

	for ( i=0; i<nRow; i++ ){
		y = i*BLOCK_SIZE + BLOCK_SIZE/2;
		for ( j=0; j<nCol; j++ ){
			x = j*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( image_buffer3[y*cxDIB+x] >= 120 ) pFPEx->Block.Data[i*nCol+j] = 0xFF;
			else pFPEx->Block.Data[i*nCol+j] = image_buffer3[y*cxDIB+x];
		}
	}
	pFPEx->Block.nCol = nCol; pFPEx->Block.nRow = nRow;
}

void get_singular_block(BYTE *OrntImg,int cxDIB,int cyDIB, int bSize, int *nNum, int *pList, int *typeList)
{	
	BYTE *BlockDir;
	int i, j, x, y, k, kk, dif, delta, poincare, id1, id2, nNumTh = 64;
	int nCol = cxDIB/bSize, nRow = cyDIB/bSize;
	int tmpList[64], tmpType[64];
	char mX25[16] = { 2,2,2,1,0,-1,-2,-2,-2,-2,-2,-1,0,1,2,2 };
	char mY25[16] = { 0,-1,-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1 };
	char mX49[8] = { 3,3,0,-3,-3,-3,0,3 };
	char mY49[8] = { 0,-3,-3,-3,0,3,3,3 };
	char mX9[8] = { 1,1,0,-1,-1,-1,0,1 };
	char mY9[8] = { 0,-1,-1,-1,0,1,1,1 };

	BlockDir = (BYTE*)malloc(sizeof(BYTE)*nRow*nCol);
	if ( BlockDir == NULL ) return;
	memset( BlockDir, 0xFF, sizeof(BYTE)*nRow*nCol);
	for ( i=0; i<nRow; i++ ){
		y = (i*bSize + bSize/2)*cxDIB;
		for ( j=0; j<nCol; j++ ){
			x = j*bSize + bSize/2; BlockDir[i*nCol+j] = OrntImg[y+x];
		}
	}
	j = 0;
	for ( i=0; i<nRow*nCol; i++ ){
		if ( (BlockDir[i] & 0x80 ) != 0 ) continue;
		x = i % nCol; y = i / nCol; poincare = 0; 
		for ( k=0; k<16; k++ ){
			if ( k == 15 ) kk = 0; 
			else kk = k + 1;
			id1 = x + mX25[k]; id2 = y + mY25[k]; 
			if ( id1 < 0 || id1 >= nCol || id2 < 0 || id2 >= nRow ) break;
			id1 = x + mX25[kk]; id2 = y + mY25[kk]; 
			if ( id1 < 0 || id1 >= nCol || id2 < 0 || id2 >= nRow ) break;
			id1 = i+mX25[k] + mY25[k]*nCol;
			id2 = i+mX25[kk] + mY25[kk]*nCol;
			if ( (BlockDir[id1] & 0x80 ) != 0 || (BlockDir[id2] & 0x80) != 0 ) continue;
			dif = BlockDir[id1] - BlockDir[id2];
			if ( abs(dif) < 60 ) delta = dif;
			if ( dif <= -60 )    delta = 120 + dif;
			if ( dif >= 60 )     delta = -120 + dif;
			poincare += delta;
		}
		if ( k < 16 ) continue;
		if ( poincare == 120 || poincare == -120 ){
			if ( poincare == 120 ) tmpType[j] = 1;
			else tmpType[j] = 0;
			tmpList[j++] = i; 
			if ( j >= nNumTh ) break;
		}
	}

	for ( i=0; i<j; i++ ){
		if ( *nNum >= nNumTh ) break; 
		pList[*nNum] = tmpList[i]; 
		typeList[*nNum] = tmpType[i]; (*nNum)++;
		x = tmpList[i] % nCol; y = tmpList[i] / nCol; 
		for ( k=0; k<8; k++ ){
			id1 = x + mX49[k]; id2 = y + mY49[k];
			if ( id1 < 0 || id1 >= nCol || id2 < 0 || id2 >= nRow ){ 
				id2 = tmpList[i] + mX9[k] + mY9[k]*nCol;
				if ( *nNum >= nNumTh ){ i = 1000; break; }
				pList[*nNum] = id2; 
				typeList[*nNum] = tmpType[i]; (*nNum)++;
				break;
			}
			id1 = tmpList[i] + mX49[k] + mY49[k]*nCol;
			if ( (BlockDir[id1] & 0x80 ) != 0 ){
				id2 = tmpList[i] + mX9[k] + mY9[k]*nCol;
				if ( *nNum >= nNumTh ){ i = 1000; break; }
				pList[*nNum] = id2; 
				typeList[*nNum] = tmpType[i]; (*nNum)++; break;
			}
		}
	}
	free(BlockDir);
}

void check_core_cand(int x,int y,int nDev,int *xList,int *yList,int *dList,int *nCandNum)
{
	int i, disTH = 10*10, dx, dy;

	for ( i=0; i<*nCandNum; i++ ){
		dx = x - xList[i]; dy = y - yList[i];
		if ( dx*dx+dy*dy < disTH ) break;
	}
	if ( i < *nCandNum ){
		if ( dList[i] < nDev ){
			xList[i] = x; yList[i] = y; dList[i] = nDev; return;
		}
	}
	else{
		xList[*nCandNum] = x; yList[*nCandNum] = y;
		dList[*nCandNum] = nDev; (*nCandNum)++;
	}
}

int get_deviation(int xx,int yy,int rr,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, j, k, nDev = 0, sx, ex, sy, ey, num;
	BYTE p, p1;

	sx = (xx > rr)? xx - rr : 0;
	ex = ( xx + rr >= cxDIB )? cxDIB : xx+rr;
	sy = (yy > rr)? yy - rr : 0;
	ey = ( yy + rr >= cyDIB )? cyDIB : yy+rr;

	p = OrntImg[yy*cxDIB+xx];
	for ( i=sy; i<ey; i++ ){
		k = i*cxDIB + sx;
		for ( j=sx; j<ex; j++,k++ ){
			p1 = abs(p - OrntImg[k]);
			if ( p1 > 60 ) p1 = 120 - p1;
			nDev += p1;
		}
	}
	num = (ex - sx)*(ey - sy);
	if ( num == 0 ) return (-1);
	nDev /= num;
	return (nDev);
}

BOOL check_outof_point(int xx,int yy,BYTE *OrntImg,int cxDIB,int cyDIB,int Radius)
{
	int ix, iy, jx, jy;

	for ( iy=-Radius; iy<=Radius; iy++ ){
		jy = iy + yy;
		if ( jy < 0 || jy >= cyDIB ) return (TRUE);
		for ( ix=-Radius; ix<=Radius; ix++ ){
			jx = ix + xx;
			if ( jx < 0 || jx >= cxDIB ) return (TRUE);
			if ( OrntImg[jy*cxDIB+jx] >= 0x7F ) return (TRUE);
		}
	}
	return (FALSE);
}

int correct_orient_core(int coreDir,int *dSum,int maxV,int minV,int type,int minP2)
{
	int i, i1, i2, j1, dir1, e1, e2;

	if(type == 1 || minP2 < 0 || minP2 >= 240){
		dir1 = 80;
		i1 = i2 = e1 = e2 = 0;
		for(i = 1;i < 120;i++){
			if(i1 == 0){
				j1 = coreDir-i;
				if(j1 < 0) j1 += 240;
				if(dSum[j1] > dir1){ i1 = 1; e1 = i; }
			}
			if(i2 == 0){
				j1 = coreDir+i;
				if(j1 >= 240) j1 -= 240;
				if(dSum[j1] > dir1){ i2 = 1; e2 = i; }
			}
			if(i1 == 1 && i2 == 1) break;
		}
		if(e1 != 0 && e2 != 0){
			i1 = (e2 - e1)/2;
			coreDir += i1;
			if(coreDir < 0) coreDir += 240;
			if(coreDir >= 240) coreDir -= 240;
		}
		return (coreDir);
	}
	i1 = i2 = abs(coreDir - minP2);
	if(i2 >= 120) i2 = 240 - i2;
	if(i2 > 63 && i2 < 116) return (-1);
	coreDir = (coreDir + minP2)/2;
	if(i1 > 120){
		coreDir += 120;
		if(coreDir >= 240) coreDir -= 240;
	}
	return (coreDir);
}

int get_orient_core(int ix,int iy,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, i1, i2, iDir, ix1, iy1, ix2, iy2;
	int maxV, minV, minP, sinV, cosV, flag, midV;
	int idx = 2, idx1 = idx + 120;
	int dDirs[250], dSum[240];

	for ( iDir=0; iDir<120; iDir++,idx++,idx1++ ){
		cosV = (30 * (int)_table_03[iDir]) >> 14;
		sinV = (30 * (int)_table_04[iDir]) >> 14;
		ix1 = ix + cosV; iy1 = iy + sinV;
		flag = 0;
		if ( ix1 < 0 || ix1 >= cxDIB || iy1 < 0 || iy1 >= cyDIB ) flag = 1;
		else if ( OrntImg[iy1*cxDIB+ix1] >= 0x7F ) flag = 1;
		if ( flag == 1 ){
			i1 = op_func_01(cxDIB/2,cyDIB/2,ix,iy);
			i2 = op_func_01(ix1,iy1,ix,iy);
			i1 = abs(i1-i2);
			if ( i1 >= 120 ) i1 = 240-i1;
			ix1 = ix + cosV/3; iy1 = iy + sinV/3;
			flag = 0;
			if ( ix1 < 0 || ix1 >= cxDIB || iy1 < 0 || iy1 >= cyDIB ) flag = 1;
			else if ( OrntImg[iy1*cxDIB+ix1] >= 0x7F ) flag = 1;
			if ( flag == 1 ) return(-1);
		}
		dDirs[idx] = abs(OrntImg[iy1*cxDIB+ix1] - iDir);
		if ( dDirs[idx] >= 60 ) dDirs[idx] = 120 - dDirs[idx];
		ix2 = ix - cosV; iy2 = iy - sinV;
		if ( ix2 < 0 || ix2 >= cxDIB || iy2 < 0 || iy2 >= cyDIB) flag = 1;
		else if ( OrntImg[iy2*cxDIB+ix2] >= 0x7F ) flag = 1;
		if(flag == 1){
			i1 = op_func_01(cxDIB/2,cyDIB/2,ix,iy);
			i2 = op_func_01(ix2,iy2,ix,iy);
			i1 = abs(i1 - i2);
			if(i1 >= 120) i1 = 240 - i1;
			ix2 = ix - cosV/3; iy2 = iy - sinV/3;
			flag = 0;
			if ( ix2 < 0 || ix2 >= cxDIB || iy2 < 0 || iy2 >= cyDIB) flag = 1;
			else if ( OrntImg[iy2*cxDIB+ix2] >= 0x7F ) flag = 1;
			if ( flag == 1 ) return (-1);
		}
		dDirs[idx1] = abs(iDir - OrntImg[iy2*cxDIB+ix2]);
		if(dDirs[idx1] >= 60) dDirs[idx1] = 120 - dDirs[idx1];
	}
	for ( i=0; i<2; i++ ){
		dDirs[i] = dDirs[239+2-i];
		dDirs[240+2+i] = dDirs[i+2];
	}
	maxV = 0; minV = 200*5;
	for ( iDir=0; iDir<240; iDir++ ){
		dSum[iDir] = 0;
		for ( i=0; i<5; i++ ){
			dSum[iDir] += dDirs[iDir+i];
		}
		if(maxV < dSum[iDir]){ maxV = dSum[iDir]; }
		if(minV > dSum[iDir]){ minV = dSum[iDir]; minP = iDir; }
	}
	midV = (maxV - minV)/3;
	if ( midV < 10 ) return (-1);
	i = 0; flag = i1 = i2 = -1;
	for ( iDir=0; iDir<240; iDir++ ){
		ix1 = minP + iDir;
		if(ix1 >= 240) ix1 -= 240;
		if(flag == -1){
			if(dSum[ix1] > maxV-midV){
				i++;
				flag = 1;
				if(i1 >= 0 && i2 == -1){ i2 = ix1; }
			}
			continue;
		}
		if(dSum[ix1] < minV+midV){
			if(flag == 1 && i2 == -1){ i1 = ix1; }
			flag = -1;
		}
	}
	ix2 = -1;
	if(i == 2){
		if(i1 >= 0 && i1 < 240 && i2 >= 0 && i2 < 240 && i1 != i2){
			iy2 = dSum[i1];
			ix2 = i1;
			for(iDir = 0;(iDir < 240 && i1 != i2);iDir++){
				i1++;
				if(i1 >= 240) i1 = 0;
				if(iy2 > dSum[i1]){ iy2 = dSum[i1]; ix2 = i1; }
			}
		}
		else i = 1;
	}
	if(minV > 200 && i == 3) i = 1;
	if(i == 1 || i == 2){
		minP = correct_orient_core(minP,dSum,maxV,minV,i,ix2);
		i1 = minP+120;
		if ( i1 >= 240 ) i1 -= 240;
		if ( dDirs[i1+2] < 30 ) minP = -1;
	}
	else if(i == 3){
		iy2 = i2-i1;
		if(iy2 < 0) iy2 += 240;
		if(iy2 > 60) minP = -1;
		else minP = -2;
	}
	else minP = -1;
	return (minP);
}

void get_core_points(SINGULAR* SingularData,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int nBlkNum = 0, pList[64], typeList[64];
	int i, k, h, kk, sx, sy, ex, ey, id, p0, p1, p2, p3, p4;
	int bSize = 8, dirTH = 6;
	int nCol = cxDIB/bSize, nRow = cyDIB/bSize;
	int nDev, num, e1, orn, flag;
	int xlist[20], ylist[20], dlist[20], nCandNum = 0;
	int coreNum, delNum;
	BYTE *Block, val[9];
	TEMPCORE Cores[10], Deltas[10], tCore;

	get_singular_block(OrntImg,cxDIB,cyDIB,bSize,&nBlkNum,pList,typeList);
	Block = (BYTE*)malloc(sizeof(BYTE)*bSize*bSize);
	if ( Block == NULL ) return;
	for ( i=0; i<nBlkNum; i++ ){
		sy = (pList[i] / nCol)*bSize; ey = sy + bSize;
		sx = (pList[i] % nCol)*bSize; ex = sx + bSize;
		memset(Block,255,sizeof(BYTE)*bSize*bSize);
		for ( k=sy; k<ey; k++ ){
			for ( h=sx; h<ex; h++ ){
				p0 = OrntImg[k*cxDIB+h]; 
				if ( p0 >= 120 ) continue; 
				if ( p0 > dirTH && p0 < 120-dirTH ) continue;
				p1 = OrntImg[k*cxDIB+h-1];
				p2 = OrntImg[k*cxDIB+h+1];
				p3 = OrntImg[(k-1)*cxDIB+h];
				p4 = OrntImg[(k+1)*cxDIB+h];
				if ( p1 >= 120 || p2 >= 120 ) continue;
				if ( p3 >= 120 || p4 >= 120 ) continue;
				if ( p1 <= 60 && p2 <= 60 ){
					if ( (p3 <= 60 && p4 <= 60) || (p3 >= 60 && p4 >= 60) ) continue;
				}
				if ( p1 >= 60 && p2 >= 60 ){
					if ( (p3 <= 60 && p4 <= 60) || (p3 >= 60 && p4 >= 60) ) continue;
				}
				id = (k-sy)*bSize + h-sx;
				Block[id] = 0;
				h++;
			}
		}
		for ( k=0; k<bSize; k++ ){
			for ( h=0; h<bSize; h++){
				if ( Block[k*bSize+h] != 0 ) continue;
				if ( h+1 >= bSize ){ val[0] = val[8] = 255; }
				else{ val[0] = val[8] = Block[k*bSize+h+1]; }
				if ( h+1 >= bSize || k+1 >= bSize ) val[1] = 255;
				else val[1] = Block[(k+1)*bSize+h+1];
				if ( k+1 >= bSize ) val[2] = 255;
				else val[2] = Block[(k+1)*bSize+h];
				if ( h-1 < 0 || k+1 >= bSize ) val[3] = 255;
				else val[3] = Block[(k+1)*bSize+h-1];
				if ( h-1 < 0 ) val[4] = 255;
				else val[4] = Block[k*bSize+h-1];
				if ( h-1 < 0 || k-1 < 0 ) val[5] = 255;
				else val[5] = Block[(k-1)*bSize+h-1];
				if ( k-1 < 0 ) val[6] = 255;
				else val[6] = Block[(k-1)*bSize+h];
				if ( h+1 >= bSize || k-1 < 0 ) val[7] = 255;
				else val[7] = Block[(k-1)*bSize+h+1];
				num = 0;
				for ( kk=0; kk<8; kk++ ){
					if ( val[kk] == 0 && val[kk+1] == 255 ) num++;
				}
				if ( num > 1 ) continue;
				nDev = get_deviation(sx+h,sy+k,12,OrntImg,cxDIB,cyDIB);
				if ( nDev < 10 ) continue;
				check_core_cand(sx+h,sy+k,nDev,xlist,ylist,dlist,&nCandNum);
				if ( nCandNum >= 10 ){ k = 100; break; }
			}
		}
		if ( k >= 100 ) break;
	}
	free(Block);
	coreNum = delNum = 0;
	for ( i=0; i<nCandNum; i++ ){
		e1 = get_orient_core(xlist[i],ylist[i],OrntImg,cxDIB,cyDIB);
		if ( e1 >= 0 || e1 == -2 ){
			k = xlist[i] / bSize; h = ylist[i] / bSize;
			for ( kk=0; kk<nBlkNum; kk++ ){
				if ( pList[kk] == nCol*h+k ){ id = kk; break;}
			}
		}
		if ( e1 < 0 ){
			if ( e1 == -2 && typeList[id] == 0){
				Deltas[delNum].x = xlist[i];
				Deltas[delNum].y = ylist[i];
				Deltas[delNum++].val = dlist[i];
			}
			continue;
		}
		if ( typeList[id] == 0 ) continue;
		Cores[coreNum].x = xlist[i]; Cores[coreNum].y = ylist[i];
		Cores[coreNum].dir = e1; Cores[coreNum++].val = dlist[i];
	}
	for ( i=0; i<coreNum-1; i++ ){
		for ( k=i+1; k<coreNum; k++ ){
			if ( Cores[k].val == 0 ) continue;
			if ( Cores[k].val > Cores[i].val ){
				tCore = Cores[i]; Cores[i] = Cores[k]; Cores[k] = tCore;
			}
			sx = Cores[k].x - Cores[i].x; sy = Cores[k].y - Cores[i].y;
			if ( sx*sx+sy*sy < 100 ) Cores[k].val = 0;
		}
	}
	e1 = 0;
	for ( i=0; i<coreNum; i++ ){
		if ( Cores[i].val == 0 ) break;
		Cores[e1++] = Cores[i];
		if ( e1 == 2 ) break;
	}
	coreNum = e1;
	if ( coreNum == 2 ){
		flag = 0;
		e1 = Cores[0].dir - Cores[1].dir;
		if ( e1 < 0 ) e1 += 240;
		if ( e1 < 65 || e1 > 175 ) flag = 1;
		else{
			orn = op_func_01(Cores[0].x,Cores[0].y,Cores[1].x,Cores[1].y);
			e1 = abs(Cores[0].dir - orn);
			if ( e1 > 120 ) e1 = 240 - e1;
			if ( e1 < 60 ) flag = 1;
		}
		if ( flag == 1 ){
			coreNum = 1;
			if ( Cores[1].val > 30 ){
				if ( check_outof_point(Cores[1].x,Cores[1].y,OrntImg,cxDIB,cyDIB,15) == FALSE ){
					if ( check_outof_point(Cores[0].x,Cores[0].y,OrntImg,cxDIB,cyDIB,15) == TRUE ){
						Cores[0] = Cores[1];
					}
				}
			}
		}
		else{
			e1 = (Cores[0].x - Cores[1].x) * (Cores[0].x - Cores[1].x);
			e1 += (Cores[0].y - Cores[1].y) * (Cores[0].y - Cores[1].y);
			if( e1 < 30*30 ){
				orn = op_func_01(Cores[1].x,Cores[1].y , Cores[0].x,Cores[0].y);
				Cores[0].dir = orn;
				orn += 120;
				if(orn >= 240)	orn -= 240;
				Cores[1].dir = orn;
			}
		}
	}
	for ( i=0; i<coreNum; i++ ){
		SingularData->nX[i] = Cores[i].x;
		SingularData->nY[i] = Cores[i].y;
		SingularData->nOrient[i] = Cores[i].dir;
		SingularData->nType[i] = 1;
	}
	SingularData->nNumber = coreNum;

	for ( i=0; i<delNum; i++ ){
		for ( k=i+1; k<delNum; k++ ){
			if ( Deltas[k].val == 0 ) continue;
			if ( Deltas[k].val > Deltas[i].val ){
				tCore = Deltas[i]; Deltas[i] = Deltas[k]; Deltas[k] = tCore;
			}
			sx = Deltas[i].x - Deltas[k].x; sy = Deltas[i].y - Deltas[k].y;
			if ( sx*sx+sy*sy < 80*80 ) Deltas[k].val = 0;
		}
	}
	e1 = 0;
	for ( i=0; i<delNum; i++ ){
		if ( Deltas[i].val == 0 ) continue;
		Deltas[e1++] = Deltas[i];
		if ( e1 == 2 ) break;
	}
	delNum = e1;
	for ( i=0; i<delNum; i++ ){
		SingularData->nX[coreNum+i] = Deltas[i].x;
		SingularData->nY[coreNum+i] = Deltas[i].y;
		SingularData->nOrient[coreNum+i] = 255;
		SingularData->nType[coreNum+i] = DELTA;
	}
	SingularData->nNumber += delNum;
}

BOOL check_whorl(SINGULAR* SData,MAINLINE* mLine)
{
	int i, n, n1, n2, cx, cy, dir, dir1, dx ,dy,
		x0, y0, x1, x2, y1, y2, d1, d2, dd, minV;
	int cross_p[4], index_p[4];

	x1 = SData->nX[0]; x2 = SData->nX[1];
	y1 = SData->nY[0]; y2 = SData->nY[1];
	cx = (x1 + x2) / 2; cy = (y1 + y2) / 2;

	dx = x2 - x1; dy = y2 - y1;
	d2 = dx*dx+dy*dy;
	d1 = op_func_02(d2);
	if(d1 == 0) return(TRUE);

	d2 = dx*cx+dy*cy;
	dir = op_func_01(SData->nX[1],SData->nY[1],SData->nX[0],SData->nY[0]);
	dir += 60;
	if ( dir >= 240 ) dir -= 240;
	for ( n=0; n<4; n++ ){
		cross_p[n] = -1000;
		index_p[n] = 0;
		n1 = mLine->nNumbers[n];
		minV = 10;
		for(i = 0;i < n1;i++){
			x0 = mLine->points_x[n][i];
			y0 = mLine->points_y[n][i];
			dd = abs(dx*x0+dy*y0-d2)/d1;
			if(minV > dd){
				minV = dd;
				dir1 = op_func_01(x0,y0,cx,cy);
				dir1 = abs(dir1-dir);
				if(dir1 >= 120) dir1 = 240-dir1;
				n2 = (x0-cx)*(x0-cx)+(y0-cy)*(y0-cy);
				dd = op_func_02(n2);
				if(dir1 > 60) dd = -dd;
				if(index_p[n] != 0 && i-index_p[n] > 3){
					break;
				}
				cross_p[n] = dd;
				index_p[n] = i;
				if(minV == 0) break;
			}
		}
	}
	d2 = 10;
	if(cross_p[0] > -900 && cross_p[2] > -900 && 
		cross_p[1] > -900 && cross_p[3] > -900){
		if(min(cross_p[0],cross_p[1]) > max(cross_p[2],cross_p[3]) || 
			min(cross_p[2],cross_p[3]) > max(cross_p[0],cross_p[1])) return(FALSE);
		if(d1 > 80){
			d2 = min(abs(cross_p[0]-cross_p[1]),abs(cross_p[2]-cross_p[3]))/2;
			dd = abs(cross_p[0]-cross_p[2]);
			if(dd < d2) return(FALSE);
			dd = abs(cross_p[1]-cross_p[3]);
			if(dd < d2) return(FALSE);
		}
		d2 = min(abs(cross_p[0]-cross_p[1]),abs(cross_p[2]-cross_p[3]))/3;
	}
	if(cross_p[0] > -900 && cross_p[2] > -900){
		dd = abs(cross_p[0]-cross_p[2]);
		if(dd < d2) return(FALSE);
	}
	if(cross_p[1] > -900 && cross_p[3] > -900){
		dd = abs(cross_p[1]-cross_p[3]);
		if(dd < d2) return(FALSE);
	}
	return(TRUE);
}

int get_distance_to_line(MAINLINE *mLine,int x,int y,int num)
{
	int i, dx, dy, len, maxlen = 100000;

	for ( i=0; i<mLine->nNumbers[num]; i++ ){
		dx = x - mLine->points_x[num][i];
		dy = y - mLine->points_y[num][i];
		len = dx*dx + dy*dy;
		if ( len < maxlen ) maxlen = len;
	}
	maxlen = op_func_02(maxlen);
	return (maxlen);
}

int check_near_line(MAINLINE *mLine,int num,int disTh)
{
	int i, x, y, dd, nn = num+1, number = 1000;

	for ( i=10; i<mLine->nNumbers[nn]; i++ ){
		x = mLine->points_x[nn][i];
		y = mLine->points_y[nn][i];
		dd = get_distance_to_line(mLine,x,y,num);
		if ( dd < disTh ){ number = i; break; }
	}
	return (number);
}

int check_line_lr(int Line_x0,int Line_y0,int Line_x1,int Line_y1,
						short *points_x, short *points_y, int nNumber)
{
	int i, dx, dy, dd, numPlus = 0, numMis = 0;

	dx = Line_x1 - Line_x0; dy = Line_y1 - Line_y0;
	for ( i=0; i<nNumber; i++ ){
		dd = dy*(points_x[i]-Line_x0)-(points_y[i]-Line_y0)*dx;
		if ( dd < 0 ) numMis++;
		if ( dd > 0 ) numPlus++;
	}
	if ( numMis > 0 && numPlus == 0 ) return (-1);
	if ( numMis == 0 && numPlus > 0 ) return (1);
	return (0);
}

BOOL check_arch(MAINLINE *mLine,COREITEMEX *Core,COREITEMEX *Delta)
{
	int i, j, ix, iy, ix0, iy0, d1, d2, dd, coreD, disTH = 15*15;

	ix = Core->x - Delta->x; iy = Core->y - Delta->y;
	coreD = ix*ix+iy*iy;

	ix0 = mLine->points_x[0][mLine->nNumbers[0]-1] - Core->x;
	iy0 = mLine->points_y[0][mLine->nNumbers[0]-1] - Core->y;
	ix = mLine->points_x[1][mLine->nNumbers[1]-1] - Core->x;
	iy = mLine->points_y[1][mLine->nNumbers[1]-1] - Core->y;
	d1 = ix0*ix0 + iy0*iy0; d2 = ix*ix + iy*iy;

	dd = d1;
	if ( dd > d2 ) dd = d2;

	if ( coreD*coreD < dd*dd && mLine->nNumbers[0] > 20 && mLine->nNumbers[1] > 20 ){
		d1 = check_line_lr(
									Core->x,Core->y,Delta->x,Delta->y,
									&mLine->points_x[0][10],
									&mLine->points_y[0][10],
									mLine->nNumbers[0]-10
								   );
		d2 = check_line_lr(
									Core->x,Core->y,Delta->x,Delta->y,
									&mLine->points_x[1][10],
									&mLine->points_y[1][10],
									mLine->nNumbers[1]-10
								   );
		if ( d1*d2 < 0 ) return TRUE;
	}
	for ( i=0; i<2; i++ ){
		for ( j=0; j<mLine->nNumbers[i]; j++ ){
			ix = mLine->points_x[i][j] - Delta->x;
			iy = mLine->points_y[i][j] - Delta->y;
			dd = ix*ix + iy*iy;
			if ( dd < disTH-10 ) return TRUE;
		}
	}
	return FALSE;
}

int get_type_line(LPFPVECTEX pVectEx,SINGULAR* SData,BYTE *OrntImg,
				  BYTE *Img,int cxDIB,int cyDIB)
{
	int i, i1, i2, j, x, y, dir, ix, iy, idir, ix0, iy0;
	int k, kk, dd, d1, d2, th, dd0 = 1;
	int directX[8], directY[8], dth = 4;
	int dx, dy, xx, yy, diff01, diff1, diff2, sDir[2];
	int co_x, co_y, co_dir, flag, circleFlag;
	int nCoreNumber = 0, nDeltaNumber = 0;
	int ccx, ccy, ddx, ddy;
	COREITEMEX Core[10], Delta[10];
	MAINLINE mLine;

	for(i = 0;i < 4;i++){
		pVectEx->MainLine.points_x[i] = 0;
		pVectEx->MainLine.points_y[i] = 0;
	}

	for ( i=0; i<SData->nNumber; i++ ){
		if ( SData->nType[i] == 1 ){
			Core[nCoreNumber].x = SData->nX[i];
			Core[nCoreNumber].y = SData->nY[i];
			Core[nCoreNumber++].dir = (BYTE)SData->nOrient[i];
		}
		else{
			Delta[nDeltaNumber].x = SData->nX[i];
			Delta[nDeltaNumber].y = SData->nY[i];
			Delta[nDeltaNumber++].dir = (BYTE)SData->nOrient[i];
		}
	}
	if ( nCoreNumber == 0 ) return (10);
	if(nCoreNumber == 2){
		ix = Core[0].x - Core[1].x; iy = Core[0].y - Core[1].y;
		if ( ix*ix+iy*iy < 50*50 ) return (1);
	}
	if(nCoreNumber == 1 && nDeltaNumber == 1){
		ix = Core[0].x - Delta[0].x; iy = Core[0].y - Delta[0].y;
		if ( ix*ix+iy*iy < 50*50 ) dth = 8;
	}

	mLine.nNumbers[0] = mLine.nNumbers[1] = 0;
	mLine.nNumbers[2] = mLine.nNumbers[3] = 0;
	circleFlag = 0;
	memset( Img, 0, sizeof(BYTE)*cxDIB*cyDIB );

	for ( i=0; i<nCoreNumber; i++ ){
		co_dir = SData->nOrient[i];
		co_x = SData->nX[i]; co_y = SData->nY[i];
		x = co_x - dth*_table_03[co_dir]/16384;
		y = co_y - dth*_table_04[co_dir]/16384;

		sDir[1] = co_dir + 60;
		if ( sDir[1] >= 240 ) sDir[1] -= 240;
		if ( sDir[1] == 240 ) sDir[1] = 0;
		sDir[0] = co_dir - 60;
		if ( sDir[0] < 0 ) sDir[0] += 240;
		if ( sDir[0] == 240 ) sDir[0] = 0;
		kk = i*2;
		for ( j=0; j<2; j++ ){
			idir = sDir[j];
			ix = x + 5*_table_03[idir]/16384;
			iy = y + 5*_table_04[idir]/16384;
			for ( k=0; k<8; k++ ){ directX[k] = x; directY[k] = y; }
			xx = 64*ix; yy = 64*iy;
			k = flag = 0;
			while (1){
				dd = dd0;
				while (1){
					dx = 64*dd*_table_03[idir]/16384;
					dy = 64*dd*_table_04[idir]/16384;
					if (dx>=64 || dy>=64 || dx<=-64 || dy<=-64) break;
					dd++;
				}
				ix0 = ix; iy0 = iy; 
				xx += dx; yy += dy;
				if ( xx < 0 || yy < 0 ) break;
				ix = xx/64; iy = yy/64;
				for ( i1=0; i1<7; i1++ ){
					directX[i1] = directX[i1+1];
					directY[i1] = directY[i1+1];
				}
				directX[7] = ix; directY[7] = iy;
				if ( ix<0 || ix>=cxDIB || iy<0 || iy>=cyDIB) break;
				if ( OrntImg[iy*cxDIB+ix] >= 120 ) break;
				if ( Img[iy*cxDIB+ix] == kk+j+1 ){
					if ( k > 10 && nCoreNumber == 1 ) circleFlag = 1;
					break;
				}
				mLine.points_x[kk+j][k] = ix;
				mLine.points_y[kk+j][k] = iy;
				k++;
				if ( k > 95 ) break;

				idir = abs(op_func_01(ix,iy,co_x,co_y) - co_dir);
				if ( idir >= 120 ) idir = 240 - idir;
				dx = co_x - ix; dy = co_y - iy;
				dd = dx*dx + dy*dy;
				if ( flag == 0 ){
					if ( idir < 60 && dd > 25*25 ) flag = 1;
				}
				else if ( idir > 110 || dd < 20*20 ){
					circleFlag = 1;
					if ( k > 10 ){
						ix0 = mLine.points_x[kk+j][k-3];
						iy0 = mLine.points_y[kk+j][k-3];
						idir = op_func_01(ix0,iy0,co_x,co_y) - co_dir;
						if ( idir < 0 ) idir += 240;
						if ( idir < 120 ) i1 = 1;
						else i1 = 0;
						if ( (idir>60 && idir<180) && j==i1 ) circleFlag = 2;
					}
					break;
				}
				Img[iy*cxDIB+ix] = kk+j+1;
				i1 = abs(OrntImg[iy*cxDIB+ix]-OrntImg[iy0*cxDIB+ix0]);
				if ( i1 >= 60 ) i1 = 120-i1;
				if ( k > 20 && i1 > 45 ){ circleFlag = 3; break; }
				idir = op_func_01(directX[7],directY[7],directX[0],directY[0]);
				i1 = abs(idir-OrntImg[iy*cxDIB+ix]);
				idir = OrntImg[iy*cxDIB+ix];
				if ( i1 >= 120 ) i1 = 240 - i1;
				if ( k > 10 && i1 > 45 && i1 < 75 ){
					circleFlag = 3; break;
				}
				if ( i1 > 60 ) idir += 120;
			}
			mLine.nNumbers[kk+j] = k;
		}
		if ( mLine.nNumbers[kk] > 30 && mLine.nNumbers[kk+1] > 30 ){
			ix = (mLine.points_x[kk][30]+mLine.points_x[kk+1][30])/2;
			iy = (mLine.points_y[kk][30]+mLine.points_y[kk+1][30])/2;
			idir = op_func_01(ix,iy,x,y);
			SData->nOrient[i] = idir; Core[i].dir = idir;
			co_dir = idir;
		}
		else{
			k = mLine.nNumbers[kk];
			if ( k > mLine.nNumbers[kk+1] ) k = mLine.nNumbers[kk+1];
			if ( k >= 15 ){
				ix = (mLine.points_x[kk][k-1]+mLine.points_x[kk+1][k-1])/2;
				iy = (mLine.points_y[kk][k-1]+mLine.points_y[kk+1][k-1])/2;
				idir = op_func_01(ix,iy,x,y);
				SData->nOrient[i] = idir; Core[i].dir = idir;
				co_dir = idir;
			}
		}
	}

	for(i = 0;i < 4;i++){
		pVectEx->MainLine.points_x[i] = 0;
		pVectEx->MainLine.points_y[i] = 0;
		if(mLine.nNumbers[i] > 0){
			pVectEx->MainLine.points_x[i] = mLine.points_x[i][mLine.nNumbers[i]-1];
			pVectEx->MainLine.points_y[i] = mLine.points_y[i][mLine.nNumbers[i]-1];
		}
	}

	if ( nCoreNumber == 2 ){
		if ( mLine.nNumbers[0] < 10 || mLine.nNumbers[1] < 10 ||
			 mLine.nNumbers[2] < 10 || mLine.nNumbers[3] < 10 ) return(3);
		dir = op_func_01(Core[1].x,Core[1].y,Core[0].x,Core[0].y);
		diff1 = abs(dir - SData->nOrient[0]);
		if ( diff1 >= 120 ) diff1 -= 120;
		if ( diff1 >= 60 ) diff1 = 120 - diff1;
		diff2 = abs(dir - SData->nOrient[1]);
		if ( diff2 >= 120 ) diff2 -= 120;
		if ( diff2 >= 60 ) diff2 = 120 - diff2;
		dir = (diff1+diff2)/2;
		if ( dir >= 17 ) return (2);
		if ( check_whorl(SData,&mLine) == TRUE ) return (0);
		else return (2);
	}

	if ( mLine.nNumbers[0] < 40 || mLine.nNumbers[1] < 40 ){
		if ( circleFlag > 0 ){
			if ( mLine.nNumbers[0] < 40 && mLine.nNumbers[1] < 40 )	return(3);
			if ( circleFlag == 1 ) return (0);
			if ( circleFlag == 2 ) return (1);
			return(3);
		}
		return(8);
	}
	if ( circleFlag > 0 ){
		kk = 0;
		for ( i=0; i<2; i++ ){
			for ( j=0; j<mLine.nNumbers[i]; j++ ){
				dx = mLine.points_x[i][j] - co_x;
				dy = mLine.points_y[i][j] - co_y;
				diff1 = dx*dx + dy*dy;
				if ( kk < diff1 ) kk = diff1;
			}
		}
		if ( kk < 50*50 ) return (1);
		if ( circleFlag == 1 ) return (0);
		if ( circleFlag == 3 ) return (3);
		return (2);
	}
	if ( nDeltaNumber == 2 ) return (0);
	ix0 = mLine.points_x[0][mLine.nNumbers[0]-1];
	iy0 = mLine.points_y[0][mLine.nNumbers[0]-1];
	ix = mLine.points_x[1][mLine.nNumbers[1]-1];
	iy = mLine.points_y[1][mLine.nNumbers[1]-1];
	dx = ix0-ix;
	dy = iy0-iy;
	diff01 = dx*dx+dy*dy;
	i1 = op_func_01(ix0,iy0,co_x,co_y);
	i2 = op_func_01(ix,iy,co_x,co_y);
	diff2 = abs(i1-i2);
	if(diff2 >= 120) diff2 = 240-diff2;

	if ( mLine.nNumbers[0] > 20 && mLine.nNumbers[1] > 20 && diff01*diff01 > 80*80 && diff2 > 50 ){
		diff1 = co_dir - i1;
		if ( diff1 < 0 ) diff1 += 240;
		if ( diff1 > 20 && diff1 < 60 ){
			diff1 = i2 - co_dir;
			if ( diff1 < 0 ) diff1 += 240;
			if ( diff1 > 20 && diff1 < 60 ) return(7);
		}
	}
	dx = ix0 - co_x; dy = iy0 - co_y;
	d1 = op_func_02(dx*dx+dy*dy);
	dx = ix - co_x; dy = iy - co_y;
	d2 = op_func_02(dx*dx+dy*dy);
	flag = 0;
	ccx = Core[0].x; ccy = Core[0].y;
	ddx = Delta[0].x; ddy = Delta[0].y;
	if ( nDeltaNumber == 1 ){
		dir = op_func_01(ddx,ddy,ccx,ccy);
		idir = abs(dir - co_dir);
		if ( idir >= 120 ) idir = 240-idir;
		if ( idir > 80 ){
			dd = d1;
			if ( dd < d2 ) dd = d2;
			if ( dd < 120 ) return (0);
			nDeltaNumber = 0;
		}
		if ( nDeltaNumber == 1 ){
			if ( idir < 10 ) return (7);
			idir = op_func_02((co_x-ddx)*(co_x-ddx)+(co_y-ddy)*(co_y-ddy));
			if ( idir >= 210 ) return (0);
			if ( nCoreNumber == 1 ){
				if ( check_arch(&mLine,&Core[0],&Delta[0]) == TRUE ) return (7);
			}
			dir -= co_dir;
			if ( dir < 0 ) dir += 240;
			if ( dir >= 120 ){
				if ( dir >= 130 ) flag = 4;
			}
			else if ( dir < 110 ) flag = 5;
		}
	}
	if ( d1 < 50 || d2 < 50 ) return (8);
	dd = check_near_line(&mLine,0,15);
	if ( dd < 1000 ){
		dx = mLine.points_x[1][dd] - co_x;
		dy = mLine.points_y[1][dd] - co_y;
		idir = dx*dx + dy*dy;
		dx = op_func_02(idir);
		dy = d1;
		if ( dy > d2 ) dy = d2;
		dy = op_func_02(dy);
		if ( dx < dy-20 ) return (0);
	}
	dd = d1;
	if ( dd < d2 ) dd = d2;
	if ( dd < 100 && diff01 < 16 ) return (8);
	th = (200 - dd)/30 + 5;
	if ( th > 10 ) th = 10;
	else if ( th < 5 ) th = 5;
	idir = i1 - i2;
	if ( idir < 0 ) idir += 240;
	if ( idir > th && idir < 120 ) return (0);
	if ( dd < 100 ){
		idir = abs(i1-i2);
		if ( idir >= 120 ) idir = 240 - idir;
		if ( idir < th/2 ) return (8);
	}
	dir = co_dir;
	if ( flag != 5 ){
		idir = i1 - dir;
		if ( idir < 0 ) idir += 240;
		if ( idir < 120 ){
			if ( idir > 10 && idir < 120 ){
				idir = abs(dir-i2);
				if ( idir >= 120 ) idir = 240 - idir;
				if ( idir > 70 ) return (2);
			}
			idir = abs(dir-i2);
			if ( idir >= 120 ) idir = 240 - idir;
			if ( idir < 4 ) return (0);
			if ( idir < 15 && d2 <= 80 ) return (8);
			return (4);
		}
		if(flag == 4){
			dx = (co_x-ix0)*(ddx-ix0)+(co_y-iy0)*(ddy-iy0);
			dy = (co_x-ix)*(ddx-ix)+(co_y-iy)*(ddy-iy);
			if ( dx > 0 && dy > 0) return (4);
		}
	}
	if ( flag != 4 ){
		idir = dir - i2;
		if ( idir < 0 ) idir += 240;
		if ( idir < 120 ){
			if ( idir > 10 && idir < 120 ){
				idir = abs(dir-i1);
				if ( idir >= 120 ) idir = 240 - idir;
				if ( idir > 70 ) return (2);
			}
			idir = abs(i1-dir);
			if ( idir >= 120 ) idir = 240-idir;
			if ( idir < 4 ) return (0);
			if ( idir < 15 && d1 <= 80 ) return (8);
			return(5);
		}
		if ( flag == 5 ){
			dx = (co_x-ix0)*(ddx-ix0)+(co_y-iy0)*(ddy-iy0);
			dy = (co_x-ix)*(ddx-ix)+(co_y-iy)*(ddy-iy);
			if ( dx > 0 && dy > 0 ) return (5);
		}
	}
	return (8);
}

void copy_core(SINGULAR* SingularData,LPFPVECTEX FPEx)
{
	int i, num = 0;

	for ( i=0 ; i<SingularData->nNumber; i++ ){
		FPEx->Core.item[num].x = SingularData->nX[i];
		FPEx->Core.item[num].y = SingularData->nY[i];
		FPEx->Core.item[num].dir = (BYTE)SingularData->nOrient[i];
		FPEx->Core.item[num].kind = (BYTE)SingularData->nType[i];
		num++;
	}
	FPEx->Core.nNumber = num;
}

void image_proc_04(BYTE *Img,int cxDIB,int cyDIB)
{
	int i, j, nSum, *pSum;
	BYTE *pTmp1, *pTmp2, *pTmp3, *pTmp;
	pTmp1 = (BYTE*)malloc(cxDIB);
	if(pTmp1 == NULL) return;
	pTmp2 = (BYTE*)malloc(cxDIB);
	if(pTmp2 == NULL) { free(pTmp1); return; }
	pTmp3 = (BYTE*)malloc(cxDIB);
	if(pTmp3 == NULL) {
		free(pTmp1); free(pTmp2); return;
	}
	pSum = (int*)calloc(cxDIB, sizeof(int));
	if(pSum == NULL) {
		free(pTmp1); free(pTmp2); free(pTmp3); return;
	}
	for(i=0; i<cyDIB+1; i++) {
		if(i >= 3) {
			for(j=0; j<cxDIB; j++) pSum[j] -= pTmp1[j];
		}
		if(i < cyDIB) {
			memcpy(pTmp1, &Img[i*cxDIB], cxDIB);
			for(j=0; j<cxDIB; j++) pSum[j] += pTmp1[j];
		}
		pTmp = pTmp1; pTmp1 = pTmp2; pTmp2 = pTmp3; pTmp3 = pTmp;
		if(i < 2) continue;
		nSum = 0;
		for(j=0; j<cxDIB+1; j++) {
			if(j >= 3) nSum -= pSum[j-3];
			if(j < cxDIB) nSum += pSum[j];
			if(j < 2) continue;
			Img[(i-1)*cxDIB+j-1] = (nSum < 1152) ? 0x00:0xFF;
		}
	}
	free(pTmp1); free(pTmp2); free(pTmp3); free(pSum);
	for(i=0; i<cxDIB; i++) {
		Img[0*cxDIB+i] = 0xFF;	Img[(cyDIB-1)*cxDIB+i] = 0xFF;
	}
	for(i=0; i<cyDIB; i++) {
		Img[i*cxDIB+0] = 0xFF; Img[i*cxDIB+cxDIB-1] = 0xFF;
	}
}

void remove_hole(BYTE* OrntImg,BYTE* Img,int cxDIB,int cyDIB)
{
	int i, j, k, ix, iy, wNum, whiteX[50], whiteY[50];
	int numTH = 32;
	BOOL flag;

	for(i = 0;i < cyDIB;i++){
		for(j = 0;j < cxDIB;j++){
			if((OrntImg[i*cxDIB+j] & 0x80) != 0) continue;
			if(Img[i*cxDIB+j] != 255) continue;
			wNum = 1;
			whiteX[0] = j;
			whiteY[0] = i;
			Img[i*cxDIB+j] = 0;
			flag = TRUE;
			for(k = 0;k < wNum;k++){
				ix = whiteX[k];
				iy = whiteY[k];
				if(ix < 1 || ix > cxDIB-2 || iy < 1 || iy > cyDIB-2){
					flag = FALSE;
					break;
				}
				if(Img[(iy-1)*cxDIB+ix] == 255){
					if(iy <= i){
						flag = FALSE;
						break;
					}
					Img[(iy-1)*cxDIB+ix] = 0;
					whiteX[wNum] = ix;
					whiteY[wNum] = iy-1;
					wNum++;
				}
				if(Img[iy*cxDIB+ix+1] == 255){
					Img[iy*cxDIB+ix+1] = 0;
					whiteX[wNum] = ix+1;
					whiteY[wNum] = iy;
					wNum++;
				}
				if(Img[(iy+1)*cxDIB+ix] == 255){
					Img[(iy+1)*cxDIB+ix] = 0;
					whiteX[wNum] = ix;
					whiteY[wNum] = iy+1;
					wNum++;
				}
				if(Img[iy*cxDIB+ix-1] == 255){
					Img[iy*cxDIB+ix-1] = 0;
					whiteX[wNum] = ix-1;
					whiteY[wNum] = iy;
					wNum++;
				}
				if(wNum > numTH){
					flag = FALSE;
					break;
				}
			}
			if(flag == FALSE){
				for(k = 0;k < wNum;k++){
					ix = whiteX[k];
					iy = whiteY[k];
					Img[iy*cxDIB+ix] = 255;
				}
			}
		}
	}
}

BYTE get_density(SINGULAR *pSingular,BYTE* OrntImg,
					 int nW,BYTE* Img,int cxDIB,int cyDIB)
{
	int nTotal=0, nLocal=0, nCoreNumber=0;
	int i, j, k, h, sum, cx, cy, nNum=0, sx, sy, ex, ey;

	for ( i=0; i<pSingular->nNumber; i++ ){
		if ( pSingular->nType[i] == -1 ) continue;
		nCoreNumber++;
	}
	if ( nCoreNumber == 0 ){
		cx = cxDIB/2; cy = cxDIB/2;
	}
	else{
		cx = cy = 0;
		for ( i=0; i<pSingular->nNumber; i++ ){
			if ( pSingular->nType[i] == -1 ) continue;
			cx += pSingular->nX[i];
			cy += pSingular->nY[i];
		}
		cx /= nCoreNumber; cy /= nCoreNumber;
	}
	sx = ( cx > nW )? cx-nW+1 : 1; 
	ex = ( cx+nW >= cxDIB )? cxDIB-2 : cx+nW-1; 
	sy = ( cy > nW )? cy-nW+1 : 1; 
	ey = ( cy+nW >= cyDIB )? cyDIB-2 : cy+nW-1;
	for ( i=sx; i<ex; i++ ){
		for ( j=sy; j<ey; j++ ){
			if ( (OrntImg[j*cxDIB+i]&0x80) != 0 ) continue;
			if ( (OrntImg[j*cxDIB+i]&0x7F) == 0x7F ) continue;
			nTotal++; sum=0;
			for ( k=-1; k<=1; k++){
				for ( h=-1; h<=1; h++){
					sum += Img[(j+k)*cxDIB+(i+h)];
				}
			}
			if ( sum<9 || sum>2286 ) continue;
			nLocal++;
		}
	}
	if ( nTotal <= 0 ) return(0);
	return ((BYTE)(255*nLocal/nTotal));
}

int get_frequency_sub(int point_x,int point_y,BYTE* Img,
						BYTE* image_buffer3,int cxDIB,int cyDIB)
{
	int i1,x0,x01,x1,x2,y0,y01,y1,y2;
	int fre = 0,dir,cosV,sinV;
	int color,color0,num,num1,numth = 3;

	dir = image_buffer3[point_y*cxDIB+point_x] + 60;
	if(dir >= 120) dir -= 120;
	cosV = _table_03[dir]; sinV = _table_04[dir];
	color = Img[point_y*cxDIB+point_x];
	x1 = point_x; y1 = point_y;
	i1 = 1;
	while(1){
		x2 = x1+i1*cosV/16384;
		y2 = y1+i1*sinV/16384;
		if(x2 < 0 || x2 >= cxDIB || y2 < 0 || y2 >= cyDIB) return(0);
		if(Img[y2*cxDIB+x2] != color) break;
		i1++;
	}
	color0 = color = Img[y2*cxDIB+x2];
	x0 = x1 = x2;	y0 = y1 = y2;
	num = 0; i1 = 1;
	while(1){
		x2 = x1+i1*cosV/16384;
		y2 = y1+i1*sinV/16384;
		if(x2 < 0 || x2 >= cxDIB || y2 < 0 || y2 >= cyDIB) break;
		if(Img[y2*cxDIB+x2] != color){
			color = Img[y2*cxDIB+x2];
			if(color == color0){
				num++;
				x0 = x2;	y0 = y2;
				if(num == numth) break;
			}
		}
		i1++;
	}
	num1 = num;
	x01 = x1;	y01 = y1;
	x1 = point_x; y1 = point_y;
	num = 0;
	color = color0 = Img[y1*cxDIB+x1];
	i1 = 1;
	while(1){
		x2 = x1-i1*cosV/16384;
		y2 = y1-i1*sinV/16384;
		if(x2 < 0 || x2 >= cxDIB || y2 < 0 || y2 >= cyDIB) break;
		if(Img[y2*cxDIB+x2] != color){
			color = Img[y2*cxDIB+x2];
			if(color == color0){
				num++;
				x01 = x2;	y01 = y2;
				if(num == numth) break;
			}
		}
		i1++;
	}
	num += num1;
	if(num == 0) return(0);

	x1 = x0-x01; x1 = x1*x1;
	y1 = y0-y01; y1 = y1*y1;
	x2 = op_func_02(x1+y1);
	fre = x2*6/num;
	return(fre);
}

BYTE get_frequency(BYTE* Img,BYTE* image_buffer3,SINGULAR *pCore,int cxDIB,int cyDIB)
{
	int i,i1,j,j1,cx,cy,hh,x1,x2,y1,y2;
	int ih = 15;
	int frequency = 0,fre[9];
	BOOL flag;

	cx = cxDIB/2; cy = cyDIB/2; hh = cxDIB/4;
	while(frequency == 0){
		for(i = 0;i < 9;i++){
			x1 = cx; y1 = cy;
			if(i < 3) x1 -= hh;
			if(i > 5) x1 += hh;
			if(i == 0 || i == 3 || i == 6) y1 -= hh;
			if(i == 2 || i == 5 || i == 8) y1 += hh;
			fre[i] = 0;
			flag = TRUE;
			for(j = -1;j < 2;j++){
				x2 = x1+j*ih;
				if(x2 < 0 || x2 >= cxDIB){
					flag = FALSE;
					break;
				}
				for(j1 = -1;j1 < 2;j1++){
					y2 = y1+j1*ih;
					if(y2 < 0 || y2 >= cyDIB){
						flag = FALSE;
						break;
					}
					if((image_buffer3[y2*cxDIB+x2]&0x80) != 0){
						flag = FALSE;
						break;
					}
				}
				if(flag == FALSE) break;
			}
			if(flag == FALSE) continue;
			for(j = 0;j < pCore->nNumber;j++){
				if(x1-ih < pCore->nX[j] && x1+ih > pCore->nX[j] &&
					y1-ih < pCore->nY[j] && y1+ih > pCore->nY[j]){
					flag = FALSE;
					break;
				}
			}
			if(flag == FALSE) continue;
			fre[i] = get_frequency_sub(x1,y1,Img,image_buffer3,cxDIB,cyDIB);
		}
		
		i1 = 0;
		for(i = 0;i < 9;i++){
			if(fre[i] == 0) continue;
			frequency += fre[i];
			i1++;
		}
		if(i1 > 0) frequency /= i1;
		cy += 30;
		if(cy+hh >= cyDIB) break;
	}
	if ( frequency >= 120 ) return (120);
	if ( frequency < 0 ) return (0);
	return ((BYTE)frequency);
}

void remove_bkgrnd(BYTE *Img,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, j;

	for ( i=0; i<cyDIB; i++ ){
		for ( j=0; j<cxDIB; j++ ){
			if ( (OrntImg[i*cxDIB+j] & 0x80) == 0 ) continue;
			Img[i*cxDIB+j] = 0xFF;
		}
	}
}

int get_connectivity(int x,int y,BYTE* pImg,int cxDIB) 
{
	BYTE a[9];
	int i, connec = 0, id0, id1, id2;

	id0 = y*cxDIB+x;
	if ( pImg[id0] != 0 ) return(-1);
	id1 = (y-1)*cxDIB+x;
	id2 = (y+1)*cxDIB+x;
	a[0] = pImg[id1];
	a[1] = pImg[id1+1];
	a[2] = pImg[id0+1];
	a[3] = pImg[id2+1];
	a[4] = pImg[id2];
	a[5] = pImg[id2-1];
	a[6] = pImg[id0-1];
	a[7] = pImg[id1-1];
	a[8] = pImg[id1];
	for ( i=0; i<8; i++ ){
		if ( a[i]==255 && a[i+1]==0 ) connec++;
	}
	return(connec);
}

void correct_thinned_image(BYTE *ppImg, int cxDIB, int cyDIB)
{
	int i, j;

	for(i=1; i<cyDIB-1; i++) {
		for(j=1; j<cxDIB-1; j++) {
			if(ppImg[i*cxDIB+j] != 0x00) continue;
			if(ppImg[i*cxDIB+j-1] == 0x00) {
				if(ppImg[(i-1)*cxDIB+j+1]==0xFF && ppImg[i*cxDIB+j+1]==0xFF && ppImg[(i+1)*cxDIB+j+1]==0xFF
						&& (ppImg[(i-1)*cxDIB+j]==0x00 || ppImg[(i+1)*cxDIB+j]==0x00)) {
					ppImg[i*cxDIB+j] = 0xFF; continue;
				}
			}
			if(ppImg[i*cxDIB+j+1] == 0x00) {
				if(ppImg[(i-1)*cxDIB+j-1]==0xFF && ppImg[i*cxDIB+j-1]==0xFF && ppImg[(i+1)*cxDIB+j-1]==0xFF
						&& (ppImg[(i-1)*cxDIB+j]==0x00 || ppImg[(i+1)*cxDIB+j]==0x00)) {
					ppImg[i*cxDIB+j] = 0xFF; continue;
				}
			}
			if(ppImg[(i+1)*cxDIB+j] == 0x00) {
				if(ppImg[(i-1)*cxDIB+j-1]==0xFF && ppImg[(i-1)*cxDIB+j]==0xFF && ppImg[(i-1)*cxDIB+j+1]==0xFF
						&& (ppImg[i*cxDIB+j-1]==0x00 || ppImg[i*cxDIB+j+1]==0x00)) {
					ppImg[i*cxDIB+j] = 0xFF; continue;
				}
			}
			if(ppImg[(i-1)*cxDIB+j] == 0x00) {
				if(ppImg[(i+1)*cxDIB+j-1]==0xFF && ppImg[(i+1)*cxDIB+j]==0xFF && ppImg[(i+1)*cxDIB+j+1]==0xFF
						&& (ppImg[i*cxDIB+j-1]==0x00 || ppImg[i*cxDIB+j+1]==0x00)) {
					ppImg[i*cxDIB+j] = 0xFF; continue;
				}
			}
			if(ppImg[(i-1)*cxDIB+j]==0x00 && ppImg[i*cxDIB+j-1]==0x00
					&& (ppImg[i*cxDIB+j+1]==0x00 || ppImg[(i+1)*cxDIB+j]==0x00)) {
				ppImg[i*cxDIB+j] = 0xFF; continue;
			}
			if(ppImg[i*cxDIB+j+1]==0x00 && ppImg[(i+1)*cxDIB+j]==0x00
					&& (ppImg[(i-1)*cxDIB+j]==0x00 || ppImg[i*cxDIB+j-1]==0x00)) {
				ppImg[i*cxDIB+j] = 0xFF; continue;
			}
		}
	}
}

void get_thinned_image(BYTE *Img, int cxDIB, int cyDIB)
{
	int nLoopCount = 12;
	int i, j, k, nNum, id;
	BYTE *pMark, *pMark0, *pTmp1, *pTmp2, *pTmp3, *pTmp, nBitSum;

	pTmp1 = (BYTE*)malloc(cxDIB);
	if(pTmp1 == NULL) return;
	pTmp2 = (BYTE*)malloc(cxDIB);
	if(pTmp2 == NULL) { free(pTmp1); return; }
	pTmp3 = (BYTE*)malloc(cxDIB);
	if(pTmp3 == NULL) { free(pTmp1); free(pTmp2); return; }
	pMark = (BYTE*)calloc(cyDIB,sizeof(BYTE));
	if(pMark == NULL) { free(pTmp1); free(pTmp2); free(pTmp3); return; }
	pMark0 = (BYTE*)calloc(cyDIB,sizeof(BYTE));
	if(pMark0 == NULL) { free(pTmp1); free(pTmp2); free(pTmp3); free(pMark);return; }

	for(k=0; k<nLoopCount; k++) {
		nNum = 0;
		for(i=0; i<cyDIB+1; i++) {
			id = (i-1)*cxDIB;
			if(i < cyDIB) memcpy(pTmp1, &Img[i*cxDIB], cxDIB);
			pTmp = pTmp1; pTmp1 = pTmp2; pTmp2 = pTmp3; pTmp3 = pTmp;
			if(i < 2) continue;
			if(pMark[i-1] != 0) continue;
			pMark[i-1] = 1;
			for(j=1; j<cxDIB-1; j++) {
				if(Img[id+j] != 0 ) continue;
				nBitSum = 0;
				if(pTmp1[j-1] == 0) nBitSum |= 0x80;
				if(pTmp2[j-1] == 0) nBitSum |= 0x40;
				if(pTmp3[j-1] == 0) nBitSum |= 0x20;
				if(pTmp3[j+0] == 0) nBitSum |= 0x10;
				if(pTmp3[j+1] == 0) nBitSum |= 0x08;
				if(pTmp2[j+1] == 0) nBitSum |= 0x04;
				if(pTmp1[j+1] == 0) nBitSum |= 0x02;
				if(pTmp1[j+0] == 0) nBitSum |= 0x01;

				if( _table_02[nBitSum] == 0) continue;
				Img[id+j] = 0xFF;
				nNum++;
				pMark[i-1] = 0;
			}
		}
		for(i=0; i<cyDIB+1; i++) {
			id = (i-1)*cxDIB;
			if(i < cyDIB) memcpy(pTmp1, &Img[i*cxDIB], cxDIB);
			pTmp = pTmp1; pTmp1 = pTmp2; pTmp2 = pTmp3; pTmp3 = pTmp;
			if(i < 2) continue;
			if(pMark0[i-1] != 0) continue;
			pMark0[i-1] = 1;
			for(j=1; j<cxDIB-1; j++) {
				if(Img[id+j] != 0 ) continue;
				nBitSum = 0;
				if(pTmp3[j+1] == 0) nBitSum |= 0x80;
				if(pTmp2[j+1] == 0) nBitSum |= 0x40;
				if(pTmp1[j+1] == 0) nBitSum |= 0x20;
				if(pTmp1[j+0] == 0) nBitSum |= 0x10;
				if(pTmp1[j-1] == 0) nBitSum |= 0x08;
				if(pTmp2[j-1] == 0) nBitSum |= 0x04;
				if(pTmp3[j-1] == 0) nBitSum |= 0x02;
				if(pTmp3[j+0] == 0) nBitSum |= 0x01;

				if( _table_02[nBitSum] == 0) continue;
				Img[id+j] = 0xFF;
				nNum++;
				pMark0[i-1] = 0;
			}
		}
		if(nNum == 0) break;
	}
	free(pTmp1); free(pTmp2); free(pTmp3); free(pMark); free(pMark0);
	correct_thinned_image(Img,cxDIB,cyDIB);
}

short get_orient_type1(int x,int y,BYTE* pImg,int cxDIB,int cyDIB)
{
	BYTE loop[9];
	int y_list[20], x_list[20];
	int xval = x, yval = y, id = 0;
	int nIters, num, h, val;
	int id0, id1, id2;

	for ( nIters=0; nIters<20; nIters++ ) {
		if ( yval < 1 || yval >= cyDIB-1 ) break;
		if ( xval < 1 || xval >= cxDIB-1 ) break;
		id0 = (yval-1)*cxDIB+xval;
		id1 = yval*cxDIB+xval;
		id2 = (yval+1)*cxDIB+xval;
		loop[0] = pImg[id0];
		loop[1] = pImg[id0+1];
		loop[2] = pImg[id1+1];
		loop[3] = pImg[id2+1];
		loop[4] = pImg[id2];
		loop[5] = pImg[id2-1];
		loop[6] = pImg[id1-1];
		loop[7] = pImg[id0-1];
		loop[8] = pImg[id0];
		num = 0;
		for ( h=0; h<8; h++ ){
			if ( loop[h] != 255 ) continue;
			if ( loop[h+1] != 0 ) continue;
			num++;
		}
		if ( num != 1 )	break;
		pImg[id1] = 255;
		y_list[id] = yval;
		x_list[id++] = xval;
		if ( loop[0] == 0 ){ yval--; continue; }
		if ( loop[2] == 0 ){ xval++; continue; }
		if ( loop[4] == 0 ){ yval++; continue; }
		if ( loop[6] == 0 ){ xval--; continue; }
		if ( loop[1] == 0 ){ yval--; xval++; continue; }
		if ( loop[3] == 0 ){ yval++; xval++; continue; }
		if ( loop[5] == 0 ){ yval++; xval--; continue; }
		if ( loop[7] == 0 ){ yval--; xval--; continue; }
	}
	for ( h=0; h<id; h++ ){
		pImg[y_list[h]*cxDIB+x_list[h]] = 0;
	}
	if ( nIters >= 10 ){
		val = op_func_01(xval, yval, x, y);
	}
	else{
		val = -20;
	}
	return (short)val;
}

short get_orient_type2(int x,int y,BYTE* pImg,int cxDIB,int cyDIB)
{
	BYTE loop[9], l_val, r_val, u_val, d_val;
	int id = 0, h, nIters, xval, yval, xflag, yflag;
	int y_list[60], x_list[60], val_list[3];
	int k, num, min_start, min_end, maxdiff, mindiff;
	int diff, diff10, diff20, diff21;
	int id0, id1, id2, yid0, yid1, yid2;

	for ( nIters=0; nIters<3; nIters++ ){
		yid0 = y*cxDIB+x;
		yid1 = (y-1)*cxDIB+x;
		yid2 = (y+1)*cxDIB+x;
		yval = y; xval = x; yflag = 0; xflag = 0;
		l_val = pImg[yid0-1];
		r_val = pImg[yid0+1];
		u_val = pImg[yid1];
		d_val = pImg[yid2];
		for ( k=0; k<20; k++ ){
			if ( yval < 1 || yval >= cyDIB-1 ) break;
			if ( xval < 1 || xval >= cxDIB-1 ) break;
			id0 = (yval-1)*cxDIB+xval;
			id1 = yval*cxDIB+xval;
			id2 = (yval+1)*cxDIB+xval;
			loop[0] = pImg[id0];
			loop[1] = pImg[id0+1];
			loop[2] = pImg[id1+1];
			loop[3] = pImg[id2+1];
			loop[4] = pImg[id2];
			loop[5] = pImg[id2-1];
			loop[6] = pImg[id1-1];
			loop[7] = pImg[id0-1];
			loop[8] = pImg[id0];
			num = 0;
			for ( h=0; h<8; h++ ){
				if ( loop[h] != 255 ) continue;
				if ( loop[h+1] != 0 ) continue;
				num++;
			}
			if ( num != 1 && k > 1 ) break;
			pImg[id1] = 255;
			y_list[id] = yval;
			x_list[id++] = xval;
			if ( loop[0] == 0 ) yval--;
			else if ( loop[2] == 0 ) xval++;
			else if ( loop[4] == 0 ) yval++;
			else if ( loop[6] == 0 ) xval--;
			else if ( loop[1] == 0 ){ yval--; xval++; }
			else if ( loop[3] == 0 ){ yval++; xval++; }
			else if ( loop[5] == 0 ){ yval++; xval--; }
			else if ( loop[7] == 0 ){ yval--; xval--; }

			if ( k == 0 ){
				if ( xval == x ) xflag = 1;
				if ( yval == y ) yflag = 1;
				if ( xflag != 0 ){
					pImg[yid0-1] = 255;
					pImg[yid0+1] = 255;
				}
				if ( yflag == 0 ) continue;
				pImg[yid1] = 255;
				pImg[yid2] = 255;
			}
			else{
				if ( k != 1 ) continue;
				if ( xflag != 0 ){
					pImg[yid0-1] = l_val;
					pImg[yid0+1] = r_val;
				}
				if ( yflag == 0 ) continue;
				pImg[yid1] = u_val;
				pImg[yid2] = d_val;
			}
		}
		if ( k < 10 ){
			val_list[nIters] = -10; break;
		}
		val_list[nIters] = op_func_01(xval, yval, x, y);
	}
	for ( h=0; h<id; h++ ){
		pImg[y_list[h]*cxDIB+x_list[h]] = 0;
	}
	for ( h=0; h<3; h++ ){
		if ( val_list[h] < 0 ) return(-666);
	}
	min_start = val_list[2], min_end = val_list[1];
	diff21 = abs( val_list[2]-val_list[1] );
	if ( diff21 > 120 ) diff21 = 240 - diff21;
	maxdiff = diff21; mindiff = diff21;

	diff10 = abs( val_list[1]-val_list[0] );
	if ( diff10 > 120 ) diff10 = 240 - diff10;
	if(diff10 > maxdiff) maxdiff = diff10;
	if ( diff10 < mindiff ){
		mindiff = diff10;	
		min_start = val_list[1]; 
		min_end = val_list[0];
	}

	diff20 = abs( val_list[2]-val_list[0] );
	if ( diff20 > 120 ) diff20 = 240 - diff20;
	if(diff20 > maxdiff) maxdiff = diff20;
	if ( diff20 < mindiff ){
		mindiff = diff20;	
		min_start = val_list[2]; 
		min_end = val_list[0];
	}

	if ( mindiff == 0 ) return(-20);
	if ( maxdiff < 60 ) return(-20);

	diff = abs( min_start - min_end );
	if ( diff > 120 ){
		diff = (240 - diff) / 2;
		if ( min_start > 120 ) diff += min_start;
		else diff += min_end;
		if ( diff >= 240 ) diff -= 240;
		return (short)diff;
	}
	else{
		return (short)((min_start+min_end)/2);
	}
}

BYTE get_point_curve(int nW, int x, int y, BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, count = 0, sum = 0, id;
	BYTE p0 = OrntImg[y*cxDIB+x] & 0x7F;
	BYTE p1, p;

	if ( p0 == 0x7F ) return(0);
	for ( i=(y<nW)?0:y-nW; i<=y+nW; i++ ){
		if ( i >= cyDIB ) break;
		id = i*cxDIB;
		for ( j=(x<nW)?0:x-nW; j<=x+nW; j++ ){
			if ( j >= cxDIB ) break;
			p1 = OrntImg[id+j] & 0x7F;
			if ( p1 == 0x7F ) continue;
			p = abs(p0-p1);
			if ( p > 60 ) p = abs(120 - p);
			sum += p;
			count++;
		}
	}
	return (BYTE)(255*sum/(60*count));
}

void get_mp_points(LPREALPVECT pVect,BYTE* Img,BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, val, id;
	short dir;
	BYTE kind;
	
	pVect->nNumber = 0;
	for ( i=8; i<cyDIB-8; i++ ){
		id = i*cxDIB;
		for ( j=8; j<cxDIB-8; j++ ){
			if ( OrntImg[id+j] >= 0x80) continue;
			val = get_connectivity(j,i,Img,cxDIB);
			if ( val != 1 && val != 3 ) continue;
			if ( val == 1 ){
				dir = get_orient_type1(j,i,Img,cxDIB,cyDIB);
				kind = END_MINUTIA; 
			}
			if ( val == 3 ){
				dir = get_orient_type2(j,i,Img,cxDIB,cyDIB);
				kind = BIF_MINUTIA; 
			}
			if ( pVect->nNumber >= MAXMINUTIATAGNUM ) return;
			pVect->item[pVect->nNumber].x = (short)j;
			pVect->item[pVect->nNumber].y = (short)i;
			pVect->item[pVect->nNumber].dir = dir;
			pVect->item[pVect->nNumber++].kind = kind;
		}
	}
}

BOOL check_false_mp(int x1,int y1,int dir1,int x2,int y2,int dir2)
{
	int dir, diff, tmp, dx, dy;

	dir = op_func_01(x2,y2, x1,y1);
	diff = abs(dir - dir1);
	if ( diff > 120 ) diff = 240 - diff;
	dx = abs(x1 - x2); dy = abs(y1 - y2);
	if ( dx > 13 || dy > 13 ) return FALSE;
	if ( dx<7 && dy<7 ){
		if ( diff < 97 ) return FALSE;
	}
	else{
		if ( diff < 100 ) return FALSE;
	}
	tmp = dir + 120;
	if ( tmp >= 240 ) tmp -= 240;
	diff = abs(dir2 - tmp);
	if ( diff > 120 ) diff = 240 - diff;
	if ( dx<7 && dy<7 ){
		if ( diff < 97 ) return FALSE;
	}
	else{
		if ( diff < 100 ) return FALSE;
	}
	return TRUE;
}

void filter_mp_points(LPREALPVECT pVect,SINGULAR* pSingular,BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, m, n, count, dx, dy;
	int x, y, startx, endx, starty, endy;
	int dTH = 16;

	for ( i=0; i<pVect->nNumber; i++ ){
		if ( pVect->item[i].kind != END_MINUTIA ) continue;
		if ( pVect->item[i].dir < 0 ) continue;
		for ( j=0; j<pVect->nNumber; j++ ){
			if ( i == j ) continue;
			if ( pVect->item[j].dir < 0 ) continue;
			if ( pVect->item[j].kind != END_MINUTIA ) continue;
			if ( check_false_mp(pVect->item[i].x,pVect->item[i].y,pVect->item[i].dir,
				pVect->item[j].x,pVect->item[j].y,pVect->item[j].dir) == TRUE ){
				pVect->item[i].dir = -1;
				pVect->item[j].dir = -1; break;
			}
		}
	}
	for ( i=0; i<pVect->nNumber; i++ ){
		count = 0;
		for ( j=0; j<pVect->nNumber; j++ ){
			if ( i == j ) continue;
			dy = pVect->item[i].y - pVect->item[j].y;
			dx = pVect->item[i].x - pVect->item[j].x;
			dx = dy*dy + dx*dx;
			if ( dx < dTH*dTH ) count++;
		}
		if ( count > 5 ) pVect->item[i].dir = -15;
	}
	for ( i=0; i<pVect->nNumber; i++ ){
		for ( j=0; j<pVect->nNumber; j++ ){
			if ( i == j ) continue;
			dy = abs(pVect->item[i].y - pVect->item[j].y);
			dx = abs(pVect->item[i].x - pVect->item[j].x);
			if ( dx*dx+dy*dy <= 4*4 ){
				pVect->item[i].dir = pVect->item[j].dir = -1;
				break;
			}
		}
	}
	for ( i=0; i<pVect->nNumber; i++ ){
		if ( pVect->item[i].dir < 0 ) continue;
		x = pVect->item[i].x; y = pVect->item[i].y;
		startx = x - 10; endx = x + 10; 
		if ( startx < 0 || endx >= cxDIB ){
			pVect->item[i].dir = -25; continue;
		}
		for ( n=startx; n<=endx; n++ ){
			starty = y - 10; endy = y + 10;
			if ( starty < 0 || endy >= cyDIB ){
				pVect->item[i].dir = -25; continue;
			}
			for ( m=starty; m<endy; m++ ){
				if ( OrntImg[m*cxDIB+n] < 0x7F ) continue;
				pVect->item[i].dir = -25; n = 10000; break;
			}
		}
	}
	for ( i=0; i<pVect->nNumber; i++ ){
		if ( pVect->item[i].dir < 0 ) continue;
		x = pVect->item[i].x; y = pVect->item[i].y;
		for ( j=0; j<pSingular->nNumber; j++ ){
			if ( pSingular->nType[j] != 1 ) continue;
			dx = x - pSingular->nX[j]; dy = y - pSingular->nY[j];
			if ( dx*dx+dy*dy < dTH*dTH ) break;
		}
		if ( j < pSingular->nNumber ){
			pVect->item[i].dir = -1;
		}
	}

	for ( i=0,j=0; i<pVect->nNumber; i++ ){
		if (pVect->item[i].dir < 0)	continue;
		pVect->item[j++] = pVect->item[i];
	}
	pVect->nNumber = j;
}

void filter_mp_points2(LPREALPVECT pVect)
{
	int i, j, dx, dy, len, diff;

	for ( i=0; i<pVect->nNumber; i++ ){
		if ( pVect->item[i].dir == 0xFF ) continue;
		if ( pVect->item[i].score >= 30 ) continue;
		if ( pVect->item[i].kind != END_MINUTIA ) continue;
		for ( j=0; j<pVect->nNumber; j++ ){
			if ( pVect->item[j].dir == 0xFF ) continue;
			if ( pVect->item[j].score >= 30 ) continue;
			if ( pVect->item[j].kind != END_MINUTIA ) continue;
			if ( i == j ) continue;
			dy = pVect->item[i].y - pVect->item[j].y;
			dx = pVect->item[i].x - pVect->item[j].x;
			len = dy*dy + dx*dx;
			if ( len >= 8*8 ) continue;
			diff = abs(pVect->item[i].dir - pVect->item[j].dir);
			if(diff > 120) diff = 240 - diff;
			if ( abs(120-diff) >= 20 ) continue;
			pVect->item[i].dir = 0xFF; pVect->item[j].dir = 0xFF;
			break;
		}
	}
	for ( i=0,j=0; i<pVect->nNumber; i++ ){
		if (pVect->item[i].dir == 0xFF )	continue;
		pVect->item[j++] = pVect->item[i];
	}
	pVect->nNumber = j;
}

void get_point_value(LPREALPVECT pVect,BYTE* Img,int cxDIB,int cyDIB)
{
	int	i, m, n, x, y, sum, num, rr = 10;

	for(i=0;i<pVect->nNumber;i++){
		x = pVect->item[i].x;
		y = pVect->item[i].y;
		sum = num = 0;
		for(m = max(0,y-rr);(m <= y+rr && m < cyDIB);m++){
			for(n = max(0,x-rr);(n <= x+rr && n < cxDIB);n++){
				if((n-x)*(n-x)+(m-y)*(m-y) > rr*rr) continue;
				sum += Img[m*cxDIB+n];
				num++;
			}
		}
		if(num == 0) sum = 0;
		else sum /= num;
		pVect->item[i].score = sum;
	}
}

void arrange_mp(LPREALPVECT pRealVect,LPMPVECTEX pVect,BYTE* OrntImg,int cxDIB,int cyDIB,COREVECTEX* pCore)
{
	int i, j, nMinIdx, nMinValue;
	TAG_POINT tmp;
	int cx = 0, cy = 0, val, score = 0;

	for ( i=0; i<pCore->nNumber; i++ ){
		cx += pCore->item[i].x; cy += pCore->item[i].y;
	}
	for ( i=0; i<pRealVect->nNumber; i++ ){
		score += pRealVect->item[i].score;
	}
	if ( pRealVect->nNumber > 0 ) score /= pRealVect->nNumber; 
	if ( pCore->nNumber > 0 && score > 45 ){
		cx /= pCore->nNumber; cy /= pCore->nNumber;
		for ( i=0; i<pRealVect->nNumber-1; i++ ){
			nMinIdx = i; nMinValue = (pRealVect->item[i].x-cx)*(pRealVect->item[i].x-cx) +
									(pRealVect->item[i].y-cy)*(pRealVect->item[i].y-cy);
			for ( j=i+1; j<pRealVect->nNumber; j++ ){
				val = (pRealVect->item[j].x-cx)*(pRealVect->item[j].x-cx) +
						(pRealVect->item[j].y-cy)*(pRealVect->item[j].y-cy);
				if ( val >= nMinValue ) continue;
				nMinIdx = j; nMinValue = val;
			}
			if(nMinIdx == i) continue;
			tmp = pRealVect->item[i];
			pRealVect->item[i] = pRealVect->item[nMinIdx];
			pRealVect->item[nMinIdx] = tmp;
		}
	}
	else {
		for ( i=0; i<pRealVect->nNumber-1; i++ ){
			nMinIdx = i; nMinValue = pRealVect->item[i].score;
			for ( j=i+1; j<pRealVect->nNumber; j++ ){
				if ( pRealVect->item[j].score <= nMinValue ) continue;
				nMinIdx = j; nMinValue = pRealVect->item[j].score;
			}
			if(nMinIdx == i) continue;
			tmp = pRealVect->item[i];
			pRealVect->item[i] = pRealVect->item[nMinIdx];
			pRealVect->item[nMinIdx] = tmp;
		}
	}
	for ( i=0,j=0; i<pRealVect->nNumber; i++ ){
		pVect->item[j].x = pRealVect->item[i].x;
		pVect->item[j].y = pRealVect->item[i].y;
		pVect->item[j].dir = (BYTE)pRealVect->item[i].dir;
		pVect->item[j].kind = pRealVect->item[i].kind;
		pVect->item[j].curv = get_point_curve(10,pVect->item[j].x,pVect->item[j].y,OrntImg,cxDIB,cyDIB);
		pVect->item[j++].score = pRealVect->item[i].score;
		if ( j >= MAX_MINUTIA_NUMBER ) break;
	}
	pVect->nNumber = j;
}

void comp_typeline(TYPELINE *pLine,BYTE *pData)
{
	pData[0] = (BYTE)(((pLine->points_x[0] << 7) & 0xFF00) >> 8);
	pData[1] = (BYTE)((pLine->points_x[0] & 0x0001) << 7);
	pData[1] |= (BYTE)(((pLine->points_x[1] << 6) & 0xFF00) >> 8);
	pData[2] = (BYTE)((pLine->points_x[1] & 0x0003) << 6);
	pData[2] |= (BYTE)(((pLine->points_x[2] << 5) & 0xFF00) >> 8);
	pData[3] = (BYTE)((pLine->points_x[2] & 0x0007) << 5);
	pData[3] |= (BYTE)(((pLine->points_x[3] << 4) & 0xFF00) >> 8);
	pData[4] = (BYTE)((pLine->points_x[3] & 0x000F) << 4);
	pData[4] |= (BYTE)(((pLine->points_y[0] << 3) & 0xFF00) >> 8);
	pData[5] = (BYTE)((pLine->points_y[0] & 0x001F) << 3);
	pData[5] |= (BYTE)(((pLine->points_y[1] << 2) & 0xFF00) >> 8);
	pData[6] = (BYTE)((pLine->points_y[1] & 0x003F) << 2);
	pData[6] |= (BYTE)(((pLine->points_y[2] << 1) & 0xFF00) >> 8);
	pData[7] = (BYTE)((pLine->points_y[2] & 0x007F) << 1);
	pData[8] |= (BYTE)(((pLine->points_y[3] << 0) & 0xFF00) >> 8);
	pData[8] = (BYTE)((pLine->points_y[3] & 0x00FF) << 0);
}

void comp_block(BLOCKVECT *pBlock,BYTE *pData)
{
	int i, k, n = pBlock->nCol*pBlock->nRow;
	BYTE val;

	for ( i=0,k=0; i<n; i+=2,k++ ){
		if ( pBlock->Data[i] == 0xFF ) val = 15;
		else val = pBlock->Data[i] / 8;
		pData[k] = (val & 0x0F) << 4;
		if ( i+1 == n ) pData[k] |= (15 & 0x0F); 
		else{
			if ( pBlock->Data[i+1] == 0xFF ) val = 15;
			else val = pBlock->Data[i+1] / 8;
			pData[k] |= (val & 0x0F);
		}
	}
}

void comp_core(COREVECTEX *pCore,BYTE *pData)
{
	int i, k;

	pData[0] = pCore->nNumber;
	for ( i=0,k=1; i<pCore->nNumber; i++,k+=4 ){
		pData[k+0] = (BYTE)(((pCore->item[i].x << 7) & 0xFF00) >> 8);
		pData[k+1] = (BYTE)((pCore->item[i].x & 0x0001) << 7);
		pData[k+1] |= (BYTE)((pCore->item[i].y & 0x0100) >> 8 );
		pData[k+2] = (BYTE)(pCore->item[i].y & 0x00FF);
		pData[k+3] = pCore->item[i].dir;
	}
}

int comp_mp(MPVECTEX *pVect,BYTE *pData)
{
	int i, k;
	BYTE nByte;

	for ( i=0,k=0; i<pVect->nNumber; i++,k+=MP_SIZE ){
		pData[k+0] = (BYTE)(((pVect->item[i].x << 7) & 0xFF00) >> 8);
		pData[k+1] = (BYTE)((pVect->item[i].x & 0x0001) << 7);
		pData[k+1] |= (BYTE)((pVect->item[i].y & 0x0100) >> 8 );
		pData[k+2] = (BYTE)(pVect->item[i].y & 0x00FF);
		pData[k+3] = pVect->item[i].dir;
		pData[k+4] = pVect->item[i].curv;
		nByte = pVect->item[i].score & 0x7F;
		if ( pVect->item[i].kind == 1 ) nByte = nByte | 0x80;
		pData[k+5] = nByte;
	}
	return (k);
}

void get_byte_template(LPFPVECTEX pFPEx,LPFPFEATURE pFeature)
{
	int n = 2, n1, n2;


	pFeature->Data[n++] = pFPEx->nHeader.nVersion;
	pFeature->Data[n++] = pFPEx->nHeader.nHeaderLength;
	pFeature->Data[n++] = (BYTE)((pFPEx->nHeader.nImageX & 0xFF00) >> 8 );
	pFeature->Data[n++] = (BYTE)(pFPEx->nHeader.nImageX & 0xFF);
	pFeature->Data[n++] = (BYTE)((pFPEx->nHeader.nImageY & 0xFF00) >> 8 );
	pFeature->Data[n++] = (BYTE)(pFPEx->nHeader.nImageY & 0xFF);

	pFeature->Data[n++] = 'A';
	pFeature->Data[n++] = 3;
	pFeature->Data[n++] = pFPEx->nType;
	pFeature->Data[n++] = pFPEx->nRidgeDensity;
	pFeature->Data[n++] = pFPEx->nFrequency;

	pFeature->Data[n++] = 'B';
	pFeature->Data[n++] = 9;
	comp_typeline(&(pFPEx->MainLine),&pFeature->Data[n]);
	n += 9;

	pFeature->Data[n++] = 'C';
	if ( pFPEx->Block.nRow*pFPEx->Block.nCol % 2 == 0 )
		n1 = 3 + (pFPEx->Block.nRow*pFPEx->Block.nCol)/2;
	else
		n1 = 4 + (pFPEx->Block.nRow*pFPEx->Block.nCol)/2;
	pFeature->Data[n++] = (BYTE)((n1 & 0xFF00) >> 8 );
	pFeature->Data[n] = (BYTE)(n1 & 0xFF);
	pFeature->Data[n+1] = pFPEx->Block.nRow;
	pFeature->Data[n+2] = pFPEx->Block.nCol;
	comp_block(&(pFPEx->Block),&pFeature->Data[n+3]);
	n += n1;

	pFeature->Data[n++] = 'D';
	pFeature->Data[n++] = 4*pFPEx->Core.nNumber+1;
	comp_core(&(pFPEx->Core),&pFeature->Data[n]);
	n += 4*pFPEx->Core.nNumber+1;

	pFeature->Data[n++] = 'E';
	n1 = comp_mp(&(pFPEx->Mp),&pFeature->Data[n+2]);
	n2 = n1+2;
	n2 /= 2;
	pFeature->Data[n++] = (BYTE)(n2);
	pFeature->Data[n++] = pFPEx->Mp.nNumber;
	n += n1+1;

	pFeature->Data[0] = (BYTE)((n & 0xFF00) >> 8 );
	pFeature->Data[1] = (BYTE)(n & 0xFF);
}

int ext_main(BYTE* pRawImage,int width,int height,LPFPFEATURE pFeatureVect)
{
	BYTE *image_buffer1, *image_buffer2, *image_buffer3, *image_buffer4;
	BYTE *buffer1, *buffer2;
	REALPVECT tempVect;
	SINGULAR SingularData, tmpSP;
	FPVECTEX FPEx;
	int nCol = width/BLOCK_SIZE, nRow = height/BLOCK_SIZE;


	memset( pFeatureVect, 0, MAX_FEATUREVECT_LEN );
	memset( &FPEx, 0, sizeof(FPVECTEX) );

	image_buffer1 = (BYTE*)malloc(sizeof(BYTE)*width*height);
	if ( image_buffer1 == NULL ) return ERR_CAN_NOT_ALLOC_MEMORY;
	memcpy(image_buffer1, pRawImage, sizeof(BYTE)*width*height);

	image_buffer2 = (BYTE*)malloc(sizeof(BYTE)*width*height);
	if ( image_buffer2 == NULL ){ free(image_buffer1); return ERR_CAN_NOT_ALLOC_MEMORY; }
	image_buffer3 = (BYTE*)calloc(width*height,sizeof(BYTE));
	if ( image_buffer3 == NULL ){ free(image_buffer2); free(image_buffer1); return ERR_CAN_NOT_ALLOC_MEMORY; }
	image_buffer4 = (BYTE*)calloc(width*height,sizeof(BYTE));
	if ( image_buffer4 == NULL ){ free(image_buffer3); free(image_buffer2); free(image_buffer1); return ERR_CAN_NOT_ALLOC_MEMORY; }
	buffer1 = (BYTE*)calloc(nCol*nRow,sizeof(BYTE));
	if ( buffer1 == NULL ){ free(image_buffer4); free(image_buffer3); free(image_buffer2); free(image_buffer1); return ERR_CAN_NOT_ALLOC_MEMORY; }
	buffer2 = (BYTE*)calloc(nCol*nRow,sizeof(BYTE));
	if ( buffer2 == NULL ){ free(buffer1); free(image_buffer4); free(image_buffer3); free(image_buffer2); free(image_buffer1); return ERR_CAN_NOT_ALLOC_MEMORY; }
	
	get_smoothed_image(image_buffer1,width,height);
	memcpy(image_buffer2, image_buffer1, sizeof(BYTE)*width*height);
	get_smoothed_image2(image_buffer2,width,height,4);
	get_sharpend_image(image_buffer2,image_buffer1,width,height,64);
	get_smoothed_image(image_buffer2,width,height);
	get_smoothed_image(image_buffer2,width,height);

	get_block_image(buffer1,buffer2,nCol,nRow,image_buffer2,width,height,17);
	filter_block_image(buffer1,buffer2,nCol,nRow,207);
	get_sp_cand(&tmpSP,buffer2,nCol,nRow);
	free(buffer2); free(buffer1);

	get_orient_image(image_buffer2,image_buffer3,width,height,12,19,image_buffer4);
	get_singular_points(&SingularData,&tmpSP,image_buffer3,width,height);
	get_bkgrnd(image_buffer2,image_buffer3,width,height,&SingularData,16,196);

	image_proc_01(image_buffer3,image_buffer1,width,height);
	get_smoothed_image(image_buffer2,width,height);
	get_binary_image2(image_buffer3,image_buffer2,image_buffer1,width,height,3,7);
	image_proc_05(image_buffer2,image_buffer3,width,height,16);
	get_smoothed_image2(image_buffer2,width,height,16);
	get_binary_image1(image_buffer3,image_buffer1,image_buffer2,width,height,130,115,3,5,7);

	image_proc_02(image_buffer1,image_buffer3,width,height,12,2);
	re_get_orient_image(image_buffer1,image_buffer2,image_buffer3,width,height,&SingularData,12,30,20,0,4);

	image_proc_03(image_buffer2,image_buffer3,image_buffer1,width,height);
	get_binary_image1(image_buffer3,image_buffer1,image_buffer2,width,height,130,115,3,5,7);

	get_block_data(&FPEx,nCol,nRow,image_buffer3,width,height);
	get_core_points(&SingularData,image_buffer3,width,height);
	FPEx.nType = (BYTE)get_type_line(&FPEx,&SingularData,image_buffer3,image_buffer2,width,height);
	copy_core(&SingularData,&FPEx);

	image_proc_04(image_buffer1,width,height);
	remove_hole(image_buffer3,image_buffer1,width,height);
	FPEx.nRidgeDensity = get_density(&SingularData,image_buffer3,64,image_buffer1,width,height);
	FPEx.nFrequency = get_frequency(image_buffer1,image_buffer3,&SingularData,width,height);

	remove_bkgrnd(image_buffer1,image_buffer3,width,height);
	get_thinned_image(image_buffer1,width,height);

	get_mp_points(&tempVect,image_buffer1,image_buffer3,width,height);
	filter_mp_points(&tempVect,&SingularData,image_buffer3,width,height);
	get_point_value(&tempVect,image_buffer4,width,height);
	filter_mp_points2( &tempVect );
	
	arrange_mp(&tempVect,&FPEx.Mp,image_buffer3,width,height,&FPEx.Core);

	free(image_buffer4); free(image_buffer3); free(image_buffer2); free(image_buffer1);
	if ( FPEx.Mp.nNumber < MIN_MINUTIA_NUMBER && FPEx.Core.nNumber == 0 ) return ERR_VECT_FAILED;

	FPEx.nHeader.nVersion = VERSION_NUMBER;
	FPEx.nHeader.nHeaderLength = MAX_HEADER_SIZE;
	FPEx.nHeader.nImageX = width;
	FPEx.nHeader.nImageY = height;
	get_byte_template(&FPEx,pFeatureVect);

ViewFPEx[0] = FPEx;
	return ERR_OK;
}


BOOL mch_sub_func_03(LPFPVECTEX pFPEx)
{
	if (pFPEx->nHeader.nVersion != VERSION_NUMBER) return FALSE;
	if (pFPEx->Core.nNumber > MAX_CORE_NUMBER) return FALSE;
	if (pFPEx->Mp.nNumber > MAX_MINUTIA_NUMBER) return FALSE;
	if (pFPEx->Mp.nNumber < MIN_MINUTIA_NUMBER) return FALSE;
	return TRUE;
}

void decomp_typeline(BYTE *pData,TYPELINE *pLine)
{
	pLine->points_x[0] = pData[0] << 1;
	pLine->points_x[0] |= (pData[1] & 0x80) >> 7;
	pLine->points_x[1] = (pData[1] & 0x7F) << 2;
	pLine->points_x[1] |= (pData[2] & 0xC0) >> 6;
	pLine->points_x[2] = (pData[2] & 0x3F) << 3;
	pLine->points_x[2] |= (pData[3] & 0xE0) >> 5;
	pLine->points_x[3] = (pData[3] & 0x1F) << 4;
	pLine->points_x[3] |= (pData[4] & 0xF0) >> 4;
	pLine->points_y[0] = (pData[4] & 0x0F) << 5;
	pLine->points_y[0] |= (pData[5] & 0xF8) >> 3;
	pLine->points_y[1] = (pData[5] & 0x07) << 6;
	pLine->points_y[1] |= (pData[6] & 0xFC) >> 2;
	pLine->points_y[2] = (pData[6] & 0x03) << 7;
	pLine->points_y[2] |= (pData[7] & 0xFE) >> 1;
	pLine->points_y[3] = (pData[7] & 0x01) << 8;
	pLine->points_y[3] |= pData[8];
}

void decomp_block(BYTE *pData,BLOCKVECT *pBlock)
{
	int i, k, n = pBlock->nCol*pBlock->nRow;
	BYTE val;

	if ( n % 2 == 0 ) n /= 2;
	else n = n/2 + 1;
	
	for ( i=0,k=0; i<n; i++ ){
		val = (pData[i] & 0xF0) >> 4;
		if ( val == 15 ) pBlock->Data[k++] = 0xFF;
		else pBlock->Data[k++] = val*8;
		val = (pData[i] & 0x0F);
		if ( val == 15 ) pBlock->Data[k++] = 0xFF;
		else pBlock->Data[k++] = val*8;
	}
}

void decomp_core(BYTE *pData,COREVECTEX *pCore)
{
	int i, k;

	pCore->nNumber = pData[0];
	for ( i=0,k=1; i<pCore->nNumber; i++,k+=4 ){
		pCore->item[i].x = pData[k+0] << 1;
		pCore->item[i].x |= (pData[k+1] & 0x80) >> 7;
		pCore->item[i].y = (pData[k+1] & 0x01) << 8;
		pCore->item[i].y |= pData[k+2];
		pCore->item[i].dir = pData[k+3];
		if ( pCore->item[i].dir == 255 ) pCore->item[i].kind = DELTA;
		else pCore->item[i].kind = CORE;
	}
}

void decomp_mp(BYTE *pData,MPVECTEX *pVect)
{
	int i, n = 1;
	BYTE nByte;

	pVect->nNumber = pData[0];
	for ( i=0; i<pVect->nNumber; i++ ){
		pVect->item[i].x = pData[n++] << 1;
		pVect->item[i].x |= (pData[n] & 0x80) >> 7;
		pVect->item[i].y = (pData[n++] & 0x01) << 8;
		pVect->item[i].y |= pData[n++];
		pVect->item[i].dir = pData[n++];
		pVect->item[i].curv = pData[n++];
		nByte = pData[n++];
		pVect->item[i].score = nByte & 0x7F;
		if ( (nByte&0x80) != 0 ) pVect->item[i].kind = 1;
		else pVect->item[i].kind = 0;
	}
}

void mch_sub_func_02(LPFPFEATURE pFeature,LPFPVECTEX pFPEx)
{
    volatile int n, nNum, nSize;
	BYTE nMark;

	memset(pFPEx,0,sizeof(FPVECTEX));

	nSize = ((pFeature->Data[0] << 8) | pFeature->Data[1]);	

	pFPEx->nHeader.nVersion = pFeature->Data[2];
	pFPEx->nHeader.nHeaderLength = pFeature->Data[3];
	pFPEx->nHeader.nImageX = ((pFeature->Data[4] << 8) | pFeature->Data[5]);
	pFPEx->nHeader.nImageY = ((pFeature->Data[6] << 8) | pFeature->Data[7]);
	n = pFPEx->nHeader.nHeaderLength;
	while (1){
		if ( n >= nSize ) break;
		nMark = pFeature->Data[n++];
		nNum = pFeature->Data[n++];
		if ( nMark != 'C' && nNum == 0 ) continue;
		
		if ( n+nNum >= nSize ) break;

		if(nMark == 'A'){
			pFPEx->nType = pFeature->Data[n];
			pFPEx->nRidgeDensity = pFeature->Data[n+1];
			pFPEx->nFrequency = pFeature->Data[n+2];
		}

		if(nMark == 'B'){
			decomp_typeline(&pFeature->Data[n],&pFPEx->MainLine);
		}

		if(nMark == 'C'){
			nNum = ((nNum << 8) | pFeature->Data[n]);
			pFPEx->Block.nRow = pFeature->Data[n+1];
			pFPEx->Block.nCol = pFeature->Data[n+2];
			decomp_block(&pFeature->Data[n+3],&pFPEx->Block);
		}

		if(nMark == 'D'){
			decomp_core(&pFeature->Data[n],&pFPEx->Core);
		}

		if(nMark == 'E'){
			nNum *= 2;
			if(nNum > MP_SIZE+1){
				decomp_mp(&pFeature->Data[n],&(pFPEx->Mp));
			}
		}
		n += nNum;
	}
}

void get_tag_item(LPMPVECTEX pVect,MATCH_TAG* pBar)
{
	int dir = op_func_01( pVect->item[pBar->nID2].x,
						  pVect->item[pBar->nID2].y,
						  pVect->item[pBar->nID1].x,
						  pVect->item[pBar->nID1].y);
	pBar->nSlope = dir;
	if ( pBar->nSlope >= 120 ) pBar->nSlope -= 120;

	pBar->nDiff1 = dir - pVect->item[pBar->nID1].dir;
	if ( pBar->nDiff1 < 0 ) pBar->nDiff1 += 240;
	if ( pBar->nDiff1 >= 240 ) pBar->nDiff1 -= 240;

	dir += 120;
	if ( dir >= 240 ) dir -= 240;

	pBar->nDiff2 = dir - pVect->item[pBar->nID2].dir;
	if ( pBar->nDiff2 < 0 ) pBar->nDiff2 += 240;
	if ( pBar->nDiff2 >= 240 ) pBar->nDiff2 -= 240;
}

int get_search_tag_sub(int nMinLen,int nMaxLen,int nMaxNum,
					 LPMPVECTEX pS,int *nMaxSearchBarLen,
					 MATCH_TAG* pBarItem,int* SDiffField,
					 int (*SArrangBarPtr)[20]
					)
{
	BARTEMP pTemp[MAX_SEARCH_BAR_NUM];
	int i, j, dx, dy, len, sum, id, nNum = 0, pSum[300];
	int nBarNum, pTmpPtr[MAX_SEARCH_BAR_NUM];
	int	nMax = nMaxLen*nMaxLen, nMin = nMinLen*nMinLen;
	int score = 0, numTh;
	
	if ( nMaxLen > 300 ) return (0);
	*nMaxSearchBarLen = 0;

	for ( i=0; i<pS->nNumber; i++ ){
		score += pS->item[i].score;
	}
	if ( pS->nNumber > 0 ) score /= pS->nNumber;
	else score = 0;
	if ( score > 45 ) numTh = MAX_SEARCH_BAR_NUM;
	else numTh = 350;

	for ( i=0; i<pS->nNumber; i++ ){
		for ( j=i+1; j<pS->nNumber; j++ ){
			if ( j == i ) continue;
			dx = pS->item[i].x - pS->item[j].x;
			dy = pS->item[i].y - pS->item[j].y;
			len = dx*dx+dy*dy;
			if ( len < nMin || len >= nMax ) continue;
			pTemp[nNum].nLen = op_func_02(len);
			pTemp[nNum].nID1 = i;
			pTemp[nNum++].nID2 = j;
			if ( nNum >= numTh ){ i = 1000; break; }
		}
	}
	if ( 2*nNum >= nMaxNum ){
		memset( pSum, 0, sizeof(int)*300 );
		for ( i=0; i<nNum; i++ ) pSum[pTemp[i].nLen]++;
		for ( i=1; i<300; i++ ) pSum[i] += pSum[i-1];
		for ( i=1; i<300; i++ ) pSum[i]--;
		for ( i=0; i<nNum; i++ ){
			len = pTemp[i].nLen;
			sum = pSum[len];
			pSum[len] = sum - 1;
			pTmpPtr[sum] = i;
		}
		nNum = nMaxNum / 2;
	}
	else{
		for ( i=0; i<nNum; i++ ) pTmpPtr[i] = i;
	}
	memset( SDiffField, 0, sizeof(int)*240 );
	nBarNum = 0;
	for ( i=0; i<nNum; i++ ){
		len = pTemp[pTmpPtr[i]].nLen;
		pBarItem[nBarNum].nLen = len;
		if ( *nMaxSearchBarLen < len ) *nMaxSearchBarLen = len + 1;
		pBarItem[nBarNum].nID1 = pTemp[pTmpPtr[i]].nID1;
		pBarItem[nBarNum].nID2 = pTemp[pTmpPtr[i]].nID2;

		get_tag_item( pS, &pBarItem[nBarNum] );

		id = pBarItem[nBarNum].nDiff1;
		SArrangBarPtr[id][SDiffField[id]] = nBarNum;
		if ( ++SDiffField[id] == 20 ){ SDiffField[id]--; }
		nBarNum++;
		pBarItem[nBarNum].nLen = pTemp[pTmpPtr[i]].nLen;
		pBarItem[nBarNum].nID1 = pTemp[pTmpPtr[i]].nID2;
		pBarItem[nBarNum].nID2 = pTemp[pTmpPtr[i]].nID1;
		pBarItem[nBarNum].nDiff1 = pBarItem[nBarNum-1].nDiff2;
		pBarItem[nBarNum].nDiff2 = pBarItem[nBarNum-1].nDiff1;
		pBarItem[nBarNum].nSlope = pBarItem[nBarNum-1].nSlope;

		id = pBarItem[nBarNum].nDiff1;
		SArrangBarPtr[id][SDiffField[id]] = nBarNum;
		if ( ++SDiffField[id] == 20 ){ SDiffField[id]--; }
		nBarNum++;
	}
	return(nBarNum);
}

void get_search_tag(LPFPVECTEX pSearch,BARVECT* pSearchBar,
				  int *nMaxSearchBarLen, int* SDiffField,
				  int (*SArrangBarPtr)[20],int nMinLen,int nMaxLen
				 )
{
	if ( nMinLen > nMaxLen ) return;
	pSearchBar->nNumber = get_search_tag_sub(nMinLen,nMaxLen,300, 
										&pSearch->Mp,
										nMaxSearchBarLen,
										&pSearchBar->item[0],
										SDiffField,
										SArrangBarPtr
										  ); 
}

int get_file_tag_sub(int nMinLen,int nMaxLen,int nMaxNum,
				   LPMPVECTEX pFT,MATCH_TAG* pBarItem,
				   int* FDiffField,int (*FArrangBarPtr)[20]
				  )
{
	int i, j, len, dx, dy, id, nBarNum = 0;
	int nMin,nMax;
	memset( FDiffField, 0, sizeof(int)*240 );
	nMin = nMinLen*nMinLen;
	nMax = nMaxLen*nMaxLen;
	for ( i=0; i<pFT->nNumber; i++ ){
		for ( j=i+1; j<pFT->nNumber; j++ ){
			dx = pFT->item[i].x - pFT->item[j].x;
			dy = pFT->item[i].y - pFT->item[j].y;
			len = dx*dx+dy*dy;
			if ( len <= nMin || len >= nMax ) continue;
			pBarItem[nBarNum].nLen = op_func_02(len);
			pBarItem[nBarNum].nID1 = i;
			pBarItem[nBarNum].nID2 = j;
			get_tag_item( pFT, &pBarItem[nBarNum] );
			id = pBarItem[nBarNum].nDiff1;
			FArrangBarPtr[id][FDiffField[id]] = nBarNum;
			if ( ++FDiffField[id] == 20 ){ FDiffField[id]--; }
			if ( ++nBarNum >= nMaxNum ) return(nBarNum);
		}
	}
	return(nBarNum);
}

void get_file_tag(int nMaxLen,LPFPVECTEX pFile,
				BARVECT* pFileBar,int* FDiffField,
				int (*FArrangBarPtr)[20],
				int* nFileCX,int* nFileCY,int nMinLen
			   )
{
	int i, nXMax = 0, nYMax = 0, nXMin = 10000, nYMin = 10000;

	pFileBar->nNumber = get_file_tag_sub(
								nMinLen-LENGTH,nMaxLen+LENGTH,600,
								&pFile->Mp,
								&pFileBar->item[0],
								FDiffField,
								FArrangBarPtr 
									  );
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		if ( nXMin > pFile->Mp.item[i].x ) nXMin = pFile->Mp.item[i].x;
		if ( nXMax < pFile->Mp.item[i].x ) nXMax = pFile->Mp.item[i].x;
		if ( nYMin > pFile->Mp.item[i].y ) nYMin = pFile->Mp.item[i].y;
		if ( nYMax < pFile->Mp.item[i].y ) nYMax = pFile->Mp.item[i].y;
	}
	*nFileCX = (nXMin+nXMax) / 2; *nFileCY = (nYMin+nYMax) / 2;
}

int mch_sub_func_01(LPCOREVECTEX pSingular,COREITEMEX *pCore,COREITEMEX *pDelta,int *nNumDelta)
{
	int i, nCoreNum = 0, nDeltaNum = 0;

	for ( i=0; i<pSingular->nNumber; i++ ){
		if ( pSingular->item[i].kind == CORE ){
			pCore[nCoreNum++] = pSingular->item[i];
			if ( nCoreNum >= 2 ) break;
		}
		else {
			if ( pDelta == NULL ) continue;
			pDelta[nDeltaNum++] = pSingular->item[i];
			if ( nDeltaNum >= 2 ) break;
		}
	}
	if ( nNumDelta != NULL ) *nNumDelta = nDeltaNum;
	return nCoreNum;
}

void transform_core(LPCOREVECTEX pVect,int cx,int cy,int nAngle,int nDiffX,int nDiffY)
{
	int i, x, y, angle, rot, nCos, nSin;

	rot = 240 - nAngle;
	if ( rot >= 240 ) rot -= 240;
	nCos = _table_03[rot]; nSin=_table_04[rot];

	for ( i=0; i<pVect->nNumber; i++ ){
		x = (pVect->item[i].x-cx)*nCos + (pVect->item[i].y-cy)*nSin;
		if ( x > 0 ) x += 8192;
		else x -= 8192;
		x = x >> 14;
		y = (pVect->item[i].y-cy)*nCos - (pVect->item[i].x-cx)*nSin;
		if ( y > 0 ) y += 8192;
		else y -= 8192;
		y = y >> 14;

		pVect->item[i].x = x + cx + nDiffX;
		pVect->item[i].y = y + cy + nDiffY;

		if ( pVect->item[i].kind != CORE ) continue;
		angle = pVect->item[i].dir + nAngle;
		if ( angle >= 240 ) angle -= 240;
		if ( angle < 0 ) angle += 240;
		pVect->item[i].dir = angle;
	}
}

void transform_mp(LPMPVECTEX pVect,int cx,int cy,int nAngle,int nDiffX,int nDiffY)
{
	int i, x, y, angle, rot, nCos, nSin;

	rot = 240 - nAngle;
	if ( rot >= 240 ) rot -= 240;
	nCos = _table_03[rot]; nSin=_table_04[rot];

	for ( i=0; i<pVect->nNumber; i++ ){
		x = (pVect->item[i].x-cx)*nCos + (pVect->item[i].y-cy)*nSin;
		if ( x > 0 ) x += 8192;
		else x -= 8192;
		x = x >> 14;
		y = (pVect->item[i].y-cy)*nCos - (pVect->item[i].x-cx)*nSin;
		if ( y > 0 ) y += 8192;
		else y -= 8192;
		y = y >> 14;

		pVect->item[i].x = x + cx + nDiffX;
		pVect->item[i].y = y + cy + nDiffY;

		angle = pVect->item[i].dir + nAngle;
		if ( angle >= 240 ) angle -= 240;
		if ( angle < 0 ) angle += 240;
		pVect->item[i].dir = angle;
	}
}

void transform_points(LPFPVECTEX pVect,int xc,int yc,
					  int nAngle,int nDiffX,int nDiffY)
{
	int i, newx, newy, angle, dx, dy;
	int	nCos = _table_03[nAngle], nSin = _table_04[nAngle];

	for(i=0; i<pVect->Mp.nNumber; i++) {
		dx = (nCos*(pVect->Mp.item[i].x-xc)) >> 14;
		dy = (nSin*(pVect->Mp.item[i].y-yc)) >> 14;
		newx = xc + dx - dy;
		dx = (nSin*(pVect->Mp.item[i].x-xc)) >> 14;
		dy = (nCos*(pVect->Mp.item[i].y-yc)) >> 14;
		newy = yc + dx + dy;
		pVect->Mp.item[i].x = newx + nDiffX;
		pVect->Mp.item[i].y = newy + nDiffY;
		angle = pVect->Mp.item[i].dir + nAngle;
		if(angle >= 240) angle -= 240;
		if(angle < 0) angle += 240;
		pVect->Mp.item[i].dir = angle;
	}
	for(i=0; i<pVect->Core.nNumber; i++) {
		dx = (nCos*(pVect->Core.item[i].x-xc)) >> 14;
		dy = (nSin*(pVect->Core.item[i].y-yc)) >> 14;
		newx = xc + dx - dy;
		dx = (nSin*(pVect->Core.item[i].x-xc)) >> 14;
		dy = (nCos*(pVect->Core.item[i].y-yc)) >> 14;
		newy = yc + dx + dy;
		pVect->Core.item[i].x = newx + nDiffX;
		pVect->Core.item[i].y = newy + nDiffY;

		if ( pVect->Core.item[i].kind != CORE ) continue;
		angle = pVect->Core.item[i].dir + nAngle;
		if(angle >= 240) angle -= 240;
		if(angle < 0) angle += 240;
		pVect->Core.item[i].dir = angle;
	}
}

int get_matched_mp_num(int nLenTh,int nAngTh,LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	int i, i1, i2, i3, j, num = 0;
	char temp[MAX_MINUTIA_NUMBER];
	BOOL flag;

	memset( temp, 0, sizeof(char)*pVect2->nNumber );
	for ( i=0; i<pVect1->nNumber; i++ ){
		flag = FALSE;
		for ( j=0; j<pVect2->nNumber; j++ ){
			i1 = pVect1->item[i].x - pVect2->item[j].x;
			i2 = pVect1->item[i].y - pVect2->item[j].y;
			i3 = i1*i1 + i2*i2;
			if ( i3 > nLenTh*nLenTh ) continue;
			i1 = abs(pVect1->item[i].dir - pVect2->item[j].dir);
			if ( i1 > 120 ) i1 = 240 - i1;
			if ( i1 > nAngTh ) continue;
			temp[j] = 1; flag = TRUE;
		}
		if( flag == TRUE ) num++;
	}
	i1 = 0;
	for ( i=0; i<pVect2->nNumber; i++ ){
		if ( temp[i] == 1 ) i1++;
	}
	if ( num > i1 ) num = i1;
	return (num);
}

int get_min_points_number(LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	POLYGON pol1, pol2;
	int i, num1, num2;

	if ( get_polygon_points(pVect1,&pol1) == FALSE ) return (0);
	if ( get_polygon_points(pVect2,&pol2) == FALSE ) return (0);
	num1 = num2 = 0;
	for ( i=0; i<pVect1->nNumber; i++ ){
		if ( check_in_polygon(pVect1->item[i].x,pVect1->item[i].y,&pol2,16) == TRUE ) num1++;
	}
	for ( i=0; i<pVect2->nNumber; i++ ){
		if ( check_in_polygon(pVect2->item[i].x,pVect2->item[i].y,&pol1,16) == TRUE ) num2++;
	}
	if ( num1 > num2 ) num1 = num2;
	return (num1);
}

int rotate_points(int cx,int cy,int* pAngle,BARVECT* pBar,LPFPVECTEX pVect)
{
	int i, j, sum, nMax, nMaxId, sumN1, sumN2, temp[300];
	int rotAngle, rot, nCos, nSin, x, y, angle;

	for ( i=0; i<240; i++ ){ 
		sum = 0;
		for ( j=i-4; j<=i+4; j++ ){
			rot = j;
			if(rot < 0) rot += 240;
			else if(rot >= 240) rot -= 240;
			sum += pAngle[rot];
		}
		temp[i] = sum;
	}
	memcpy( pAngle, temp, sizeof(int)*240 );
	nMax = 0; nMaxId = 0;
	for ( i=0; i<240; i++ ){
		if ( pAngle[i] > nMax ) {
			nMax = pAngle[i]; nMaxId = i;
		}
	}
	for ( i=0; i<10; i++ ) temp[i] = pAngle[230+i];
	for ( i=0; i<240; i++ ) temp[i+10] = pAngle[i];
	for ( i=0; i<10; i++ ) temp[i+250] = pAngle[i];

	nMax /= 2;
	sumN1 = sumN2 = 0;
	for ( i=nMaxId; i<nMaxId+20; i++) {
		if ( temp[i] <= nMax ) continue;
		if ( temp[i] <= 20 ) continue;
		sumN1 += (temp[i] - nMax) * i;
		sumN2 += (temp[i] - nMax);
	}
	if(sumN2 == 0) rotAngle = 0;
	else{
		rotAngle = (100*sumN1/sumN2 + 50) / 100;
	}
	rotAngle -= 10;
	if(rotAngle < 0) rotAngle += 240;
	if(rotAngle >= 240) rotAngle -= 240;

	rot = 240 - rotAngle;
	if ( rot >= 240 ) rot -= 240;
	nCos = _table_03[rot]; nSin = _table_04[rot];
	for(i=0; i<pVect->Mp.nNumber; i++) {
		x = (pVect->Mp.item[i].x-cx)*nCos + (pVect->Mp.item[i].y-cy)*nSin;
		if ( x > 0 ) x += 8192;
		else x -= 8192;
		x = x >> 14;
		y = (pVect->Mp.item[i].y-cy)*nCos - (pVect->Mp.item[i].x-cx)*nSin;
		if ( y > 0 ) y += 8192;
		else y -= 8192;
		y = y >> 14;

		pVect->Mp.item[i].x = x + cx;
		pVect->Mp.item[i].y = y + cy;

		angle = pVect->Mp.item[i].dir + rotAngle;
		if ( angle >= 240 ) angle -= 240;
		if ( angle < 0 ) angle += 240;
		pVect->Mp.item[i].dir = angle;
	}
	for ( i=0; i<pBar->nNumber; i++ ) {
		angle = pBar->item[i].nSlope + rotAngle;
		if ( angle >= 240 ) angle -= 240;
		if ( angle < 0 ) angle += 240;
		if ( angle >= 120 ) angle -= 120;
		pBar->item[i].nSlope = angle;
	}

	for ( i=0; i<pVect->Core.nNumber; i++ ){
		x = (pVect->Core.item[i].x-cx)*nCos + (pVect->Core.item[i].y-cy)*nSin;
		if ( x > 0 ) x += 8192;
		else x -= 8192;
		x = x >> 14;
		y = (pVect->Core.item[i].y-cy)*nCos - (pVect->Core.item[i].x-cx)*nSin;
		if ( y > 0 ) y += 8192;
		else y -= 8192;
		y = y >> 14;

		pVect->Core.item[i].x = x + cx;
		pVect->Core.item[i].y = y + cy;

		angle = pVect->Core.item[i].dir + rotAngle;
		if ( angle >= 240 ) angle -= 240;
		if ( angle < 0 ) angle += 240;
		pVect->Core.item[i].dir = angle;
	}
	return(rotAngle);
}

void get_shift_param(int nTH,int nScore,MATCH_TAG* pFBar,MATCH_TAG* pSBar, 
				   int* XField,int* YField,LPMPVECTEX pFile,LPMPVECTEX pSearch)
{
	int nSid1, nSid2, nFid1, nFid2, dx1, dx2, dy1, dy2, dx, dy;
	if ( nScore == 0 ) return;

	nSid1 = pSBar->nID1; nSid2 = pSBar->nID2;
	nFid1 = pFBar->nID1; nFid2 = pFBar->nID2;
	dx1 = pSearch->item[nSid1].x - pFile->item[nFid1].x;
	dx2 = pSearch->item[nSid2].x - pFile->item[nFid2].x;
	dy1 = pSearch->item[nSid1].y - pFile->item[nFid1].y;
	dy2 = pSearch->item[nSid2].y - pFile->item[nFid2].y;

	dx = abs( dx1 - dx2 );
	if ( dx >= nTH ) return;
	dy = abs( dy1 - dy2 );
	if ( dy >= nTH ) return;
	if ( abs(dx2) >= MAX_SHIFT_X_SIZE ) return;
	if ( abs(dy2) >= MAX_SHIFT_Y_SIZE ) return;
	if ( abs(dx1) >= MAX_SHIFT_X_SIZE ) return;
	if ( abs(dy1) >= MAX_SHIFT_Y_SIZE ) return;

	dx = dx1 + dx2;
	if ( dx < 0 ) dx++;
	dx /= 2;
	XField[dx+MAX_SHIFT_X_SIZE/2] += nScore;
	dy = dy1 + dy2;
	if ( dy < 0 ) dy++;
	dy /= 2;
	YField[dy+MAX_SHIFT_Y_SIZE/2] += nScore;
}

void shift_points(int* nXoffset,int* nYoffset,LPFPVECTEX pVect,int* XField,int* YField)
{
	int pTmp[MAX_SHIFT_SIZE];
	int i, j, starti, endi, sum, nMax, nMaxId, nSum1, nSum2;

	memset( pTmp, 0, sizeof(int)*MAX_SHIFT_SIZE );
	for( i=5; i<MAX_SHIFT_X_SIZE-5; i++ ){
		sum = 0;
		for ( j=i-5; j<i+5; j++ ) sum += XField[j];
		pTmp[i] = sum;
	}
	memcpy( XField, pTmp, sizeof(int)*MAX_SHIFT_X_SIZE );
	nMax = nMaxId = 0;
	for( i=0; i<MAX_SHIFT_X_SIZE; i++ ){
		if ( XField[i] <= nMax ) continue;
		nMax = XField[i]; nMaxId = i;
	}
	nSum1 = nSum2 = 0;
	starti = ( nMaxId-10 < 0 ) ? 0 : nMaxId-10;
	endi = ( nMaxId+10 > MAX_SHIFT_X_SIZE-1 ) ? MAX_SHIFT_X_SIZE-1 : nMaxId+10;
	nMax = (nMax*2) / 3;
	for( i=starti; i<endi; i++ ){
		if ( XField[i] <= nMax ) continue;
		nSum1 += XField[i] * i;
		nSum2 += XField[i];
	}
	*nXoffset = ( nSum2 == 0 ) ? 0 : (nSum1/nSum2) - MAX_SHIFT_X_SIZE/2;
	
	for( i=5; i<MAX_SHIFT_Y_SIZE-5; i++ ){
		sum = 0;
		for ( j=i-5; j<i+5; j++ ) sum += YField[j];
		pTmp[i] = sum;
	}
	memcpy( YField, pTmp, sizeof(int)*MAX_SHIFT_Y_SIZE );
	nMax = nMaxId = 0;
	for( i=0; i<MAX_SHIFT_Y_SIZE; i++ ){
		if ( YField[i] <= nMax ) continue;
		nMax = YField[i]; nMaxId = i;
	}
	nSum1 = nSum2 = 0;
	starti = ( nMaxId-10 < 0 ) ? 0 : nMaxId-10;
	endi = ( nMaxId+10 > MAX_SHIFT_Y_SIZE-1 ) ? MAX_SHIFT_Y_SIZE-1 : nMaxId+10;
	nMax = (nMax*2) / 3;
	for( i=starti; i<endi; i++ ){
		if ( YField[i] <= nMax ) continue;
		nSum1 += YField[i] * i;
		nSum2 += YField[i];
	}
	*nYoffset = ( nSum2 == 0 ) ? 0 : (nSum1/nSum2) - MAX_SHIFT_Y_SIZE/2;

	for( i=0; i<pVect->Mp.nNumber; i++ ){
		pVect->Mp.item[i].x += *nXoffset; pVect->Mp.item[i].y += *nYoffset;
	}
	for( i=0; i<pVect->Core.nNumber; i++ ){
		pVect->Core.item[i].x += *nXoffset; pVect->Core.item[i].y += *nYoffset;
	}
}

BOOL re_arrange_point(PAIRBAR* PairList,int* LastList,int nPairNum,int *nLastNum,
					   LPFPVECTEX pFile,LPFPVECTEX pSearch,
					   BARVECT* pFBar,BARVECT* pSBar)
{
	int i, j, k, x, y, bx, by, nNum, nRow, nCol, nDiv;
	int nList[MAX_BLOCK_NUMBER];
	int nNewId[MAX_BAR_NUM], nNewBarId[MAX_BAR_NUM];
	int nEmptyPair[MAX_BAR_NUM], nEmptyPairNum;

	nNum = 0;
	for ( i=0; i<pFile->Block.nCol*pFile->Block.nRow; i++ ){
		if ( pFile->Block.Data[i] == 0xFF ) continue;
		if ( pSearch->Block.Data[i] == 0xFF ) continue;
		nList[nNum++] = i;
	}
	if ( nNum == 0 ) return FALSE;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ) nNewId[i] = -1;
	nDiv = pFile->Block.nCol;
	for ( i=0,k=0; i<pFile->Mp.nNumber; i++ ){
		x = pFile->Mp.item[i].x; y = pFile->Mp.item[i].y;
		for ( j=0; j<nNum; j++ ){
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( abs(x-bx)>8 || abs(y-by)>8 ) continue;
			nNewId[i] = k;
			pFile->Mp.item[k++] = pFile->Mp.item[i];
			break;
		}
	}
	if ( k == 0 ) return FALSE;
	pFile->Mp.nNumber = k;

	for ( i=0; i<MAX_BAR_NUM; i++ ) nNewBarId[i] = -1;
	for ( i=0,j=0; i<pFBar->nNumber; i++ ){
		if ( nNewId[pFBar->item[i].nID1] == -1 ) continue;
		if ( nNewId[pFBar->item[i].nID2] == -1 ) continue;
		nNewBarId[i] = j;
		pFBar->item[j] = pFBar->item[i];
		pFBar->item[j].nID1 = nNewId[pFBar->item[i].nID1];
		pFBar->item[j++].nID2 = nNewId[pFBar->item[i].nID2];
	}
	if ( j == 0 ) return FALSE;
	pFBar->nNumber = j;

	nEmptyPairNum = 0;
	for ( i=0; i<nPairNum; i++ ){
		j = nNewBarId[PairList[i].fid];
		if ( j == -1 ) nEmptyPair[nEmptyPairNum++] = i;
		else PairList[i].fid = j;
		if(nEmptyPairNum >= MAX_BAR_NUM) break;
	}
	
	for ( i=0,k=0; i<*nLastNum; i++ ){
		for ( j=0; j<nEmptyPairNum; j++ ){
			if ( LastList[i] == nEmptyPair[j] ) break;
		}
		if ( j >= nEmptyPairNum ) nNewId[k++] = LastList[i];
	}
	for ( i=0; i<k; i++ ) LastList[i] = nNewId[i];
	*nLastNum = k;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ) nNewId[i] = -1;
	nDiv = pSearch->Block.nCol;
	for ( i=0,k=0; i<pSearch->Mp.nNumber; i++ ){
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		for ( j=0; j<nNum; j++ ){
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( abs(x-bx)>8 || abs(y-by)>8 ) continue;
			nNewId[i] = k;
			pSearch->Mp.item[k++] = pSearch->Mp.item[i];
			break;
		}
	}
	if ( k == 0 ) return FALSE;
	pSearch->Mp.nNumber = k;
	
	for ( i=0; i<MAX_BAR_NUM; i++ ) nNewBarId[i] = -1;
	for ( i=0,j=0; i<pSBar->nNumber; i++ ){
		if ( nNewId[pSBar->item[i].nID1] == -1 ) continue;
		if ( nNewId[pSBar->item[i].nID2] == -1 ) continue;
		nNewBarId[i] = j;
		pSBar->item[j] = pSBar->item[i];
		pSBar->item[j].nID1 = nNewId[pSBar->item[i].nID1];
		pSBar->item[j++].nID2 = nNewId[pSBar->item[i].nID2];
	}
	if ( j == 0 ) return FALSE;
	pSBar->nNumber = j;

	nEmptyPairNum = 0;
	for ( i=0; i<nPairNum; i++ ){
		j = nNewBarId[PairList[i].sid];
		if ( j == -1 ) nEmptyPair[nEmptyPairNum++] = i;
		else PairList[i].sid = j;
		if(nEmptyPairNum >= MAX_BAR_NUM) break;
	}
	
	for ( i=0,k=0; i<*nLastNum; i++ ){
		for ( j=0; j<nEmptyPairNum; j++ ){
			if ( LastList[i] == nEmptyPair[j] ) break;
		}
		if ( j >= nEmptyPairNum ) nNewId[k++] = LastList[i];
	}
	for ( i=0; i<k; i++ ) LastList[i] = nNewId[i];
	*nLastNum = k;
	return TRUE;
}

BOOL check_limit(int nTH,MATCH_TAG* pFBar,MATCH_TAG* pSBar,LPFPVECTEX pF,LPFPVECTEX pS,int cx,int cy)
{
	int nSid = pSBar->nID1, nFid = pFBar->nID1;
	int dx, dy, th, len, nRidge, nVal;

	nRidge = max(pF->nRidgeDensity,pS->nRidgeDensity);
	if ( nRidge > 200 ) nVal = 9;
	else  nVal = 10;

	dx = pF->Mp.item[nFid].x - cx; dy = pF->Mp.item[nFid].y - cy;
	len = op_func_02(dx*dx+dy*dy);
	th = ( len <= 150 ) ? nVal+len/50 : nVal+3;

	if ( abs(pS->Mp.item[nSid].x - pF->Mp.item[nFid].x) >= th ) return FALSE;
	if ( abs(pS->Mp.item[nSid].y - pF->Mp.item[nFid].y) >= th ) return FALSE;

	nSid = pSBar->nID2; nFid = pFBar->nID2;
	dx = pF->Mp.item[nFid].x - cx; dy = pF->Mp.item[nFid].y - cy;
	len = op_func_02(dx*dx+dy*dy);
	th = ( len <= 150 ) ? nVal+len/50 : nVal+3;

	if ( abs(pS->Mp.item[nSid].x - pF->Mp.item[nFid].x) >= th ) return FALSE;
	if ( abs(pS->Mp.item[nSid].y - pF->Mp.item[nFid].y) >= th ) return FALSE;
	return(TRUE);
}

BOOL check_exist(int x,int y,int dir,int nID,int nLenTh,int nAngTh,
			 LPMPVECTEX pVect,PAIRVECT *pPair,
			 BOOL nPairFlag,BOOL nSimple,BOOL nForS
			)
{
	int i, j, id, len, diff, th = nLenTh*nLenTh;

	for ( i=0; i<pVect->nNumber; i++ ){
		if ( i == nID ) continue;
		if ( nPairFlag == TRUE ){
			for ( j=0; j<pPair->nNumber; j++ ){
				id = (nForS==FALSE)?pPair->nFileID[j]:pPair->nSearchID[j];
				if ( i == id ) break;
			}
			if ( j < pPair->nNumber ) continue;
		}
		len = (x-pVect->item[i].x)*(x-pVect->item[i].x) + (y-pVect->item[i].y)*(y-pVect->item[i].y);
		diff = abs(dir - pVect->item[i].dir);
		if ( diff >= 120 ) diff = 240 - diff;
		if ( len < th ){
			if ( nSimple == FALSE ){
				if ( diff < nAngTh ) return TRUE;
				continue;
			}
			return TRUE;
		}
	}
	return FALSE;
}

int dec_func_01(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int nRetScore = score, nFileScore = 0, nSearchScore = 0;
	int i, maxlen, diff1, len, diflen, diff, difdir, dx, dy;
	BOOL nScoreFlag = TRUE;
	int nRot = pPair->nRot;
	int nFCoreNum, nSCoreNum, nFDeltaNum = 0, nSDeltaNum = 0;
	COREITEMEX FileCore[2], SearchCore[2];
	COREITEMEX FileDelta[2], SearchDelta[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,FileDelta,&nFDeltaNum);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,SearchDelta,&nSDeltaNum);

	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		nFileScore += pFile->Mp.item[i].score;
	}
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		nSearchScore += pSearch->Mp.item[i].score;
	}
	if ( pFile->Mp.nNumber==0 || pSearch->Mp.nNumber==0 ) return 0;
	nFileScore /= pFile->Mp.nNumber;
	nSearchScore /= pSearch->Mp.nNumber;
	if ( nFileScore < 30 || nSearchScore < 30 ) nScoreFlag = FALSE;

	if ( nFCoreNum == 1 && nSCoreNum == 1 ){
		dx = FileCore[0].x - SearchCore[0].x;
		dy = FileCore[0].y - SearchCore[0].y;
		len = op_func_02(dx*dx + dy*dy);
		diff = abs(FileCore[0].dir - SearchCore[0].dir);
		if ( diff >= 120 ) diff = 240 - diff;
		if ( len < LENGTH  && diff >= ANGLE*4 ){
			if ( nScoreFlag == TRUE ) nRetScore = nRetScore/2;
			else nRetScore = nRetScore*2/3;
		}
		else if ( len >= LENGTH*3 && diff < ANGLE*2){
			if ( nScoreFlag == TRUE ) nRetScore = nRetScore/2;
			else {
				nRetScore = nRetScore*68/100;
			}
		}
		else if ( 2*len >= LENGTH*5 && diff < ANGLE) nRetScore = nRetScore*2/3;
		else if ( 4*len >= LENGTH*5 && diff < ANGLE) nRetScore = nRetScore*84/100; 
		else if ( len > LENGTH && diff > 3*ANGLE) nRetScore = nRetScore*80/100; 
		if ( nRot >= 90 && nRot < 150 ){
			if ( len >= LENGTH*3 && diff >= ANGLE*3 && nScoreFlag ){
				nRetScore = nRetScore/2;
			}
		}
		diff = abs(120-diff);
		if ( len >= LENGTH*7 && diff >= ANGLE*3 ){
			nRetScore = nRetScore/2;
		}
		else{
			if ( len >= 120 && nScoreFlag){
				if ( (pFile->nType == 8 && (pSearch->nType == 4 || pSearch->nType == 5)) || 
					 (pSearch->nType == 8 && (pFile->nType == 4 || pFile->nType == 5)) ) nRetScore = nRetScore*3/5;
				else nRetScore = nRetScore/3;
			}
			else if ( len > 100 ) nRetScore = nRetScore*2/3;
		}
	}
	if ( nFCoreNum == 2 && nSCoreNum == 2 ){
		dx = FileCore[0].x - FileCore[1].x;
		dy = FileCore[0].y - FileCore[1].y;
		len = op_func_02(dx*dx + dy*dy);
		dx = SearchCore[0].x - SearchCore[1].x;
		dy = SearchCore[0].y - SearchCore[1].y;
		maxlen = op_func_02(dx*dx + dy*dy);
		len = abs(len-maxlen);
		diff = op_func_01(FileCore[0].x,FileCore[0].y,
						  FileCore[1].x,FileCore[1].y);
		if ( diff >= 120 ) diff -= 120;
		diff1 = op_func_01(SearchCore[0].x,SearchCore[0].y,
						   SearchCore[1].x,SearchCore[1].y);
		if ( diff1 >= 120 ) diff1 -= 120;
		diff = abs(diff-diff1);
		if ( diff >= 60 ) diff = 120 - diff;
		diff = op_func_02(len*len + diff*diff);
		if ( diff >= 50 ){
			nRetScore = nRetScore/2;
		}
	}
	if ( nFCoreNum != nSCoreNum && nFCoreNum > 0 && nSCoreNum > 0 ){
		if ( nFCoreNum == 1 ){
			dx = FileCore[0].x - SearchCore[0].x;
			dy = FileCore[0].y - SearchCore[0].y;
			len = dx*dx + dy*dy;
			dx = FileCore[0].x - SearchCore[1].x;
			dy = FileCore[0].y - SearchCore[1].y;
			maxlen = dx*dx + dy*dy;
			if ( maxlen > len ) maxlen = len;
		}
		else{
			dx = FileCore[0].x - SearchCore[0].x;
			dy = FileCore[0].y - SearchCore[0].y;
			len = dx*dx + dy*dy;
			dx = FileCore[1].x - SearchCore[0].x;
			dy = FileCore[1].y - SearchCore[0].y;
			maxlen = dx*dx + dy*dy;
			if ( maxlen > len ) maxlen = len;
		}
		if ( maxlen >= LENGTH*LENGTH*100 ) nRetScore = 0;
		if ( maxlen >= LENGTH*LENGTH*64 ) nRetScore = nRetScore/2;
	}

	if ( nFDeltaNum == 0 || nSDeltaNum == 0 ) return nRetScore;
	if ( nFDeltaNum == 2 && nSDeltaNum == 2 ){
		dx = FileDelta[0].x-FileDelta[1].x;
		dy = FileDelta[0].y-FileDelta[1].y;
		len = op_func_02(dx*dx+dy*dy);
		dx = SearchDelta[0].x-SearchDelta[1].x;
		dy = SearchDelta[0].y-SearchDelta[1].y;
		maxlen = op_func_02(dx*dx+dy*dy);
		diflen = abs(len-maxlen);
		diff = op_func_01(FileDelta[0].x,FileDelta[0].y,
						  FileDelta[1].x,FileDelta[1].y);
		if ( diff >= 120 ) diff -= 120;
		diff1 = op_func_01(SearchDelta[0].x,SearchDelta[0].y,
						   SearchDelta[1].x,SearchDelta[1].y);
		if ( diff1 >= 120 ) diff1 -= 120;
		difdir = abs(diff-diff1);
		if ( difdir >= 60 ) difdir = 120 - difdir;
		diff = op_func_02(diflen*diflen + difdir*difdir);
		if ( diff >= 200 || (nScoreFlag == TRUE && diff >= 100) ) nRetScore = 0;
		if ( diff >= 60 || difdir > 4*ANGLE ){
			nRetScore = nRetScore - nRetScore*diff/200;
		}
		if ( diflen < LENGTH && difdir < ANGLE ) nRetScore = nRetScore*6/5;
	}

	if ( nFCoreNum == 1 && nFDeltaNum == 1 && nSCoreNum == 1 && nSDeltaNum == 1 ){
		if ( pFile->nType == pSearch->nType && ( pFile->nType == 4 || pFile->nType == 5 || pFile->nType == 7 )){
			dx = FileCore[0].x-FileDelta[0].x;
			dy = FileCore[0].y-FileDelta[0].y;
			len = op_func_02(dx*dx+dy*dy);
			dx = SearchCore[0].x-SearchDelta[0].x;
			dy = SearchCore[0].y-SearchDelta[0].y;
			maxlen = op_func_02(dx*dx+dy*dy);
			diflen = abs(len-maxlen);
			diff = op_func_01(FileCore[0].x,FileCore[0].y,
							  FileDelta[0].x,FileDelta[0].y);
			if ( diff >= 120 ) diff -= 120;
			diff1 = op_func_01(SearchCore[0].x,SearchCore[0].y,
							   SearchDelta[0].x,SearchDelta[0].y);
			if ( diff1 >= 120 ) diff1 -= 120;
			difdir = abs(diff-diff1);
			if ( difdir >= 60 ) difdir = 120 - difdir;
			if ( nScoreFlag == TRUE && ( diflen >= 10*LENGTH )) nRetScore = 0;
			if ( diflen > 4*LENGTH || difdir > 3*ANGLE ){
				if ( nScoreFlag == TRUE ) nRetScore = nRetScore/3;
				else nRetScore = nRetScore*2/3;
			}
			if ( diflen <= 10 && difdir <= 7 && score > 50 ) nRetScore = nRetScore*6/5;
		}
	}
	return nRetScore;
}

int dec_func_02(int score,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int nRetScore = score, nDiffNum = 0;
	int i, j, k, fid, sid, dx, dy, nFNum, nSNum;
	int nLenTh = 35*35;

	for ( i=0; i<pPair->nNumber; i++ ){
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		if ( pFile->item[fid].score >= 45 && pSearch->item[sid].score >= 30 ){
			nFNum = 0;
			for ( j=0; j<pFile->nNumber; j++ ){
				if ( pFile->item[j].score < 45 ) continue;
				for ( k=0; k<pPair->nNumber; k++ ){
					if ( j == pPair->nFileID[k] ) break;
				}
				if ( k < pPair->nNumber ) continue;
				dx = pFile->item[fid].x - pFile->item[j].x;
				dy = pFile->item[fid].y - pFile->item[j].y;
				if ( dx*dx+dy*dy < nLenTh ) nFNum++;
			}
			nSNum = 0;
			for ( j=0; j<pSearch->nNumber; j++ ){
				for ( k=0; k<pPair->nNumber; k++ ){
					if ( j == pPair->nSearchID[k] ) break;
				}
				if ( k < pPair->nNumber ) continue;
				dx = pSearch->item[sid].x - pSearch->item[j].x;
				dy = pSearch->item[sid].y - pSearch->item[j].y;
				if ( dx*dx+dy*dy < nLenTh ) nSNum++;
			}
			if ( nSNum == 0 && nFNum >= 3 ) nDiffNum++;
		}
		else if ( pSearch->item[sid].score >= 45 && pFile->item[fid].score >= 30 ){
			nFNum = 0;
			for ( j=0; j<pFile->nNumber; j++ ){
				for ( k=0; k<pPair->nNumber; k++ ){
					if ( j == pPair->nFileID[k] ) break;
				}
				if ( k < pPair->nNumber ) continue;
				dx = pFile->item[fid].x - pFile->item[j].x;
				dy = pFile->item[fid].y - pFile->item[j].y;
				if ( dx*dx+dy*dy < nLenTh ) nFNum++;
			}
			nSNum = 0;
			for ( j=0; j<pSearch->nNumber; j++ ){
				if ( pSearch->item[j].score < 45 ) continue;
				for ( k=0; k<pPair->nNumber; k++ ){
					if ( j == pPair->nSearchID[k] ) break;
				}
				if ( k < pPair->nNumber ) continue;
				dx = pSearch->item[sid].x - pSearch->item[j].x;
				dy = pSearch->item[sid].y - pSearch->item[j].y;
				if ( dx*dx+dy*dy < nLenTh ) nSNum++;
			}
			if ( nFNum == 0 && nSNum >= 3 ) nDiffNum++;
		}
	}
	if ( nDiffNum >= 2 ) nRetScore = nRetScore*3/4;
	return nRetScore;
}

int dec_func_03(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
					  int nRidge,int nBlockScore)
{
	int retScore = score, th = (LENGTH+5)*(LENGTH+5);
	int i, j, len, fx, fy, sx, sy, nFlag = 0, nRTh = 220;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);
	if ( nFCoreNum==0 || nSCoreNum==0 ) return retScore;

	for ( i=0; i<nFCoreNum; i++ ){
		fx = FileCore[i].x; fy = FileCore[i].y;
		for ( j=0; j<nSCoreNum; j++ ){
			sx = SearchCore[j].x; sy = SearchCore[j].y;
			len = (fx-sx)*(fx-sx) + (fy-sy)*(fy-sy);
			if ( len < th ){ nFlag = 1; i = 3; break; }
		}
	}
	if ( nFlag == 0 ) return retScore;
	if ( nRidge <= 235 && nBlockScore < 87 ) retScore = retScore/2;
	else if ( nRidge <= 222 && nBlockScore <= 96 ) retScore = 2*retScore/3;
	else if ( nRidge < 245 && nBlockScore < 85 ) retScore = 2*retScore/3;
	else if ( nRidge <= 235 && nBlockScore <= 94 ) retScore = 9*retScore/10;
	else if ( nRidge <= 240 && nBlockScore <= 93 ) retScore = 9*retScore/10;
	return retScore;
}

int dec_func_06(int score,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int nRetScore = score;
	int i, fid, sid, nNum = 0, maxscore, minscore;

	if ( pPair->nNumber == 0 ) return 0;
	for ( i=0; i<pPair->nNumber; i++ ){
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		maxscore = max(pFile->item[fid].score,pSearch->item[sid].score);
		minscore = min(pFile->item[fid].score,pSearch->item[sid].score);
		if ( maxscore < 40 && minscore < 25 )  continue;
		if ( pFile->item[fid].kind != pSearch->item[sid].kind ) nNum++;
	}
	if ( pPair->nNumber >= 5 ){
		if ( nNum*10 >= pPair->nNumber*7 ) nRetScore = nRetScore/2;
		else if ( nNum*3 >= pPair->nNumber*2 ) nRetScore = nRetScore*6/8;
		else if ( nNum*2 >= pPair->nNumber ) nRetScore = nRetScore*7/8;
	}
	return nRetScore;
}

int dec_func_07(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
						PAIRVECT *pPair,int nBlockScore,int nRidge, int mScore)
{
	int nRetScore = score, th = (LENGTH+5)*(LENGTH+5);
	int i, j, k, dx, dy, diff, num0, num1, num2, dScore = 4, sTh = 3;
	int nNum = pPair->nNumber, num;
	int nMinNum = min(pFile->Mp.nNumber,pSearch->Mp.nNumber);

	if ( nNum > 9 && nNum*100 > nMinNum*30 ){
		num = 0;
		for ( i=0; i<pFile->Mp.nNumber; i++ ){
			if ( pFile->Mp.item[i].score < 35 ) continue;
			for ( j=0; j<nNum; j++ ){
				if ( i == pPair->nFileID[j] ) break;
			}
			if ( j < nNum ) continue;
			if ( check_exist(pFile->Mp.item[i].x,pFile->Mp.item[i].y,pFile->Mp.item[i].dir,
						-1,20,ANGLE,&pSearch->Mp,pPair,FALSE,FALSE,TRUE) == TRUE ){
				num++;
			}
		}
		num = nNum + num;
		if ( 100*num >= nMinNum*63 && nNum >= 12 && nRidge+nBlockScore >= 345 ) return nRetScore;
		if ( 100*num > nMinNum*60 && nNum >= 13 && nBlockScore > 88 && nRidge > 237 ) return nRetScore;
		if ( 100*num > nMinNum*80 && nNum >= 11 && nBlockScore > 90 ) return nRetScore;
	}
	num = 0;
	if ( nNum*100 >= nMinNum*40 ){
		if ( nBlockScore >= 96 && nNum > 6 && nRidge > 245 ) dScore = 2;
		if ( nBlockScore >= 92 && nNum >= 10 && nRidge >= 250 ) dScore = 2;
		else if ( nBlockScore >= 92 && nNum > 10 && nRidge > 245 ) dScore = 3;
	}
	if ( nNum <= 6 ){
		if ( nNum <= 5 && nBlockScore <= 92 ) dScore = 5;
		if ( nBlockScore <= 89 || (nBlockScore <= 91 && nRidge < 245) ) dScore = 5;
	}
	if ( mScore > 1280 && nNum >= 12 && nNum*100 >= nMinNum*33 ) dScore = 3;
	if ( mScore < 350 && nNum <= 7 ) dScore = 5; 
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		if ( pFile->Mp.item[i].score <= 25 ) continue;
		num++;
		for ( j=0; j<pPair->nNumber; j++ ){
			if ( i == pPair->nFileID[j] ) break;
		}
		if ( j < pPair->nNumber ) continue;
		num0 = num1 = num2 = 0;
		for ( j=0; j<pSearch->Mp.nNumber; j++ ){
			if ( pSearch->Mp.item[j].score <= 25 ) continue;
			for ( k=0; k<pPair->nNumber; k++ ){
				if ( j == pPair->nSearchID[k] ) break;
			}
			if ( k < pPair->nNumber ) continue;
			dx = pFile->Mp.item[i].x - pSearch->Mp.item[j].x;
			dy = pFile->Mp.item[i].y - pSearch->Mp.item[j].y;
			dx = dx*dx + dy*dy;
			if ( dx >= th ) continue;
			num0++;
	 		diff = abs(pFile->Mp.item[i].dir-pSearch->Mp.item[j].dir);
			if ( diff >= 120 ) diff = 240 - diff;
			if ( diff >= 8*ANGLE ){ num2++; nRetScore -= 4; continue;}
			if ( diff >= 4*ANGLE ){ num1++; nRetScore -= 2; continue;}
			if ( diff >= 2*ANGLE ){ nRetScore -= 1; continue;}
		}
		if ( num0 == 0 ) nRetScore -= dScore;
	}

	if ( num == 0 ) nRetScore /= 2;
	if ( nRetScore < 0 ) nRetScore = 0;

	return nRetScore;
}

int dec_func_08(int score,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	POLYGON polyF, polyS;
	int retScore = score, nTh = 40*40, nPairNum = pPair->nNumber;
	int *nFileID = &(pPair->nFileID[0]);
	int *nSearchID = &(pPair->nSearchID[0]);
	int i, j, k, h, fid1, fid2, sid1, sid2, flen, slen, len;
	int fx1, fy1, fx2, fy2, sx1, sy1, sx2, sy2, fx, fy, sx, sy;
	int x, y, dir, nFalseNum = 0, fnum = 0, snum = 0;
	int flist[MAX_MINUTIA_NUMBER], slist[MAX_MINUTIA_NUMBER];

	if ( get_polygon_points(pSearch,&polyS) == FALSE ) return retScore;
	if ( get_polygon_points(pFile,&polyF) == FALSE ) return retScore;

	for ( i=0; i<nPairNum-1; i++ ){
		fid1 = nFileID[i]; sid1 = nSearchID[i];
		if ( pFile->item[fid1].score < 45 ) continue;
		if ( pSearch->item[sid1].score < 45 ) continue;
		fx1 = pFile->item[fid1].x; fy1 = pFile->item[fid1].y;
		sx1 = pSearch->item[sid1].x; sy1 = pSearch->item[sid1].y;
		for ( j=i+1; j<nPairNum; j++ ){
			fid2 = nFileID[j]; sid2 = nSearchID[j];
			if ( pFile->item[fid2].score < 45 ) continue;
			if ( pSearch->item[sid2].score < 45 ) continue;
			fx2 = pFile->item[fid2].x; fy2 = pFile->item[fid2].y;
			sx2 = pSearch->item[sid2].x; sy2 = pSearch->item[sid2].y;
			flen = (fx1-fx2)*(fx1-fx2) + (fy1-fy2)*(fy1-fy2);
			slen = (sx1-sx2)*(sx1-sx2) + (sy1-sy2)*(sy1-sy2);
			fx = (fx1 + fx2)/2; fy = (fy1 + fy2)/2;
			sx = (sx1 + sx2)/2; sy = (sy1 + sy2)/2;
			if ( max(flen,slen) >= nTh ) continue;
			for ( k=0; k<pFile->nNumber; k++ ){
				if ( k == fid1 || k == fid2 ) continue;
				if ( pFile->item[k].score < 45 ) continue;
				for ( h=0; h<fnum; h++ ){
					if ( flist[h] == k ) break;
				}
				if ( h < fnum ) continue;
				for ( h=0; h<nPairNum; h++ ){
					if ( k == nFileID[h] ) break;
				}
				if ( h < nPairNum ) continue;
				len = (fx-pFile->item[k].x)*(fx-pFile->item[k].x)
					+ (fy-pFile->item[k].y)*(fy-pFile->item[k].y);
				if ( len > nTh ) continue;
				x = sx - (fx - pFile->item[k].x);
				y = sy - (fy - pFile->item[k].y);
				dir = pFile->item[k].dir;
				if ( check_in_polygon(x,y,&polyS,-1) == FALSE ) continue;
				if ( check_exist(x,y,dir,-1,LENGTH+10,25,pSearch,pPair,TRUE,FALSE,TRUE) == FALSE ){
					retScore = retScore*9/10;
					flist[fnum++] = k;
				}
			}
			for ( k=0; k<pSearch->nNumber; k++ ){
				if ( k == sid1 || k == sid2 ) continue;
				if ( pSearch->item[k].score < 45 ) continue;
				for ( h=0; h<snum; h++ ){
					if ( slist[h] == k ) break;
				}
				if ( h < snum ) continue;
				for ( h=0; h<nPairNum; h++ ){
					if ( k == nSearchID[h] ) break;
				}
				if ( h < nPairNum ) continue;
				len = (sx-pSearch->item[k].x)*(sx-pSearch->item[k].x)
					+ (sy-pSearch->item[k].y)*(sy-pSearch->item[k].y);
				if ( len > nTh ) continue;
				x = fx - (sx - pSearch->item[k].x);
				y = fy - (sy - pSearch->item[k].y);
				dir = pSearch->item[k].dir;
				if ( check_in_polygon(x,y,&polyF,0) == FALSE ) continue;
				if ( check_exist(x,y,dir,-1,LENGTH+10,25,pFile,pPair,TRUE,FALSE,FALSE) == FALSE ){
					retScore = retScore*9/10;
					slist[snum++] = k;
				}
			}
		}
	}
	return retScore;
}

int dec_func_09(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int retScore = score;
	int fx, fy, fdir, sx, sy, sdir, len, diff, nFlag = 0;
	int fx1, fy1, fdir1, sx1, sy1, sdir1;
	int *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber, num;
	POLYGON polyF, polyS;
	int i, j, k, th = (LENGTH+5)*(LENGTH+5);
	int nFCoreNum, nSCoreNum;
	MPVECTEX *pFMp, *pSMp;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&pFile->Core,FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&pSearch->Core,SearchCore,NULL,NULL);

	pFMp = &(pFile->Mp); pSMp = &(pSearch->Mp);

	if ( nFCoreNum==0 || nSCoreNum==0 ) return retScore;
	for ( i=0; i<nFCoreNum; i++ ){
		fx = FileCore[i].x; fy = FileCore[i].y;
		fdir = FileCore[i].dir;
		for ( j=0; j<nSCoreNum; j++ ){
			sx = SearchCore[j].x; sy = SearchCore[j].y;
			sdir = SearchCore[j].dir;
			diff = abs(fdir - sdir);
			if ( diff >= 120 ) diff = 240 - diff;
			if ( diff > ANGLE ) continue;
			len = (fx-sx)*(fx-sx) + (fy-sy)*(fy-sy);
			if ( len < LENGTH*LENGTH ){ nFlag = 1; i = 3; break; }
		}
	}
	if ( nFlag == 0 ) return retScore;
	if ( get_polygon_points(pFMp,&polyF) == FALSE ) return retScore;
	if ( get_polygon_points(pSMp,&polyS) == FALSE ) return retScore;
	for ( i=0; i<pFMp->nNumber; i++ ){
		if ( pFMp->item[i].score < 20 ) continue;
		for ( j=0; j<nPairNum; j++ ){ if (i == nFileID[j]) break; }
		if ( j < nPairNum ) continue;
		fx1 = pFMp->item[i].x; fy1 = pFMp->item[i].y;
		fdir1 = pFMp->item[i].dir;
		len = (fx-fx1)*(fx-fx1) + (fy-fy1)*(fy-fy1);
		if ( len < LENGTH*LENGTH ) continue;
		if ( len >= 50*50 ) continue;
		if ( check_in_polygon(fx1,fy1,&polyS,0) == FALSE ) continue;
		num = 0;
		for ( j=0; j<pSMp->nNumber; j++ ){
			for ( k=0; k<nPairNum; k++ ){ if (j == nSearchID[k]) break; }
			if ( k < nPairNum ) continue;
			sx1 = pSMp->item[j].x; sy1 = pSMp->item[j].y;
			sdir1 = pSMp->item[j].dir;
			len = (sx-sx1)*(sx-sx1) + (sy-sy1)*(sy-sy1);
			if ( len < LENGTH*LENGTH ) continue;
			if ( len >= 50*50 ) continue;
			len = (fx1-sx1)*(fx1-sx1) + (fy1-sy1)*(fy1-sy1);
			if ( len >= th ) continue;
			num++;
	 		diff = abs(fdir1 - sdir1);
			if ( diff >= 120 ) diff = 240 - diff;
			if ( diff >= 4*ANGLE ){ retScore -= 5; continue; }
			if ( 2*diff >= 5*ANGLE ){ retScore -= 2; continue; }
		}
		if ( num == 0 ){
			if ( pFMp->item[i].score < 40 ) continue;
			retScore -= 7;
		}
	}
	if ( retScore < 0 ) retScore = 0;
	return retScore;
}

BOOL check_overlap(LPCOREVECTEX pFile,LPCOREVECTEX pSearch)
{
	int i, j, len, diff, fx, fy, fdir, sx, sy, sdir;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(pFile,FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(pSearch,SearchCore,NULL,NULL);

	if ( nFCoreNum==0 || nSCoreNum==0 ) return FALSE;
	for ( i=0; i<nFCoreNum; i++ ){
		fx = FileCore[i].x; fy = FileCore[i].y;
		fdir = FileCore[i].dir;
		for ( j=0; j<nSCoreNum; j++ ){
			sx = SearchCore[j].x; sy = SearchCore[j].y;
			sdir = SearchCore[j].dir;
			len = (fx-sx)*(fx-sx) + (fy-sy)*(fy-sy);
			diff = abs(fdir-sdir);
			if ( diff >= 120 ) diff = 240 - diff;
			if ( len < 16*16 && diff < 7 ) return TRUE;
		}
	}
	return FALSE;
}

int dec_func_04(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int retScore = score;
	int *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber;
	BOOL nCoreFlag = TRUE, nFalse = FALSE;
	int sid, fid, fx, fy, sx, sy, flen, slen, sid1, fid1, sid2, fid2;
	int fdir1, fdir2, fdir, sdir1, sdir2, sdir, diff;
	int i, j, k;

	if ( nPairNum < 3 || nPairNum > 10 ) return retScore;
	if ( pFile->Core.nNumber==0 || pSearch->Core.nNumber==0 ){
		nCoreFlag = FALSE;
	}
	else{
		if ( check_overlap(&pFile->Core,&pSearch->Core) == FALSE ){
			nCoreFlag = FALSE;
		}
	}
	for ( i=0; i<nPairNum; i++ ){
		sid = nSearchID[i]; fid = nFileID[i];
		fx = pFile->Mp.item[fid].x; fy = pFile->Mp.item[fid].y;
		sx = pSearch->Mp.item[sid].x; sy = pSearch->Mp.item[sid].y;
		if ( pFile->Mp.item[fid].score < 30 ) continue;
		if ( pSearch->Mp.item[sid].score < 30 ) continue;
		for ( j=0; j<nPairNum; j++ ){
			if ( i == j ) continue;
			sid1 = nSearchID[j]; fid1 = nFileID[j];
			if ( pFile->Mp.item[fid1].score < 20 ) continue;
			if ( pSearch->Mp.item[sid1].score < 20 ) continue;
			flen = (fx-pFile->Mp.item[fid1].x)*(fx-pFile->Mp.item[fid1].x) 
				 + (fy-pFile->Mp.item[fid1].y)*(fy-pFile->Mp.item[fid1].y);
			slen = (sx-pSearch->Mp.item[sid1].x)*(sx-pSearch->Mp.item[sid1].x) 
				 + (sy-pSearch->Mp.item[sid1].y)*(sy-pSearch->Mp.item[sid1].y);
			if ( flen >= 9*LENGTH*LENGTH || slen >= 9*LENGTH*LENGTH ) continue;
			fdir1 = op_func_01(fx,fy,pFile->Mp.item[fid1].x,pFile->Mp.item[fid1].y);
			sdir1 = op_func_01(sx,sy,pSearch->Mp.item[sid1].x,pSearch->Mp.item[sid1].y);
			for ( k=0; k<nPairNum; k++ ){
				if ( k == j || k == i ) continue;
				sid2 = nSearchID[k]; fid2 = nFileID[k];
				if ( pFile->Mp.item[fid2].score < 20 ) continue;
				if ( pSearch->Mp.item[sid2].score < 20 ) continue;
				flen = (fx-pFile->Mp.item[fid2].x)*(fx-pFile->Mp.item[fid2].x) 
					 + (fy-pFile->Mp.item[fid2].y)*(fy-pFile->Mp.item[fid2].y);
				slen = (sx-pSearch->Mp.item[sid2].x)*(sx-pSearch->Mp.item[sid2].x) 
					 + (sy-pSearch->Mp.item[sid2].y)*(sy-pSearch->Mp.item[sid2].y);
				if ( flen >= 9*LENGTH*LENGTH || slen >= 9*LENGTH*LENGTH ) continue;
				fdir2 = op_func_01(fx,fy,pFile->Mp.item[fid2].x,pFile->Mp.item[fid2].y);
				sdir2 = op_func_01(sx,sy,pSearch->Mp.item[sid2].x,pSearch->Mp.item[sid2].y);
				fdir = abs(fdir1 - fdir2);
				if ( fdir >= 120 ) fdir = 240 - fdir;
				sdir = abs(sdir1 - sdir2);
				if ( sdir >= 120 ) sdir = 240 - sdir;
				diff = abs(fdir - sdir);
				if ( diff >= 120 ) diff = 240 - diff;
				if ( diff >= ANGLE*3 ){
					j = 100; i = 100; nFalse = TRUE; break;
				}
			}
		}
	}
	if ( nFalse == TRUE ){
		if ( nCoreFlag == FALSE ) retScore = retScore*9/10;
		else retScore = retScore*7/10;
	}
	return retScore;
}

void get_neighbor(int cx,int cy,LPMPVECTEX pVect,int *pPairID,int nPairNum,
				 BOOL nPairFlag,int nLenTh,BOOL nScoreFlag,int nScoreTh,
				 BOOL nNumFlag,int nNumTh,LPMPVECTEX pNewVect)
{
	int i, x, y, len, th = nLenTh*nLenTh;
	int j, nNum = 0, list[100], lenlist[100], nMinIdx, nMinValue, tmp;

	pNewVect->nNumber = 0;
	for (i=0; i<pVect->nNumber; i++ ){
		if ( nScoreFlag == TRUE ){
			if ( pVect->item[i].score < nScoreTh ) continue;
		}
		x = pVect->item[i].x; y = pVect->item[i].y;
		if ( x == cx && y == cy ) continue;
		if ( nPairFlag == TRUE ){
			for ( j=0; j<nPairNum; j++ ){
				if ( i == pPairID[j] ) break;
			}
			if ( j < nPairNum ) continue;
		}
		len = (x-cx)*(x-cx) + (y-cy)*(y-cy);
		if ( len >= th ) continue;
		list[nNum] = i; lenlist[nNum++] = len;
		pNewVect->item[pNewVect->nNumber++] = pVect->item[i];
	}
	if ( nNumFlag == TRUE ){
		if ( nNum > nNumTh ){
			for ( i=0; i<nNum-1; i++ ){
				nMinIdx = i;
				nMinValue = lenlist[i];
				for ( j=i+1; j<nNum; j++ ){
					if ( lenlist[j] >= nMinValue) continue;
					nMinIdx = j; nMinValue = lenlist[j];
				}
				if ( nMinIdx == i ) continue;
				tmp = list[i]; list[i] = list[nMinIdx]; list[nMinIdx] = tmp;
				tmp = lenlist[i]; lenlist[i] = lenlist[nMinIdx]; lenlist[nMinIdx] =tmp;
			}
			pNewVect->nNumber = 0;
			for ( i=0; i<nNumTh; i++ ){
				pNewVect->item[pNewVect->nNumber++] = pVect->item[list[i]];
			}
		}
	}
}

BOOL check_neighbor(int nFid,int nSid,LPMPVECTEX tmpF,LPMPVECTEX tmpS,
				   LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	POLYGON polyF, polyS;
	int i, x, y, dir, nFNum = 0, nSNum = 0;

	if ( get_polygon_points(&pFile->Mp,&polyF) == FALSE ) return TRUE;
	if ( get_polygon_points(&pSearch->Mp,&polyS) == FALSE ) return TRUE;
	for ( i=0; i<tmpF->nNumber; i++ ){
		x = tmpF->item[i].x; y = tmpF->item[i].y;
		dir = tmpF->item[i].dir;
		if ( check_in_polygon(x,y,&polyS,0) == FALSE ) continue;
		if ( check_exist(x,y,dir,nSid,20,20,&pSearch->Mp,NULL,FALSE,FALSE,FALSE) == FALSE ){
			nFNum++;
		}
	}
	for ( i=0; i<tmpS->nNumber; i++ ){
		x = tmpS->item[i].x; y = tmpS->item[i].y;
		dir = tmpS->item[i].dir;
		if ( check_in_polygon(x,y,&polyF,0) == FALSE ) continue;
		if ( check_exist(x,y,dir,nFid,20,20,&pFile->Mp,NULL,FALSE,FALSE,FALSE) == FALSE ){
			nSNum++;
		}
	}
	if ( nSNum > 0 && nSNum == tmpS->nNumber ) return FALSE;
	if ( nFNum > 0 && nFNum == tmpF->nNumber ) return FALSE;
	return TRUE;
}

int dec_func_10(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair,int nBlockScore,int nRidge, int mScore )
{
	int retScore = score;
	int *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber;
	BOOL nCoreFlag, nScoreFlag = TRUE;
	MPVECTEX tmpF, tmpS;
	int scoreF = 0, scoreS = 0, nTotalNum, nFalseNum;
	int i, fid, sid, fx, fy, sx, sy;

	if ( nPairNum > 10 ) return retScore;
	if ( nPairNum >= 10 && nBlockScore >= 95 && nRidge > 250 && mScore > 900 ) return retScore;

	if ( pFile->Mp.nNumber==0 || pSearch->Mp.nNumber==0 ) return 0;
	nCoreFlag = check_overlap(&pFile->Core,&pSearch->Core);
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		scoreF += pFile->Mp.item[i].score;
	}
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		scoreS += pSearch->Mp.item[i].score;
	}
	scoreF /= pFile->Mp.nNumber; scoreS /= pSearch->Mp.nNumber;
	if ( scoreF < 35 || scoreS < 35 ) nScoreFlag = FALSE;

	nTotalNum = 0; nFalseNum = 0;
	for ( i=0; i<nPairNum; i++ ){
		fid = nFileID[i]; sid = nSearchID[i];
		if ( pFile->Mp.item[fid].score < 30 ) continue;
		if ( pSearch->Mp.item[sid].score < 30 ) continue;
		fx = pFile->Mp.item[fid].x; fy = pFile->Mp.item[fid].y;
		sx = pSearch->Mp.item[sid].x; sy = pSearch->Mp.item[sid].y;
		nTotalNum++;
		get_neighbor(fx,fy,&pFile->Mp,nFileID,nPairNum,TRUE,50,TRUE,30,FALSE,0,&tmpF);
		get_neighbor(sx,sy,&pSearch->Mp,nSearchID,nPairNum,TRUE,50,TRUE,30,FALSE,0,&tmpS);
		if ( check_neighbor(fid,sid,&tmpF,&tmpS,pFile,pSearch) == FALSE ){
			nFalseNum++;
		}
	}
	if ( nTotalNum > 0 ){
		if ( nFalseNum >= 5 ){
			if ( nCoreFlag == TRUE && nScoreFlag == TRUE )
				retScore = retScore/2;
			else
				retScore = retScore*33/50;
		}
		else if ( nFalseNum >= 4 ){
			if ( nScoreFlag == TRUE ) retScore -= nFalseNum*6;
			else retScore -= nFalseNum*5;
		}
		else if ( nFalseNum >= 3 ){
			if ( nScoreFlag == TRUE || nFalseNum*4 >= nTotalNum*3 ) retScore -= nFalseNum*5;
			else retScore -= 10*nFalseNum/3;
		}
		else if ( nFalseNum >= 2 ) retScore -= 3;
		else if ( nFalseNum >= 1 ) retScore -= 1;
	}
	else{
		retScore = retScore*84/100;
	}
	return retScore;
}

BOOL check_core(LPCOREVECTEX pFile,LPCOREVECTEX pSearch,int nLenTh,int nAngTh)
{
	int len, flen, slen, diff, fx, fy, fdir, sx, sy, sdir;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(pFile,FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(pSearch,SearchCore,NULL,NULL);

	if ( nFCoreNum==0 || nSCoreNum==0 ) return FALSE;

	if ( nFCoreNum != nSCoreNum ) return FALSE;

	if ( nFCoreNum == 1 ){
		fx = FileCore[0].x; fy = FileCore[0].y;
		sx = SearchCore[0].x; sy = SearchCore[0].y;
		len = op_func_02((fx-sx)*(fx-sx)+(fy-sy)*(fy-sy));
		fdir = FileCore[0].dir; sdir = SearchCore[0].dir;
	}
	else{
		flen = op_func_02((FileCore[0].x-FileCore[1].x)*(FileCore[0].x-FileCore[1].x)
					+ (FileCore[0].y-FileCore[1].y)*(FileCore[0].y-FileCore[1].y));
		slen = op_func_02((SearchCore[0].x-SearchCore[1].x)*(SearchCore[0].x-SearchCore[1].x)
					+ (SearchCore[0].y-SearchCore[1].y)*(SearchCore[0].y-SearchCore[1].y));
		len = abs(flen-slen);
		fdir = op_func_01(FileCore[0].x,FileCore[0].y,
						  FileCore[1].x,FileCore[1].y);
		if ( fdir >= 120 ) fdir -= 120;
		sdir = op_func_01(SearchCore[0].x,SearchCore[0].y,
						  SearchCore[1].x,SearchCore[1].y);
		if ( sdir >= 120 ) sdir -= 120;
	}
	diff = abs(fdir-sdir);
	if ( diff >= 120 ) diff = 240 - diff;
	if ( len < nLenTh && diff < nAngTh ) return TRUE;
	return FALSE;
}

int dec_func_05(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
						PAIRVECT *pPair,int nMatchScoreTh,int nRidge,
						int nBlockScore,int score_old)
{
	int retScore = score;
	MPVECTEX tmpF, tmpS;
	POLYGON polyF, polyS;
	int *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber, nFScore = 0, nSScore = 0;
	int i, j, x, y, dir, nFNum = 0, nSNum = 0;
	BOOL nCoreFlag, nScoreFlag = TRUE;

	if ( pFile->Mp.nNumber==0 || pSearch->Mp.nNumber==0 ) return 0;
	if ( nPairNum == 0 ) return 0;
	if ( nPairNum > 10 || score >= nMatchScoreTh*3 || score_old >= nMatchScoreTh*4 ) return retScore;
	if ( nRidge >= 245 && nBlockScore >= 95 ) return retScore;

	nCoreFlag = check_core(&pFile->Core,&pSearch->Core,16,7);
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		nFScore += pFile->Mp.item[i].score;
	}
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		nSScore += pSearch->Mp.item[i].score;
	}
	nFScore /= pFile->Mp.nNumber; nSScore /= pSearch->Mp.nNumber;
	if ( nFScore < 35 || nSScore < 35 ) nScoreFlag = FALSE;

	tmpF.nNumber = tmpS.nNumber = nPairNum;
	for ( i=0; i<nPairNum; i++ ){
		tmpF.item[i] = pFile->Mp.item[nFileID[i]];
		tmpS.item[i] = pSearch->Mp.item[nSearchID[i]];
	}
	if ( get_polygon_points(&tmpF,&polyF) == FALSE ) return retScore;
	if ( get_polygon_points(&tmpS,&polyS) == FALSE ) return retScore;

	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		if ( pFile->Mp.item[i].score < 30 ) continue;
		for ( j=0; j<nPairNum; j++ ){ if (i==nFileID[j]) break; }
		if ( j < nPairNum ) continue;
		x = pFile->Mp.item[i].x; y = pFile->Mp.item[i].y;
		dir = pFile->Mp.item[i].dir;
		if ( check_in_polygon(x,y,&polyF,0) == FALSE ) continue;
		if ( check_exist(x,y,dir,-1,20,20,&pSearch->Mp,NULL,FALSE,FALSE,TRUE) == FALSE ){
			nFNum++;
		}
	}
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		if ( pSearch->Mp.item[i].score < 30 ) continue;
		for ( j=0; j<nPairNum; j++ ){ if (i==nSearchID[j]) break; }
		if ( j < nPairNum ) continue;
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		dir = pSearch->Mp.item[i].dir;
		if ( check_in_polygon(x,y,&polyS,0) == FALSE ) continue;
		if ( check_exist(x,y,dir,-1,20,20,&pFile->Mp,NULL,FALSE,FALSE,FALSE) == FALSE ){
			nSNum++;
		}
	}
	nFNum = nFNum + nSNum;
	if ( nFNum >= 5 ){
		if ( nCoreFlag == TRUE && nScoreFlag == TRUE )
			retScore = retScore/2;
		else
			retScore = retScore*7/10;
	}
	else if ( nFNum >= 3 ){
		if ( nScoreFlag == TRUE ) retScore -= nFNum*3;
		else retScore -= nFNum*2;
	}
	return retScore;
}

int dec_func_11(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
					   PAIRVECT *pPair,int nMatchScoreTh,int nRidge,int nBlockScore)
{
	int retScore = score;
	int i, dx, dy, len, fid, sid, nFalseNum = 0;

	if ( pPair->nNumber >= 10 ) return retScore;
	if ( nRidge >= 250 && nBlockScore >= 90 ) return retScore;

	for ( i=0; i<pPair->nNumber; i++ ){
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		if ( pFile->Mp.item[fid].score < 30 ) continue;
		if ( pSearch->Mp.item[sid].score < 30 ) continue;
		dx = pFile->Mp.item[fid].x-pSearch->Mp.item[sid].x;
		dy = pFile->Mp.item[fid].y-pSearch->Mp.item[sid].y;
		len = op_func_02(dx*dx+dy*dy);
		if ( len >= 8  && (pFile->Mp.item[fid].kind != pSearch->Mp.item[sid].kind)) nFalseNum++;
	}
	if ( nFalseNum >= 5 ){
		retScore = retScore/2;
	}
	else if ( nFalseNum >= 3 ){
		if ( nFalseNum >= 4 ) retScore -= nFalseNum*5;
		else retScore -= nFalseNum*3;
	}
	else if ( nFalseNum >= 2 ) retScore -= 2;
	else if ( nFalseNum >= 1 ) retScore -= 1;

	return retScore;
}

int dec_func_12(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
					   PAIRVECT *pPair,int nMatchScoreTh,int nBlockScore)
{
	int retScore = score;
	int i, diff, nFalseNum = 0;

	if ( pPair->nNumber >= 15 || score >= nMatchScoreTh*3/2 ) return retScore;

	for ( i=0; i<pPair->nNumber; i++ ){
		diff = abs(pFile->Mp.item[pPair->nFileID[i]].score - 
				   pSearch->Mp.item[pPair->nSearchID[i]].score);
		if ( diff >= 20 ) nFalseNum++;
	}
	if ( nFalseNum >= 5 && nFalseNum*2 >= pPair->nNumber ){
		retScore = retScore/2;
	}
	else if ( nFalseNum >= 3 ){
		if ( nBlockScore >= 97 ) retScore -= 2*nFalseNum;
		else if ( nFalseNum*7 >= pPair->nNumber*3 && nBlockScore < 83 ) retScore -= 13*nFalseNum/3;
		else retScore -= 10*nFalseNum/3;
	}
	else if ( nFalseNum >= 2 ){
		retScore -= nFalseNum*2;
	}
	else if ( nFalseNum >= 1 ){
		retScore -= 1;
	}
	return retScore;
}

int dec_func_13(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int nRetScore = score;
	int i, j, k, dx, dy, fx, fy, sx, sy, fid, sid;
	int num = 0, num1 = 0, fnum, snum;
	
	for ( i=0; i<pPair->nNumber; i++ ){
		fid = pPair->nFileID[i];
		fx = pFile->Mp.item[fid].x; fy = pFile->Mp.item[fid].y;
		fnum = 0;
		sid = pPair->nSearchID[i];
		sx = pSearch->Mp.item[sid].x; sy = pSearch->Mp.item[sid].y;
		snum = 0;
		if ( pSearch->Mp.item[sid].score < 40 ) continue;
		for ( k=0; k<pFile->Mp.nNumber; k++ ){
			for ( j=0; j<pPair->nNumber; j++ ){
				if ( k == pPair->nFileID[j] ) break;
			}
			if ( j < pPair->nNumber ) continue;
			dx = fx - pFile->Mp.item[k].x; dy = fy - pFile->Mp.item[k].y;
			if ( dx*dx+dy*dy < 45*45 ) fnum++;
		}
		for ( k=0; k<pSearch->Mp.nNumber; k++ ){
			for ( j=0; j<pPair->nNumber; j++ ){
				if ( k == pPair->nSearchID[j] ) break;
			}
			if ( j < pPair->nNumber ) continue;
			dx = sx - pSearch->Mp.item[k].x; dy = sy - pSearch->Mp.item[k].y;
			if ( dx*dx+dy*dy < 45*45 ) snum++;
		}
		if ( fnum >= 3 && snum == 0 ) num++;
		if ( fnum >= 5 && snum == 1 ) num1++;
	}
	if ( num >= 3 ) nRetScore /= 2;
	if ( num == 2 ) nRetScore = nRetScore*7/10;
	if ( num1 > 0 ) nRetScore = nRetScore*9/10;
	return nRetScore;
}

void transform_block(int nRot,int xOff,int yOff,int cx,int cy,BLOCKVECT* pBlock)
{
	BYTE pTmp[MAX_BLOCK_NUMBER];
	int c, nCos, nSin, d1, d2, i, j, x, y, bx, by;

	memcpy(pTmp,pBlock->Data,sizeof(BYTE)*MAX_BLOCK_NUMBER);
	memset(pBlock->Data,0xFF,sizeof(BYTE)*MAX_BLOCK_NUMBER);

	nCos = _table_03[nRot]; nSin = _table_04[nRot];

	for( i=0; i<pBlock->nRow; i++ ){
		for ( j=0; j<pBlock->nCol; j++ ){
			bx = (j*BLOCK_SIZE + BLOCK_SIZE/2) - xOff; 
			by = (i*BLOCK_SIZE + BLOCK_SIZE/2) - yOff;
			x = (bx-cx)*nCos + (by-cy)*nSin;
			if ( x > 0 ) x += 8192;
			else x -= 8192;
			x = x >> 14;
			y = (by-cy)*nCos - (bx-cx)*nSin;
			if ( y > 0 ) y += 8192;
			else y -= 8192;
			y = y >> 14;

			d1 = x + cx; d2 = y + cy;
			if ( d1 < 0 ) d1 -= BLOCK_SIZE-1;
			if ( d2 < 0 ) d2 -= BLOCK_SIZE-1;
			d1 /= BLOCK_SIZE; d2 /= BLOCK_SIZE;
			if ( d1 < 0 || d1 >= pBlock->nCol ) continue;
			if ( d2 < 0 || d2 >= pBlock->nRow ) continue;

			c = pTmp[d2*pBlock->nCol+d1];
			if ( c == 0xFF ) pBlock->Data[i*pBlock->nCol+j] = 0xFF;
			else{
				c += nRot;
				if ( c >= 240 ) c -= 240;
				if ( c < 0 ) c += 240;
				if ( c >= 120 ) c -= 120;
				pBlock->Data[i*pBlock->nCol+j] = (BYTE)c;
			}
		}
	}
}

int check_block(int Th,int nNumTh,BLOCKVECT* pFBlock,BLOCKVECT* pSBlock)
{
	int i, retVal = 0, num = 0, nDiv = 0;
	BYTE c;

	for ( i=0; i<pFBlock->nCol*pFBlock->nRow; i++ ){
		if ( pFBlock->Data[i] == 0xFF ) continue;
		if ( pSBlock->Data[i] == 0xFF ) continue;
		c = abs(pFBlock->Data[i] - pSBlock->Data[i]);
		if ( c > 60 ) c = 120 - c;
		if ( c < 5 ) c = 0;
		else if ( c > Th ) c = 60;
		retVal += 60 - c; nDiv += 60; num++;
	}
	if ( nNumTh*num < pFBlock->nCol*pFBlock->nRow ) return (0);
	if ( nDiv == 0 ) return (0);
	retVal = 100*retVal/nDiv;
	return (retVal);
}

BOOL arrange_points_sub(int cx,int cy,int nRot,int nXDiff,int nYDiff,
						 LPMPVECTEX pFile,BLOCKVECT* tmpFBlock,LPFPVECTEX pSearch,
						 BARVECT *pSBar,int *nMaxSBarLen,
						 int (*SBarPtr)[20],int *SDiffField
						)
{
	MPVECTEX tmpFMp = *pFile;
	int i, j, k, x, y, bx, by, nNum, nRow, nCol, nDiv;
	int nNewId[MAX_MINUTIA_NUMBER];
	int nList[MAX_BLOCK_NUMBER];
  

	transform_mp(&tmpFMp,cx,cy,nRot,nXDiff,nYDiff);

	nNum = 0;
	for ( i=0; i<tmpFBlock->nCol*tmpFBlock->nRow; i++ ){
		if ( tmpFBlock->Data[i] == 0xFF ) continue;
		if ( pSearch->Block.Data[i] == 0xFF ) continue;
		nList[nNum++] = i;
	}
	if ( nNum == 0 ) return FALSE;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ) nNewId[i] = -1;
	nDiv = tmpFBlock->nCol;
	for ( i=0,k=0; i<tmpFMp.nNumber; i++ ){
		x = tmpFMp.item[i].x; y = tmpFMp.item[i].y;
		for ( j=0; j<nNum; j++ ){
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( abs(x-bx) > LENGTH+8 || abs(y-by) > LENGTH+8 ) continue;
			nNewId[i] = k;
			pFile->item[k++] = pFile->item[i];
			break;
		}
	}
	pFile->nNumber = k;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ) nNewId[i] = -1;
	nDiv = pSearch->Block.nCol;
	for ( i=0,k=0; i<pSearch->Mp.nNumber; i++ ){
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		for ( j=0; j<nNum; j++ ){
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( abs(x-bx) > LENGTH+8 || abs(y-by) > LENGTH+8 ) continue;
			nNewId[i] = k;
			pSearch->Mp.item[k++] = pSearch->Mp.item[i];
			break;
		}
	}
	pSearch->Mp.nNumber = k;

	for ( i=0,j=0; i<pSBar->nNumber; i++ ){
		if ( nNewId[pSBar->item[i].nID1] == -1 ) continue;
		if ( nNewId[pSBar->item[i].nID2] == -1 ) continue;
		pSBar->item[j] = pSBar->item[i];
		pSBar->item[j].nID1 = nNewId[pSBar->item[i].nID1];
		pSBar->item[j++].nID2 = nNewId[pSBar->item[i].nID2];
	}
	if ( j == 0 ) return FALSE;
	pSBar->nNumber = j;

	memset( SDiffField, 0, sizeof(int)*240 );
	*nMaxSBarLen = 0;
	for ( i=0; i<pSBar->nNumber; i++ ){
		if ( *nMaxSBarLen < pSBar->item[i].nLen ){
			*nMaxSBarLen = pSBar->item[i].nLen + 1;
		}
		k = pSBar->item[i].nDiff1;
		SBarPtr[k][SDiffField[k]] = i;
		if ( ++SDiffField[k] == 20 ){ SDiffField[k]--; }
	}
	return TRUE;
}

int arrange_points(LPFPVECTEX pFile,LPFPVECTEX pSearch,BARVECT *pSearchBar,int *nMaxSBarLen,
				   int (*SBarPtr)[20],int *SDiffField,BOOL nCommonFlag)
{
	MPVECTEX tmpMP;
	BLOCKVECT tmpBlk;
	int i, k1, k2, nBlkScore;
	int fcx, fcy, scx, scy, dx, dy, fdir, sdir, rot;
	int maxrot, maxdx, maxdy, cx, cy;
	COREITEMEX FileCore[2], SearchCore[2];
	int nFCoreNum, nSCoreNum;
	BOOL nTypeFlag = TRUE;

	if ( pFile->Mp.nNumber == 0 || pSearch->Mp.nNumber == 0 ) return(-1);

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);
	if ( nFCoreNum == 0 || nSCoreNum == 0 ) return (0);

	if ( nFCoreNum == nSCoreNum ){
		fcx = fcy = scx = scy = 0;
		for ( i=0; i<nFCoreNum; i++ ){
			fcx += FileCore[i].x; fcy += FileCore[i].y;
			scx += SearchCore[i].x; scy += SearchCore[i].y;
		}
		fcx /= nFCoreNum; fcy /= nFCoreNum;
		scx /= nFCoreNum; scy /= nFCoreNum;
		dx = scx - fcx; dy = scy - fcy;
		if ( nFCoreNum == 1 ){
			rot = SearchCore[0].dir - FileCore[0].dir;
			if ( rot < 0 ) rot += 240;
		}
		else{
			fdir = op_func_01(FileCore[0].x,FileCore[0].y,FileCore[1].x,FileCore[1].y);
			sdir = op_func_01(SearchCore[0].x,SearchCore[0].y,SearchCore[1].x,SearchCore[1].y);
			rot = sdir - fdir;
			if ( rot < 0 ) rot += 240;
			if ( rot >= 120 ) rot -= 120;
			tmpMP = pFile->Mp;
			transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
			k1 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
			maxrot = rot;
			rot = 120 + rot;
			if ( rot >= 240 ) rot = rot - 240;
			tmpMP = pFile->Mp;
			transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
			k2 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
			if ( k1 > k2 ) rot = maxrot;
		}
		cx = fcx; cy = fcy;
	}
	else{
		if ( nFCoreNum == 1 ){
			fcx = FileCore[0].x; fcy = FileCore[0].y;
			fdir = FileCore[0].dir;
			k2 = 0;
			for ( i=0; i<2; i++ ){
				dx = SearchCore[i].x - fcx; dy = SearchCore[i].y - fcy;
				rot = SearchCore[i].dir - fdir;
				if ( rot < 0 ) rot += 240;
				tmpMP = pFile->Mp;
				transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
				k1 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
				if ( k2 < k1 ){
					k2 = k1; maxrot = rot; maxdx = dx; maxdy = dy;
				}
			}
			cx = fcx; cy = fcy;
		}
		else{
			scx = SearchCore[0].x; scy = SearchCore[0].y;
			sdir = SearchCore[0].dir;
			k2 = 0;
			for ( i=0; i<2; i++ ){
				fcx = FileCore[i].x; fcy = FileCore[i].y;
				dx = scx - fcx; dy = scy - fcy;
				rot = sdir - FileCore[i].dir;
				if ( rot < 0 ) rot += 240;
				tmpMP = pFile->Mp;
				transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
				k1 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
				if ( k2 < k1 ){
					k2 = k1; maxrot = rot; maxdx = dx; maxdy = dy;
					cx = fcx; cy = fcy;
				}
			}
		}
		if ( k2 == 0 ) return (-1);
		rot = maxrot; dx = maxdx; dy = maxdy;
	}
	if ( nCommonFlag == TRUE ){
		tmpBlk = pFile->Block;
		transform_block(rot,dx,dy,cx,cy,&tmpBlk);
		nBlkScore = check_block(30,5,&tmpBlk,&pSearch->Block);
		if ( nBlkScore > 80 ){
			arrange_points_sub(cx,cy,rot,dx,dy,&pFile->Mp,&tmpBlk,pSearch,pSearchBar,nMaxSBarLen,SBarPtr,SDiffField);
		}
	}
	return(0);
}

int get_score_sub(LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	int i, j, x, y, dir, dx, dy, diff, val, minval, score = 0;

	for ( i=0; i<pVect1->nNumber; i++ ){
		x = pVect1->item[i].x; y = pVect1->item[i].y;
		dir = pVect1->item[i].dir;
		minval = 10000;
		for ( j=0; j<pVect2->nNumber; j++ ){
			dx = abs(pVect2->item[j].x - x);
			if ( dx > LENGTH+3 ) continue;
			dy = abs(pVect2->item[j].y - y);
			if ( dy > LENGTH+3 ) continue;
			diff = abs(pVect2->item[j].dir - dir);
			if ( diff >= 120 ) diff = 240 - diff;
			if ( diff > ANGLE ) continue;
			val = dx + dy + diff;
			if ( minval > val ) minval = val;
			if ( minval < 20) break;
		}
		if ( minval < 35 ) score += (35 - minval);
	}
	val = (pVect1->nNumber + pVect2->nNumber) / 2;
	if ( val == 0 ) return (0);
	score = (score * 100 ) / val;
	return (score);
}

int get_point_score(LPFPVECTEX pFile,LPFPVECTEX pSearch) 
{
	int i, j, k, cx, cy, dx, dy, rot;
	int retScore = 0, score, maxscore;
	int nNumF = 0, nMaxF = 0, nMaxS = 0, coreflag = 0;
	int fcx, fcy, scx, scy, flen, slen, fdir, fdir0, fdir1, sdir, sdir0, sdir1, difdir, difdir1;
	int scoreF = 0, scoreS = 0, nTh = 18, AngTh = 30;
	MP_POINT tmpF[7];
	MPVECTEX tmpMp = pFile->Mp;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	if ( pFile->Mp.nNumber < MIN_MINUTIA_NUMBER || pSearch->Mp.nNumber < MIN_MINUTIA_NUMBER ) return (0);

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);

	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		if ( nMaxF < pFile->Mp.item[i].curv ) nMaxF = pFile->Mp.item[i].curv;
		scoreF += pFile->Mp.item[i].score;
	}
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		if ( nMaxS < pSearch->Mp.item[i].curv ) nMaxS = pSearch->Mp.item[i].curv;
		scoreS += pSearch->Mp.item[i].score;
	}
	if ( nMaxF > nMaxS ) nMaxF = nMaxS;
	scoreF /= pFile->Mp.nNumber; scoreS /= pSearch->Mp.nNumber;
	if ( scoreF < 35 || scoreS < 35 ){ nTh = 30; AngTh = 40; }

	for ( i=0,j=0; i<pFile->Mp.nNumber; i++ ){
		if ( pFile->Mp.item[i].curv >= nMaxF+8 ) continue;
		for ( k=0; k<j; k++ ){
			dx = pFile->Mp.item[i].x - tmpF[k].x;
			dy = pFile->Mp.item[i].y - tmpF[k].y;
			if ( dx*dx+dy*dy <= 30*30 ) break;
		}
		if ( k < j ) continue;
		tmpF[j++] = pFile->Mp.item[i];
		if ( j >= 7 ) break;
	}
	nNumF = j;
	pFile->Mp = tmpMp;

	if ( nFCoreNum > 0 ){
		fcx = fcy = 0;
		for ( i=0; i<nFCoreNum; i++ ){
			fcx += FileCore[i].x; fcy += FileCore[i].y;
		}
		fcx /= nFCoreNum; fcy /= nFCoreNum;
		if ( nFCoreNum == 2 ){
			fdir = op_func_01(FileCore[0].x,FileCore[0].y,FileCore[1].x,FileCore[1].y);
		}
	}
	if ( nSCoreNum > 0 ){
		scx = scy = 0;
		for ( i=0; i<nSCoreNum; i++ ){
			scx += SearchCore[i].x; scy += SearchCore[i].y;
		}
		scx /= nSCoreNum; scy /= nSCoreNum;
		if ( nSCoreNum == 2 ){
			sdir = op_func_01(SearchCore[0].x,SearchCore[0].y,SearchCore[1].x,SearchCore[1].y);
		}
	}
	if ( nFCoreNum == nSCoreNum && nFCoreNum > 0 ) coreflag = 1;

	for ( i=0; i<nNumF; i++ ){
		cx = tmpF[i].x; cy = tmpF[i].y;
		if ( coreflag == 1 ) flen = op_func_02((cx-fcx)*(cx-fcx)+(cy-fcy)*(cy-fcy));
		if ( nFCoreNum == 1 ){
			fdir0 = tmpF[i].dir - FileCore[0].dir;
			if ( fdir0 < 0 ) fdir0 += 240;
		}
		if ( nFCoreNum == 2 ){
			fdir0 = tmpF[i].dir - fdir;
			if ( fdir0 < 0 ) fdir0 += 240;
			fdir1 = fdir0 + 120;
			if ( fdir1 >= 240 ) fdir1 -= 240; 
		}
		maxscore = 0;
		for ( j=0; j<pSearch->Mp.nNumber; j++ ){
			if ( abs(tmpF[i].curv-pSearch->Mp.item[j].curv) >= 6 ) continue;
			rot = pSearch->Mp.item[j].dir - tmpF[i].dir;
			if ( rot < 0 ) rot += 240;

			if ( nSCoreNum == 1 ){
				sdir0 = pSearch->Mp.item[j].dir - SearchCore[0].dir;
				if ( sdir0 < 0 ) sdir0 += 240;
				if ( nFCoreNum == 1 ){
					difdir = abs(fdir0 - sdir0);
					if ( difdir >= 120 ) difdir = 240 - difdir;
					if ( difdir > AngTh ) continue;
				}
				if ( nFCoreNum == 2 ){
					difdir = abs(fdir0 - sdir0);
					if ( difdir >= 120 ) difdir = 240 - difdir;
					difdir1 = abs(fdir1 - sdir0);
					if ( difdir1 >= 120 ) difdir1 = 240 - difdir1;
					if ( difdir > AngTh && difdir1 > AngTh ) continue;
				}
			}
			if ( nSCoreNum == 2 ){
				sdir0 = pSearch->Mp.item[j].dir - sdir;
				if ( sdir0 < 0 ) sdir0 += 240;
				sdir1 = sdir0 + 120;
				if ( sdir1 >= 240 ) sdir1 -= 240;
				if ( nFCoreNum == 1 ){
					difdir = abs(fdir0 - sdir0);
					if ( difdir >= 120 ) difdir = 240 - difdir;
					difdir1 = abs(fdir0 - sdir1);
					if ( difdir1 >= 120 ) difdir1 = 240 - difdir1;
					if ( difdir > AngTh && difdir1 > AngTh) continue;
				}
				if ( nFCoreNum == 2 ){
					difdir = abs(fdir0 - sdir0);
					if ( difdir >= 120 ) difdir = 240 - difdir;
					difdir1 = abs(fdir0 - sdir1);
					if ( difdir1 >= 120 ) difdir1 = 240 - difdir1;
					if ( difdir > AngTh && difdir1 > AngTh) continue;
					difdir = abs(fdir1 - sdir0);
					if ( difdir >= 120 ) difdir = 240 - difdir;
					difdir1 = abs(fdir1 - sdir1);
					if ( difdir1 >= 120 ) difdir1 = 240 - difdir1;
					if ( difdir > AngTh && difdir1 > AngTh) continue;
				}
			}
			if ( coreflag == 1 ){
				slen = op_func_02((pSearch->Mp.item[j].x-scx)*(pSearch->Mp.item[j].x-scx) +
						(pSearch->Mp.item[j].y-scy)*(pSearch->Mp.item[j].y-scy));
				if ( abs(flen - slen) > nTh ) continue;
			}

			dx = pSearch->Mp.item[j].x - cx; 
			dy = pSearch->Mp.item[j].y - cy;
			transform_mp(&pFile->Mp,cx,cy,rot,dx,dy);
			score = get_score_sub(&pFile->Mp,&pSearch->Mp);
			if ( maxscore < score ) maxscore = score;
			if ( maxscore > 1700 ) return (retScore);
			pFile->Mp = tmpMp;
		}
		if ( retScore < maxscore ) retScore = maxscore;
	}
	return (retScore);
}

BOOL check_point_kind(LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int i, fid, sid, nNum = 0;

	if ( pPair->nNumber == 0 ) return FALSE;
	for ( i=0; i<pPair->nNumber; i++ ){
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		if ( pFile->item[fid].kind == pSearch->item[sid].kind ) nNum++;
	}
	if ( nNum == pPair->nNumber ) return TRUE;
	return FALSE;
}

int point_matching(LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair,
					BOOL nCommonFlag,int nScoreTh,BOOL nSortFlag,int *nRetScore) 
{
	BARVECT FileBar, SearchBar;
	MATCH_TAG FBar, SBar;
	PAIRBAR PairList[1000], PairTemp[256];
	int AngleField[240], SDiffField[240], FDiffField[240];
	int XField[MAX_SHIFT_X_SIZE], YField[MAX_SHIFT_Y_SIZE], LastList[500];
	int SArrangBarPtr[240][20], FArrangBarPtr[240][20];
	int nPairNum, nLastNum, nMaxSearchBarLen;
	int nSid1, nSid2, nFid1, nFid2;
	int nScurv1, nScurv2, nFcurv1, nFcurv2, curv1, curv2;
	int nFileCX = 0, nFileCY = 0, nRot, nXoffset, nYoffset, nRidge;
	int score_sum, score, score_max, nBlockScore, score_old;
	int val, nSdir, nFdir, dslope, nCheck;
	int i, j, h, k, nUpper, nLower, id, diff1, diff2, lendiff;
	int nNum, snum, fnum, nSid, nFid, id1, id2;
	int newCX, newCY, nMatchScoreTh = nScoreTh;
	BOOL nRepeat, hflag, sflag, fflag, flag;
	BOOL nReCommonFlag = FALSE;
	int  nTh = ((ANGLE+LENGTH+3)*600)/1000 ;
	FPVECTEX saveFile = *pFile, saveSearch = *pSearch;

	get_search_tag(pSearch,&SearchBar,&nMaxSearchBarLen,SDiffField,SArrangBarPtr,20,200);
	if ( SearchBar.nNumber <= 0 ) return (0);

	arrange_points(pFile,pSearch,&SearchBar,&nMaxSearchBarLen,SArrangBarPtr,SDiffField,nCommonFlag);

	get_file_tag(nMaxSearchBarLen,pFile,&FileBar,FDiffField,FArrangBarPtr,&nFileCX,&nFileCY,20);

	score_sum = 0; nPairNum = 0;

	memset( AngleField, 0, sizeof(int)*240 );

	for ( i=0; i<240; i++ ){
		nUpper = i + ANGLE; nLower = i - ANGLE;
		for ( j=0; j<SDiffField[i]; j++ ){
			score = 0; score_max = 0;
			SBar = SearchBar.item[SArrangBarPtr[i][j]];
			id1 = SArrangBarPtr[i][j];
			nSid1 = SBar.nID1; nSid2 = SBar.nID2;
			for ( h=nLower; h<nUpper; h++ ){
				id = h - 240;
				if ( id < -240 ) id += 480;
				else{
					if ( id < 0 ) id += 240;
				}
				for ( k=0; k<FDiffField[id]; k++ ){
					FBar = FileBar.item[FArrangBarPtr[id][k]];
					id2 = FArrangBarPtr[id][k];
					nFid1 = FBar.nID1; nFid2 = FBar.nID2;
					diff1 = abs(FBar.nDiff1 - SBar.nDiff1);
					if ( diff1 >= 120 ) diff1 = 240 - diff1;
					if ( diff1 >= ANGLE ) continue;
					lendiff = abs(FBar.nLen - SBar.nLen);
					if ( lendiff >= LENGTH ) continue;
					diff2 = abs(FBar.nDiff2 - SBar.nDiff2);
					if ( diff2 >= 120 ) diff2 = 240 - diff2;
					if ( diff2 >= ANGLE ) continue;

					val = op_func_02(diff1*diff1+lendiff*lendiff+diff2*diff2);
					score = nTh - val;
					if ( score < 0 ) continue;
					nScurv1 = pSearch->Mp.item[nSid1].curv;
					nScurv2 = pSearch->Mp.item[nSid2].curv;
					nFcurv1 = pFile->Mp.item[nFid1].curv;
					nFcurv2 = pFile->Mp.item[nFid2].curv;
					curv1 = 30 - abs(nScurv1-nFcurv1);
					curv2 = 30 - abs(nScurv2-nFcurv2);
					if ( curv1 < 0 ) continue;
					if ( curv2 < 0 ) continue;
					score = (score*curv1*curv2)/900;
					if ( pFile->Mp.item[nFid1].kind != pSearch->Mp.item[nSid1].kind ) score = score*9/10;
					if ( pFile->Mp.item[nFid2].kind != pSearch->Mp.item[nSid2].kind ) score = score*9/10;
					if ( score > 0 ){
						PairList[nPairNum].score = score;
						PairList[nPairNum].fid = id2;
						PairList[nPairNum].sid = id1;
						if ( ++nPairNum == 1000 ){ nPairNum--; }

						nSdir = pSearch->Mp.item[nSid1].dir;
						nFdir = pFile->Mp.item[nFid1].dir;
						diff1 = nSdir - nFdir;
						if ( diff1 < 0 ) diff1 += 240;
						if ( diff1 >= 240 ) diff1 -= 240;

						nSdir = pSearch->Mp.item[nSid2].dir;
						nFdir = pFile->Mp.item[nFid2].dir;
						diff2 = nSdir - nFdir;
						if ( diff2 < 0 ) diff2 += 240;
						if ( diff2 >= 240 ) diff2 -= 240;

						AngleField[diff1] += score;
						AngleField[diff2] += score;
					}
					if ( score > score_max ) score_max = score;
				}
			}
			score_sum += score_max;
		}
	}
	if ( nCommonFlag == TRUE && 2*score_sum < nPairNum*5 && score_sum < nMatchScoreTh*3 ) return (0);
	if (SearchBar.nNumber > 100) 
		score_sum = (score_sum * 200) / SearchBar.nNumber;
	else 
		score_sum = score_sum * 2;
	score_sum = (score_sum*929 + 1137 + 500) / 1000;
	if ( score_sum < nMatchScoreTh ) return (0);

	nRot = rotate_points(nFileCX,nFileCY,AngleField,&FileBar,pFile);

	memset( XField, 0, sizeof(int)*MAX_SHIFT_X_SIZE );
	memset( YField, 0, sizeof(int)*MAX_SHIFT_Y_SIZE );
	nLastNum = 0; score_sum = score_max = 0;
	id = PairList[0].sid;
	for ( i=0; i<nPairNum; i++ ){
		score = PairList[i].score;
		FBar = FileBar.item[PairList[i].fid];
		SBar = SearchBar.item[PairList[i].sid];
		dslope = abs(FBar.nSlope - SBar.nSlope);
		if ( dslope >= 60 ) dslope = 120 - dslope;
		if ( dslope >= ANGLE ) continue;
		get_shift_param(LENGTH,score,&FBar,&SBar,XField,YField,&pFile->Mp,&pSearch->Mp);
		LastList[nLastNum++] = i;
		if ( nLastNum == 500 ){ break; }
		if ( id != PairList[i].sid ){
			score_sum += score_max; score_max = 0;
			id = PairList[i].sid;
		}
		if ( score_max < PairList[i].score ) score_max = PairList[i].score;
	}
	score_sum += score_max;
	if (SearchBar.nNumber > 100) 
		score_sum = (score_sum * 200) / SearchBar.nNumber;
	else 
		score_sum = score_sum * 2;
	score_sum = (score_sum*929 + 1137 + 500) / 1000;
	if ( score_sum < nMatchScoreTh ) return (0);

	shift_points(&nXoffset,&nYoffset,pFile,XField,YField);
	transform_block(nRot,nXoffset,nYoffset,nFileCX,nFileCY,&pFile->Block);
	nBlockScore = check_block(30,6,&pFile->Block,&pSearch->Block);

	if ( re_arrange_point(PairList,LastList,nPairNum,&nLastNum,pFile,pSearch,&FileBar,&SearchBar) == FALSE ){
		return (0);
	}

	newCX = 0; newCY = 0;
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		newCX += pFile->Mp.item[i].x; newCY += pFile->Mp.item[i].y;
	}
	newCX /= pFile->Mp.nNumber; newCY /= pFile->Mp.nNumber;

	score_sum = 0;
	for ( i=0,j=0; i<nLastNum; i++ ){
		FBar = FileBar.item[PairList[LastList[i]].fid];
		SBar = SearchBar.item[PairList[LastList[i]].sid];
		nCheck = check_limit(LENGTH,&FBar,&SBar,pFile,pSearch,newCX,newCY);
		if ( nCheck == FALSE ) continue;
		LastList[j++] = LastList[i];
	}
	if(j <= 0) return (0);

	nLastNum = j;
	for ( i=0; i<nLastNum; i++ ){
		if ( LastList[i] == -1 ) continue;
		score_max = PairList[LastList[i]].score;
		PairTemp[0] = PairList[LastList[i]];
		nNum = 1;
		do {
			nRepeat = FALSE;
			for ( j=0; j<nLastNum; j++ ){
				if ( i == j ) continue;
				if ( LastList[j] == -1 ) continue;
				nSid = PairList[LastList[j]].sid;
				nFid = PairList[LastList[j]].fid;
				for ( k=0; k<nNum; k++ ){
					if ( PairTemp[k].sid == nSid ) break;
					if ( PairTemp[k].fid == nFid ) break;
				}
				if ( k >= nNum ) continue;
				nRepeat = TRUE;
				PairTemp[nNum++] = PairList[LastList[j]];
				if ( nNum >= 255 ){ nRepeat = FALSE; break; }
				if ( score_max < PairList[LastList[j]].score )
					score_max = PairList[LastList[j]].score;
				LastList[j] = -1;
			}
		} while ( nRepeat == TRUE );
		if ( nNum <= 3 ){ score_sum += score_max; continue; }
		else{
			snum = fnum = 0;
			for ( j=0; j<nNum; j++ ){
				sflag = fflag = TRUE;
				nSid = PairTemp[j].sid; nFid = PairTemp[j].fid;
				for ( k=0; k<j; k++ ){
					if ( PairTemp[k].sid == nSid ) sflag = FALSE;
					if ( PairTemp[k].fid == nFid ) fflag = FALSE;
				}
				if ( sflag == TRUE ) snum++;
				if ( fflag == TRUE ) fnum++;
			}
			if ( snum > fnum ){ flag = TRUE; snum = fnum; }
			else flag = FALSE;
			for ( j=0; j<snum; j++ ){
				hflag = TRUE;
				for ( k=0; k<nNum; k++ ){
					if ( PairTemp[k].score == -1 ) continue;
					if ( hflag == TRUE ){
						if ( flag == TRUE ){
							id1 = PairTemp[k].fid;
							score_max = PairTemp[k].score; hflag = FALSE;
						}
						else{
							id1 = PairTemp[k].sid;
							score_max = PairTemp[k].score; hflag = FALSE;
						}
					}
					else{
						if ( flag == TRUE ) id2 = PairTemp[k].fid;
						else id2 = PairTemp[k].sid;
						if ( id1 != id2 ) continue;
						if ( score_max < PairTemp[k].score ) 
							score_max = PairTemp[k].score;
					}
					PairTemp[k].score = -1;
				}
				score_sum += score_max;
			}
		}
	}

	if ( pPair != NULL ){
		id = 0;
		for ( i=0; i<nLastNum; i++ ){
			if ( LastList[i] == -1 ) continue;
			FBar = FileBar.item[PairList[LastList[i]].fid];
			nFid = FBar.nID1;
			SBar = SearchBar.item[PairList[LastList[i]].sid];
			nSid = SBar.nID1;
			for ( j=0; j<id; j++ ){
				if ( nFid == pPair->nFileID[j] ) break;
				if ( nSid == pPair->nSearchID[j] ) break;
			}
			if ( j >= id ){
				pPair->nFileID[id] = nFid; 
				pPair->nSearchID[id++] = nSid;
			}
			nFid = FBar.nID2; nSid = SBar.nID2;
			for ( j=0; j<id; j++ ){
				if ( nFid == pPair->nFileID[j] ) break;
				if ( nSid == pPair->nSearchID[j] ) break;
			}
			if ( j >= id ){
				pPair->nFileID[id] = nFid; 
				pPair->nSearchID[id++] = nSid;
			}
		}
		pPair->nNumber = id; pPair->nRot = nRot; 
		pPair->nXOffset = nXoffset; pPair->nYOffset = nYoffset;
		pPair->nXC = nFileCX; pPair->nYC = nFileCY;
	}

	if ( score_sum < nMatchScoreTh/2 ) return (score_sum);

	nNum = pFile->Mp.nNumber;
	if ( nNum > pSearch->Mp.nNumber ) nNum = pSearch->Mp.nNumber;
	nRidge = 255 - abs(pSearch->nRidgeDensity - pFile->nRidgeDensity);
	dslope = 0;
	if ( pFile->nFrequency != 0 && pSearch->nFrequency != 0 ){
		fnum = snum = pFile->nFrequency;
		if ( pSearch->nFrequency > fnum ) fnum = pSearch->nFrequency;
		if ( pSearch->nFrequency < snum ) snum = pSearch->nFrequency;
		dslope = fnum - snum;
	}

	if ( *nRetScore > 0 ) score = *nRetScore;
	else {
		score = get_point_score(&saveFile,&saveSearch);
		*nRetScore = score;
		if ( score > 1700 ) return (nMatchScoreTh*2);
	}

	score_old = score_sum;
	if ( SearchBar.nNumber < 100 ){
		if ( pPair->nNumber < 5 && pPair->nNumber*3 <= nNum*2 && score < 600 && nBlockScore < 92 )
			score_sum = (score_sum * 15) / 10;
		else if ( nBlockScore < 85 && pPair->nNumber <= 6 || score < 270 || ( nBlockScore <= 90 && nRidge < 230 && score < 400 ) ) 
			score_sum = (score_sum * 11) / 10;
		else if ( nBlockScore < 93 && pPair->nNumber <= 6 && score < 850 && SearchBar.nNumber < 85 ) 
			score_sum = (score_sum * 15) / 10;
		else if ( nBlockScore < 95 && pPair->nNumber <= 5 && score < 500 ) 
			score_sum = (score_sum * 15) / 10;
		else
			score_sum = score_sum * 2;
	}
	else {
		if ( score > 1210 && SearchBar.nNumber >= 250 ) score_sum = (score_sum * 280) / SearchBar.nNumber;
		else if ( score > 1000 && SearchBar.nNumber > 250 && nRidge >= 230 ) score_sum = (score_sum * 240) / SearchBar.nNumber;
		else if ( score < 350 && SearchBar.nNumber < 120 && nBlockScore < 85) score_sum = (score_sum * 130) / SearchBar.nNumber;
		else if ( (score < 700 && SearchBar.nNumber < 180) || nBlockScore <= 80 ) score_sum = (score_sum * 150) / SearchBar.nNumber;
		else score_sum = (score_sum * 200) / SearchBar.nNumber;
	}
	score_sum = (score_sum*929 + 1137 + 500) / 1000;

	score_sum = (score_sum*nRidge) / 255;
	if ( score < 1175 && nRidge+nBlockScore<334 && nBlockScore<95 && (pPair->nNumber<10 || (score_sum<pPair->nNumber*10 && pPair->nNumber<13))){
		score_sum = (score_sum*nRidge)/255;
	}

	if ( dslope > 21 && score < 1200 ){
		score_sum = score_sum*snum/fnum;
		if ( (pFile->nType == 7 && (pSearch->nType==4 || pSearch->nType==5))
		  || (pSearch->nType == 7 && (pFile->nType==4 || pFile->nType==5)) ){
			score_sum = score_sum*snum/fnum;
		}
	}

	if ( pPair->nNumber >= 17 && pPair->nNumber*100 >= 40*nNum && (nBlockScore > 88 || nRidge >= 235) ){
		return (score_old);
	}
	if ( pPair->nNumber >= 18 && pPair->nNumber*100 >= 40*nNum && nBlockScore >= 85 ){
		return (score_old);
	}
	if ( score_sum > nMatchScoreTh && nBlockScore >= 83 ){
		if ( score > 1050 && pPair->nNumber >= 11 && pPair->nNumber*100 > nNum*50 )
			score_sum = score_sum*2;
		if ( nRidge >= 250 && score > 600 ){
			if ( nBlockScore > 90 && pPair->nNumber > 11 && pPair->nNumber*100 >= nNum*55 ){
				return (score_sum);
			}
			if ( nBlockScore > 86 && pPair->nNumber >= 9 && pPair->nNumber*100 > nNum*55 ){
				if ( check_core(&pFile->Core,&pSearch->Core,8,10) == TRUE )
					return (score_sum);
			}
			if ( pPair->nNumber >= 13 && pPair->nNumber*100 > nNum*40 && score > 1000){
				if ( check_core(&pFile->Core,&pSearch->Core,10,10) == TRUE )	return (score_sum);
			}
			if ( pPair->nNumber >= 10 && pPair->nNumber*100 >= nNum*45 && score >= 750 ){
				if ( check_point_kind(&pFile->Mp,&pSearch->Mp,pPair) == TRUE ){
					score_sum = score_old*2;
				}
			}
		}
		if ( pPair->nNumber > 13 && pPair->nNumber*100 >= 56*nNum && nRidge > 241 ){
			return (score_sum);
		}
		if ( pPair->nNumber >= 15 && pPair->nNumber*100 >= 65*nNum && score > 1150 ){
			return (score_sum);
		}
		if ( score_old*5 >= nMatchScoreTh*12 && pPair->nNumber >= 10 && pPair->nNumber*100 > 50*nNum ){
			if ( (nRidge >= 243 && nBlockScore > 88) && dslope < 8 ){
				if ( pPair->nNumber > 14 ) return (score_sum);
				if ( (check_core(&pFile->Core,&pSearch->Core,10,12) == TRUE) && score > 1000 ){
					if ( (pFile->nType < 4 && pSearch->nType < 4) || (pFile->nType >= 4 && pSearch->nType >= 4) ) 
						return (score_sum);
				}
			}
		}
		if ( pPair->nNumber >= 13 && pPair->nNumber*100 >= 45*nNum && nRidge >= 250 ){
			if ( pPair->nNumber >= 15 && nBlockScore > 90 ) return (score_sum);
			if ( nBlockScore > 91 ){
				if ( ((pFile->nType <= 1 && (pSearch->nType==4 || pSearch->nType==5))
				  || (pSearch->nType <= 1 && (pFile->nType==4 || pFile->nType==5)))
				  && (pFile->Core.nNumber != pSearch->Core.nNumber) )
					score_sum = (score_sum * nRidge)/255;
				else
					return (score_sum);
			}
		}
		if ( pPair->nNumber >= 10 && pPair->nNumber*100 >= nNum*40 && score >= 1000 ){
			if ( check_point_kind(&pFile->Mp,&pSearch->Mp,pPair) == TRUE ){
				score_sum = score_old*2;
			}
		}
	}
	score_sum = score_sum*nBlockScore*nBlockScore/10000;

	if ( score_old >= 100 && SearchBar.nNumber > 280 && score >= 1100 ){
		if ( nRidge > 245 && nBlockScore >= 91 && pPair->nNumber >= 10 ){
			if ( check_core(&pFile->Core,&pSearch->Core,9,7) == TRUE ){
				score_sum = score_old;
			}
		}
	}

	sflag = ((TRUE == check_core(&pFile->Core,&pSearch->Core,16,7)) && nRidge+nBlockScore < 330 && score < 1250);
	fflag = ( score < 1100 || pPair->nNumber <= 10 || pPair->nNumber*100 < nNum*40 );
	hflag = ( score < 1165 || pPair->nNumber <= 10 || pPair->nNumber*100 < nNum*38 );

	if ( score_sum < nMatchScoreTh*4 ){
		score_sum = dec_func_01(score_sum,pFile,pSearch,pPair);
		score_sum = dec_func_02(score_sum,&pFile->Mp,&pSearch->Mp,pPair);
		if ( fflag || sflag )
		score_sum = dec_func_03(score_sum,pFile,pSearch,nRidge,nBlockScore);
		score_sum = dec_func_04(score_sum,pFile,pSearch,pPair);
		score_sum = dec_func_05(score_sum,pFile,pSearch,pPair,
										nMatchScoreTh,nRidge,nBlockScore,score_old);
	}
	if ( score_sum < 10*nMatchScoreTh/5 ){
		score_sum = dec_func_06(score_sum,&pFile->Mp,&pSearch->Mp,pPair);
		if ( (hflag || sflag) && score < 1500 )
		score_sum = dec_func_07(score_sum,pFile,pSearch,pPair,nBlockScore,nRidge,score);
		score_sum = dec_func_08(score_sum,&pFile->Mp,&pSearch->Mp,pPair);
		score_sum = dec_func_09(score_sum,pFile,pSearch,pPair);
		if ( fflag && score < 1400 )
		score_sum = dec_func_10(score_sum,pFile,pSearch,pPair,nBlockScore,nRidge,score);
		score_sum = dec_func_11(score_sum,pFile,pSearch,pPair,nMatchScoreTh,nRidge,nBlockScore);
		score_sum = dec_func_12(score_sum,pFile,pSearch,pPair,nMatchScoreTh,nBlockScore);
		if ( SearchBar.nNumber < 160 && pPair->nNumber <= 7 && score_sum >= nMatchScoreTh && score_sum < nMatchScoreTh*16/10 ){
			score_sum = dec_func_13(score_sum,pFile,pSearch,pPair);
		}
	}
	return score_sum;
}

void get_matched_points_number(LPMPVECTEX pVect1,LPMPVECTEX pVect2,int *nNum1,int *nNum2)
{
	int i, j, dx, dy, ang, num1 = 0, num2 = 0;
	char temp1[MAX_MINUTIA_NUMBER], temp2[MAX_MINUTIA_NUMBER];
	BOOL flag1, flag2;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ){ temp1[i] = temp2[i] = 0; }
	*nNum1 = *nNum2 = 0;
	for ( i=0; i<pVect1->nNumber; i++ ){
		flag1 = flag2 = FALSE;
		for ( j=0; j<pVect2->nNumber; j++ ){
			dx = pVect1->item[i].x - pVect2->item[j].x;
			dy = pVect1->item[i].y - pVect2->item[j].y;
			dx = dx*dx + dy*dy;
			if ( dx > 12*12 ) continue;
			ang = abs(pVect1->item[i].dir - pVect2->item[j].dir);
			if ( ang > 120 ) ang = 240 - ang;
			if ( ang <= 7 ){ temp1[j] = 1; flag1 = TRUE; }
			if ( ang <= ANGLE ){ temp2[j] = 1; flag2 = TRUE; }
		}
		if( flag1 == TRUE ) num1++;
		if( flag2 == TRUE ) num2++;
	}
	dx = dy = 0;
	for ( i=0; i<pVect2->nNumber; i++ ){
		if ( temp1[i] == 1 ) dx++;
		if ( temp2[i] == 1 ) dy++;
	}
	if ( num1 > dx ) num1 = dx;
	if ( num2 > dy ) num2 = dy;
	*nNum1 = num1; *nNum2 = num2;
}

int coarse_matching(LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	FPVECTEX pVect;
	int i, j, cx, cy, dx, dy, rot, nNum1 = 0, nNum2 = 0, nTotalNum = 0, nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);

	if ( nFCoreNum == 0 || nSCoreNum == 0 ) return (0);
	if ( pFile->Mp.nNumber == 0 || pSearch->Mp.nNumber == 0 ) return (-1);

	for ( i=0; i<nFCoreNum; i++ ){
		cx = FileCore[i].x; cy = FileCore[i].y;
		for ( j=0; j<nSCoreNum; j++ ){
			dx = SearchCore[j].x - cx;
			dy = SearchCore[j].y - cy;
			rot = SearchCore[j].dir - FileCore[i].dir;
			if ( rot < 0 ) rot += 240;
			if ( rot >= 240 ) rot -= 240;
			pVect = *pFile;
			transform_points(&pVect,cx,cy,rot,dx,dy);
			nTotalNum = get_min_points_number(&pVect.Mp,&pSearch->Mp);
			get_matched_points_number(&pVect.Mp,&pSearch->Mp,&nNum1,&nNum2);
			if ( nNum1 > 6 && nNum1*100 > 80*nTotalNum && nTotalNum >= 9) return (1);
			if ( nNum2 > 12 && nNum2*100 > 70*nTotalNum ) return (1);
		}
	}
	if ( nFCoreNum == nSCoreNum && nFCoreNum == 1 ){
		if ( nNum1 > 6 && nNum2 >= 10 && nNum2*100 > 40*nTotalNum ){
			transform_block(rot,dx,dy,cx,cy,&pVect.Block);
			nNum1 = check_block(30,5,&pVect.Block,&pSearch->Block);
			if ( nNum1 >= 90 ) return (2);
		}
	}
	return (0);
}

int get_distance(LPFPVECTEX pVect,int CoreNumber,int MinMaxFlag)
{
	int x, y, dx, dy, nn, len1, len2;
	COREITEMEX FileCore[2];

	mch_sub_func_01(&(pVect->Core),FileCore,NULL,NULL);
	
	if ( CoreNumber < 0 || CoreNumber > 1 ) return (0);
	nn = CoreNumber*2;
	x = FileCore[CoreNumber].x;
	y = FileCore[CoreNumber].y;
	dx = x - pVect->MainLine.points_x[nn];
	dy = y - pVect->MainLine.points_y[nn];
	len1 = dx*dx + dy*dy;
	nn++;
	dx = x - pVect->MainLine.points_x[nn];
	dy = y - pVect->MainLine.points_y[nn];
	len2 = dx*dx + dy*dy;

	if ( MinMaxFlag == 0 ){ 
		nn = len1;
		if ( nn > len2 ) nn = len2;
	}
	else{
		nn = len1;
		if ( nn < len2 ) nn = len2;
	}
	nn = op_func_02(nn);
	return (nn);
}

BOOL type_matching(LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	int i, t1, t2, i1, i2, dx ,dy, n1, pScore1, pScore2;
	int nFCoreNum, nSCoreNum, nFDeltaNum, nSDeltaNum;
	COREITEMEX FileCore[2], SearchCore[2];
	COREITEMEX FileDelta[2], SearchDelta[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,FileDelta,&nFDeltaNum);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,SearchDelta,&nSDeltaNum);

	if ( pFile->nType < pSearch->nType ){
		t1 = pFile->nType; t2 = pSearch->nType;
	}
	else{
		t1 = pSearch->nType; t2 = pFile->nType;
	}
	if ( t1 == t2 ) return TRUE;
	if ( pFile->Mp.nNumber==0 || pSearch->Mp.nNumber==0 ) return FALSE;
	pScore1 = pScore2 = 0;
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		pScore1 += pFile->Mp.item[i].score;
	}
	pScore1 /= pFile->Mp.nNumber;
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		pScore2 += pSearch->Mp.item[i].score;
	}
	pScore2 /= pSearch->Mp.nNumber;

	i1 = i2 = pScore1;
	if ( i1 > pScore2 ) i1 = pScore2;
	if ( i2 < pScore2 ) i2 = pScore2;
	if ( i1 <= 35 && i2 <= 41 ) return TRUE;

	if ( t1 < 3 ){
		if ( t2 == 3 ) return TRUE;
		if ( t2 < 3 ){
			i1 = i2 = 0;
			if ( nFCoreNum == 2 ){
				dx = FileCore[0].x - FileCore[1].x;
				dy = FileCore[0].y - FileCore[1].y;
				i1 = op_func_02(dx*dx+dy*dy);
			}
			if ( nSCoreNum == 2 ){
				dx = SearchCore[0].x - SearchCore[1].x;
				dy = SearchCore[0].y - SearchCore[1].y;
				i2 = op_func_02(dx*dx+dy*dy);
			}
			if ( abs(i1-i2) <= 55 ) return TRUE;
			if ( pFile->nType == 1 && nSCoreNum == 1) return TRUE;
			if ( pSearch->nType == 1 && nFCoreNum == 1) return TRUE;
			return FALSE ;
		}
		if ( t1 == 2 && (t2 == 4 || t2 == 5) ) return TRUE;
		if ( t2 > 3 && t2 < 8 ){
			if ( t1 == 0 && (t2 == 4 || t2 == 5) ){
				if ( pFile->nType == 0 ){
					dx = pSearch->MainLine.points_x[0] - pSearch->MainLine.points_x[1];
					dy = pSearch->MainLine.points_y[0] - pSearch->MainLine.points_y[1];
				}
				else{
					dx = pFile->MainLine.points_x[0] - pFile->MainLine.points_x[1];
					dy = pFile->MainLine.points_y[0] - pFile->MainLine.points_y[1];
				}
				if ( dx*dx+dy*dy <= 40*40 ) return TRUE;
			}
			if(t1 == 1){
				if(pFile->nType == 1) i2 = get_distance(pSearch,0,0);
				else i2 = get_distance(pFile,0,0);
				if ( i2 > 70 ) return FALSE;
				return TRUE;
			}
			return FALSE;
		}
		if ( t1 != 2 && t2 == 8 ){
			if ( pScore2 > 35 && nFCoreNum == 2 && nSCoreNum == 1 ){
				dx = FileCore[0].x - FileCore[1].x;
				dy = FileCore[0].y - FileCore[1].y;
				i1 = op_func_02(dx*dx+dy*dy);
				i2 = get_distance(pSearch,0,0);
				if ( t1 == 1 ){
					if ( i1 < i2-40 ) return FALSE;
					return TRUE;
				}
				if ( i1 < i2-25 ) return FALSE;
				if ( i1 < i2-5 ){
					dx = pSearch->MainLine.points_x[0] - pSearch->MainLine.points_x[1];
					dy = pSearch->MainLine.points_y[0] - pSearch->MainLine.points_y[1];
					if ( dx*dx+dy*dy > 40*40 ) return FALSE;
				}
			}
			if ( pScore1 > 35 && nFCoreNum == 1 && nSCoreNum == 2 ){
				dx = SearchCore[0].x - SearchCore[1].x;
				dy = SearchCore[0].y - SearchCore[1].y;
				i1 = op_func_02(dx*dx+dy*dy);
				i2 = get_distance(pFile,0,0);
				if ( t1 == 1 ){
					if ( i1 < i2-40 ) return FALSE;
					return TRUE;
				}
				if ( i1 < i2-25 ) return FALSE;
				if ( i1 < i2-5 ){
					dx = pFile->MainLine.points_x[0] - pFile->MainLine.points_x[1];
					dy = pFile->MainLine.points_y[0] - pFile->MainLine.points_y[1];
					if ( dx*dx+dy*dy > 40*40 ) return FALSE;
				}
			}
		}
		return TRUE;
	}
	if ( t1 == 3 ){
		if ( t2 > 8 || pScore1 <= 35 || pScore2 <= 35 ) return TRUE;
		if ( nFCoreNum == 1 && nSCoreNum == 1 ){
			i1 = get_distance(pFile,0,0);
			i2 = get_distance(pSearch,0,0);
			if ( pSearch->nType == 3 ){
				n1 = i1; i1 = i2; i2 = n1;
			}
			if ( i1 > i2-30 ) return TRUE;
		}
	}
	if ( t1 == 4 ){
		if ( t2 == 5 || t2 == 6 ) return FALSE;
	}
	if ( t1 == 5 && t2 == 6 ) return FALSE;
	if ( t1 == 6 && t2 == 7 ) return FALSE;
	if ( t1 == 4 || t1 == 5 ){
		if ( t2 == 7 ){
			if ( pFile->nType == 7 && nFDeltaNum == 1 ){
				i = op_func_01(FileDelta[0].x,FileDelta[0].y,FileCore[0].x,FileCore[0].y);
				i1 = i - FileCore[0].dir;
				if ( i1 > 0 ){ 
					if ( pSearch->nType == 4 ) return FALSE;
				}
				else if ( i1 < 0  ){
					if ( pSearch->nType == 5 ) return FALSE;
				}
			}
			else if ( pSearch->nType == 7 && nSDeltaNum == 1 ){
				i = op_func_01(SearchDelta[0].x,SearchDelta[0].y,SearchCore[0].x,SearchCore[0].y);
				i2 = i - SearchCore[0].dir;
				if ( i2 > 0 ){ 
					if ( pFile->nType == 4 ) return FALSE;
				}
				else if ( i2 < 0 ){
					if ( pFile->nType == 5 ) return FALSE;
				}
			}
		}
	}
	if ( t1 == 7 && t2 == 7 && nFDeltaNum == 1 && nSDeltaNum == 1 ){
			i = op_func_01(FileDelta[0].x,FileDelta[0].y,FileCore[0].x,FileCore[0].y);
			i1 = i - FileCore[0].dir;
			i = op_func_01(SearchDelta[0].x,SearchDelta[0].y,SearchCore[0].x,SearchCore[0].y);
			i2 = i - SearchCore[0].dir;
			if ( i1 > 0 && i2 < 0 ) return FALSE;
			if ( i1 < 0 && i2 > 0 ) return FALSE;
	}
	return TRUE;
}

BOOL reget_points_sub(int nRot,int nXDiff,int nYDiff,LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	FPVECTEX tmpF = *pFile;
	int i, j, k, x, y, fx = 0, fy = 0;
	int nList[MAX_BLOCK_NUMBER], nNum, nRow, nCol, bx, by;
	int nNewId[MAX_MINUTIA_NUMBER];

	if ( pFile->Mp.nNumber == 0 || pSearch->Mp.nNumber == 0 ) return FALSE;
	for ( i=0; i<tmpF.Mp.nNumber; i++ ){
		fx += tmpF.Mp.item[i].x; fy += tmpF.Mp.item[i].y;
	}
	fx /= tmpF.Mp.nNumber; fy /= tmpF.Mp.nNumber;

	transform_block(nRot,nXDiff,nYDiff,fx,fy,&tmpF.Block);
	transform_mp(&tmpF.Mp,fx,fy,nRot,nXDiff,nYDiff);

	nNum = 0;
	for ( i=0; i<pFile->Block.nCol*pFile->Block.nRow; i++ ){
		if ( tmpF.Block.Data[i] == 0xFF ) continue;
		if ( pSearch->Block.Data[i] == 0xFF ) continue;
		nList[nNum++] = i;
	}
	if ( nNum == 0 ) return FALSE;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ) nNewId[i] = -1;
	for ( i=0,k=0; i<tmpF.Mp.nNumber; i++ ){
		x = tmpF.Mp.item[i].x; y = tmpF.Mp.item[i].y;
		for ( j=0; j<nNum; j++ ){
			nRow = nList[j] / MAX_BLOCK_COL; 
			nCol = nList[j] % MAX_BLOCK_COL;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2; 
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( abs(x-bx)>LENGTH || abs(y-by)>LENGTH ) continue;
			nNewId[i] = k;
			pFile->Mp.item[k++] = pFile->Mp.item[i];
			break;
		}
	}
	pFile->Mp.nNumber = k;

	for ( i=0; i<MAX_MINUTIA_NUMBER; i++ ) nNewId[i] = -1;
	for ( i=0,k=0; i<pSearch->Mp.nNumber; i++ ){
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		for ( j=0; j<nNum; j++ ){
			nRow = nList[j] / MAX_BLOCK_COL; 
			nCol = nList[j] % MAX_BLOCK_COL;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2; 
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if ( abs(x-bx)>LENGTH || abs(y-by)>LENGTH ) continue;
			nNewId[i] = k;
			pSearch->Mp.item[k++] = pSearch->Mp.item[i];
			break;
		}
	}
	pSearch->Mp.nNumber = k;
	return TRUE;
}

BOOL get_parameter(int *nRot,int *nXDiff,int *nYDiff,int *Field,UINT *CoorField,int num)
{
	int i, j, dir1, dx1, dy1, dir2, dx2, dy2, ddir, dx, dy, maxV, maxP;
	int angTH = 4, lenTH = 5;
	int tempArray[2000];
	UINT u1;

	memset(tempArray, 0, sizeof(int)*2000 );
	for ( i=0; i<num; i++ ){
		tempArray[i] += Field[i];
		u1 = CoorField[i];
		dy1 = u1 % 1000; u1 /= 1000;
		dx1 = u1 % 1000; dir1 = u1 / 1000;
		for ( j=i+1; j<num; j++ ){
			u1 = CoorField[j];
			dy2 = u1 % 1000; u1 /= 1000;
			dx2 = u1 % 1000; dir2 = u1 / 1000;
			ddir = abs(dir1 - dir2);
			if ( ddir >= 120 ) ddir = 240 - ddir;
			if ( ddir > angTH ) continue;
			dx = abs(dx1 - dx2);
			if ( dx > lenTH ) continue;
			dy = abs(dy1 - dy2);
			if ( dy > lenTH ) continue;
			tempArray[i] += Field[j];
			tempArray[j] += Field[i];
		}
	}
	maxV = 0;
	for ( i=0; i<num; i++ ){
		if ( maxV < tempArray[i] ){ maxV = tempArray[i]; maxP = i; }
	}
	if ( maxV == 0 ) return FALSE;
	u1 = CoorField[maxP];
	dy1 = u1 % 1000; u1 /= 1000;
	dx1 = u1 % 1000; dir1 = u1 / 1000;
	*nRot = dir1; *nXDiff = dx1 - 500; *nYDiff = dy1 - 500;
	return TRUE;
}

BOOL reget_points(LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	BLOCKVECT BlockF;
	UINT CoorField[2000];
	int Field[2000], nn = 0, nRot, nXDiff, nYDiff, nf, ns;
	int i, j, cx, cy, nCos, nSin, x, y, fx, fy, fdir, nBlockScore;
	int dcurv;

	if ( pFile->Mp.nNumber == 0 || pSearch->Mp.nNumber == 0 ) return FALSE;
	cx = cy = 0;
	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		cx += pFile->Mp.item[i].x; cy += pFile->Mp.item[i].y;
	}
	cx /= pFile->Mp.nNumber; cy /= pFile->Mp.nNumber;

	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		fx = pFile->Mp.item[i].x; fy = pFile->Mp.item[i].y;
		fdir = pFile->Mp.item[i].dir;
		for ( j=0; j<pSearch->Mp.nNumber; j++ ){
			nRot = pSearch->Mp.item[j].dir - fdir;
			if ( nRot < 0 ) nRot += 240;
			if ( nRot > 40 && nRot < 200 ) continue;
			dcurv = abs(pFile->Mp.item[i].curv-pSearch->Mp.item[j].curv);
			if ( dcurv > 10 ) continue;

			nCos = _table_03[nRot]; nSin = _table_04[nRot];
			x = (fx-cx)*nCos - (fy-cy)*nSin;
			x = x >> 14;
			y = (fy-cy)*nCos + (fx-cx)*nSin;
			y = y >> 14;
			nXDiff = pSearch->Mp.item[j].x - x - cx;
			nYDiff = pSearch->Mp.item[j].y - y - cy;
			if ( nXDiff >= 500 ) continue;
			if ( nYDiff >= 500 ) continue;

			BlockF = pFile->Block;
			transform_block(nRot,nXDiff,nYDiff,cx,cy,&BlockF);
			nBlockScore = check_block(30,6,&BlockF,&pSearch->Block);
			if ( nBlockScore < 70 ) continue;

			CoorField[nn] = (nRot*1000+nXDiff+500)*1000+nYDiff+500;
			Field[nn++] = nBlockScore;
			if ( nn >= 2000 ){ i = 1000; break; }
		}
	}
	if ( get_parameter(&nRot,&nXDiff,&nYDiff,Field,CoorField,nn) == FALSE ) return FALSE;
	if ( nRot > 40 && nRot < 200 ) return FALSE;
	nf = pFile->Mp.nNumber; ns = pSearch->Mp.nNumber;

	if ( reget_points_sub(nRot,nXDiff,nYDiff,pFile,pSearch) == FALSE ) return FALSE;

	fx = nf - pFile->Mp.nNumber; fy = ns - pSearch->Mp.nNumber;
	if ( (nf > 40 || ns > 40) && fx < 4 && fy < 4 ) return FALSE;
	if ( nf > 10 && ns > 10 && fx < 3 && fy < 3 ) return FALSE;
	return TRUE;
}

int matching_main(LPFPFEATURE pFeatureVect1,LPFPFEATURE pFeatureVect2,int securitylevel)
{
	FPVECTEX Vect1, Vect2, TmpVect1, TmpVect2;
	PAIRVECT pPair;
	int score = 0, nTh = MATCH_TH_MEDIUM;
	int nCoarse = 0, nGlobalScore = 0;

	if ( securitylevel == HIGH_LEVEL )	nTh = MATCH_TH_HIGH;
	if ( securitylevel == LOW_LEVEL )	nTh = MATCH_TH_LOW;

	mch_sub_func_02(pFeatureVect1,&Vect1);
	if ( mch_sub_func_03(&Vect1) == FALSE ) return (0);
	TmpVect1 = Vect1;
	mch_sub_func_02(pFeatureVect2,&Vect2);
	if ( mch_sub_func_03(&Vect2) == FALSE ) return (0);
	TmpVect2 = Vect2;

    ViewFPEx[0] = Vect1; ViewFPEx[1] = Vect2;

	nCoarse = coarse_matching(&Vect1,&Vect2);
	if ( nCoarse == -1 ) return (0);
	if ( nCoarse == 1 ) return (1000);

	if ( type_matching(&Vect1,&Vect2) == FALSE ) return (0);

	score = point_matching(&Vect1,&Vect2,&pPair,TRUE,nTh,FALSE,&nGlobalScore);
	if ( score < nTh ){
		if ( nCoarse == 2 && score >= nTh/2 ) return (2*score);
		Vect1 = TmpVect1; Vect2 = TmpVect2;
		nGlobalScore = 0;
		score = point_matching(&Vect2,&Vect1,&pPair,TRUE,nTh,FALSE,&nGlobalScore);
		if ( score < nTh ){
			if ( nCoarse == 2 && score >= nTh/2 ) return (2*score);
			if ( nGlobalScore > 900 ){
				Vect1 = TmpVect1; Vect2 = TmpVect2;
				if ( reget_points(&Vect1,&Vect2) == TRUE ){
					score = point_matching(&Vect1,&Vect2,&pPair,FALSE,nTh,FALSE,&nGlobalScore);
					if ( score >= nTh ) return (score);
				}
			}
			return (score);
		}
		return (score);
	}
	return (score);
}

void sch_sub_func_03(LPMPVECTEX pVect,int cx,int cy,int nAngle,int nDiffX,int nDiffY)
{
	int i, x, y, angle, rot, nCos, nSin, dx, dy;
	int nX = cx + nDiffX, nY = cy + nDiffY;

	rot = 240 - nAngle;
	if ( rot >= 240 ) rot -= 240;
	nCos = _table_03[rot]; nSin=_table_04[rot];

	for ( i=0; i<pVect->nNumber; i++ ){
		dx = pVect->item[i].x-cx; dy = pVect->item[i].y-cy;
		x = dx*nCos + dy*nSin;
		x = x >> 14;
		y = dy*nCos - dx*nSin;
		y = y >> 14;
		pVect->item[i].x = x + nX;
		pVect->item[i].y = y + nY;

		angle = pVect->item[i].dir + nAngle;
		if ( angle >= 240 ) angle -= 240;
		if ( angle < 0 ) angle += 240;
		pVect->item[i].dir = angle;
	}
}

int sch_sub_func_04(LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	int i, j, x, y, dir, dx, dy, diff, val, minval, score = 0;

	for ( i=0; i<pVect1->nNumber; i++ ){
		x = pVect1->item[i].x; y = pVect1->item[i].y;
		dir = pVect1->item[i].dir;
		minval = 10000;
		for ( j=0; j<pVect2->nNumber; j++ ){
			dx = abs(pVect2->item[j].x - x);
			if ( dx > LENGTH+3 ) continue;
			dy = abs(pVect2->item[j].y - y);
			if ( dy > LENGTH+3 ) continue;
			diff = abs(pVect2->item[j].dir - dir);
			if ( diff >= 120 ) diff = 240 - diff;
			if ( diff > ANGLE ) continue;
			val = dx + dy + diff;
			if ( minval > val ) minval = val;
			if ( minval < 20) break;
		}
		if ( minval < 35 ) score += (35 - minval);
	}
	val = (pVect1->nNumber + pVect2->nNumber) / 2;
	if ( val == 0 ) return (0);
	score = (score * 100 ) / val;
	return (score);
}

int sch_sub_func_05(LPFPVECTEX pFPEx0,LPFPVECTEX pFPEx1)
{
	MPVECTEX tmpMPEx0, tmpMPEx1;
	int i, j, k, nn, cx0, cy0, dir0, cx1, cy1, dir1, dx, dy, rot;
	int score = 0, radTh, nNumTh = 5;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFPEx0->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pFPEx1->Core),SearchCore,NULL,NULL);
	
	for ( i=0; i<nFCoreNum; i++ ){
		cx0 = FileCore[i].x; cy0 = FileCore[i].y;
		dir0 = FileCore[i].dir;
		radTh = 100; nn = 0;
		while (1){
			if ( radTh > 200 ) break;
			for ( k=0; k<pFPEx0->Mp.nNumber; k++ ){
				dx = pFPEx0->Mp.item[k].x - cx0;
				dy = pFPEx0->Mp.item[k].y - cy0;
				if ( dx*dx+dy*dy < radTh*radTh ){
					tmpMPEx0.item[nn] = pFPEx0->Mp.item[k];
					nn++;
				}
			}
			if ( nn < nNumTh ){ radTh += 20; continue; }
			tmpMPEx0.nNumber = nn; break;
		}
		for ( j=0; j<nSCoreNum; j++ ){
			cx1 = SearchCore[j].x; cy1 = SearchCore[j].y;
			dir1 = SearchCore[j].dir;

			rot = dir0 - dir1;
			if ( rot < 0 ) rot += 240;

			radTh = 100; nn = 0;
			while (1){
				if ( radTh > 200 ) break;
				for ( k=0; k<pFPEx1->Mp.nNumber; k++ ){
					dx = pFPEx1->Mp.item[k].x - cx1;
					dy = pFPEx1->Mp.item[k].y - cy1;
					if ( dx*dx+dy*dy < radTh*radTh ){
						tmpMPEx1.item[nn] = pFPEx1->Mp.item[k];
						nn++;
					}
				}
				if ( nn < nNumTh ){ radTh += 20; continue; }
				tmpMPEx1.nNumber = nn; break;
			}
			dx = cx0 - cx1; dy = cy0 - cy1;
			sch_sub_func_03(&tmpMPEx1,cx1,cy1,rot,dx,dy);
			nn = sch_sub_func_04(&tmpMPEx0,&tmpMPEx1);
			if ( score < nn ) score = nn;
		}
	}
	return(score);
}

int sch_sub_func_01(LPFPVECTEX pFile,LPFPVECTEX pSearch) 
{
	int i, j, k, cx, cy, dx, dy, rot;
	int retScore = 0, score, maxscore, nNumF = 0, nMaxF = 0, nMaxS = 0;
	MP_POINT tmpF[5];
	MPVECTEX tmpVect = pFile->Mp;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);

	if ( pFile->Mp.nNumber < MIN_MINUTIA_NUMBER || pSearch->Mp.nNumber < MIN_MINUTIA_NUMBER ) return (0);

	if ( nFCoreNum != 0 && nSCoreNum != 0 ){
		retScore = sch_sub_func_05(pSearch,pFile);
		return (retScore);
	}

	for ( i=0; i<pFile->Mp.nNumber; i++ ){
		if ( nMaxF < pFile->Mp.item[i].curv ) nMaxF = pFile->Mp.item[i].curv;
	}
	for ( i=0; i<pSearch->Mp.nNumber; i++ ){
		if ( nMaxS < pSearch->Mp.item[i].curv ) nMaxS = pSearch->Mp.item[i].curv;
	}
	if ( nMaxF > nMaxS ) nMaxF = nMaxS;

	for ( i=0,j=0; i<pFile->Mp.nNumber; i++ ){
		if ( pFile->Mp.item[i].curv >= nMaxF+8 ) continue;
		for ( k=0; k<j; k++ ){
			dx = pFile->Mp.item[i].x - tmpF[k].x;
			dy = pFile->Mp.item[i].y - tmpF[k].y;
			if ( dx*dx+dy*dy <= 40*40 ) break;
		}
		if ( k < j ) continue;
		tmpF[j++] = pFile->Mp.item[i];
		if ( j >= 5 ) break;
	}
	nNumF = j;
	pFile->Mp = tmpVect;
	for ( i=0; i<nNumF; i++ ){
		cx = tmpF[i].x; cy = tmpF[i].y;
		maxscore = 0;
		for ( j=0; j<pSearch->Mp.nNumber; j++ ){
			if ( abs(tmpF[i].curv-pSearch->Mp.item[j].curv) >= 6 ) continue;
			rot = pSearch->Mp.item[j].dir - tmpF[i].dir;
			if ( rot < 0 ) rot += 240;
			dx = pSearch->Mp.item[j].x - cx; dy = pSearch->Mp.item[j].y - cy;
			sch_sub_func_03(&pFile->Mp,cx,cy,rot,dx,dy);
			score = sch_sub_func_04(&pFile->Mp,&pSearch->Mp);
			if ( maxscore < score ) maxscore = score;
			pFile->Mp = tmpVect;
		}
		if ( retScore < maxscore ) retScore = maxscore;
	}
	return (retScore);
}

void sch_sub_func_02(int *pScore,int nSize,int *pIndex)
{
	int i, j, tmp, nMin;

    for ( i=0; i<nSize; i++ ) {
        pIndex[i] = i;
    }
	nMin = nSize - 1;
    if ( nMin > N_TEMPLATES_SCORED ) {
        nMin = N_TEMPLATES_SCORED;
    }
	for ( i=0; i<nMin; i++ ){
		for ( j=i+1; j<nSize; j++ ){
			if ( pScore[pIndex[i]] < pScore[pIndex[j]] ){
				tmp = pIndex[i]; pIndex[i] = pIndex[j]; pIndex[j] = tmp;
			}
		}
	}
}

int sch_sub_func(LPFPFEATURE pFeatureVect,LPFPFEATURE pDataBaseVects,int nDataBaseSize,int* pIndex)
{
	int i, *pScore;
	FPVECTEX FileFPEx, SearchFPEx;
	
	pScore = (int*)malloc(sizeof(int)*nDataBaseSize);
	if ( pScore == NULL ) return ERR_GENERAL_ERROR;

	mch_sub_func_02(pFeatureVect,&SearchFPEx);

	for ( i=0; i<nDataBaseSize; i++ ){
		mch_sub_func_02(&(pDataBaseVects[i]),&FileFPEx);
		pScore[i] = sch_sub_func_01(&FileFPEx,&SearchFPEx);
	}

	sch_sub_func_02(pScore,nDataBaseSize,pIndex);

	free(pScore);
	return ERR_OK;
}


// exported function //

/*
 *	fingerprint matching function.
 *  parameter ;
 *		pFeature1 : pointer to first fingerprint feature template buffer
 *		pFeature2 : pointer to second fingerprint feature template buffer
 *		securitylevel : value of security level (default - MEDIUM_LEVEL)
 *  return value ;
 *		the simility value (0~1000).
 *		if this value is greater than matching threshold value which is determined by
 *		security level, two fingerprint template is matched.
 *		if this value is less than matching threshold value, matching is failed.
 */
int finger_match(BYTE* pFeature1,BYTE* pFeature2,int securitylevel)
{
	int res;
	int nTh = MATCH_TH_MEDIUM;
	LPFPFEATURE pVect1 = (LPFPFEATURE)(pFeature1);
	LPFPFEATURE pVect2 = (LPFPFEATURE)(pFeature2);

	if ( securitylevel == HIGH_LEVEL )	nTh = MATCH_TH_HIGH;
	if ( securitylevel == LOW_LEVEL )	nTh = MATCH_TH_LOW;

	if ( pFeature1 == NULL || pFeature2 == NULL ) return (0);
	
	res = matching_main(pVect1,pVect2,securitylevel);
	if (res < 0) res = 0;
	if (res > 1000) res = 1000;

    /* Envia res diretamente (score)
	// Padroniza retorno em 1 (match) ou <0 (no match)
	if ( res >= nTh ) return ERR_OK;
	return ERR_MATCH_FAILED;
    */
    return res;
}

// Modificação: insere score como parâmetro de saída.
/*
 *	fingerprint searching function.
 *  parameter ;
 *		pFeature : pointer to inputed fingerprint feature template buffer
 *		pDBFeature : pointer to buffer of fingerprint features(size is 512 * nDBSize bytes)
 *		nDBSize : number of fingerprint feature in database
 *		securitylevel : value of security level (default - MEDIUM_LEVEL)
 *              pscore : score level of matched template (out parameter)
 *  return value ;
 *		if success, return the index value of matched fingerprint feature (0 ~ nDBSize).
 *		if failed, return error code ( < 0 ).
 */
int finger_search(BYTE* pFeature,BYTE* pDBFeature,int nDBSize,int securitylevel, int* pscore)
{
	int i, idx, res, nRealIdentNum;
	int *pIndexArray;
	LPFPFEATURE pVect = (LPFPFEATURE)(pFeature);
	LPFPFEATURE pDBVect = (LPFPFEATURE)(pDBFeature);
	int nTh = MATCH_TH_MEDIUM;

	if ( securitylevel == HIGH_LEVEL )	nTh = MATCH_TH_HIGH;
	if ( securitylevel == LOW_LEVEL )	nTh = MATCH_TH_LOW;

	if ( nDBSize < 1 ) return ERR_GENERAL_ERROR;
	pIndexArray = (int*)malloc(sizeof(int)*nDBSize);
	if ( pIndexArray == NULL ) return ERR_CAN_NOT_ALLOC_MEMORY;

	if ( nDBSize == 1 ){ pIndexArray[0] = 0; }
	else{
		if ( ERR_OK != sch_sub_func(pVect,pDBVect,nDBSize,pIndexArray) ){
			free(pIndexArray); return ERR_GENERAL_ERROR;
		}
	}

        // Não percorre o banco inteiro de templates,
        // apenas uma fração, para otimizar tempo de busca.
        nRealIdentNum = N_TEMPLATES_SCORED /*nDBSize*/;

        for ( i=0; i<nRealIdentNum; i++ ) {
            idx = pIndexArray[i];
            if ( idx < 0 || idx >= nDBSize ) continue;
            res = matching_main(&pDBVect[idx], pVect, securitylevel);
            if ( res >= nTh ) {
                free(pIndexArray);
                *pscore = res;
                return (idx);
            }
        }

        free(pIndexArray);
        *pscore = res;
        return (ERR_MATCH_FAILED);
}

/*
 *	feature extraction function.
 *  parameter ;
 *		pImage : pointer to fingerprint image_buffer1 buffer
 *		nWidth : width of fingerprint image_buffer1.
 *		nHeight : height of fingerprint image_buffer1.
 *		pFeature : pointer to fingerprint feature template buffer
 *  return value ;
 *		if success, return 1, otherwise return other value.
 */
int create_template(BYTE* pImage,int nWidth,int nHeight,BYTE* pFeature)
{
	BYTE *pTmpImg;
	int width, height;
	LPFPFEATURE pVect = (LPFPFEATURE)(pFeature);

	if ( pImage == NULL || pFeature == NULL ) return (ERR_GENERAL_ERROR);
	if ( nWidth < 0 || nWidth > MAX_IMG_WIDTH || nHeight < 0 || nHeight > MAX_IMG_HEIGHT ) return ERR_INVALID_IMAGESIZE;

	width = nWidth; height = nHeight;
	pTmpImg = (BYTE*)malloc(sizeof(BYTE)*width*height);
	if (pTmpImg == NULL) return ERR_CAN_NOT_ALLOC_MEMORY;

	memcpy(pTmpImg,pImage,sizeof(BYTE)*nWidth*nHeight);

	int ret = ext_main(pTmpImg,width,height,pVect);
	
	free(pTmpImg);
	return ret;
}
