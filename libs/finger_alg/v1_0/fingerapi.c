#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <fcntl.h>
#include "consttable.h"
#include "fingerapi.h"

int op_func_01(int ex,int ey,int sx,int sy)
{
	int nDiffX = abs(sx-ex), nDiffY = abs(sy-ey);
	int nAngle;

	while (nDiffX>=50 || nDiffY>=50) {
		nDiffX >>= 1; nDiffY >>= 1;
	}
	nAngle = _table_01[nDiffY*50+nDiffX];
	if (ex > sx) {
		if(ey < sy) nAngle = 240 - nAngle;
	} 
	else {
		if(ey > sy) nAngle = 120 - nAngle;
		else nAngle += 120;
	}
	nAngle = ANGLE_0_240(nAngle);

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
		if (sqr_val == tmp1) return sqrt_val;
		if (sqr_val > tmp1) sqrt_val += tmp2;
		else sqrt_val -= tmp2;
	} while (tmp2 > 1);
	tmp1 = sqrt_val*sqrt_val;
	if (tmp1 == sqr_val) return sqrt_val;
	if (tmp1 > sqr_val) {
		tmp2 = tmp1 - sqrt_val;
		if (tmp2 >= sqr_val) return sqrt_val-1;
	} 
	else {
		tmp2 = tmp1 + sqrt_val;
		if (tmp2 < sqr_val) return sqrt_val+1;
	}
	return sqrt_val;
}

BOOL check_in_polygon(int point_x,int point_y,POLYGON *polygon,int limit)
{
	int i, i1, i2, num, x1, x2, x3, y1, y2, y3, lim2, nSign, mx0, mx1, my0, my1;

	if (polygon->orderNum < 3) return FALSE;
	num = polygon->orderNum; nSign = 1;
	if (limit < 0) nSign = -1;
	lim2 = limit*limit;
	x1 = polygon->x[0]; y1 = polygon->y[0];
	x3 = polygon->x[num-1]; y3 = polygon->y[num-1];

	for (i = 0; i < num; i++) {
		if (i+1 < num) {
			x2 = polygon->x[i+1]; y2 = polygon->y[i+1];
		}
		else {
			x2 = polygon->x[0]; y2 = polygon->y[0];
		}
		i1 = (y1-y2)*(point_x-x2) - (x1-x2)*(point_y-y2);
		if (i1 == 0) {
			if (x1 < x2) {	mx0 = x1; mx1 =	x2;	}
			else { mx0 = x2; mx1 = x1; }
			if (y1 < y2) {	my0 = y1; my1 =	y2;	}
			else { my0 = y2; my1 = y1; }
			if (mx0 <= point_x && point_x <= mx1 && my0 <= point_y && point_y <= my1) {
				if (nSign > 0) return TRUE;
				return FALSE;
			}
		}
		i2 = (y1-y2)*(x3-x2) - (x1-x2)*(y3-y2);
		if ((i1>0 && i2<0) || (i1<0 && i2>0)) {
			if (limit == 0 || nSign < 0) return FALSE;
			if (x1 == x2 ) i1 = (point_x-x2)*(point_x-x2);
			else {
				i2 = op_func_02((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)) * 100;
				if (i2 == 0) return FALSE;
				i2 = (i1 * 100)/i2;
				i1 = i2*i2;
			}
			if (i1 > lim2) return FALSE;
		}
		if (nSign < 0) {
			if (x1 == x2) i1 = (point_x-x2)*(point_x-x2);
			else {
				i2 = op_func_02((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)) * 100;
				if (i2 == 0) return FALSE;
				i2 = (i1 * 100)/i2;
				i1 = i2*i2;
			}
			if (i1 < lim2) return FALSE;
		}
		x3 = x1; y3 = y1; x1 = x2; y1 = y2;
	}
	return TRUE;
}

BOOL get_polygon_points_sub(int *point_x,int *point_y,int pNum,POLYGON *polygon)
{
	int i, j, i1, i2, i3, j1, j2, k1, k2;
	int boundNum = 1;
	int endFlag = 0,repeatFlag;

	if (pNum < 3) return FALSE;

	i1 = 1000;
	for (i = 0; i < pNum; i++) {
		if (point_x[i] < i1) {
			i1 = point_x[i];
			polygon->x[0] = i1;
			polygon->y[0] = point_y[i];
		}
		else if (point_x[i] == i1 && polygon->y[0] > point_y[i]) {
			i1 = point_x[i];
			polygon->x[0] = i1;
			polygon->y[0] = point_y[i];
		}
	}
	while(1) {
		repeatFlag = 0;
		for (i = 0;i < pNum;i++)	{
			if(point_x[i] == polygon->x[boundNum-1] && 
				point_y[i] == polygon->y[boundNum-1]) continue;
			if(boundNum > 1 && point_x[i] == polygon->x[boundNum-2] && 
				point_y[i] == polygon->y[boundNum-2]) continue;

			j1 = point_x[i];
			j2 = point_y[i];

			i1 = i2 = 0;
			for (j = 0; j < pNum; j++) {
				if (i == j) continue;
				if (point_x[j] == polygon->x[boundNum-1] && 
					point_y[j] == polygon->y[boundNum-1]) continue;
				i3 = (point_y[j] - polygon->y[boundNum-1])*(j1-polygon->x[boundNum-1]);
				i3 -= (point_x[j] - polygon->x[boundNum-1])*(j2-polygon->y[boundNum-1]);
				if (i3 < 0) i1++;
				if (i3 > 0) i2++;
				if (i3 == 0) {
					k1 = point_x[j] - j1;
					k2 = point_x[j] - polygon->x[boundNum-1];
					if (k1*k2 > 0 && abs(k1) < abs(k2))	{
						i1++;	i2++;
					}
					k1 = point_y[j] - j2;
					k2 = point_y[j] - polygon->y[boundNum-1];
					if (k1*k2 > 0 && abs(k1) < abs(k2)) {
						i1++;	i2++;
					}
					if (i1 == 0 || i2 == 0)	{
						if (point_x[j] == polygon->x[0] && point_y[j] == polygon->y[0])	{
							endFlag = 1; break;
						}
					}
				}
				if (i1 > 0 && i2 > 0) break;
			}
			if (i1 > 0 && i2 > 0) continue;
			if (polygon->x[0] == j1 && polygon->y[0] == j2) endFlag = 1;
			if (endFlag == 1) break;
			polygon->x[boundNum] = j1;
			polygon->y[boundNum] = j2;
			boundNum++;
			repeatFlag = 1;
			break;
		}
		if (endFlag == 1 || repeatFlag == 0) break;
	}

	if (endFlag == 0) {
		polygon->orderNum = 0;
		return FALSE;
	}
	polygon->orderNum = boundNum;
	if (boundNum < 3) return FALSE;
	return TRUE;
}

BOOL get_polygon_points(LPMPVECTEX pVect,POLYGON *polygon)
{
	int i, xx[MAX_MINUTIA_NUMBER], yy[MAX_MINUTIA_NUMBER], num;

	num = pVect->nNumber;
	for (i = 0; i < num; i++) {
		xx[i] = pVect->item[i].x;
		yy[i] = pVect->item[i].y;
	}
	return(get_polygon_points_sub(xx,yy,num,polygon));
}

void get_smoothed_image(BYTE* img,int cxDIB,int cyDIB)
{
	int i, j, k, nSum, rownum = 0;
	int st = cxDIB - 1, dptr = 0;
	int pSum[MAX_IMG_WIDTH];
	BYTE pTmp[3*MAX_IMG_WIDTH], *row, *ptr = img - cxDIB + 1;
	
	memset( pSum, 0, sizeof(int)*cxDIB );
	
	for (i = 0; i < cyDIB+1; i++, dptr+=cxDIB) {
		k = i % 3;
		if (i >= 3) {
			row = pTmp + k * cxDIB;
			for (j = 0; j < cxDIB; j++) pSum[j] -= row[j];
			rownum--;
		}
		if (i < cyDIB) {
			row = pTmp + k * cxDIB;
			memcpy(row, img + dptr, sizeof(BYTE)*cxDIB);
			for (j = 0; j < cxDIB; j++) pSum[j] += row[j];
			rownum++;
		}
		if (i < 1) continue;
		row = ptr + dptr; j = st;
		nSum = pSum[j];
		j--;
		nSum += pSum[j];
		if (rownum == 2) {
			row[j] =  (BYTE)(nSum >> 2);
			j--;
			nSum += pSum[j];
			row[j] =  div6_table[nSum];
			j--;
			for (; j>=0; j--) {
				nSum += pSum[j]; nSum -= pSum[j+3];
				row[j] = div6_table[nSum];
			}
			nSum -= pSum[2];
			row[-1] = (BYTE)(nSum >> 2);
			continue;
		}
		row[j] =  div6_table[nSum];
		j--;
		nSum += pSum[j];
		row[j] =  div9_table[nSum];
		j--;
		for (; j >= 0; j--) {
			nSum += pSum[j]; nSum -= pSum[j+3];
			row[j] = div9_table[nSum];
		}
		nSum -= pSum[2];
		row[-1] = div6_table[nSum];
	}
}

void histogram_smooth(int *Histo,int nLen,int nSize)
{
	int i, k = 0, n = 0, sum = 0, nWnd = nSize*2 + 1;
	int pSrc[MAX_IMG_WIDTH];
	
	for (i = 0; i < nLen; i++) { pSrc[i] = Histo[i]; }
	
	for (i = 0; i < nLen+nSize; i++) {
		if (i < nLen) { sum += pSrc[i]; n++; }
		if (i < nSize) continue;
		if (i >= nWnd) { sum -= pSrc[i-nWnd]; n--; }
		Histo[k++] = (int)(sum / n);
	}
}

void get_segmentation(BYTE *pImg,BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, k, h = 0, m, k1, k2 , mean, th, szBlk = 16, wid = cxDIB/4, hei = cyDIB/4;
	int *varHist, xHist[MAX_IMG_WIDTH/4], yHist[MAX_IMG_HEIGHT/4];
	int lf = 0, rf = cxDIB-1, tf = 0, bf = cyDIB-1, nl, nr;
	int rows = cyDIB / BLOCK_SIZE;
	BYTE *ptr;

	if (cyDIB % 16 != 0) rows++;

	varHist = (int*)calloc((cxDIB/4)*rows,sizeof(int));
	if (varHist == NULL) return;

	memset(xHist,0,sizeof(int)*cxDIB/4);
	memset(yHist,0,sizeof(int)*cyDIB/4);

	for (i = 0; i < cyDIB; i += szBlk, h++) {
		if (i+14 >= cyDIB) break;
		ptr = pImg + i*cxDIB;
		k = 0;
		for (j = 0; j < cxDIB; j += 4) {
			mean = (ptr[j] + ptr[2*cxDIB+j] + ptr[4*cxDIB+j] + ptr[6*cxDIB+j]);
			mean += (ptr[8*cxDIB+j] + ptr[10*cxDIB+j] + ptr[12*cxDIB+j] + ptr[14*cxDIB+j]);
			mean = mean / 8;

			varHist[h*wid+k] = (abs(mean-ptr[j]) + abs(mean-ptr[2*cxDIB+j]) + abs(mean-ptr[4*cxDIB+j]) + abs(mean-ptr[6*cxDIB+j]));
			varHist[h*wid+k] += (abs(mean-ptr[8*cxDIB+j]) + abs(mean-ptr[10*cxDIB+j]) + abs(mean-ptr[12*cxDIB+j]) + abs(mean-ptr[14*cxDIB+j]));

			xHist[k] += varHist[h*wid+k];

			k = k + 1;
		}

		histogram_smooth(&varHist[h*wid], wid, 2);
	}
	mean = 0;
	for (j = 0; j < wid; j++)
		if (mean < xHist[j]) mean = xHist[j];
	
		if (mean > 2000) {
			th = mean / 5;
			if (th < 1000) th = 1000;
		}
		else th = mean / 5;

	for (j = 0; j < wid; j++) {
		if (xHist[j] > th) break;
	} 
	lf = 4*j;
	if (lf > cxDIB) lf = cxDIB;

	for (j = wid-1; j >= 0; j--) {
		if (xHist[j] > th) break; 
	}
	rf = 4*j;
	if (j == wid-1) rf += 3;
	if (rf < 0) rf = 0;

	for (i = lf; i <= rf; i += 8) {
		if (i+8 >= cxDIB) break;
		k = 0;
		for (j = 0; j < cyDIB; j += 4) {
			ptr = pImg + j*cxDIB; h = j / 16;
			mean = (ptr[i] + ptr[i+2] + ptr[i+4] + ptr[i+6] ) / 4;
			yHist[k] += abs(mean-ptr[i]) + abs(mean-ptr[i+2]) + abs(mean-ptr[i+4]) + abs(mean-ptr[i+6]);
			yHist[k++] += varHist[h*wid+i/4]/4;
		}
	}

	mean = 0;
	for (j = 0; j < hei; j++)
		if (mean < yHist[j]) mean = yHist[j];
	
	th = mean / 7;

	if (th > 1000) th = 1000;
	for (j = 0; j < hei; j++) {
		if (yHist[j] > th) break; 
	} 
	tf = 4*j;
	if (tf > cyDIB) tf = cyDIB;

	if (mean > 3000) {
		th = mean / 5;
		if (th < 1600) th = 1600;
	}
	else th = mean / 5;
	for (j = hei-1; j >= 0; j--) {
		if (yHist[j] > th) {
			if (j-4 > 0) {
				if (yHist[j-4] > th) break;
				continue;
			}
			break;
		}
	}
	bf = 4*j;
	if (j == hei-1) bf += 3;
	if (bf < 0) bf = 0;

	ptr = OrntImg;
	for (i = 0; i < cyDIB; i++) {
		for (j = 0; j < lf; j++) { ptr[j] = 255; }
		for (j = cxDIB-1; j >= rf; j--) { ptr[j] = 255; }
		ptr += cxDIB;
	}
	ptr = OrntImg;
	for (i = 0; i < tf; i++) {
		for (j = lf; j < rf; j++) { ptr[j] = 255; }
		ptr += cxDIB;
	}
	ptr = OrntImg + (cyDIB-1)*cxDIB;
	for (i = cyDIB-1; i >= bf; i--) {
		for (j = lf; j < rf; j++){ ptr[j] = 255; }
		ptr -= cxDIB;
	}

	h = 100;
	if (tf % 16 != 0) h = tf / 16;
	k1 = lf / 4; k2 = rf / 4;
	for (i = tf; i < bf; i++) {
		k = i / 16;
		if (k == h) szBlk = (tf/16 + 1)*16 - tf;
		else szBlk = 16;
		if (i+szBlk >= bf) szBlk = bf-i+1;
		mean = 0; 
		for (j = k1; j <= k2; j++)
			if (mean < varHist[k*wid+j]) mean = varHist[k*wid+j];
		
			if (mean > 350) th = mean / 5;
			else {
				th = mean / 5;
				if (th < 70) th = 70;
			}
		
		for (j = k1; j <= k2; j++) {
			if (varHist[k*wid+j] > th) break;
		}
		nl = 4*j;
		if (nl >= cxDIB) nl = cxDIB-1;
		for (j = k2; j >= k1; j--) {
			if (varHist[k*wid+j] > th) break;
		}
		if (j < 0) j = 0; 
		nr = 4*j;
		ptr = OrntImg + i*cxDIB;
		for (m = 0; m < szBlk; m++) {
			for (j = lf; j < nl; j++) { ptr[j] = 255; }
			for (j = rf; j > nr; j--) { ptr[j] = 255; }
			ptr += cxDIB;
		}
		i += szBlk-1;			
	}

	free(varHist);
}

void get_smoothed_image4(BYTE* img,int cxDIB,int cyDIB)
{
	int i, j, k, nSum, rownum = 0;
	int *pSum;
	BYTE *pTmp, *row, *ptr;
	int dptr = 0;
	

	pTmp = (BYTE*)malloc(sizeof(BYTE)*9*cxDIB);
	if (pTmp == NULL) return;

	pSum = (int*)calloc(cxDIB,sizeof(int));
	if (pSum == NULL) { free(pTmp); return; }

	ptr = img - 4*cxDIB + 4;
	dptr = 0;

	for (i = 0; i < cyDIB+4; i++, dptr+=cxDIB) {
		k = i % 9;
		if (i >= 9) {
			row = pTmp + k * cxDIB;
			for (j = 0; j < cxDIB; j++) pSum[j] -= row[j];
			rownum--;
		}
		if (i < cyDIB) {
			row = pTmp + k * cxDIB;
			memcpy(row, img + dptr, sizeof(BYTE)*cxDIB);
			for (j = 0; j < cxDIB; j++) pSum[j] += row[j];
			rownum++;
		}
		if (i < 4) continue;

		row = ptr + dptr;

		j = cxDIB - 1;
		nSum = pSum[j];
		j--;
		nSum += pSum[j];
		j--;
		nSum += pSum[j];
		j--;
		nSum += pSum[j];
		j--;
		nSum += pSum[j];

		if (rownum == 5) {
			row[j] =  (BYTE)(nSum / 25);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 30);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 35);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 40);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)( (91*nSum) >> 12 );
			j--;
			for (; j >= 0; j--) {
				nSum += pSum[j];  nSum -= pSum[j+9];
				row[j] = (BYTE)( (91*nSum) >> 12 );
			}
			nSum -= pSum[8];
			row[-1] = (BYTE)(nSum / 40);
			nSum -= pSum[7];
			row[-2] = (BYTE)(nSum / 35);
			nSum -= pSum[6];
			row[-3] = (BYTE)(nSum / 30);
			nSum -= pSum[5];
			row[-4] = (BYTE)(nSum / 25);
			continue;
		}

		if (rownum == 6) {
			row[j] =  (BYTE)(nSum / 30);
			j--;
			nSum += pSum[j];
			 row[j] =  (BYTE)(nSum / 36);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 42);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 48);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)( (19*nSum) >> 10 );
			j--;
			for (; j >= 0; j--) {
				nSum += pSum[j];  nSum -= pSum[j+9];
				row[j] = (BYTE)( (19*nSum) >> 10 );
			}
			nSum -= pSum[8];
			row[-1] = (BYTE)(nSum / 48);
			nSum -= pSum[7];
			row[-2] = (BYTE)(nSum / 42);
			nSum -= pSum[6];
			row[-3] = (BYTE)(nSum / 36);
			nSum -= pSum[5];
			row[-4] = (BYTE)(nSum / 30);
			continue;
		}

		if (rownum == 7) {
			row[j] =  (BYTE)(nSum / 35);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 42);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 49);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 56);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)( (65*nSum) >> 12 );
			j--;
			for (; j >= 0; j--) {
				nSum += pSum[j];  nSum -= pSum[j+9];
				row[j] = (BYTE)( (65*nSum) >> 12 );
			}
			nSum -= pSum[8];
			row[-1] = (BYTE)(nSum / 56);
			nSum -= pSum[7];
			row[-2] = (BYTE)(nSum / 49);
			nSum -= pSum[6];
			row[-3] = (BYTE)(nSum / 42);
			nSum -= pSum[5];
			row[-4] = (BYTE)(nSum / 35);
			continue;
		}

		if (rownum == 8) {
			row[j] =  (BYTE)(nSum / 40);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 48);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)(nSum / 56);
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)( nSum >> 6 );
			j--;
			nSum += pSum[j];
			row[j] =  (BYTE)( (455*nSum) >> 15 );
			j--;
			for (; j >= 0; j--) {
				nSum += pSum[j];  nSum -= pSum[j+9];
				row[j] = (BYTE)( (455*nSum) >> 15 );
			}
			nSum -= pSum[8];
			row[-1] = (BYTE)( nSum >> 6 );
			nSum -= pSum[7];
			row[-2] = (BYTE)(nSum / 56);
			nSum -= pSum[6];
			row[-3] = (BYTE)(nSum / 48);
			nSum -= pSum[5];
			row[-4] = (BYTE)(nSum / 40);
			continue;
		}

		row[j] =  (BYTE)( (91*nSum) >> 12 );
		j--;
		nSum += pSum[j];
		row[j] =  (BYTE)( (19*nSum) >> 10 );
		j--;
		nSum += pSum[j];
		row[j] =  (BYTE)( (65*nSum) >> 12 );
		j--;
		nSum += pSum[j];
		row[j] =  (BYTE)( (455*nSum) >> 15 );
		j--;
		nSum += pSum[j];
		row[j] =  (BYTE)( (809*nSum) >> 16 );
		j--;
		for (; j >= 0; j--) {
			nSum += pSum[j];  nSum -= pSum[j+9];
			row[j] = (BYTE)( (809*nSum) >> 16 );
		}
		nSum -= pSum[8];
		row[-1] = (BYTE)( (455*nSum) >> 15 );
		nSum -= pSum[7];
		row[-2] = (BYTE)( (65*nSum) >> 12 );
		nSum -= pSum[6];
		row[-3] = (BYTE)( (19*nSum) >> 10 );
		nSum -= pSum[5];
		row[-4] = (BYTE)( (91*nSum) >> 12 );
		
	}

	free(pSum);  free(pTmp);
}

void get_sharpend_image(BYTE *Img,BYTE *RefImg,BYTE *OrntImg,int cxDIB,int cyDIB,int nStep)
{
	int nWindow = 2*nStep + 1,  cx2 = cxDIB / 2,  nW = nStep + 1;
	int i, j, k, ii, rownum = 0, colnum = 0, nSum, nMean, nImgVal, nRefVal;
	int id = -1, idd = -cx2, idx, n;
	int *pSum;
	BYTE *pAbs,  *row,  *refptr,  *imgptr;
	int offset[4];

	offset[0] = 0; offset[1] = 1; offset[2] = cxDIB; offset[3] = cxDIB + 1; 

	pSum = (int*)calloc(cxDIB, sizeof(int));
	if (pSum == NULL)  return;
	
	pAbs = (BYTE*)malloc(sizeof(BYTE) * nW * cx2);
	if (pAbs == NULL) { free(pSum);  return; }
	
	for (i = 0; i < cyDIB+nStep; i += 2) {
		if (id >= nStep) { id = idd = 0; }
		else { id++; idd += cx2; } 

		if (i >= nWindow) {
			row = pAbs + idd;
			for (j = 0; j < cxDIB; j += 2) { 
				pSum[j] -= row[j >> 1];
			}
			rownum--;
		}

		if (i < cyDIB) {
			row = pAbs + idd;  imgptr = Img + i*cxDIB;  refptr = RefImg + i*cxDIB;
			for (j = 0; j < cxDIB; j += 2) {
				row[j >> 1] = abs(imgptr[j] - refptr[j]);
				pSum[j] += row[j >> 1];
			}
			rownum++;
		}

		if (i < nStep)  continue;

		nSum = colnum = 0;

		k = (i - nStep) * cxDIB - nStep;

		for (j = 0; j < cxDIB+nStep; j += 2, k += 2) {
			if (j < cxDIB) {
				colnum++;  nSum += pSum[j];
			}

			if (j < nStep)  continue;

			if (j >= nWindow) {
				colnum--;  nSum -= pSum[j-nWindow-1];
			}

			nMean = (nSum * divX_table1[rownum]) >> divX_table2[rownum] ;
			nMean = (nMean * divX_table1[colnum]) >> divX_table2[colnum];

			for (ii = 3; ii >= 0; ii--)	{
				idx = k + offset[ii];

				if (OrntImg[idx] == 255)  continue;

				if (nMean == 0)	{
					Img[idx] = RefImg[idx];  continue;
				}

				nImgVal = Img[idx];  nRefVal = RefImg[idx];
				
				if (nRefVal <= nImgVal-nMean) {
					Img[idx] = (BYTE)(nRefVal >> 1);
					continue;
				}

				if (nRefVal >= nImgVal+nMean) { 
					Img[idx] = (BYTE)((255 + nRefVal) >> 1);
					continue;
				}
			
				if (nImgVal-nMean <= 0)	{
					if (nImgVal+nMean >= 255) {
						Img[idx] = nRefVal;  continue;
					}

					n = nImgVal + nMean;
					nImgVal =  (nRefVal * 255 * divX_table1[n] ) >> divX_table2[n];
					Img[idx] = (BYTE)((nImgVal + nRefVal) >> 1);
					continue;
				}

				if (nImgVal+nMean >= 255) {
					n = 255 - nImgVal + nMean;
					nImgVal = 255 * (nRefVal - nImgVal + nMean) * divX_table1[n] >> divX_table2[n];
					Img[idx] = (BYTE)((nImgVal + nRefVal) >> 1);
					continue;
				}

				n = 2 * nMean;
				nImgVal = 255 * (nRefVal - nImgVal + nMean) * divX_table1[n] >> divX_table2[n];
				Img[idx] = (BYTE)((nImgVal + nRefVal) >> 1);
			}
		}
 	}

	free(pAbs);  free(pSum);
}

void get_orient_image(BYTE* OrntImg,BYTE* Img,int cxDIB,int cyDIB,BYTE *image_buffer4)
{
	int dir0, dir1, dir2, dir3, sum0, sum1, sum2, sum3, var0, var1, var2, var3;
	int nMean, nVar, nMax, nAvg, nVal, s1, s2;
	int i, j, k, kk, idx1 = 0, idx2 = 0, idx3 = 0;
	int cx2 = cxDIB/2 - 1, kkk = 0;
	int *pSum, *sumptr0, *sumptr1, *sumptr2, *sumptr3;
	unsigned short *pAbs[19], *absptr0, *absptr1, *absptr2, *absptr3;
	BYTE *pTmp, *tmpptr0, *tmpptr3, *tmpptr6, *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8;
	BYTE *imgptr, *ptr2, *dstptr1, *dstptr2, *orntptr;

	pTmp = (BYTE*)malloc(sizeof(BYTE)*38*cxDIB);
	if (pTmp == NULL) return;
	pSum = (int*)calloc(2*cxDIB, sizeof(int));
	if (pSum == NULL) { free(pTmp); return; }

	for (i = 0; i < 19; i++) {
		pAbs[i] = (unsigned short*)malloc(sizeof(unsigned short)*2*cxDIB);
		if (pAbs[i] == NULL) { 
			for (j = 0; j < i; j++) free(pAbs[j]); 
			free(pSum); free(pTmp); return;
		}
	}

	tmpptr0 = pTmp + 2;  tmpptr3 = tmpptr0 - cxDIB;  tmpptr6 = tmpptr0 + cxDIB;
	sumptr0 = pSum; sumptr1 = pSum +1; sumptr2 = pSum + 2; sumptr3 = pSum + 3;

	imgptr = Img;  orntptr = OrntImg - 18*cxDIB;

	for (i = 0; i < cyDIB+18; i += 2) {
		k = i % 38;  dstptr1 = pTmp + k*cxDIB;  dstptr2 = dstptr1 + cxDIB;
		if (i < cyDIB) {
			ptr2 = imgptr + cxDIB;
			for (j = 0; j < cxDIB; j++) {
				dstptr1[j] = imgptr[j]; dstptr2[j] = ptr2[j];
			}
			imgptr += 2*cxDIB; 
		}
		if (i < 4) continue;
		if (i < cyDIB) {
			if (idx3 >= 19) { idx3 = 1; tmpptr3 = pTmp + cxDIB + 2; }
			else { idx3++; tmpptr3 += 2*cxDIB; }
			if (idx1 >= 18) { idx1 = 0; tmpptr0 = pTmp + 2; tmpptr6 = tmpptr0 + cxDIB; }
			else { idx1++; tmpptr0 += 2*cxDIB; tmpptr6 += 2*cxDIB; }
			absptr0 = pAbs[idx1]; absptr1 = absptr0 + 1; absptr2 = absptr1 + 1; absptr3 = absptr2 + 1;
			p0 = tmpptr0; p1 = p0 - 1;  p2 = p0 + 1;
			p3 = tmpptr3; p4 = p3 - 1;  p5 = p3 + 1;
			p6 = tmpptr6; p7 = p6 - 1;  p8 = p6 + 1;
			for (k = 1; k < cx2; k++) {
				absptr0[k*4] = abs(*p0 - *p1) + abs(*p0 - *p2);
				absptr1[k*4] = abs(*p0 - *p4) + abs(*p0 - *p8);
				absptr2[k*4] = abs(*p0 - *p3) + abs(*p0 - *p6);
				absptr3[k*4] = abs(*p0 - *p7) + abs(*p0 - *p5);
				sumptr0[k*4] += absptr0[k*4];
				sumptr1[k*4] += absptr1[k*4];
				sumptr2[k*4] += absptr2[k*4];
				sumptr3[k*4] += absptr3[k*4];
				p0 += 2; p1 += 2; p2 += 2; p3 += 2; p4 += 2; p5 += 2; p6 += 2; p7 += 2; p8 += 2;
			}
		}
		if (i < 18) continue;
		if (i >= 36) {
			if (idx2 >= 18) idx2 = 0;  
			else idx2++;	
			absptr0 = pAbs[idx2]; absptr1 = absptr0 + 1; absptr2 = absptr1 + 1; absptr3 = absptr2 + 1;
			for (k = 1; k < cx2; k++) {
				sumptr0[k*4] -= absptr0[k*4];
				sumptr1[k*4] -= absptr1[k*4];
				sumptr2[k*4] -= absptr2[k*4];
				sumptr3[k*4] -= absptr3[k*4];
			}
		}
	
		ptr2 = orntptr + i*cxDIB;
		k = 0;  kk = 0;
		sum0 = sum1 = sum2 = sum3 = 0;
		for (j = 0; j < cxDIB+18; j += 2) {
			if (j < cxDIB) {
				sum0 += sumptr0[k];
				sum1 += sumptr1[k];
				sum2 += sumptr2[k];
				sum3 += sumptr3[k];
				k += 4;
			}
			if (j < 18) continue;
			if (j >= 36) {
				sum0 -= sumptr0[kk];
				sum1 -= sumptr1[kk];
				sum2 -= sumptr2[kk];
				sum3 -= sumptr3[kk];
				kk += 4;
			}

			if (*ptr2 == 255) { ptr2 += 2; image_buffer4[kkk++] = 0; continue; }

			nMean = var0 = dir0 = sum0;
			var1 = dir1 = (91 * sum1) >> 7;
			nMean += var1;
			var2 = dir2 = sum2;
			nMean += var2;
			var3 = dir3 = (91 * sum3) >> 7;
			nMean += var3;
			nMean /= 4;
			
			if (nMean != 0) {
				nVar = abs(nMean - var0);
				nVar += abs(nMean - var1);
				nVar += abs(nMean - var2);
				nVar += abs(nMean - var3);
				nVar = nVar*40/nMean;
				image_buffer4[kkk++]= nVar;
			}
			else
				image_buffer4[kkk++]= 0;
			
			nMax = dir0;
			if (nMax < dir1) nMax = dir1;
			if (nMax < dir2) nMax = dir2;
			if (nMax < dir3) nMax = dir3;

			nAvg = (nMean - nMax) * 4;
			if (nAvg == 0) { *ptr2 = 0x7F; ptr2 += 2; continue; }

			nVal = 45; s1 = var2 + var3; s2 = var2 + var1;
			if (s1 < s2) {
				s2 = s1; nVal = 75;
				dir0 = var1;	dir1 = var2;
				dir2 = var3;	dir3 = var0;
			}
			s1 = var0 + var3;
			if (s1 < s2) {
				s2 = s1; nVal = 105;
				dir0 = var2;	dir1 = var3;
				dir2 = var0;	dir3 = var1;
			}
			s1 = var0 + var1;
			if (s1 < s2) {
				nVal = 15;
				dir0 = var3;	dir1 = var0;
				dir2 = var1;	dir3 = var2;
			}
			nVal += 15 * (3*(dir3 - dir0) - dir1 + dir2) / nAvg;
			if (nVal == 120) nVal = 0;
			*ptr2 = (BYTE)(nVal);  ptr2 += 2;
		}
	}
	p0 = OrntImg; p2 = p0 + cxDIB;
	for (i = 0; i < cyDIB; i += 2) {
		p1 = p0 + 1; p3 = p2 + 1;
		for (j = 0; j < cxDIB; j += 2) {
			p1[j] = p2[j] = p3[j] = p0[j];
		}
		p0 += 2*cxDIB; p2 = p0 + cxDIB;
	}

	for (i = 0; i < 19; i++) free(pAbs[i]);
	free(pSum); free(pTmp);
}

void image_proc_01(BYTE* Img,BYTE* OrntImg,BYTE* pTmpImg,int cxDIB,int cyDIB)
{
	int i, j, k, n = 0, nSum, idx;
	BYTE nOrnt;
	short _tmptbl[120*18];

	memcpy( pTmpImg, Img, sizeof(BYTE)*cxDIB*cyDIB );

	for (i = 0; i < 120*18; i++) {
		_tmptbl[i] = _table3[i]*cxDIB + _table4[i];
	}

	for (i = 0; i < cyDIB-12; i += 2) {
		n = (i+6)*cxDIB + 6;
		for (j = cxDIB-14; j >= 0; j -= 2, n += 2) {
			nOrnt = OrntImg[n];
			if (nOrnt >= 120) continue; 
			nSum = pTmpImg[n] * _table1[nOrnt];
			idx = nOrnt * 18;
			for (k = _table2[nOrnt]-1; k >= 0; k--, idx++) {
				nSum += (pTmpImg[n + _tmptbl[idx]]+pTmpImg[n - _tmptbl[idx]]) * _table5[idx];
			}
			Img[n] = (BYTE)(nSum >> 14);
		}
	}

	for (i = 0; i < cyDIB-12; i += 2) {
		for (j = 0; j < cxDIB-14; j += 2) {
			n = (i+6)*cxDIB + j+6;
			Img[n+1] = ( Img[n] + Img[n+2] ) / 2;
		}
		Img[n+1] = Img[n];
	}

	for (i = 0; i < cyDIB-15; i += 2) {
		for (j = 0; j < cxDIB-12; j++) {
			n = (i+6+1)*cxDIB + j+6;
			Img[n] = (Img[n-cxDIB] + Img[n+cxDIB]) / 2;
		}
	}
	memcpy( pTmpImg, Img, sizeof(BYTE)*cxDIB*cyDIB );
}

void get_binary_image2(BYTE *OrntImg,BYTE *BinImg,BYTE *RefImg,int cxDIB,int cyDIB,int nStep1,int nStep2)
{
	int i, j, nTh, nWindow1 = 2*nStep1+1, nWindow2 = 2*nStep2+1;
	int *pSum1, *pSum2, nSum, nNum, nNum1 = 0, nNum2 = 0;

	pSum1 = (int*)calloc(cxDIB, sizeof(int));
	if (pSum1 == NULL) return;
	pSum2 = (int*)calloc(cxDIB, sizeof(int));
	if (pSum2 == NULL) { free(pSum1); return; }

	for (i = 0; i < cyDIB+nStep2; i++) {
		if (i < cyDIB) {
			for (j = 0; j < cxDIB; j++) {
				pSum1[j] += RefImg[i*cxDIB+j];
				pSum2[j] += RefImg[i*cxDIB+j];
			}
			nNum1++; nNum2++;
		}
		if (i-nStep1 >= 0 && i-nStep1 < cyDIB) {
			if (i >= nWindow1) {
				for (j = 0; j < cxDIB; j++) {
					pSum1[j] -= RefImg[(i-nWindow1)*cxDIB+j];
				}
				nNum1--;
			}
			nSum =0; nNum = 0;
			for (j = 0; j < cxDIB+nStep1; j++) {
				if (j < cxDIB) {
					nSum += pSum1[j]; nNum += nNum1;
				}
				if (j < nStep1) continue;
				if (j >= nWindow1) {
					nSum -= pSum1[j-nWindow1]; nNum -= nNum1;
				}
				if ((OrntImg[(i-nStep1)*cxDIB+(j-nStep1)]&0x80) == 0) {
					BinImg[(i-nStep1)*cxDIB+(j-nStep1)] = (BYTE)(nSum/nNum);
				}
			}
		}
		if (i < nStep2) continue;
		if (i >= nWindow2) {
			for (j = 0; j < cxDIB; j++) {
				pSum2[j] -= RefImg[(i-nWindow2)*cxDIB+j];
			}
			nNum2--;
		}
		nSum =0; nNum = 0;
		for (j = 0; j < cxDIB+nStep2; j++) {
			if (j < cxDIB) {
				nSum += pSum2[j]; nNum += nNum2;
			}
			if (j < nStep2) continue;
			if (j >= nWindow2) {
				nSum -= pSum2[j-nWindow2]; nNum -= nNum2;
			}
			if ((OrntImg[(i-nStep2)*cxDIB+(j-nStep2)] & 0x80) == 0) {
				nTh = (nSum/nNum + BinImg[(i-nStep2)*cxDIB+(j-nStep2)]) / 2;
				if (RefImg[(i-nStep2)*cxDIB+(j-nStep2)] >= nTh)
					BinImg[(i-nStep2)*cxDIB+(j-nStep2)] = 0xFF;
				else 
					BinImg[(i-nStep2)*cxDIB+(j-nStep2)] = 0x00;
			}
		}
	}
	free(pSum1); free(pSum2);
}

void image_proc_04(BYTE *Img,int cxDIB,int cyDIB)
{
	int i, j, nSum, *pSum;
	BYTE *pTmp1, *pTmp2, *pTmp3, *pTmp;

	pTmp1 = (BYTE*)malloc(cxDIB);
	if (pTmp1 == NULL) return;

	pTmp2 = (BYTE*)malloc(cxDIB);
	if (pTmp2 == NULL) { free(pTmp1); return; }

	pTmp3 = (BYTE*)malloc(cxDIB);
	if (pTmp3 == NULL) {
		free(pTmp1); free(pTmp2); return;
	}

	pSum = (int*)calloc(cxDIB, sizeof(int));
	if (pSum == NULL) {
		free(pTmp1); free(pTmp2); free(pTmp3); return;
	}

	for (i = 0; i < cyDIB+1; i++) {
		if (i >= 3) {
			for (j = 0; j < cxDIB; j++) pSum[j] -= pTmp1[j];
		}
		if (i < cyDIB) {
			memcpy(pTmp1, &Img[i*cxDIB], cxDIB);
			for (j = 0; j < cxDIB; j++) pSum[j] += pTmp1[j];
		}
		pTmp = pTmp1; pTmp1 = pTmp2; pTmp2 = pTmp3; pTmp3 = pTmp;
		if (i < 2) continue;
		nSum = 0;
		for (j = 0; j < cxDIB+1; j++) {
			if (j >= 3) nSum -= pSum[j-3];
			if (j < cxDIB) nSum += pSum[j];
			if (j < 2) continue;
			Img[(i-1)*cxDIB+j-1] = (nSum < 1152) ? 0x00:0xFF;
		}
	}

	free(pTmp1); free(pTmp2); free(pTmp3); free(pSum);
}

void remove_hole(BYTE* OrntImg,BYTE* Img,int cxDIB,int cyDIB)
{
	int i, j, k, ix, iy, wNum, whiteX[50], whiteY[50];
	int numTH = 20;
	BOOL flag;

	for (i = 0; i < cyDIB; i++) {
		for (j = 0; j < cxDIB; j++) {
			if ((OrntImg[i*cxDIB+j] & 0x80) != 0) continue;
			if (Img[i*cxDIB+j] != 255) continue;
			wNum = 1;
			whiteX[0] = j;
			whiteY[0] = i;
			Img[i*cxDIB+j] = 0;
			flag = TRUE;
			for (k = 0; k < wNum; k++) {
				ix = whiteX[k];
				iy = whiteY[k];
				if (ix < 1 || ix > cxDIB-2 || iy < 1 || iy > cyDIB-2) {
					flag = FALSE;
					break;
				}
				if (Img[(iy-1)*cxDIB+ix] == 255) {
					if (iy <= i) {
						flag = FALSE;
						break;
					}
					Img[(iy-1)*cxDIB+ix] = 0;
					whiteX[wNum] = ix;
					whiteY[wNum] = iy-1;
					wNum++;
				}
				if (Img[iy*cxDIB+ix+1] == 255) {
					Img[iy*cxDIB+ix+1] = 0;
					whiteX[wNum] = ix+1;
					whiteY[wNum] = iy;
					wNum++;
				}
				if (Img[(iy+1)*cxDIB+ix] == 255) {
					Img[(iy+1)*cxDIB+ix] = 0;
					whiteX[wNum] = ix;
					whiteY[wNum] = iy+1;
					wNum++;
				}
				if (Img[iy*cxDIB+ix-1] == 255) {
					Img[iy*cxDIB+ix-1] = 0;
					whiteX[wNum] = ix-1;
					whiteY[wNum] = iy;
					wNum++;
				}
				if (wNum > numTH) {
					flag = FALSE;
					break;
				}
			}
			if (flag == FALSE) {
				for (k = 0; k < wNum; k++) {
					ix = whiteX[k];
					iy = whiteY[k];
					Img[iy*cxDIB+ix] = 255;
				}
			}
		}
	}
}

void re_get_orient_image(BYTE* OrntImg,BYTE* Img,int cxDIB,int cyDIB)
{
	int dir0, dir1, dir2, dir3, sum0, sum1, sum2, sum3, var0, var1, var2, var3;
	int nMean, nMax, nAvg, nVal, s1, s2;
	int i, j, k, kk, idx1 = 0, idx2 = 0, idx3 = 0;
	int cx2 = cxDIB/2 - 1;
	int *pSum, *sumptr0, *sumptr1, *sumptr2, *sumptr3;
	unsigned short *pAbs[19], *absptr0, *absptr1, *absptr2, *absptr3;
	BYTE *pTmp, *tmpptr0, *tmpptr3, *tmpptr6, *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8;
	BYTE *imgptr, *ptr2, *dstptr1, *dstptr2, *orntptr;

	pTmp = (BYTE*)malloc(sizeof(BYTE)*38*cxDIB);
	if (pTmp == NULL) return;
	pSum = (int*)calloc(2*cxDIB, sizeof(int));
	if (pSum == NULL) { free(pTmp); return; }

	for (i = 0; i < 19; i++) {
		pAbs[i] = (unsigned short*)malloc(sizeof(unsigned short)*2*cxDIB);
		if (pAbs[i] == NULL) { 
			for (j = 0; j < i; j++) free(pAbs[j]); 
			free(pSum); free(pTmp); return;
		}
	}

	tmpptr0 = pTmp + 2;  tmpptr3 = tmpptr0 - cxDIB;  tmpptr6 = tmpptr0 + cxDIB;
	sumptr0 = pSum; sumptr1 = pSum +1; sumptr2 = pSum + 2; sumptr3 = pSum + 3;

	imgptr = Img;  orntptr = OrntImg - 18*cxDIB;

	for (i = 0; i < cyDIB+18; i += 2) {
		k = i % 38;  dstptr1 = pTmp + k*cxDIB;  dstptr2 = dstptr1 + cxDIB;
		if (i < cyDIB) {
			ptr2 = imgptr + cxDIB;
			for (j = 0; j < cxDIB; j++) {
				dstptr1[j] = imgptr[j]; dstptr2[j] = ptr2[j];
			}
			imgptr += 2*cxDIB; 
		}
		if (i < 4) continue;
		if (i < cyDIB) {
			if (idx3 >= 19) { idx3 = 1; tmpptr3 = pTmp + cxDIB + 2; }
			else { idx3++; tmpptr3 += 2*cxDIB; }
			if (idx1 >= 18) { idx1 = 0; tmpptr0 = pTmp + 2; tmpptr6 = tmpptr0 + cxDIB; }
			else { idx1++; tmpptr0 += 2*cxDIB; tmpptr6 += 2*cxDIB; }
			absptr0 = pAbs[idx1]; absptr1 = absptr0 + 1; absptr2 = absptr1 + 1; absptr3 = absptr2 + 1;
			p0 = tmpptr0; p1 = p0 - 1;  p2 = p0 + 1;
			p3 = tmpptr3; p4 = p3 - 1;  p5 = p3 + 1;
			p6 = tmpptr6; p7 = p6 - 1;  p8 = p6 + 1;
			for (k = 1; k < cx2; k++) {
				absptr0[k*4] = abs(*p0 - *p1) + abs(*p0 - *p2);
				absptr1[k*4] = abs(*p0 - *p4) + abs(*p0 - *p8);
				absptr2[k*4] = abs(*p0 - *p3) + abs(*p0 - *p6);
				absptr3[k*4] = abs(*p0 - *p7) + abs(*p0 - *p5);

				sumptr0[k*4] += absptr0[k*4];
				sumptr1[k*4] += absptr1[k*4];
				sumptr2[k*4] += absptr2[k*4];
				sumptr3[k*4] += absptr3[k*4];
				p0 += 2; p1 += 2; p2 += 2; p3 += 2; p4 += 2; p5 += 2; p6 += 2; p7 += 2; p8 += 2;
			}
		}
		if (i < 18) continue;
		if (i >= 36) {
			if (idx2 >= 18) idx2 = 0;
			else idx2++;	
			absptr0 = pAbs[idx2]; absptr1 = absptr0 + 1; absptr2 = absptr1 + 1; absptr3 = absptr2 + 1;
			for (k = 1; k < cx2; k++) {
				sumptr0[k*4] -= absptr0[k*4];
				sumptr1[k*4] -= absptr1[k*4];
				sumptr2[k*4] -= absptr2[k*4];
				sumptr3[k*4] -= absptr3[k*4];
			}
		}
	
		ptr2 = orntptr + i*cxDIB;
		k = 0;  kk = 0;
		sum0 = sum1 = sum2 = sum3 = 0;
		for (j = 0; j < cxDIB+18; j += 2) {
			if (j < cxDIB) {
				sum0 += sumptr0[k];
				sum1 += sumptr1[k];
				sum2 += sumptr2[k];
				sum3 += sumptr3[k];
				k += 4;
			}
			if (j < 18) continue;
			if (j >= 36) {
				sum0 -= sumptr0[kk];
				sum1 -= sumptr1[kk];
				sum2 -= sumptr2[kk];
				sum3 -= sumptr3[kk];
				kk += 4;
			}

			if (*ptr2 == 255) { ptr2 += 2;  continue; }

			nMean = var0 = dir0 = sum0;
			var1 = dir1 = (91 * sum1) >> 7;
			nMean += var1;
			var2 = dir2 = sum2;
			nMean += var2;
			var3 = dir3 = (91 * sum3) >> 7;
			nMean += var3;
			nMean /= 4;
			
			nMax = dir0;
			if (nMax < dir1) nMax = dir1;
			if (nMax < dir2) nMax = dir2;
			if (nMax < dir3) nMax = dir3;

			nAvg = (nMean - nMax) * 4;
			if (nAvg == 0) { *ptr2 = 0x7F; ptr2 += 2; continue; }

			nVal = 45; s1 = var2 + var3; s2 = var2 + var1;
			if (s1 < s2) {
				s2 = s1; nVal = 75;
				dir0 = var1;	dir1 = var2;
				dir2 = var3;	dir3 = var0;
			}
			s1 = var0 + var3;
			if (s1 < s2) {
				s2 = s1; nVal = 105;
				dir0 = var2;	dir1 = var3;
				dir2 = var0;	dir3 = var1;
			}
			s1 = var0 + var1;
			if (s1 < s2) {
				nVal = 15;
				dir0 = var3;	dir1 = var0;
				dir2 = var1;	dir3 = var2;
			}
			nVal += 15 * (3*(dir3 - dir0) - dir1 + dir2) / nAvg;
			if (nVal == 120) nVal = 0;
			*ptr2 = (BYTE)(nVal);  ptr2 += 2;
		}
	}
	p0 = OrntImg; p2 = p0 + cxDIB;
	for (i = 0; i < cyDIB; i += 2) {
		p1 = p0 + 1; p3 = p2 + 1;
		for (j = 0; j < cxDIB; j += 2) {
			p1[j] = p2[j] = p3[j] = p0[j];
		}
		p0 += 2*cxDIB; p2 = p0 + cxDIB;
	}

	for (i = 0; i < 19; i++) free(pAbs[i]);
	free(pSum); free(pTmp);
}

void get_block_data(BYTE* OrntImg,int cxDIB,int cyDIB,BLOCKVECT* pBlock,int nCol,int nRow)
{
	int i, j;
	BYTE *ptr = OrntImg + cxDIB*BLOCK_SIZE/2 + BLOCK_SIZE/2,  *blkptr = pBlock->Data;

	for (i = 0; i < nRow; i++) {
		for (j = nCol-1; j >= 0; j--) {
			blkptr[j] = ptr[j*BLOCK_SIZE];
		}
		ptr += BLOCK_SIZE * cxDIB;  blkptr += nCol;
	}
	pBlock->nCol = nCol; pBlock->nRow = nRow;
}

void get_singular_block(BYTE *OrntImg,int cxDIB,int cyDIB,int *nNum,int *pList,int *typeList)
{	
	BYTE *BlockDir;
	int i, x, y, j, k, dif, poincare, id1, id2;
	int col = cxDIB/8, row = cyDIB/8;
	int tmpList[64], tmpType[64];
	short mX49[8] = { 3,3,0,-3,-3,-3,0,3 };
	short mY49[8] = { 0,-3,-3,-3,0,3,3,3 };
	short mX9[8] = { 1,1,0,-1,-1,-1,0,1 };
	short mY9[8] = { 0,-1,-1,-1,0,1,1,1 };


	BYTE *ptr = OrntImg + cxDIB*8/2 + 8/2, *blkptr;
	BYTE loop[17];

	BlockDir = (BYTE*)malloc(sizeof(BYTE)*col*row);
	if (BlockDir == NULL)  return;

	memset(BlockDir, 255, sizeof(BYTE)*col*row);

	blkptr = BlockDir;

	for (i = 0; i < row; i++) {
		for (j = col-1; j >= 0; j--) {
			blkptr[j] = ptr[j*8];
		}
		ptr += 8*cxDIB;  blkptr += col;
	}

	for (i = 0, j = 0; i < col*row; i++) {
		if (BlockDir[i] == 255) continue;
		y = i / col; 
		if (y-2 < 0 || y+2 >= row) continue;
		x = i % col; 
		if (x-2 < 0 || x+2 >= col) continue;
		loop[0] = BlockDir[i+2]; 
		loop[1] = BlockDir[i-col+2]; 
		loop[2] = BlockDir[i-2*col+2];
		loop[3] = BlockDir[i-2*col+1]; 
		loop[4] = BlockDir[i-2*col]; 
		loop[5] = BlockDir[i-2*col-1];
		loop[6] = BlockDir[i-2*col-2];
		loop[7] = BlockDir[i-col-2];
		loop[8] = BlockDir[i-2];
		loop[9] = BlockDir[i+col-2];
		loop[10] = BlockDir[i+2*col-2];
		loop[11] = BlockDir[i+2*col-1];
		loop[12] = BlockDir[i+2*col];
		loop[13] = BlockDir[i+2*col+1];
		loop[14] = BlockDir[i+2*col+2];
		loop[15] = BlockDir[i+col+2];
		loop[16] = loop[0];
		poincare = 0; 
		for (k = 0; k < 16; k++) {
			if (loop[k] == 255 || loop[k+1] == 255) continue; 
			dif = loop[k] - loop[k+1];
			if (dif <= -60) { poincare += dif + 120;  continue; }
			if (dif >= 60) {  poincare += dif - 120;  continue; }
			poincare += dif;
		}
		if (poincare == 120) {
			tmpType[j] = 1;
			tmpList[j++] = i;
			if (j >= 64) break;
			continue;
		}
		if (poincare == -120) {
			tmpType[j] = 0;
			tmpList[j++] = i;
			if (j >= 64) break;
		}
	}

	for (i = 0; i < j; i++) {
		if (*nNum >= 64) break; 
		pList[*nNum] = tmpList[i]; 
		typeList[*nNum] = tmpType[i]; (*nNum)++;
		y = tmpList[i] / col; 
		x = tmpList[i] % col; 
		for (k = 0; k < 8; k++) {
			id1 = x + mX49[k]; id2 = y + mY49[k];
			if (id1 < 0 || id1 >= col || id2 < 0 || id2 >= row) { 
				id2 = tmpList[i] + mX9[k] + mY9[k]*col;
				if (*nNum >= 32) { i = 1000; break; }
				pList[*nNum] = id2; 
				typeList[*nNum] = tmpType[i]; (*nNum)++;
				break;
			}
			id1 = tmpList[i] + mX49[k] + mY49[k]*col;
			if (BlockDir[id1] == 255) {
				id2 = tmpList[i] + mX9[k] + mY9[k]*col;
				if (*nNum >= 32) { i = 1000; break; }
				pList[*nNum] = id2; 
				typeList[*nNum] = tmpType[i]; (*nNum)++; break;
			}
		}
	}

	free(BlockDir);
}

void check_core_cand(int x,int y,int nDev,int *xList,int *yList,int *dList,int *nCandNum)
{
	int i, dx, dy;

	for (i = 0; i < *nCandNum; i++) {
		dx = x - xList[i]; dy = y - yList[i];
		if (dx*dx+dy*dy < 15*15) break;
	}
	if (i < *nCandNum) {
		if (dList[i] < nDev) {
			xList[i] = x; yList[i] = y; dList[i] = nDev; return;
		}
	}
	else {
		xList[*nCandNum] = x; yList[*nCandNum] = y;
		dList[*nCandNum] = nDev; (*nCandNum)++;
	}
}

int get_deviation(int x,int y,int size,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, j, val, num, nDev = 0, sx, ex, sy, ey, xsize, ysize;
	BYTE p0 = OrntImg[y*cxDIB+x];
	BYTE *ptr, *p;

	sx = (x > size)? x - size : 0;
	ex = ( x + size >= cxDIB )? cxDIB : x + size;
	sy = (y > size)? y - size : 0;
	ey = ( y + size >= cyDIB )? cyDIB : y + size;

	ptr = OrntImg + sy*cxDIB + sx;  xsize = ex - sx;  ysize = ey - sy;
	for (i = 0; i < ysize; i++, ptr+=cxDIB) {
		p = ptr;
		for (j = 0; j < xsize; j++, p++) {
			val = abs(p0 - *p);
			if (val > 60) val = 120 - val;
			nDev += val; 
		}
	}
	num = xsize * ysize;
	if (num == 0) return (-1);
	return (nDev/num);
}

BOOL check_outof_point(int x,int y,int size,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, j;
	BYTE *ptr, *p;

	if (y < size || y >= cyDIB-size || x < size || x >= cxDIB-size) return (TRUE);

	ptr = OrntImg + (y-size)*cxDIB + x - size;
	for (i = 0; i <= 2*size; i++, ptr += cxDIB) {
		p = ptr;
		for (j = 0; j <= 2*size; j++, p++) {
			if (*p >= 120) return (TRUE);
		}
	}
	return (FALSE);
}

int correct_orient_core(int coreDir,int *dSum,int maxV,int minV,int type,int minP2)
{
	int i, i1, i2, j1, dir1, e1, e2;

	if (type == 1 || minP2 < 0 || minP2 >= 240) {
		dir1 = 80;
		i1 = i2 = e1 = e2 = 0;
		for (i = 1; i < 120; i++) {
			if (i1 == 0) {
				j1 = coreDir-i;
				j1 = ANGLE_0(j1);
				if (dSum[j1] > dir1) { i1 = 1; e1 = i; }
			}
			if (i2 == 0) {
				j1 = coreDir+i;
				j1 = ANGLE_240(j1);
				if (dSum[j1] > dir1) { i2 = 1; e2 = i; }
			}
			if (i1 == 1 && i2 == 1) break;
		}
		if (e1 != 0 && e2 != 0) {
			i1 = (e2 - e1)/2;
			coreDir += i1;
			coreDir = ANGLE_0_240(coreDir);
		}
		return (coreDir);
	}
	i1 = i2 = abs(coreDir - minP2);
	i2 = IANGLE_120(i2);
	if (i2 > 50) return (-1);
	coreDir = (coreDir + minP2)/2;
	if (i1 > 120) {
		coreDir += 120;
		coreDir = ANGLE_240(coreDir);
	}
	return (coreDir);
}

int get_orient_core(int ix,int iy,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int i, i1, i2, k, ix1, iy1, ix2, iy2;
	int maxV, minV, minP = 0, sinV, cosV, flag, midV;
	int dDirs[250], dSum[240];
	int *psum, *pdir;

	pdir = dDirs + 2;  psum = pdir + 120;
	for (i = 0; i < 120; i++, pdir++, psum++) {
		cosV = (30 * (int)_table_03[i]) >> 14;
		sinV = (30 * (int)_table_04[i]) >> 14;
		ix1 = ix + cosV; iy1 = iy + sinV;
		flag = 0;
		if (ix1 < 0 || ix1 >= cxDIB || iy1 < 0 || iy1 >= cyDIB) flag = 1;
		else if (OrntImg[iy1*cxDIB+ix1] >= 0x7F) flag = 1;
		if (flag == 1) {
			i1 = op_func_01(cxDIB/2,cyDIB/2,ix,iy);
			i2 = op_func_01(ix1,iy1,ix,iy);
			i1 = abs(i1 - i2);
			i1 = IANGLE_120(i1);
			if (i1 < 45) return (-1);
			ix1 = ix + cosV/4; iy1 = iy + sinV/4;
			if (ix1 < 0 || ix1 >= cxDIB || iy1 < 0 || iy1 >= cyDIB) return (-1);
			if (OrntImg[iy1*cxDIB+ix1] >= 0x7F) return (-1);
		}
		*pdir = abs(OrntImg[iy1*cxDIB+ix1] - i);
		*pdir = IANGLE_60(*pdir);
		flag = 0;
		ix2 = ix - cosV; iy2 = iy - sinV;
		if (ix2 < 0 || ix2 >= cxDIB || iy2 < 0 || iy2 >= cyDIB) flag = 1;
		else if (OrntImg[iy2*cxDIB+ix2] >= 0x7F) flag = 1;
		if (flag == 1) {
			i1 = op_func_01(cxDIB/2,cyDIB/2,ix,iy);
			i2 = op_func_01(ix2,iy2,ix,iy);
			i1 = abs(i1 - i2);
			i1 = IANGLE_120(i1);
			if (i1 < 45) return (-1);
			ix2 = ix - cosV/4; iy2 = iy - sinV/4;
			if (ix2 < 0 || ix2 >= cxDIB || iy2 < 0 || iy2 >= cyDIB) return (-1);
			if (OrntImg[iy2*cxDIB+ix2] >= 0x7F) return (-1);
		}
		*psum = abs(OrntImg[iy2*cxDIB+ix2] - i);
		*psum = IANGLE_60(*psum);
	}
	for (i = 0; i < 2; i++) {
		dDirs[i] = dDirs[239+2-i];
		dDirs[240+2+i] = dDirs[i+2];
	}
	maxV = 0;  minV = 200 * 5;
	psum = dSum;  pdir = dDirs;
	for (i = 0; i < 240; i++, psum++, pdir++) {
		*psum = pdir[0] + pdir[1] + pdir[2] + pdir[3] + pdir[4];
		if (maxV < *psum) { maxV = *psum; }
		if (minV > *psum) { minV = *psum; minP = i; }
	}
	midV = (maxV - minV)/3;
	if (midV < 10) return (-1);
	i = 0;  flag = i1 = i2 = -1;
	for (k = 0; k < 240; k++) {
		ix1 = minP + k;
		ix1 = ANGLE_240(ix1);
		if (flag == -1) {
			if (dSum[ix1] > maxV-midV) {
				i++;
				flag = 1;
				if (i1 >= 0 && i2 == -1) { i2 = ix1; }
			}
			continue;
		}
		if (dSum[ix1] < minV+midV) {
			if (flag == 1 && i2 == -1) { i1 = ix1; }
			flag = -1;
		}
	}
	ix2 = -1;
	if (i == 2) {
		if (i1 >= 0 && i1 < 240 && i2 >= 0 && i2 < 240 && i1 != i2) {
			iy2 = dSum[i1];
			ix2 = i1;
			for (k = 0; (k < 240 && i1 != i2); k++) {
				i1++;
				if (i1 >= 240) i1 = 0;
				if (iy2 > dSum[i1]) { iy2 = dSum[i1]; ix2 = i1; }
			}
		}
		else i = 1;
	}
	if (minV > 200 && i == 3) i = 1;
	if (i == 1 || i == 2) {
		minP = correct_orient_core(minP,dSum,maxV,minV,i,ix2);
		i1 = minP + 120;
		i1 = ANGLE_240(i1);
		if (dDirs[i1+2] < 20) minP = -1;
	}
	else if (i == 3) {
		iy2 = i2-i1;
		iy2 = ANGLE_0(iy2);
		if (iy2 > 60) minP = -1;
		else minP = -2;
	}
	else minP = -1;
	return (minP);
}

void get_core_points(SINGULAR* SingularData,BYTE *OrntImg,int cxDIB,int cyDIB)
{
	int nBlkNum = 0, pList[64], typeList[64];
	int i, k, h, n, sx, sy, id = 0, kk;
	int nDev, num, e1, orn, flag;
	int xlist[20], ylist[20], dlist[20], nCandNum = 0;
	int coreNum, delNum;
	BYTE Block[64], loop[9];
	TEMPCORE Cores[10], Deltas[10], tCore;

	int col = cxDIB / 8;
	BYTE *ptr, *blkptr, *p0, *p1, *p2, *p3, *p4;

	get_singular_block(OrntImg,cxDIB,cyDIB,&nBlkNum,pList,typeList);

	for (i = 0; i < nBlkNum; i++) {
		sy = (pList[i] / col) * 8;  sx = (pList[i] % col) * 8;

		memset( Block, 1, sizeof(BYTE)*64 );

		ptr = OrntImg + sy*cxDIB + sx;  blkptr = Block;
		for (k = 0; k < 8; k++, ptr += cxDIB) {
			p0 = ptr; p1 = p0 - 1; p2 = p0 + 1; p3 = p0 - cxDIB; p4 = p0 +cxDIB;
			for (h = 0; h < 8; h++, p0++, p1++, p2++, p3++, p4++, blkptr++) {
				if (*p0 >= 120) continue;
				if (*p0 > 6 && *p0 < 60-6 && *p0 > 60+6 && *p0 < 120-6) continue;
				if (*p1 >= 120 || *p2 >= 120) continue;
				if (*p3 >= 120 || *p4 >= 120) continue;
				if (*p1 <= 60 && *p2 <= 60) {
					if ((*p3 <= 60 && *p4 <= 60) || (*p3 >= 60 && *p4 >= 60)) continue;
				}
				if (*p1 >= 60 && *p2 >= 60) {
					if ((*p3 <= 60 && *p4 <= 60) || (*p3 >= 60 && *p4 >= 60)) continue;
				}
				*blkptr = 0;
			}
		}
		for (k = 0; k < 8; k++) {
			p0 = Block + k*8; p1 = p0 + 8; p2 = p0 - 8;
			for (h = 0; h < 8; h++, p0++, p1++, p2++ ) {
				if (*p0 != 0) continue;
				if (h < 7) loop[0] = *(p0 + 1);
				else loop[0] = 1;
				if (h < 7 && k < 7)  loop[1] = *(p1 + 1);
				else loop[1] = 1;
				if (k < 7) loop[2] = *(p1 + 0);
				else loop[2] = 1;
				if (k < 7 && h >= 1) loop[3] = *(p1 - 1);
				else loop[3] = 1;
				if (h >= 1) loop[4] = *(p0 - 1);
				else loop[4] = 1;
				if (h >= 1 && k >= 1) loop[5] = *(p2 - 1);
				else loop[5] = 1;
				if (k >= 1) loop[6] = *(p2 + 0);
				else loop[6] = 1;
				if (k >= 1 && h < 7) loop[7] = *(p2 + 1);
				else loop[7] = 1;
				loop[8] = loop[0];
				num = 0;
				for (kk = 0; kk < 9; kk++) {
					if (kk == 8) n = 0;
					else n = kk + 1;
					if (loop[kk] == 0 && loop[n] != 0) num++;
				}
				if (num > 1) continue;

				nDev = get_deviation(sx+h, sy+k, 12, OrntImg, cxDIB, cyDIB);
				if (nDev < 10) continue;

				check_core_cand(sx+h, sy+k, nDev, xlist, ylist, dlist, &nCandNum);
				if (nCandNum >= 10) { k = 100; break; }
			}
		}
		if (k >= 100) break;
	}

	coreNum = delNum = 0;
	for (i = 0; i < nCandNum; i++) {
		e1 = get_orient_core(xlist[i],ylist[i],OrntImg,cxDIB,cyDIB);
		if (e1 >= 0 || e1 == -2) {
			k = xlist[i] / 8; h = ylist[i] / 8;
			for (n = 0; n < nBlkNum; n++) {
				if (pList[n] == col*h+k) { id = n; break; }
			}
		}
		if (e1 < 0) {
			if (e1 == -2 && typeList[id] == 0) {
				Deltas[delNum].x = xlist[i];
				Deltas[delNum].y = ylist[i];
				Deltas[delNum++].val = dlist[i];
			}
			continue;
		}
		if (typeList[id] == 0) continue;
		Cores[coreNum].x = xlist[i]; Cores[coreNum].y = ylist[i];
		Cores[coreNum].dir = e1; Cores[coreNum++].val = dlist[i];
	}
	for (i = 0; i < coreNum-1; i++) {
		for (k = i+1; k < coreNum; k++) {
			if (Cores[k].val == 0) continue;
			if (Cores[k].val > Cores[i].val) {
				tCore = Cores[i]; Cores[i] = Cores[k]; Cores[k] = tCore;
			}
			sx = Cores[k].x - Cores[i].x; sy = Cores[k].y - Cores[i].y;
			if (sx*sx+sy*sy < 100) Cores[k].val = 0;
		}
	}
	e1 = 0;
	for (i = 0; i < coreNum; i++) {
		if (Cores[i].val == 0) break;
		Cores[e1++] = Cores[i];
		if (e1 == 2) break;
	}
	coreNum = e1;
	if (coreNum == 2) {
		flag = 0;
		e1 = Cores[0].dir - Cores[1].dir;
		e1 = ANGLE_0(e1);
		if (e1 < 65 || e1 > 175) flag = 1;
		else {
			orn = op_func_01(Cores[0].x,Cores[0].y,Cores[1].x,Cores[1].y);
			e1 = abs(Cores[0].dir - orn);
			if (e1 > 120) e1 = 240 - e1;
			if (e1 < 60) flag = 1;
		}
		if (flag == 1) {
			coreNum = 1;
			if (Cores[1].val > 30) {
				if (check_outof_point(Cores[1].x,Cores[1].y,15,OrntImg,cxDIB,cyDIB) == FALSE) {
					if (check_outof_point(Cores[0].x,Cores[0].y,15,OrntImg,cxDIB,cyDIB) == TRUE) {
						Cores[0] = Cores[1];
					}
				}
			}
		}
		else {
			e1 = (Cores[0].x - Cores[1].x) * (Cores[0].x - Cores[1].x);
			e1 += (Cores[0].y - Cores[1].y) * (Cores[0].y - Cores[1].y);
			if (e1 < 30*30) {
				orn = op_func_01(Cores[1].x,Cores[1].y , Cores[0].x,Cores[0].y);
				Cores[0].dir = orn;
				orn += 120;
				orn = ANGLE_240(orn);
				Cores[1].dir = orn;
			}
		}
	}
	for (i = 0; i < coreNum; i++) {
		SingularData->nX[i] = Cores[i].x;
		SingularData->nY[i] = Cores[i].y;
		SingularData->nOrient[i] = Cores[i].dir;
		SingularData->nType[i] = 1;
	}
	SingularData->nNumber = coreNum;
	// for Delta
	for (i = 0; i < delNum; i++) {
		for (k = i+1; k < delNum; k++) {
			if (Deltas[k].val == 0) continue;
			if (Deltas[k].val > Deltas[i].val) {
				tCore = Deltas[i]; Deltas[i] = Deltas[k]; Deltas[k] = tCore;
			}
			sx = Deltas[i].x - Deltas[k].x; sy = Deltas[i].y - Deltas[k].y;
			if (sx*sx+sy*sy < 80*80) Deltas[k].val = 0;
		}
	}
	for (i = 0; i < delNum; i++) {
		if (coreNum == 1) {
			orn = op_func_01(Deltas[i].x,Deltas[i].y , Cores[0].x,Cores[0].y);
			e1 = abs(Cores[0].dir-orn);
			if (e1 > 120) e1 = 240 - e1;
			if (e1 > 60) Deltas[i].val = 0;
		}
	}
	e1 = 0;
	for (i = 0; i < delNum; i++) {
		if (Deltas[i].val == 0) continue;
		Deltas[e1++] = Deltas[i];
		if (e1 == 2) break;
	}
	delNum = e1;
	for (i = 0; i < delNum; i++) {
		SingularData->nX[coreNum+i] = Deltas[i].x;
		SingularData->nY[coreNum+i] = Deltas[i].y;
		SingularData->nOrient[coreNum+i] = 255;
		SingularData->nType[coreNum+i] = DELTA;
	}
	SingularData->nNumber += delNum;
}

BOOL check_whorl(SINGULAR* SData,MAINLINE* mLine)
{
	int i, n, n1, n2, cx, cy, dir, dir1, dx ,dy;
	int x0, y0, x1, x2, y1, y2, d1, d2, dd, minV;
	int cross_p[4], index_p[4];
	int m1, m2, m3, m4;

	x1 = SData->nX[0]; x2 = SData->nX[1];
	y1 = SData->nY[0]; y2 = SData->nY[1];
	cx = (x1 + x2) / 2; cy = (y1 + y2) / 2;

	dx = x2 - x1; dy = y2 - y1;
	d2 = dx*dx + dy*dy;
	d1 = op_func_02(d2);
	if (d1 == 0) return(TRUE);
	d2 = dx*cx + dy*cy;
	dir = op_func_01(SData->nX[1],SData->nY[1],SData->nX[0],SData->nY[0]);
	dir += 60;
	dir = ANGLE_240(dir);
	for (n = 0; n < 4; n++) {
		cross_p[n] = -1000; index_p[n] = 0;
		n1 = mLine->nNumbers[n];
		minV = 10;
		for (i = 0; i < n1; i++) {
			x0 = mLine->points_x[n][i]; y0 = mLine->points_y[n][i];
			dd = dx*x0 + dy*y0 - d2;
			dd = abs(dd) / d1;
			if (minV > dd) {
				minV = dd;
				dir1 = op_func_01(x0,y0,cx,cy);
				dir1 = abs(dir1 - dir);
				dir1 = IANGLE_120(dir1);
				n2 = (x0-cx)*(x0-cx) + (y0-cy)*(y0-cy);
				dd = op_func_02(n2);
				if (dir1 > 60) dd = -dd;
				if (index_p[n] != 0 && i-index_p[n] > 3) break;
				cross_p[n] = dd; index_p[n] = i;
				if (minV == 0) break;
			}
		}
	}
	d2 = 10;
	if (cross_p[0]>-900 && cross_p[2]>-900 && cross_p[1]>-900 && cross_p[3]>-900) {
		m1 = m2 = cross_p[0];
		if (m1 > cross_p[1]) m1 = cross_p[1];
		if (m2 < cross_p[1]) m2 = cross_p[1];
		m3 = m4 = cross_p[2];
		if (m3 > cross_p[3]) m3 = cross_p[3];
		if (m4 < cross_p[3]) m4 = cross_p[3];
		if (m1 > m4 || m3 > m2) return(FALSE);
		if (d1 > 80) {
			m1 = abs(cross_p[0] - cross_p[1]);
			m2 = abs(cross_p[2] - cross_p[3]);
			if (m1 > m2) m1 = m2;
			d2 = m1 / 2;
			dd = abs(cross_p[0] - cross_p[2]);
			if (dd < d2) return (FALSE);
			dd = abs(cross_p[1] - cross_p[3]);
			if (dd < d2) return (FALSE);
		}
		m1 = abs(cross_p[0] - cross_p[1]);
		m2 = abs(cross_p[2] - cross_p[3]);
		if (m1 > m2) m1 = m2;
		d2 = m1 / 3;
	}
	if (cross_p[0] > -900 && cross_p[2] > -900) {
		dd = abs(cross_p[0] - cross_p[2]);
		if (dd < d2) return (FALSE);
	}
	if (cross_p[1] > -900 && cross_p[3] > -900) {
		dd = abs(cross_p[1] - cross_p[3]);
		if (dd < d2) return (FALSE);
	}
	return (TRUE);
}

int get_distance_to_line(MAINLINE *mLine,int x,int y,int num)
{
	int i, dx, dy, len, maxlen = 10000;

	for (i = 0; i < mLine->nNumbers[num]; i++) {
		dx = x - mLine->points_x[num][i];
		dy = y - mLine->points_y[num][i];
		len = dx*dx + dy*dy;
		if (len < maxlen) maxlen = len;
	}
	maxlen = op_func_02(maxlen);
	return (maxlen);
}

int check_near_line(MAINLINE *mLine,int num,int disTh)
{
	int i, x, y, dd, nn = num+1, number = 1000;

	for (i = 10; i < mLine->nNumbers[nn]; i++) {
		x = mLine->points_x[nn][i]; y = mLine->points_y[nn][i];
		dd = get_distance_to_line(mLine,x,y,num);
		if (dd < disTh) { number = i; break; }
	}
	return (number);
}

int check_line_lr(int Line_x0,int Line_y0,int Line_x1,int Line_y1,
						short *points_x, short *points_y, int nNumber)
{
	int i, dx, dy, dd, numPlus = 0, numMis = 0;

	dx = Line_x1 - Line_x0; dy = Line_y1 - Line_y0;
	for (i = 0; i < nNumber; i++) {
		dd = dy*(points_x[i]-Line_x0) - (points_y[i]-Line_y0)*dx;
		if (dd < 0) numMis++;
		if (dd > 0) numPlus++;
	}
	if (numMis > 0 && numPlus == 0) return (-1);
	if (numMis == 0 && numPlus > 0) return (1);
	return (0);
}

BOOL check_arch(MAINLINE *mLine,COREITEMEX *Core,COREITEMEX *Delta)
{
	int i, j, ix, iy, ix0, iy0, d1, d2, dd, coreD, disTH = 15*15;

	ix = Core->x - Delta->x; iy = Core->y - Delta->y;
	coreD = op_func_02(ix*ix+iy*iy);

	ix0 = mLine->points_x[0][mLine->nNumbers[0]-1] - Core->x;
	iy0 = mLine->points_y[0][mLine->nNumbers[0]-1] - Core->y;
	ix = mLine->points_x[1][mLine->nNumbers[1]-1] - Core->x;
	iy = mLine->points_y[1][mLine->nNumbers[1]-1] - Core->y;
	d1 = ix0*ix0 + iy0*iy0; d2 = ix*ix + iy*iy;
	dd = d1;
	if (dd > d2) dd = d2;
	dd = op_func_02(dd);

	if (coreD < dd && mLine->nNumbers[0] > 20 && mLine->nNumbers[1] > 20) {
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
		if (d1*d2 < 0) return TRUE;
	}
	for (i = 0; i < 2; i++) {
		for (j = 0; j < mLine->nNumbers[i]; j++) {
			ix = mLine->points_x[i][j] - Delta->x;
			iy = mLine->points_y[i][j] - Delta->y;
			dd = ix*ix + iy*iy;
			if (dd < disTH-10) return TRUE;
		}
	}
	return FALSE;
}

BYTE get_type_line(LPFPVECTEX pVectEx,SINGULAR* SData,BYTE* OrntImg,BYTE* Img,int cxDIB,int cyDIB)
{
	int i, i1, i2, j, x, y, dir, ix, iy, idir, ix0, iy0;
	int k, kk, dd, d1, d2, th, dd0 = 1;
	int directX[8], directY[8], dth = 4;
	int dx, dy, xx, yy, diff01, diff1, diff2, sDir[2];
	int co_x = 0, co_y = 0, co_dir = 0, flag, circleFlag;
	int nCoreNumber = 0, nDeltaNumber = 0;
	int ccx, ccy, ddx, ddy;
	COREITEMEX Core[10], Delta[10];
	MAINLINE mLine;
	BYTE byteP;
	
	for (i = 0; i < 4; i++) {
		pVectEx->MainLine.points_x[i] = 0;
		pVectEx->MainLine.points_y[i] = 0;
	}
	for (i = 0; i < SData->nNumber; i++) {
		if (SData->nType[i] == 1) {
			Core[nCoreNumber].x = SData->nX[i];
			Core[nCoreNumber].y = SData->nY[i];
			Core[nCoreNumber++].dir = (BYTE)SData->nOrient[i];
		}
		else {
			Delta[nDeltaNumber].x = SData->nX[i];
			Delta[nDeltaNumber].y = SData->nY[i];
			Delta[nDeltaNumber++].dir = (BYTE)SData->nOrient[i];
		}
	}
	if (nCoreNumber == 0) return (10);
	if (nCoreNumber == 2) {
		ix = (Core[0].x - Core[1].x);
		iy = (Core[0].y - Core[1].y);
		ix = ix*ix + iy*iy;
		if (ix < 50*50) return (1);
	}
	if (nCoreNumber == 1 && nDeltaNumber == 1) {
		ix = (Core[0].x - Delta[0].x);
		iy = (Core[0].y - Delta[0].y);
		ix = ix*ix + iy*iy;
		if (ix < 50*50) dth = 8;
	}
	mLine.nNumbers[0] = mLine.nNumbers[1] = 0;
	mLine.nNumbers[2] = mLine.nNumbers[3] = 0;
	circleFlag = 0;
	
	for (i = 0; i < cxDIB*cyDIB; i++) Img[i] = 0;

	for (i = 0; i < nCoreNumber; i++) {
		co_dir = SData->nOrient[i];
		co_x = SData->nX[i]; co_y = SData->nY[i];
		x = co_x - ((dth * (int)_table_03[co_dir]) >> 14);
		y = co_y - ((dth * (int)_table_04[co_dir]) >> 14);
		sDir[1] = co_dir + 60;
		sDir[1] = ANGLE_240(sDir[1]);
		if (sDir[1] == 240) sDir[1] = 0;
		sDir[0] = co_dir - 60;
		sDir[0] = ANGLE_0(sDir[0]);
		if (sDir[0] == 240) sDir[0] = 0;
		kk = i*2;
		for (j = 0; j < 2; j++) {
			idir = sDir[j];
			ix = x + _5cos_table[idir];
			iy = y + _5sin_table[idir];
			for (k = 0; k < 8; k++) { directX[k] = x; directY[k] = y; }
			xx = 64*ix; yy = 64*iy;
			k = flag = 0;
			while (1) {
				dd = dd0;
				while (1) {
					dx = (64*dd * (int)_table_03[idir]) >> 14;
					dy = (64*dd * (int)_table_04[idir]) >> 14;
					if (dx>=64 || dy>=64 || dx<=-64 || dy<=-64) break;
					dd++;
				}
				ix0 = ix; iy0 = iy; 
				xx += dx; yy += dy;
				if (xx < 0 || yy < 0) break;
				ix = xx / 64; iy = yy / 64;
				for (i1 = 0; i1 < 7; i1++) {
					directX[i1] = directX[i1+1]; directY[i1] = directY[i1+1];
				}
				directX[7] = ix; directY[7] = iy;
				if (ix < 0 || ix >= cxDIB || iy < 0 || iy >= cyDIB) break;
				if (OrntImg[iy*cxDIB+ix] >= 120) break;
				if (Img[iy*cxDIB+ix] == kk+j+1) {
					if (k > 10 && nCoreNumber == 1) circleFlag = 1;
					break;
				}
				mLine.points_x[kk+j][k] = ix;
				mLine.points_y[kk+j][k] = iy;
				k++;
				if (k > 95) break;
				idir = op_func_01(ix,iy,co_x,co_y);
				idir = abs(idir - co_dir);
				idir = IANGLE_120(idir);
				dx = co_x - ix; dy = co_y - iy;
				dd = dx*dx + dy*dy;
				if (flag == 0) {
					if (idir < 60 && dd > 25*25) flag = 1;
				}
				else if (idir > 110 || dd < 20*20) {
					circleFlag = 1;
					if (k > 10) {
						ix0 = mLine.points_x[kk+j][k-3];
						iy0 = mLine.points_y[kk+j][k-3];
						idir = op_func_01(ix0,iy0,co_x,co_y);
						idir = idir - co_dir;
						idir = ANGLE_0(idir);
						if (idir < 120) i1 = 1;
						else i1 = 0;
						if ((idir > 60 && idir < 180) && j == i1) circleFlag = 2;
					}
					break;
				}
				Img[iy*cxDIB+ix] = (BYTE)(kk+j+1);
				byteP = OrntImg[iy*cxDIB+ix];
				i1 = abs(byteP - OrntImg[iy0*cxDIB+ix0]);
				i1 = IANGLE_60(i1);
				if (k > 20 && i1 > 45) { circleFlag = 3; break; }
				idir = op_func_01(directX[7],directY[7],directX[0],directY[0]);
				i1 = abs(idir - byteP);
				idir = byteP;
				i1 = IANGLE_120(i1);
				if (k > 10 && i1 > 45 && i1 < 75) {
					circleFlag = 3; break;
				}
				if (i1 > 60) idir += 120;
			}
			mLine.nNumbers[kk+j] = k;
		}
		if (mLine.nNumbers[kk] > 30 && mLine.nNumbers[kk+1] > 30) {
			ix = (mLine.points_x[kk][30] + mLine.points_x[kk+1][30])/2;
			iy = (mLine.points_y[kk][30] + mLine.points_y[kk+1][30])/2;
			idir = op_func_01(ix,iy,x,y);
			SData->nOrient[i] = idir; Core[i].dir = idir;
			co_dir = idir;
		}
		else {
			k = mLine.nNumbers[kk];
			if (k > mLine.nNumbers[kk+1]) k = mLine.nNumbers[kk+1];
			if (k >= 15) {
				ix = (mLine.points_x[kk][k-1] + mLine.points_x[kk+1][k-1])/2;
				iy = (mLine.points_y[kk][k-1] + mLine.points_y[kk+1][k-1])/2;
				idir = op_func_01(ix,iy,x,y);
				SData->nOrient[i] = idir; Core[i].dir = idir;
				co_dir = idir;
			}
		}
	}

	for (i = 0; i < 4; i++) {
		pVectEx->MainLine.points_x[i] = 0;
		pVectEx->MainLine.points_y[i] = 0;
		if (mLine.nNumbers[i] > 0) {
			pVectEx->MainLine.points_x[i] = mLine.points_x[i][mLine.nNumbers[i]-1];
			pVectEx->MainLine.points_y[i] = mLine.points_y[i][mLine.nNumbers[i]-1];
		}
	}


	if (nCoreNumber == 2) {
		if (mLine.nNumbers[0] < 10 || mLine.nNumbers[1] < 10 ||
			 mLine.nNumbers[2] < 10 || mLine.nNumbers[3] < 10) return(3);
		dir = op_func_01(Core[1].x,Core[1].y,Core[0].x,Core[0].y);
		diff1 = abs(dir - SData->nOrient[0]);
		diff1 = ANGLE_120(diff1);
		diff1 = IANGLE_60(diff1);
		diff2 = abs(dir - SData->nOrient[1]);
		diff2 = ANGLE_120(diff2);
		diff2 = IANGLE_60(diff2);
		dir = (diff1+diff2)/2;
		if (dir >= 17) return (2);
		if (check_whorl(SData,&mLine) == TRUE) return (0);
		else return (2);
	}

	if (mLine.nNumbers[0] < 40 || mLine.nNumbers[1] < 40) {
		if (circleFlag > 0) {
			if (mLine.nNumbers[0] < 40 && mLine.nNumbers[1] < 40) return(3);
			if (circleFlag == 1) return (0);
			if (circleFlag == 2) return (1);
			return (3);
		}
		return (8);
	}
	if (circleFlag > 0) {
		kk = 0;
		for (i = 0; i < 2; i++) {
			for (j = 0; j < mLine.nNumbers[i]; j++) {
				dx = mLine.points_x[i][j] - co_x;
				dy = mLine.points_y[i][j] - co_y;
				diff1 = dx*dx + dy*dy;
				if (kk < diff1) kk = diff1;
			}
		}
		if (kk < 50*50) return (1);
		if (circleFlag == 1) return (0);
		if (circleFlag == 3) return (3);
		return (2);
	}
	if (nDeltaNumber == 2) return (0);
	ix0 = mLine.points_x[0][mLine.nNumbers[0]-1];
	iy0 = mLine.points_y[0][mLine.nNumbers[0]-1];
	ix = mLine.points_x[1][mLine.nNumbers[1]-1];
	iy = mLine.points_y[1][mLine.nNumbers[1]-1];
	dx = (ix0-ix)*(ix0-ix); dy = (iy0-iy)*(iy0-iy);
	diff01 = op_func_02(dx+dy);
	i1 = op_func_01(ix0,iy0,co_x,co_y);
	i2 = op_func_01(ix,iy,co_x,co_y);
	diff2 = abs(i1 - i2);
	diff2 = IANGLE_120(diff2);

	if (mLine.nNumbers[0] > 20 && mLine.nNumbers[1] > 20 && diff01 > 80 && diff2 > 50) {
		diff1 = co_dir - i1;
		diff1 = ANGLE_0(diff1);
		if (diff1 > 20 && diff1 < 60) {
			diff1 = i2 - co_dir;
			diff1 = ANGLE_0(diff1);
			if (diff1 > 20 && diff1 < 60) return (7);
		}
	}
	dx = (ix0 - co_x)*(ix0 - co_x); dy = (iy0 - co_y)*(iy0 - co_y);
	d1 = op_func_02(dx+dy);
	dx = (ix - co_x)*(ix - co_x); dy = (iy - co_y)*(iy - co_y);
	d2 = op_func_02(dx+dy);
	flag = 0;
	ccx = Core[0].x; ccy = Core[0].y;
	ddx = Delta[0].x; ddy = Delta[0].y;
	if (nDeltaNumber == 1) {
		dir = op_func_01(ddx,ddy,ccx,ccy);
		idir = abs(dir - co_dir);
		idir = IANGLE_120(idir);
		if (idir > 80) {
			dd = d1;
			if (dd < d2) dd = d2;
			if (dd < 120) return (0);
			nDeltaNumber = 0;
		}
		if (nDeltaNumber == 1) {
			if (idir < 10) return (7);
			idir = (co_x-ddx)*(co_x-ddx) + (co_y-ddy)*(co_y-ddy);
			idir = op_func_02(idir);
			if (idir >= 215) return (0);
			if (nCoreNumber == 1) {
				if (check_arch(&mLine,&Core[0],&Delta[0]) == TRUE) return (7);
			}
			dir -= co_dir;
			dir = ANGLE_0(dir);
			if (dir >= 120) {
				if (dir >= 130) flag = 4;
			}
			else if (dir < 110) flag = 5;
		}
	}
	if (d1 < 50 || d2 < 50) return (8);
	dd = check_near_line(&mLine,0,15);
	if (dd < 1000) {
		dx = mLine.points_x[1][dd] - co_x;
		dy = mLine.points_y[1][dd] - co_y;
		idir = dx*dx + dy*dy;
		dx = op_func_02(idir);
		dy = d1;
		if (dy > d2) dy = d2;
		dy = op_func_02(dy);
		if (dx < dy-20) return (0);
	}
	dd = d1;
	if (dd < d2) dd = d2;
	if (dd < 100 && diff01 < 16) return (8);
	th = (200 - dd)/30 + 5;
	if (th > 10) th = 10;
	else if (th < 5) th = 5;
	idir = i1 - i2;
	if (dd < 100) {
		idir = abs(i1-i2);
		idir = IANGLE_120(idir);
		if (idir < th/2) return (8);
	}
	dir = co_dir;
	if (flag != 5) {
		idir = i1 - dir;
		idir = ANGLE_0(idir);
		if (idir < 120) {
			if (idir > 10 && idir < 120) {
				idir = abs(dir-i2);
				idir = IANGLE_120(idir);
				if (idir > 70) return (2);
			}
			idir = abs(dir-i2);
			idir = IANGLE_120(idir);
			if (idir < 4) return (0);
			if (idir < 15 && d2 <= 80) return (8);
			return (4);
		}
		if (flag == 4) {
			dx = (co_x-ix0)*(ddx-ix0)+(co_y-iy0)*(ddy-iy0);
			dy = (co_x-ix)*(ddx-ix)+(co_y-iy)*(ddy-iy);
			if (dx > 0 && dy > 0) return (4);
		}
	}
	if (flag != 4) {
		idir = dir - i2;
		idir = ANGLE_0(idir);
		if (idir < 120) {
			if (idir > 10 && idir < 120) {
				idir = abs(dir-i1);
				idir = IANGLE_120(idir);
				if (idir > 70) return (2);
			}
			idir = abs(i1-dir);
			idir = IANGLE_120(idir);
			if (idir < 4) return (0);
			if (idir < 15 && d1 <= 80) return (8);
			return (5);
		}
		if (flag == 5) {
			dx = (co_x-ix0)*(ddx-ix0)+(co_y-iy0)*(ddy-iy0);
			dy = (co_x-ix)*(ddx-ix)+(co_y-iy)*(ddy-iy);
			if (dx > 0 && dy > 0) return (5);
		}
	}
	return (8);
}

void copy_core(SINGULAR* SingularData,LPFPVECTEX FPEx)
{
	int i, num = 0;

	for (i = 0 ; i < SingularData->nNumber; i++) {
		FPEx->Core.item[num].x = SingularData->nX[i];
		FPEx->Core.item[num].y = SingularData->nY[i];
		FPEx->Core.item[num].dir = (BYTE)SingularData->nOrient[i];
		FPEx->Core.item[num].kind = (BYTE)SingularData->nType[i];
		num++;
		if (num >= MAX_CORE_NUMBER) break;
	}
	FPEx->Core.nNumber = num;
}

int get_frequency_sub(int point_x,int point_y,BYTE* Img,
						BYTE* image_buffer3,int cxDIB,int cyDIB)
{
	int i1, x0, x01, x1, x2, y0, y01, y1, y2;
	int fre = 0, dir, cosV, sinV;
	int color, color0, num, num1, numth = 3;

	dir = image_buffer3[point_y*cxDIB+point_x] + 60;
	dir = ANGLE_120(dir);
	cosV = _table_03[dir]; sinV = _table_04[dir];
	color = Img[point_y*cxDIB+point_x];
	x1 = point_x; y1 = point_y;
	i1 = 1;
	while (1) {
		x2 = x1 + (i1*cosV >> 14);
		y2 = y1 + (i1*sinV >> 14);
		if (x2 < 0 || x2 >= cxDIB || y2 < 0 || y2 >= cyDIB) return(0);
		if (Img[y2*cxDIB+x2] != color) break;
		i1++;
	}
	color0 = color = Img[y2*cxDIB+x2];
	x0 = x1 = x2;	y0 = y1 = y2;
	num = 0; i1 = 1;
	while (1) {
		x2 = x1 + (i1*cosV >> 14);
		y2 = y1 + (i1*sinV >> 14);
		if (x2 < 0 || x2 >= cxDIB || y2 < 0 || y2 >= cyDIB) break;
		if (Img[y2*cxDIB+x2] != color) {
			color = Img[y2*cxDIB+x2];
			if (color == color0) {
				num++;
				x0 = x2;	y0 = y2;
				if (num == numth) break;
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
	while (1) {
		x2 = x1 - (i1*cosV >> 14);
		y2 = y1 - (i1*sinV >> 14);
		if (x2 < 0 || x2 >= cxDIB || y2 < 0 || y2 >= cyDIB) break;
		if (Img[y2*cxDIB+x2] != color) {
			color = Img[y2*cxDIB+x2];
			if (color == color0) {
				num++;
				x01 = x2;	y01 = y2;
				if (num == numth) break;
			}
		}
		i1++;
	}
	num += num1;
	if (num == 0) return(0);

	x1 = x0-x01; x1 = x1*x1;
	y1 = y0-y01; y1 = y1*y1;
	x2 = op_func_02(x1+y1);
	fre = x2*6/num;
	return(fre);
}

BYTE get_frequency(BYTE* Img,BYTE* image_buffer3,SINGULAR *pCore,int cxDIB,int cyDIB)
{
	int i, i1, j, j1, cx, cy, hh, x1, x2, y1, y2;
	int ih = 15;
	int frequency = 0, fre[9];
	BOOL flag;

	cx = cxDIB/2; cy = cyDIB/2; hh = cxDIB/4;
	while (frequency == 0) {
		for (i = 0; i < 9; i++) {
			x1 = cx; y1 = cy;
			if (i < 3) x1 -= hh;
			if (i > 5) x1 += hh;
			if (i == 0 || i == 3 || i == 6) y1 -= hh;
			if (i == 2 || i == 5 || i == 8) y1 += hh;
			fre[i] = 0;
			flag = TRUE;
			for (j = -1; j < 2; j++) {
				x2 = x1+j*ih;
				if (x2 < 0 || x2 >= cxDIB) {
					flag = FALSE;
					break;
				}
				for (j1 = -1; j1 < 2; j1++) {
					y2 = y1+j1*ih;
					if (y2 < 0 || y2 >= cyDIB) {
						flag = FALSE;
						break;
					}
					if ((image_buffer3[y2*cxDIB+x2]&0x80) != 0) {
						flag = FALSE;
						break;
					}
				}
				if (flag == FALSE) break;
			}
			if (flag == FALSE) continue;
			for (j = 0; j < pCore->nNumber; j++) {
				if (x1-ih < pCore->nX[j] && x1+ih > pCore->nX[j] &&
					y1-ih < pCore->nY[j] && y1+ih > pCore->nY[j]) {
					flag = FALSE;
					break;
				}
			}
			if (flag == FALSE) continue;
			fre[i] = get_frequency_sub(x1,y1,Img,image_buffer3,cxDIB,cyDIB);
		}
		
		i1 = 0;
		for (i = 0; i < 9; i++) {
			if (fre[i] == 0) continue;
			frequency += fre[i];
			i1++;
		}
		if (i1 > 0) frequency /= i1;
		cy += 30;
		if (cy+hh >= cyDIB) break;
	}
	if (frequency >= 120) return (120);
	if (frequency < 0) return (0);
	return ((BYTE)frequency);
}

void get_point_section(int x1,int y1,int x2,int y2,LPSECTION pSection)
{
	int	x, y, sx, sy, ex, ey, dx, dy, err, c1, c2, num = 0;
	
	sx = x1; sy = y1; ex = x2; ey = y2;
	if (sx > ex) {
		x = sx;	sx = ex; ex = x; y = sy; sy = ey; ey = y;
	}	
	dx = ex - sx;	dy = ey - sy; x = sx; y = sy;
	if (dy > 0) {
		if (dx >= dy) {
			c1 = 2*dy; err = c1-dx; c2 = err-dx;
			while (x <= ex) {
				pSection->point[num].x = x;
				pSection->point[num++].y = y;
				x++;
				if (err < 0) err = err + c1;
				else { y++; err = err + c2; }
			}
		}
		else {
			c1 = 2*dx; err = c1-dy; c2 = err-dy;
			while (y <= ey) {
				pSection->point[num].x = x;
				pSection->point[num++].y = y;
				y++;
				if (err < 0) err = err + c1;
				else { x++; err = err + c2; }
			}
		}
	}
	else {
		dy = -dy;
		if (dy > dx) {
			c1 = 2*dx; err = c1-dy; c2 = err-dy;
			while (y >= ey) {
				pSection->point[num].x = x;
				pSection->point[num++].y = y;
				y--;
				if (err < 0) err = err + c1;
				else { x++; err = err + c2; }
			}
		}
		else {
			c1 = 2*dy; err = c1-dx; c2 = err-dx;
			while (x <= ex) {
				pSection->point[num].x = x;
				pSection->point[num++].y = y;
				x++;
				if (err < 0) err = err + c1;
				else { y--; err = err + c2; }
			}
		}
	}
	pSection->nNumber = num;
}

void filter_section(LPPOINTTAG pPoint,int dir,BOOL nFlag,LPSECTION pSection,BYTE* Img,int cxDIB,int cyDIB)
{
	int	i, x1, x2, y1, y2, angle = dir - 60;
	int val[MAX_SECTION_NUMBER];
	
	angle = ANGLE_0(angle);
	
	if (nFlag != FALSE) {
		x1 = pPoint->x + _9cos_table[angle];  y1 = pPoint->y + _9sin_table[angle];
		x2 = pPoint->x - _9cos_table[angle];  y2 = pPoint->y - _9sin_table[angle];
		
		get_point_section(x1, y1, x2, y2, pSection);
		
		for (i = pSection->nNumber-1; i >= 0; i--) {
			if (pSection->point[i].x < 0 || pSection->point[i].x >= cxDIB || pSection->point[i].y < 0 || pSection->point[i].y >= cyDIB)
				val[i] = 255;
			else 
				val[i] = Img[pSection->point[i].y*cxDIB+pSection->point[i].x];
		}
		
		for (i = pSection->nNumber-3; i >= 2; i--) {
			pSection->nValue[i] = (BYTE)((val[i-2] + 4*val[i-1] + 6*val[i] + 4*val[i+1] + val[i+2]) >> 4);
		}
	}
	else {
		x1 = pPoint->x + _5cos_table[angle];  y1 = pPoint->y + _5sin_table[angle];
		x2 = pPoint->x - _5cos_table[angle];  y2 = pPoint->y - _5sin_table[angle];
		
		get_point_section(x1, y1, x2, y2, pSection);
		
		for (i = pSection->nNumber-1; i >= 0; i--) {
			if (pSection->point[i].x < 0 || pSection->point[i].x >= cxDIB || pSection->point[i].y < 0 || pSection->point[i].y >= cyDIB)
				val[i] = 255;
			else 
				val[i] = Img[pSection->point[i].y*cxDIB+pSection->point[i].x];
		}
		
		for (i = pSection->nNumber-2; i >= 1; i--) {
			pSection->nValue[i] = (BYTE)((val[i-1] + 6*val[i] + val[i+1]) >> 3);
		}
	}
}

int get_max_index(LPPOINTTAG pPoint,LPSECTION pSection,int nInterval)
{
	int i, len1, len2, maxlen = 10000, maxid = -1;
	int	localmaxid[MAX_SECTION_NUMBER], localmaxnum = 0;
	int *ptr;
	
	for (i = nInterval+1; i < pSection->nNumber-nInterval-1; i++) {
		if (pSection->nValue[i-1] < pSection->nValue[i]) continue;
		if (pSection->nValue[i+1] < pSection->nValue[i]) continue;
		localmaxid[localmaxnum++] = i;
	}
	if (localmaxnum <= 0) return (-1);
	ptr = localmaxid;
	for (i = localmaxnum-1; i >= 0; i--, ptr++) {
		len1 = (pPoint->x - pSection->point[*ptr].x)*(pPoint->x - pSection->point[*ptr].x);
		len2 = (pPoint->y - pSection->point[*ptr].y)*(pPoint->y - pSection->point[*ptr].y);
		len1 += len2;
		if (maxlen >= len1) {
			maxlen = len1; maxid = *ptr;
		}
	}
	return (maxid);
}

POINTTAG get_local_maximum(LPPOINTTAG pPoint,int dir,BOOL nInitFlag,
						 BYTE *Img,USHORT *image_buffer1,int cxDIB,int cyDIB)
{
	POINTTAG null = {NULLVALUE,NULLVALUE};
	SECTION section;
	BYTE byteP, nByte;
	BOOL nFlag = FALSE;
	int i, maxid, size, sid, eid, x, y;
	
	filter_section(pPoint,dir,TRUE,&section,Img,cxDIB,cyDIB);

	maxid = get_max_index(pPoint,&section,2);

	if (maxid != -1) {
		x = section.point[maxid].x;   y = section.point[maxid].y;
		if (x < 0 || x >= cxDIB || y < 0 || y >= cyDIB) return (null);
		if (nInitFlag != FALSE) {
			if (image_buffer1[y*cxDIB+x] == 0) return (section.point[maxid]);
		}
		size = abs(maxid - section.nNumber/2);
		if (size <= 2) return (section.point[maxid]);
		if (size == 3) {
			byteP = Img[pPoint->y*cxDIB+pPoint->x];
			sid = maxid; eid = section.nNumber/2;
			if (sid > eid) { sid = eid; eid = maxid; }
			for (i = sid+1; i < eid; i++) {
				nByte = Img[section.point[i].y*cxDIB+section.point[i].x];
				if (nByte > byteP || nByte > 40) { nFlag = TRUE; break; }
			}
			if (nFlag == FALSE) return (section.point[maxid]);
		}
	}

	filter_section(pPoint,dir,FALSE,&section,Img,cxDIB,cyDIB);

	maxid = get_max_index(pPoint,&section,1);
	if (maxid == -1) return (null);
	if (section.point[maxid].x < 0 || section.point[maxid].x >= cxDIB 
	  || section.point[maxid].y < 0 || section.point[maxid].y >= cyDIB) 
		return (null);
	return (section.point[maxid]);
}

void labelling(LPPOINTTAG pPoint1,LPPOINTTAG pPoint2,USHORT* image_buffer1,int cxDIB,int cyDIB,USHORT *nLabelNum)
{
	SECTION section;
	int i, dx, dy, idx;

	get_point_section(pPoint1->x,pPoint1->y,pPoint2->x,pPoint2->y,&section);

	dx = abs(pPoint1->x - pPoint2->x);  dy = abs(pPoint1->y - pPoint2->y);
	dx -= dy;

	for (i = section.nNumber-1; i >= 0; i--) {
		if (section.point[i].y <= 0 || section.point[i].y >= cyDIB-1)  continue;
		if (section.point[i].x <= 0 || section.point[i].x >= cxDIB-1)  continue;
		idx = section.point[i].y*cxDIB + section.point[i].x;
		image_buffer1[idx] = *nLabelNum;
		if (dx >= 0) {
			image_buffer1[idx+cxDIB] = *nLabelNum;  image_buffer1[idx-cxDIB] = *nLabelNum;
		}
		else {
			image_buffer1[idx+1] = *nLabelNum;  image_buffer1[idx-1] = *nLabelNum;
		}
	}
}

int check_stop_criteria(LPPOINTTAG pPoint,USHORT* image_buffer1,int cxDIB,USHORT *nRetLabel,LPPOINTTAG pArray,int nNum)
{
	USHORT val;
	int i;

	if (pPoint->x == NULLVALUE && pPoint->y == NULLVALUE) return (END_POINT);
	val = image_buffer1[pPoint->y*cxDIB+pPoint->x]; 
	if (val != 0) {
		*nRetLabel = val; return (BIF_POINT);
	}

	for (i = 0; i < nNum; i++) {
		if (pArray[i].x == pPoint->x && pArray[i].y == pPoint->y) return (LOOP_STOP);
	}

	return (NO_STOP);
}

int get_point_orient1(LPPOINTTAG pResult,LPPOINTTAG pStart,int cdir,USHORT nRetLabel,
			  BYTE *OrntImg,USHORT *image_buffer1,int cxDIB,int cyDIB)
{
	int cx = pResult->x, cy = pResult->y, nTh = 30;
	int r, x, y, dir, i, val, diff, rdir, rdiff;
	int n, px, py, pdir[2];

	x = cx - pStart->x; y = cy - pStart->y;
	r = op_func_02(x*x + y*y);

	n = 0; px = py = 10000;
	for (i = -nTh; i <= nTh; i++) {
		dir = cdir + i;
		dir = ANGLE_0_240(dir);
		x = cx + ((r * _table_03[dir]) >> 14);
		y = cy + ((r * _table_04[dir]) >> 14);
		if (x < 0 || x >= cxDIB || y < 0 || y >= cyDIB) continue;
		if (image_buffer1[y*cxDIB+x] != nRetLabel) continue;
		if (abs(x-px) <= 1 && abs(y-py) <= 1) {
			px = x; py = y; continue;
		}
		px = x; py = y;
		pdir[n++] = op_func_01(x,y,cx,cy);
		if (n >= 2) break;
	}
	if (n > 0) {
		if (n == 2) {
			dir = OrntImg[cy*cxDIB+cx];
			rdir = dir + 120;
			rdir = ANGLE_240(rdir);
			diff = abs(dir - pdir[0]);
			diff = IANGLE_120(diff);
			rdiff = abs(rdir - pdir[0]);
			rdiff = IANGLE_120(rdiff);
			if (rdiff < diff) diff = rdiff;
			val = abs(dir - pdir[1]);
			val = IANGLE_120(val);
			rdiff = abs(rdir - pdir[1]);
			rdiff = IANGLE_120(rdiff);
			if (rdiff < val) val = rdiff;
			dir = pdir[0];
			if (val < diff) dir = pdir[1];
		}
		else {
			dir = pdir[0];
		}
		diff = abs( cdir - dir );
		if (diff > 120) {
			diff = (240 - diff) / 2;
			if (cdir > 120) diff += cdir;
			else diff += dir;
			diff = ANGLE_240(diff);
			return (diff);
		}
		else {
			return ((cdir+dir)/2);
		}
	}
	return (-1);
}

int follow_ridge_point(LPPOINTTAG pStart,int nStartDir,BOOL nBackFlag,
				   BYTE *Img,USHORT *image_buffer1,BYTE *OrntImg,int cxDIB,int cyDIB,
				   LPPOINTTAG pResult,int *nResultDir,
				   LPPOINTTAG pRetArray,int *nRetNum,BOOL *nRetFlag,
				   USHORT *nLabelNum,USHORT *nRetLabel
				   )
{
	POINTTAG currentp, tempp, newp, backPoint;
	POINTTAG rarray[21], parray[8];
	BOOL nReTrace = FALSE, nReFlag = FALSE;
	int currentdir, prevdir = 0, result = 0, dy, dx;
	int nArrayId = 0, i, dir, num1, num2;
	int nStep0 = 3, nReId = 0;

	*nRetLabel = 0; *nRetNum = 0; *nRetFlag = FALSE;
	currentp = *pStart; currentdir = nStartDir;
	if (nBackFlag == FALSE) {
		currentdir = currentdir + 120;
		currentdir = ANGLE_240(currentdir);
	}
	prevdir = currentdir;
	parray[1] = currentp;
	
	while (1) {
		if (nReFlag == TRUE) {
			if (nReId == 21) {
				for (i = 0; i < 20; i++) rarray[i] = rarray[i+1];
				nReId = 20;
			}
			rarray[nReId++] = currentp;
		}
		if (nArrayId == 8) {
			for (i = 0; i < 7; i++) parray[i] = parray[i+1];
			nArrayId = 7;
		}
		parray[nArrayId++] = currentp;
		if (*nRetNum < 8 && nStep0 == 3) {
			pRetArray[*nRetNum] = currentp;
			(*nRetNum)++;
		}
		if (nStep0 == 3) {
			tempp.x = currentp.x + ((3 * (int)_table_03[currentdir]) >> 14);
			tempp.y = currentp.y + ((3 * (int)_table_04[currentdir]) >> 14);
		}
		else {
			if (currentdir < 15) { dx = 1; dy = 0; }
			else if (currentdir < 45) { dx = 1; dy = 1; }
			else if (currentdir < 75) { dx = 0; dy = 1; }
			else if (currentdir < 105) { dx = -1; dy = 1; }
			else if (currentdir < 135) { dx = -1; dy = 0; }
			else if (currentdir < 165) { dx = -1; dy = -1; }
			else if (currentdir < 195) { dx = 0; dy = -1; }
			else if (currentdir < 225) { dx = 1; dy = -1; }
			else { dx = 1; dy = 0; }
			tempp.x = currentp.x + dx;
			tempp.y = currentp.y + dy;
			if (abs(currentp.x-pStart->x) <= 1 && abs(currentp.y-pStart->y) <= 1) {
				tempp.x = tempp.x + dx;
				tempp.y = tempp.y + dy;
			}
		}
		if (tempp.x < 0 || tempp.x >= cxDIB || tempp.y < 0 || tempp.y >= cyDIB 
			|| OrntImg[tempp.y*cxDIB+tempp.x] == 255) {
			for (i = 0 ; i < nArrayId-1 ; i++)
				labelling(&parray[i],&parray[i+1],image_buffer1,cxDIB,cyDIB,nLabelNum);
			return (EXCEPTION_STOP);
		}

		newp = get_local_maximum(&tempp,currentdir,FALSE,Img,image_buffer1,cxDIB,cyDIB);

		if (newp.x != NULLVALUE && newp.y != NULLVALUE && (newp.x < 0 || newp.x >= cxDIB || newp.y < 0 || newp.y >= cyDIB) ) {
			return (EXCEPTION_STOP);
		}
		result = check_stop_criteria(&newp,image_buffer1,cxDIB,nRetLabel,parray,nArrayId);
		if (result == NO_STOP) {
			if (nArrayId >= 7)
				labelling(&parray[0],&parray[1],image_buffer1,cxDIB,cyDIB,nLabelNum);
			currentp = newp;
			currentdir = OrntImg[currentp.y*cxDIB+currentp.x];
			if (currentp.x < 0 || currentp.x >= cxDIB || currentp.y < 0 || currentp.y >= cyDIB || currentdir == 255)	
				return (EXCEPTION_STOP);
			dir = currentdir + 120;
			dir = ANGLE_240(dir);
			num1 = abs(prevdir - currentdir);
			num1 = IANGLE_120(num1);

			num2 = abs(prevdir - dir);
			num2 = IANGLE_120(num2);
			if (num2 < num1) currentdir = dir;
			prevdir = currentdir;
		}
		else {
			if (result == END_POINT)
				*pResult = currentp;
			else if (result == BIF_POINT)
				*pResult = newp;
			if (nArrayId < 7 || nReTrace == TRUE) {
				if (nStep0 == 3) {
					currentp = parray[0];
					currentdir = OrntImg[currentp.y*cxDIB+currentp.x];
					if (currentp.x < 0 || currentp.x >= cxDIB || currentp.y < 0 || currentp.y >= cyDIB || currentdir == 255)
						return (EXCEPTION_STOP);
					dir = currentdir + 120;
					dir = ANGLE_240(dir);
					num1 = abs(prevdir - currentdir);
					num1 = IANGLE_120(num1);
					num2 = abs(prevdir - dir);
					num2 = IANGLE_120(num2);
					if (num2 < num1) currentdir = dir;
					prevdir = currentdir;
					nStep0 = 1; nArrayId = 1;
					nReTrace = TRUE;
					nReId = 0; nReFlag = TRUE;
					result = NO_STOP;
					continue;
				}
				if (nReId >= 20) {
					dir = op_func_01(rarray[0].x,rarray[0].y,pResult->x,pResult->y);
					if (result == BIF_POINT) {
						dir = get_point_orient1(pResult,&rarray[0],dir,*nRetLabel,OrntImg,image_buffer1,cxDIB,cyDIB);
					}
				}
				else {
					dir = -1; *nRetFlag = TRUE;
				}
				*nResultDir = dir;
				for (i = 0 ; i < nArrayId-1; i++)
					labelling(&parray[i],&parray[i+1],image_buffer1,cxDIB,cyDIB,nLabelNum);
				break;
			}
			if (nStep0 == 3) {
				currentp = parray[4];
				backPoint = parray[0];
				for (i = 4; i < 7 ; i++) parray[i] = parray[3];
				parray[7] = currentp;
				labelling(&parray[0],&parray[1],image_buffer1,cxDIB,cyDIB,nLabelNum);
				currentdir = OrntImg[currentp.y*cxDIB+currentp.x];
				dir = currentdir + 120;
				dir = ANGLE_240(dir);
				num1 = abs(prevdir - currentdir);
				num1 = IANGLE_120(num1);
				num2 = abs(prevdir - dir);
				num2 = IANGLE_120(num2);
				if (num2 < num1) currentdir = dir;
				prevdir = currentdir;
				nStep0 = 1; nReId = 0; nReFlag = TRUE;
				result = NO_STOP;
				continue;
			}
			if (nReId >= 20) {
				dir = op_func_01(rarray[0].x,rarray[0].y,pResult->x,pResult->y);
				backPoint = rarray[0];
			}
			else {
				dir = op_func_01(backPoint.x,backPoint.y,pResult->x,pResult->y);
			}
			*nResultDir = dir;
			if (result == BIF_POINT) {
				*nResultDir = get_point_orient1(pResult,&backPoint,dir,*nRetLabel,OrntImg,image_buffer1,cxDIB,cyDIB);
			}
			for (i = 0; i < nArrayId-1; i++)
				labelling(&parray[i],&parray[i+1],image_buffer1,cxDIB,cyDIB,nLabelNum);
			labelling(&parray[nArrayId-1],&currentp,image_buffer1,cxDIB,cyDIB,nLabelNum);
			break;
		}
	}
	return (result);
}

BOOL trace_ridge(LPPOINTTAG pPoint,int dir,BOOL nFlag,
				  BYTE *Img,USHORT *image_buffer1,BYTE *OrntImg,int cxDIB,int cyDIB,
				  REALMINUTIA *pMinutia,LPPOINTTAG pRetArray,int *nRetNum,BOOL *nRetFlag,
				  USHORT *nLabelNum,USHORT *nRetLabel
				  )
{
	POINTTAG resultp;
	int resultdir, result;

	result = follow_ridge_point(pPoint,dir,nFlag,Img,image_buffer1,OrntImg,cxDIB,cyDIB,
							&resultp,&resultdir,pRetArray,nRetNum,nRetFlag,
							nLabelNum,nRetLabel);
	if (result == END_POINT || result == BIF_POINT) {
		pMinutia->x = resultp.x;
		pMinutia->y = resultp.y;
		pMinutia->dir = (short)resultdir;
		pMinutia->kind = (result==END_POINT)? END_MINUTIA: BIF_MINUTIA;
		return TRUE;
	}
	return FALSE;
}

BOOL start_trace(LPPOINTTAG pPoint,BYTE *Img,USHORT *image_buffer1,
					BYTE *OrntImg,int cxDIB,int cyDIB,
					LPREALPVECT pVect,USHORT *nLabelNum
					)
{
	REALMINUTIA minutiap1, minutiap2;
	POINTTAG currentp, array1[8], array2[8], temp;
	int startdir, i, j, num1, num2, dir, idx;
	BOOL flag, flag1 = FALSE, flag2 = FALSE;
	USHORT label1, label2;
	USHORT *ptr1, *ptr2;

	startdir = OrntImg[pPoint->y*cxDIB+pPoint->x];

	currentp = get_local_maximum(pPoint,startdir,TRUE,Img,image_buffer1,cxDIB,cyDIB);

	if (currentp.y == NULLVALUE && currentp.x == NULLVALUE) return FALSE;

	if (currentp.y < 3 || currentp.y >= cyDIB-3) return FALSE;
	if (currentp.x < 3 || currentp.x >= cxDIB-3) return FALSE;

	idx = currentp.y*cxDIB + currentp.x;

	if (OrntImg[idx] == 255) return FALSE;

	if (Img[idx] > 225) return FALSE;

	ptr1 = image_buffer1 + idx - 3*cxDIB;  ptr2 = image_buffer1 + idx - 3;
	for (i = 6; i >= 0; i--, ptr1 += cxDIB, ptr2++) {
		if (*ptr1 != 0) return FALSE;
		if (*ptr2 != 0) return FALSE;
	}

	if (trace_ridge(&currentp,startdir,TRUE,Img,image_buffer1,OrntImg,cxDIB,cyDIB,
					  &minutiap1,array1,&num1,&flag,nLabelNum,&label1) == TRUE) {
		if (flag == FALSE) {
			pVect->item[pVect->nNumber++] = minutiap1;
			if (pVect->nNumber >= MAXMINUTIATAGNUM) pVect->nNumber--;
		}
		else flag1 = TRUE;
	}
	if (trace_ridge(&currentp,startdir,FALSE,Img,image_buffer1,OrntImg,cxDIB,cyDIB,
					  &minutiap2,array2,&num2,&flag,nLabelNum,&label2) == TRUE) {
		if (flag == FALSE) {
			pVect->item[pVect->nNumber++] = minutiap2;
			if (pVect->nNumber >= MAXMINUTIATAGNUM) pVect->nNumber--;
		}
		else {
			if (flag1 == FALSE && num1+num2 >= 8) {
				dir = op_func_01(array1[7-num2].x,array1[7-num2].y,minutiap2.x,minutiap2.y);
				if (minutiap2.kind == BIF_MINUTIA) {
					temp.x = minutiap2.x; temp.y = minutiap2.y;
					dir = get_point_orient1(&temp,&array1[7-num2],dir,label2,OrntImg,image_buffer1,cxDIB,cyDIB);
				}
				minutiap2.dir = dir;
				pVect->item[pVect->nNumber++] = minutiap2;
				if (pVect->nNumber >= MAXMINUTIATAGNUM) pVect->nNumber--;
			}
			flag2 = TRUE;
		}
	}
	if (flag1 == TRUE && flag2 == FALSE && num1+num2 >= 8) {
		dir = op_func_01(array2[7-num1].x,array2[7-num1].y,minutiap1.x,minutiap1.y);
		if (minutiap1.kind == BIF_MINUTIA) {
			temp.x = minutiap1.x; temp.y = minutiap1.y;
			dir = get_point_orient1(&temp,&array2[7-num1],dir,label1,OrntImg,image_buffer1,cxDIB,cyDIB);
		}
		minutiap1.dir = dir;
		pVect->item[pVect->nNumber++] = minutiap1;
		if (pVect->nNumber >= MAXMINUTIATAGNUM) pVect->nNumber--;
	}
	if (flag1 == TRUE && flag2 == TRUE && num1+num2 >= 7) {
		if (num1+num2 == 7) j = 6;
		else j = 7;
		dir = op_func_01(array2[j-num1].x,array2[j-num1].y,minutiap1.x,minutiap1.y);
		if (minutiap1.kind == BIF_MINUTIA) {
			temp.x = minutiap1.x; temp.y = minutiap1.y;
			dir = get_point_orient1(&temp,&array2[j-num1],dir,label1,OrntImg,image_buffer1,cxDIB,cyDIB);
		}
		minutiap1.dir = dir;
		pVect->item[pVect->nNumber++] = minutiap1;
		if (pVect->nNumber >= MAXMINUTIATAGNUM) pVect->nNumber--;

		dir = op_func_01(array1[j-num2].x,array1[j-num2].y,minutiap2.x,minutiap2.y);
		if (minutiap2.kind == BIF_MINUTIA) {
			temp.x = minutiap2.x; temp.y = minutiap2.y;
			dir = get_point_orient1(&temp,&array1[j-num2],dir,label2,OrntImg,image_buffer1,cxDIB,cyDIB);
		}
		minutiap2.dir = dir;
		pVect->item[pVect->nNumber++] = minutiap2;
		if (pVect->nNumber >= MAXMINUTIATAGNUM) pVect->nNumber--;
	}
	(*nLabelNum)++;
	return TRUE;
}

void get_mp_points(BYTE *Img,USHORT *image_buffer1,BYTE *OrntImg,int cxDIB,int cyDIB,
				LPREALPVECT pVect)
{	
	int i, j;
	USHORT nLabelNum = 1;
	POINTTAG startp;
	BYTE* orntbase = OrntImg + 4*cxDIB + 4,  *orntptr;
	USHORT* markbase = image_buffer1 + 4*cxDIB + 4,  *markptr;

	memset(image_buffer1, 0, sizeof(USHORT)*cxDIB*cyDIB);
	pVect->nNumber = 0;  startp.y = 4;  
	for (i = 0; i < cyDIB; i += 8) {
		orntptr = orntbase;  markptr = markbase;  startp.x = 4;
		for (j = 0; j < cxDIB; j += 8, orntptr += 8, markptr += 8, startp.x += 8) {
			if (*markptr != 0) continue;
			if (*orntptr >= 120) continue;
			start_trace(&startp,Img,image_buffer1,OrntImg,cxDIB,cyDIB,pVect,&nLabelNum);
		}
		orntbase += 8 * cxDIB;  markbase += 8 * cxDIB;  startp.y += 8;
	}
}

BOOL check_false_mp(int x1,int y1,int dir1,int x2,int y2,int dir2)
{
	int dir, diff, tmp, dx, dy;
	
	dir = op_func_01(x2,y2, x1,y1);
	diff = abs(dir - dir1);
	diff = IANGLE_120(diff);
	dx = abs(x1 - x2); dy = abs(y1 - y2);
	if (dx > 13 || dy > 13) return FALSE;
	if (dx < 7 && dy < 7) {
		if (diff < 97) return FALSE;
	}
	else {
		if (diff < 100) return FALSE;
	}
	tmp = dir + 120;
	tmp = ANGLE_240(tmp);
	diff = abs(dir2 - tmp);
	diff = IANGLE_120(diff);
	if (dx < 7 && dy < 7) {
		if (diff < 97) return FALSE;
	}
	else {
		if (diff < 100) return FALSE;
	}
	return TRUE;
}

void filter_mp_points(LPREALPVECT pVect,SINGULAR* pSingular,BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, m, n, count, dx, dy, x, y;
	int dTH = 16;
	BYTE *orntbase, *orntptr;
	
	for (i = pVect->nNumber-1; i >= 0; i--) {
		if (pVect->item[i].kind != END_MINUTIA) continue;
		if (pVect->item[i].dir < 0) continue;
		for (j = pVect->nNumber-1; j >= 0; j--) {
			if (pVect->item[j].kind != END_MINUTIA) continue;
			if (pVect->item[j].dir < 0) continue;
			if (i == j) continue;
			if (FALSE != check_false_mp(pVect->item[i].x,pVect->item[i].y,pVect->item[i].dir,pVect->item[j].x,pVect->item[j].y,pVect->item[j].dir)) {
				pVect->item[i].dir = -1;  pVect->item[j].dir = -1;
				break;
			}
		}
	}
	for (i = pVect->nNumber-1; i >= 0; i--) {
		count = 0;
		for (j = pVect->nNumber-1; j >= 0; j--) {
			if (i == j) continue;
			dy = pVect->item[i].y - pVect->item[j].y;
			dx = pVect->item[i].x - pVect->item[j].x;
			dx = dy*dy + dx*dx;
			if (dx < dTH*dTH) count++;
		}
		if (count > 5) pVect->item[i].dir = -15;
	}
	for (i = pVect->nNumber-1; i >= 0; i--) {
		for (j = pVect->nNumber-1; j >= 0; j--) {
			if (i == j) continue;
			dy = pVect->item[i].y - pVect->item[j].y;
			dx = pVect->item[i].x - pVect->item[j].x;
			if (dx*dx+dy*dy <= 4*4) {
				pVect->item[i].dir = pVect->item[j].dir = -1;
				break;
			}
		}
	}
	for (i = 0; i < pVect->nNumber; i++) {
		if (pVect->item[i].dir < 0) continue;
		x = pVect->item[i].x;  y = pVect->item[i].y;
		if (x < 8 || x >= cxDIB-8 || y < 8 || y >= cyDIB-8) { 
			pVect->item[i].dir = -25;  continue;
		}
		orntbase = OrntImg + (y - 8) * cxDIB + (x - 8);
		for (m = 15; m >= 0; m -= 2, orntbase += 2*cxDIB) {
			orntptr = orntbase;
			for (n = 15; n >= 0; n -= 2, orntptr += 2) {
				if (*orntptr >= 120) {  pVect->item[i].dir = -25; break;  }
			}
			if (n >= 0) break;
		}
	}
	for (i = 0; i < pVect->nNumber; i++) {
		if (pVect->item[i].dir < 0) continue;
		x = pVect->item[i].x; y = pVect->item[i].y;
		for (j = 0; j < pSingular->nNumber; j++) {
			if (pSingular->nType[j] != 1) continue;
			dx = x - pSingular->nX[j]; dy = y - pSingular->nY[j];
			if (dx*dx+dy*dy < dTH*dTH) break;
		}
		if (j < pSingular->nNumber) pVect->item[i].dir = -1;
	}

	for (i = 0,j = 0; i < pVect->nNumber; i++) {
		if (pVect->item[i].dir < 0)	continue;
		pVect->item[j++] = pVect->item[i];
	}
	pVect->nNumber = j;
}

BYTE get_density(SINGULAR *pSingular,BYTE* OrntImg,int nW,USHORT* Img,int cxDIB,int cyDIB)
{
	int nTotal=0, nLocal=0, nCoreNumber=0;
	int i, j, sum, cx, cy, sx, sy, ex, ey;
	
	BYTE *orntbase, *orntptr;
	USHORT *imgbase, *ptr0, *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptr6, *ptr7, *ptr8;

	for (i = 0; i < pSingular->nNumber; i++) {
		if (pSingular->nType[i] == -1) continue;
		nCoreNumber++;
	}
	if (nCoreNumber == 0) {
		cx = cxDIB/2;  cy = cyDIB/2;
	}
	else {
		cx = cy = 0;
		for (i = 0; i < pSingular->nNumber; i++) {
			if (pSingular->nType[i] == -1) continue;
			cx += pSingular->nX[i];  cy += pSingular->nY[i];
		}
		cx /= nCoreNumber; cy /= nCoreNumber;
	}

	sx = (cx > nW) ? cx-nW+1 : 1; 
	ex = (cx+nW >= cxDIB) ? cxDIB-2 : cx+nW-1; 
	sy = (cy > nW) ? cy-nW+1 : 1; 
	ey = (cy+nW >= cyDIB) ? cyDIB-2 : cy+nW-1;

	orntbase = OrntImg + sy*cxDIB + sx;  imgbase = Img + sy*cxDIB + sx;

	for (j = sy; j < ey; j++, orntbase += cxDIB, imgbase += cxDIB) {
		orntptr = orntbase;  ptr0 = imgbase;
		ptr1 = ptr0 + 1;  ptr2 = ptr0 - 1;  ptr3 = ptr0 + cxDIB;  ptr4 = ptr0 - cxDIB;
		ptr5 = ptr3 + 1;  ptr6 = ptr3 - 1;  ptr7 = ptr4 + 1;  ptr8 = ptr4 - 1;
		for (i = ex-sx-1; i >= 0; i--, orntptr++, ptr0++, ptr1++, ptr2++, ptr3++, ptr4++, ptr5++, ptr6++, ptr7++, ptr8++) {
			if (*orntptr >= 120) continue;
			nTotal++;
			sum = 0;
			if (*ptr0 == 0)  sum++;
			if (*ptr1 == 0)  sum++;
			if (*ptr2 == 0)  sum++;
			if (*ptr3 == 0)  sum++;
			if (*ptr4 == 0)  sum++;
			if (*ptr5 == 0)  sum++;
			if (*ptr6 == 0)  sum++;
			if (*ptr7 == 0)  sum++;
			if (*ptr8 == 0)  sum++;
			if (sum == 0 || sum >= 9) continue;
			nLocal++;
		}
	}

	if (nTotal <= 0) return (0);
	return ((BYTE)(255*nLocal/nTotal));
}

void filter_mp_points2(LPREALPVECT pVect)
{
	int i, j, dx, dy, len, diff;

	for (i = 0; i < pVect->nNumber; i++) {
		if (pVect->item[i].score >= 35) continue;
		for (j = 0; j < pVect->nNumber; j++) {
			if (i == j) continue;
			if (pVect->item[j].score >= 35) continue;
			dy = pVect->item[i].y - pVect->item[j].y;
			dx = pVect->item[i].x - pVect->item[j].x;
			len = dy*dy + dx*dx;
			if (len >= 8*8) continue;
			diff = abs(pVect->item[i].dir - pVect->item[j].dir);
			if (diff > 120) diff = 240 - diff;
			if (120 - diff >= 20) continue;
			pVect->item[i].dir = -1; 
			pVect->item[j].dir = -1;
			break;
		}
	}

	for (i = 0,j = 0; i < pVect->nNumber; i++) {
		if (pVect->item[i].dir < 0) continue;
		pVect->item[j++] = pVect->item[i];
	}
	pVect->nNumber = j;
}

BYTE get_point_curve(int x,int y,BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, count = 0, sum = 0;
	int sx, sy, ex, ey;
	BYTE p, p0 = OrntImg[y*cxDIB+x];
	BYTE *orntbase,  *orntptr;

	sy = (y < 10) ? 0 : y-10;
	ey = (y > cyDIB-11) ? cyDIB-1 : y+10;
	sx = (x < 10) ? 0 : x-10;
	ex = (x > cxDIB-11) ? cxDIB-1 : x+10;
	orntbase = OrntImg + sy*cxDIB + sx;
	for (i = sy; i <= ey; i += 2, orntbase += 2*cxDIB) {
		orntptr = orntbase;
		for (j = ex-sx; j >= 0; j -= 2, orntptr += 2) {
			p = *orntptr;
			if (p == 0xFF) continue;
			p = abs(p0 - p);
			if (p > 60) p = 120 - p;
			sum += p; count++;
		}
	}
	if (count == 0)	return 0;
	sum = (255*sum) / (60*count);
	if (sum > 127) sum = 127;
	return ((BYTE)(sum));
}

void get_point_value(LPREALPVECT pVect,BYTE* Img,int cxDIB,int cyDIB)
{
	int i, m, n, x, y, sum, num, mstart, nstart, mend, nend, cx2 = cxDIB >> 1;
	BYTE *ptr, *ptrbase, p;

	for (i = 0; i < pVect->nNumber; i++) {
		x = pVect->item[i].x; y = pVect->item[i].y;
		sum = num = 0;
		mstart = (y > 10) ? (y-10) >> 1 : 0;
		mend = (y+10 > cyDIB-1) ? (cyDIB-1) >> 1 : (y+10) >> 1;
		nstart = (x > 10) ? (x-10) >> 1 : 0;
		nend = (x+10 > cyDIB-1) ? (cxDIB-1) >> 1 : (x+10) >> 1;
		ptrbase = Img + mstart*cx2+nstart;
		for (m = mstart; m <= mend; m++, ptrbase += cx2) {
			ptr = ptrbase;
			for (n = nstart; n <= nend; n++, ptr++) {
				p = *ptr;
				if (p > 63) sum += 63;
				else sum += p; 
				num++;
			}
		}
		if (num == 0) sum = 0;
		else sum /= num;
		pVect->item[i].score = (BYTE)(sum);
	}
}

void arrange_mp(LPREALPVECT pRealVect,LPMPVECTEX pVect,BYTE* OrntImg,int cxDIB,int cyDIB)
{
	int i, j, nMinIdx, nMinValue;
	REALMINUTIA tmp;

	for (i = 0; i < pRealVect->nNumber; i++) {
		if (pRealVect->item[i].score > 10) continue;
		pRealVect->item[i].dir = -1;
	}

	filter_mp_points2(pRealVect);

	for (i = 0; i < pRealVect->nNumber-1; i++) {
		nMinIdx = i; nMinValue = pRealVect->item[i].score;
		for (j = i+1; j < pRealVect->nNumber; j++) {
			if (pRealVect->item[j].score <= nMinValue) continue;
			nMinIdx = j; nMinValue = pRealVect->item[j].score;
		}
		if (nMinIdx == i) continue;
		tmp = pRealVect->item[i];
		pRealVect->item[i] = pRealVect->item[nMinIdx];
		pRealVect->item[nMinIdx] = tmp;
	}
	for (i = 0,j = 0; i < pRealVect->nNumber; i++) {
		pVect->item[j].x = pRealVect->item[i].x;
		pVect->item[j].y = pRealVect->item[i].y;
		pVect->item[j].dir = (BYTE)pRealVect->item[i].dir;
		pVect->item[j].kind = pRealVect->item[i].kind;
		pVect->item[j].curv = get_point_curve(pVect->item[j].x, pVect->item[j].y, OrntImg, cxDIB, cyDIB);
		pVect->item[j++].score = pRealVect->item[i].score;
		if (j >= MAX_MINUTIA_NUMBER) break;
	}
	pVect->nNumber = j;

}

void comp_typeline(DESMAINLINE *pLine,BYTE *pData)
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
	int i, k, n = pBlock->nRow*pBlock->nCol;
	BYTE val;

	for (i = 0,k = 0; i < n; i += 2,k++) {
		if (pBlock->Data[i] == 0xFF) val = 15;
		else val = pBlock->Data[i] / 8;
		pData[k] = (val & 0x0F) << 4;
		if (pBlock->Data[i+1] == 0xFF) val = 15;
		else val = pBlock->Data[i+1] / 8;
		pData[k] |= (val & 0x0F);
	}
}

void comp_core(COREVECTEX *pCore,BYTE *pData)
{
	int i, k;

	for (i = 0,k = 0; i < pCore->nNumber; i++,k += 4) {
		pData[k+0] = (BYTE)(((pCore->item[i].x << 7) & 0xFF00) >> 8);
		pData[k+1] = (BYTE)((pCore->item[i].x & 0x0001) << 7);
		pData[k+1] |= (BYTE)((pCore->item[i].y & 0x0100) >> 8 );
		pData[k+2] = (BYTE)(pCore->item[i].y & 0x00FF);
		pData[k+3] = pCore->item[i].dir;
	}
}

void get_byte_template(LPFPVECTEX pFPEx,LPFPFEATURE pFeature)
{
	int i, n, n1, n2, buf;

	pFeature->Data[0] = pFPEx->nHeader.nVersion;
	pFeature->Data[1] = pFPEx->nHeader.nHeaderLength;
	pFeature->Data[2] = (BYTE)(pFPEx->nHeader.nImageX/2);
	pFeature->Data[3] = (BYTE)(pFPEx->nHeader.nImageY/2);

	n = pFPEx->nHeader.nHeaderLength;
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
	n1 = 2 + (pFPEx->Block.nRow*pFPEx->Block.nCol)/2;
	pFeature->Data[n++] = n1;
	pFeature->Data[n] = pFPEx->Block.nRow;
	pFeature->Data[n+1] = pFPEx->Block.nCol;
	comp_block(&(pFPEx->Block),&pFeature->Data[n+2]);
	n += n1;

	pFeature->Data[n++] = 'D';
	n1 = 1+4*pFPEx->Core.nNumber;
	pFeature->Data[n++] = n1;
	pFeature->Data[n] = pFPEx->Core.nNumber;
	comp_core(&(pFPEx->Core),&pFeature->Data[n+1]);
	n += n1;

	pFeature->Data[n++] = 'E';
	n1 = n+2;
	for (i = 0; i < pFPEx->Mp.nNumber; i++) {
		if (MAX_FEATUREVECT_LEN-n1 < 5) break;
		buf = (pFPEx->Mp.item[i].x << 15) | (pFPEx->Mp.item[i].y << 6) | (pFPEx->Mp.item[i].score);
		pFeature->Data[n1++] = (BYTE)((buf >> 16) & 0xFF);
		pFeature->Data[n1++] = (BYTE)((buf >> 8) & 0xFF);
		pFeature->Data[n1++] = (BYTE)(buf & 0xFF);
		pFeature->Data[n1++] = pFPEx->Mp.item[i].dir;
		pFeature->Data[n1++] = (pFPEx->Mp.item[i].curv << 1) | pFPEx->Mp.item[i].kind; 
	}
	n2 = i*5+2;
	pFeature->Data[n++] = (BYTE)(n2);
	pFeature->Data[n] = (BYTE)(i);
}

int ext_main(unsigned char* image,int width,int height,LPFPFEATURE pFeatureVect)
{
	BYTE *image_buffer3 = NULL, *image_buffer4 = NULL, *image_buffer2 = NULL;
	USHORT *image_buffer1 = NULL;
	REALPVECT tempVect;
	SINGULAR SingularData;
	FPVECTEX FPEx;
	int cxDIB = width, cyDIB = height;
	int nCol = cxDIB / BLOCK_SIZE, nRow = cyDIB / BLOCK_SIZE;

	memset(pFeatureVect, 0, sizeof(MINUTIAVECT));
	memset(&FPEx, 0, sizeof(FPVECTEX));

	image_buffer3 = (BYTE*)calloc(cxDIB*cyDIB, sizeof(BYTE));
	if (image_buffer3 == NULL)  return ERR_CAN_NOT_ALLOC_MEMORY;

	image_buffer1 = (USHORT*)malloc(sizeof(USHORT)*cxDIB*cyDIB);
	if (image_buffer1 == NULL) { free(image_buffer3);  return ERR_CAN_NOT_ALLOC_MEMORY; }

	image_buffer2 = (BYTE*)malloc(sizeof(BYTE)*cxDIB*cyDIB);
	if (image_buffer2 == NULL) { free(image_buffer1); free(image_buffer3); return ERR_CAN_NOT_ALLOC_MEMORY; }

	image_buffer4 = (BYTE*)malloc(sizeof(BYTE)*cxDIB*cyDIB);
	if (image_buffer4 == NULL) { free(image_buffer2); free(image_buffer1); free(image_buffer3); return ERR_CAN_NOT_ALLOC_MEMORY; }

	get_smoothed_image(image,cxDIB,cyDIB);
	memcpy(image_buffer2, image, sizeof(BYTE)*cxDIB*cyDIB);
	get_segmentation(image,image_buffer3,cxDIB,cyDIB);

	get_smoothed_image4(image_buffer2,cxDIB,cyDIB);
	get_sharpend_image(image_buffer2, image, image_buffer3, cxDIB, cyDIB, 64);

	memcpy(image, image_buffer2, sizeof(BYTE)*cxDIB*cyDIB);

	get_smoothed_image(image_buffer2,cxDIB,cyDIB);

	get_orient_image(image_buffer3,image_buffer2,cxDIB,cyDIB,image_buffer4);

	image_proc_01(image,image_buffer3,image_buffer2,cxDIB,cyDIB);

	get_smoothed_image4(image,cxDIB,cyDIB);
	get_sharpend_image(image, image_buffer2, image_buffer3, cxDIB,cyDIB,64);

	get_binary_image2(image_buffer3,image_buffer2,image,cxDIB,cyDIB,3,7);
	
	image_proc_04(image_buffer2,cxDIB,cyDIB);

	remove_hole(image_buffer3,image_buffer2,cxDIB,cyDIB);

	re_get_orient_image(image_buffer3,image_buffer2,cxDIB,cyDIB);
	
	get_block_data(image_buffer3,cxDIB,cyDIB,&FPEx.Block,nCol,nRow);

	get_core_points(&SingularData,image_buffer3,cxDIB,cyDIB);

	FPEx.nType = get_type_line(&FPEx,&SingularData,image_buffer3,(BYTE*)image_buffer1,cxDIB,cyDIB);

	copy_core(&SingularData,&FPEx);
	
	FPEx.nFrequency = get_frequency(image_buffer2,image_buffer3,&SingularData,cxDIB,cyDIB);

	image_proc_01(image_buffer2,image_buffer3,image,cxDIB,cyDIB);

	get_mp_points(image,image_buffer1,image_buffer3,cxDIB,cyDIB,&tempVect);

	filter_mp_points(&tempVect,&SingularData,image_buffer3,cxDIB,cyDIB);

	FPEx.nRidgeDensity = get_density(&SingularData,image_buffer3,64,image_buffer1,cxDIB,cyDIB);

	get_point_value(&tempVect,image_buffer4,cxDIB,cyDIB);

	arrange_mp(&tempVect,&FPEx.Mp,image_buffer3,cxDIB,cyDIB);

	free(image_buffer4);  free(image_buffer2);  free(image_buffer1); free(image_buffer3); 

	if ( FPEx.Mp.nNumber < MIN_MINUTIA_NUMBER && FPEx.Core.nNumber == 0 ) return ERR_VECT_FAILED;

	FPEx.nHeader.nVersion = CURRENT_VERSION;
	FPEx.nHeader.nHeaderLength = MAX_HEADER_SIZE;
	FPEx.nHeader.nImageX = cxDIB;
	FPEx.nHeader.nImageY = cyDIB;

	get_byte_template(&FPEx,pFeatureVect);

	return ERR_OK;
}


void decomp_typeline(BYTE *pData,DESMAINLINE *pLine)
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
	int i, k, n = pBlock->nRow*pBlock->nCol/2;
	BYTE val;

	for (i = 0,k = 0; i < n; i++) {
		val = (pData[i] & 0xF0) >> 4;
		if (val == 15) pBlock->Data[k++] = 0xFF;
		else pBlock->Data[k++] = val*8;
		val = (pData[i] & 0x0F);
		if (val == 15) pBlock->Data[k++] = 0xFF;
		else pBlock->Data[k++] = val*8;
	}
}

void decomp_core(BYTE *pData,COREVECTEX *pCore)
{
	int i, k;

	for (i = 0,k = 0; i < pCore->nNumber; i++,k += 4) {
		pCore->item[i].x = pData[k+0] << 1;
		pCore->item[i].x |= (pData[k+1] & 0x80) >> 7;
		pCore->item[i].y = (pData[k+1] & 0x01) << 8;
		pCore->item[i].y |= pData[k+2];
		pCore->item[i].dir = pData[k+3];
		if (pCore->item[i].dir == 255) pCore->item[i].kind = DELTA;
		else pCore->item[i].kind = CORE;
	}
}

int mch_sub_func_04(LPMPVECTEX pVect)
{
	int i, score = 0;

	if (pVect->nNumber == 0) return (0);
	for (i = 0; i < pVect->nNumber; i++)
		score += pVect->item[i].score;
	
	score /= pVect->nNumber;

	return score;
}

void mch_sub_func_02(LPFPFEATURE pFeature,LPFPVECTEX pFPEx)
{
	int i, n, n1, nNum, buf;
	BYTE nMark;

	memset(pFPEx, 0, sizeof(FPVECTEX));

	pFPEx->nHeader.nVersion = pFeature->Data[0];
	pFPEx->nHeader.nHeaderLength = pFeature->Data[1];
	pFPEx->nHeader.nImageX = pFeature->Data[2]*2;
	pFPEx->nHeader.nImageY = pFeature->Data[3]*2;
	n = pFPEx->nHeader.nHeaderLength;
	while (1) {
		if (n > MAX_FEATUREVECT_LEN-2) break;
		nMark = pFeature->Data[n++];
		nNum = pFeature->Data[n++];
		if (nNum == 0) continue;
		if (n+nNum >= MAX_FEATUREVECT_LEN) break;

		if (nMark == 'A') {
			pFPEx->nType = pFeature->Data[n];
			pFPEx->nRidgeDensity = pFeature->Data[n+1];
			pFPEx->nFrequency = pFeature->Data[n+2];
		}

		if (nMark == 'B') {
			decomp_typeline(&pFeature->Data[n],&pFPEx->MainLine);
		}

		if (nMark == 'C') {
			pFPEx->Block.nRow = pFeature->Data[n];
			pFPEx->Block.nCol = pFeature->Data[n+1];
			decomp_block(&pFeature->Data[n+2],&pFPEx->Block);
		}

		if (nMark == 'D') {
			pFPEx->Core.nNumber = pFeature->Data[n];
			decomp_core(&pFeature->Data[n+1],&pFPEx->Core);
		}

		if (nMark == 'E') {
			if (nNum > 5+1) {
				pFPEx->Mp.nNumber = pFeature->Data[n];
				n1 = n+1;
				for (i = 0; i < pFPEx->Mp.nNumber; i++, n1 += 5) {
					if (MAX_FEATUREVECT_LEN-n1 < 5) break;
					buf = (pFeature->Data[n1] << 16) | (pFeature->Data[n1+1] << 8) | pFeature->Data[n1+2];
					pFPEx->Mp.item[i].x = (buf >> 15) & 0x1FF;
					pFPEx->Mp.item[i].y = (buf >> 6) & 0x1FF;
					pFPEx->Mp.item[i].score = buf & 0x3F;
					pFPEx->Mp.item[i].dir = pFeature->Data[n1+3];
					pFPEx->Mp.item[i].curv = (pFeature->Data[n1+4] >> 1) & 0x7F;
					pFPEx->Mp.item[i].kind = pFeature->Data[n1+4] & 0x1;
				}
				pFPEx->Mp.nNumber = i;
			}
		}
		n += nNum;
	}

	pFPEx->Mp.quality = mch_sub_func_04(&pFPEx->Mp);
}

BOOL mch_sub_func_03(LPFPVECTEX pFPEx)
{
	if ( pFPEx->nHeader.nVersion != CURRENT_VERSION ) 
		 return FALSE;
	if (pFPEx->Core.nNumber > MAX_CORE_NUMBER) return FALSE;
	if (pFPEx->nType > 10) return FALSE;
	if (pFPEx->Mp.nNumber > MAX_MINUTIA_NUMBER) return FALSE;
	if (pFPEx->Mp.nNumber < MIN_MINUTIA_NUMBER) return FALSE;
	return TRUE;
}

int mch_sub_func_01(LPCOREVECTEX pSingular,COREITEMEX *pCore,COREITEMEX *pDelta,int *nNumDelta)
{
	int i, nCoreNum = 0, nDeltaNum = 0;

	for (i = 0; i < pSingular->nNumber; i++) {
		if (pSingular->item[i].kind == CORE) {
			if (nCoreNum >= 2) break;
			pCore[nCoreNum++] = pSingular->item[i];
		}
		else {
			if (pDelta == NULL) continue;
			if (nDeltaNum >= 2) break;
			pDelta[nDeltaNum++] = pSingular->item[i];
		}
	}
	if (nNumDelta != NULL) *nNumDelta = nDeltaNum;
	return nCoreNum;
}

void transform_core(LPCOREVECTEX pVect,int cx,int cy,int nAngle,int nDiffX,int nDiffY)
{
	int i, x, y, angle, rot, nCos, nSin;

	rot = 240 - nAngle;
	rot = ANGLE_240(rot);
	nCos = _table_03[rot]; nSin=_table_04[rot];

	for (i = 0; i < pVect->nNumber; i++) {
		x = (pVect->item[i].x-cx)*nCos + (pVect->item[i].y-cy)*nSin;
		x = ROUND(x);
		x = x >> 14;
		y = (pVect->item[i].y-cy)*nCos - (pVect->item[i].x-cx)*nSin;
		y = ROUND(y);
		y = y >> 14;

		pVect->item[i].x = x + cx + nDiffX;
		pVect->item[i].y = y + cy + nDiffY;

		if (pVect->item[i].kind != CORE) continue;
		angle = pVect->item[i].dir + nAngle;
		angle = ANGLE_0_240(angle);
		pVect->item[i].dir = angle;
	}
}

void transform_points(LPFPVECTEX pVect,int xc,int yc,int nAngle,int nDiffX,int nDiffY)
{
	int i, newx, newy, angle, dx, dy;
	int	nCos = _table_03[nAngle], nSin = _table_04[nAngle];

	for (i = 0; i < pVect->Mp.nNumber; i++)	{
		dx = (nCos*(pVect->Mp.item[i].x-xc)) >> 14;
		dy = (nSin*(pVect->Mp.item[i].y-yc)) >> 14;
		newx = xc + dx - dy;
		dx = (nSin*(pVect->Mp.item[i].x-xc)) >> 14;
		dy = (nCos*(pVect->Mp.item[i].y-yc)) >> 14;
		newy = yc + dx + dy;
		pVect->Mp.item[i].x = newx + nDiffX;
		pVect->Mp.item[i].y = newy + nDiffY;
		angle = pVect->Mp.item[i].dir + nAngle;
		angle = ANGLE_0_240(angle);
		pVect->Mp.item[i].dir = angle;
	}
	for (i = 0; i < pVect->Core.nNumber; i++) {
		dx = (nCos*(pVect->Core.item[i].x-xc)) >> 14;
		dy = (nSin*(pVect->Core.item[i].y-yc)) >> 14;
		newx = xc + dx - dy;
		dx = (nSin*(pVect->Core.item[i].x-xc)) >> 14;
		dy = (nCos*(pVect->Core.item[i].y-yc)) >> 14;
		newy = yc + dx + dy;
		pVect->Core.item[i].x = newx + nDiffX;
		pVect->Core.item[i].y = newy + nDiffY;

		if (pVect->Core.item[i].kind != CORE) continue;
		angle = pVect->Core.item[i].dir + nAngle;
		angle = ANGLE_0_240(angle);
		pVect->Core.item[i].dir = angle;
	}
}

void transform_mp(LPMPVECTEX pVect,int cx,int cy,int nAngle,int nDiffX,int nDiffY)
{
	int i, x, y, dx, dy, angle, rot, nCos, nSin;

	rot = 240 - nAngle;
	rot = ANGLE_240(rot);
	nCos = _table_03[rot]; nSin=_table_04[rot];

	for (i = 0; i < pVect->nNumber; i++) {
		dx = pVect->item[i].x - cx; dy = pVect->item[i].y - cy;
		x = dx*nCos + dy*nSin;
		x = ROUND(x);
		x = x >> 14;
		y = dy*nCos - dx*nSin;
		y = ROUND(y);
		y = y >> 14;

		pVect->item[i].x = x + cx + nDiffX;
		pVect->item[i].y = y + cy + nDiffY;

		angle = pVect->item[i].dir + nAngle;
		angle = ANGLE_0_240(angle);
		pVect->item[i].dir = angle;
	}
}

int get_matched_mp_num(int nLenTh,int nAngTh,LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	int i, i1, i2, i3, j, num = 0;
	char temp[MAX_MINUTIA_NUMBER];
	BOOL flag;

	memset( temp, 0, sizeof(char)*pVect2->nNumber );
	for (i = 0; i < pVect1->nNumber; i++) {
		flag = FALSE;
		for (j = 0; j < pVect2->nNumber; j++) {
			i1 = pVect1->item[i].x - pVect2->item[j].x;
			i2 = pVect1->item[i].y - pVect2->item[j].y;
			i3 = i1*i1 + i2*i2;
			if (i3 > nLenTh*nLenTh) continue;
			i1 = abs(pVect1->item[i].dir - pVect2->item[j].dir);
			i1 = IANGLE_120(i1);
			if (i1 > nAngTh) continue;
			temp[j] = 1; flag = TRUE;
		}
		if (flag == TRUE) num++;
	}
	i1 = 0;
	for (i = 0; i < pVect2->nNumber; i++) {
		if (temp[i] == 1) i1++;
	}
	if (num > i1) num = i1;
	return (num);
}

int get_min_points_number(LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	POLYGON pol1, pol2;
	int i, num1, num2;

	if (get_polygon_points(pVect1,&pol1) == FALSE) return (0);
	if (get_polygon_points(pVect2,&pol2) == FALSE) return (0);
	num1 = num2 = 0;
	for (i = 0; i < pVect1->nNumber; i++) {
		if (TRUE == check_in_polygon(pVect1->item[i].x,pVect1->item[i].y,&pol2,16)) num1++;
	}
	for (i = 0; i < pVect2->nNumber; i++) {
		if (TRUE == check_in_polygon(pVect2->item[i].x,pVect2->item[i].y,&pol1,16)) num2++;
	}
	if (num1 > num2) num1 = num2;
	return (num1);
}

void get_matched_points_number(LPMPVECTEX pVect1,LPMPVECTEX pVect2,int *nNum1,int *nNum2)
{
	int i, j, dx, dy, ang, num1 = 0, num2 = 0;
	char temp1[MAX_MINUTIA_NUMBER], temp2[MAX_MINUTIA_NUMBER];
	BOOL flag1, flag2;

	for (i = 0; i < MAX_MINUTIA_NUMBER; i++) {
		temp1[i] = temp2[i] = 0; 
	}
	*nNum1 = *nNum2 = 0;
	for (i = 0; i < pVect1->nNumber; i++) {
		flag1 = flag2 = FALSE;
		for (j = 0; j < pVect2->nNumber; j++) {
			dx = pVect1->item[i].x - pVect2->item[j].x;
			dy = pVect1->item[i].y - pVect2->item[j].y;
			dx = dx*dx + dy*dy;
			if (dx > 12*12) continue;
			ang = abs(pVect1->item[i].dir - pVect2->item[j].dir);
			ang = IANGLE_120(ang);
			if (ang <= 7) { temp1[j] = 1; flag1 = TRUE; }
			if (ang <= ANGLE) { temp2[j] = 1; flag2 = TRUE; }
		}
		if (flag1 == TRUE) num1++;
		if (flag2 == TRUE) num2++;
	}
	dx = dy = 0;
	for (i = 0; i < pVect2->nNumber; i++) {
		if (temp1[i] == 1) dx++;
		if (temp2[i] == 1) dy++;
	}
	if (num1 > dx) num1 = dx;
	if (num2 > dy) num2 = dy;
	*nNum1 = num1; *nNum2 = num2;
}

void transform_block(int nRot,int xOff,int yOff,int cx,int cy,BLOCKVECT* pBlock)
{
	
	int c, nCos, nSin, d1, d2, dd1, dd2, i, j;
	int x, y, bx, by, d, dx, dy;
	int nCol = pBlock->nCol, nRow = pBlock->nRow, size = nCol * nRow; 
	BYTE pTmp[MAX_BLOCK_ROW*MAX_BLOCK_COL], *ptr;
	int x_sin[40], x_cos[40];

	memcpy(pTmp,pBlock->Data,sizeof(BYTE)*size);

	nCos = _table_03[nRot]; nSin = _table_04[nRot];

	d = (BLOCK_SIZE/2) - xOff - cx;
	x_cos[0] = d * nCos; x_sin[0] = d * nSin;
	dd1 = BLOCK_SIZE * nCos; dd2 = BLOCK_SIZE * nSin;

	for (i = 1; i < nCol; i++) {
		x_cos[i] = x_cos[i-1] + dd1; x_sin[i] = x_sin[i-1] + dd2;
	}

	ptr = pBlock->Data + size - 1;

	d = ((nRow - 1)*BLOCK_SIZE+BLOCK_SIZE/2) - yOff - cy;
	dx = d*nSin; dy = d*nCos;
	
	for (i = nRow-1; i >= 0; i--) {
		for (j = nCol-1; j >= 0; j--, ptr--) {
			bx = x_cos[j]; by = x_sin[j]; 
			x = bx + dx;
			x = ROUND(x);
			x = x >> 14;

			d1 = x + cx; 
			if (d1 < 0) { *ptr = 0xFF; continue; }
			d1 = d1 >> 4; 
			if (d1 >= nCol) { *ptr = 0xFF; continue; }

			y = dy - by;
			y = ROUND(y);
			y = y >> 14;

			d2 = y + cy;
			if (d2 < 0) { *ptr = 0xFF; continue; }
			d2 = d2 >> 4;
			if (d2 >= nRow) { *ptr = 0xFF; continue; }

			c = pTmp[d2*nCol+d1];
			if (c == 0xFF) *ptr = 0xFF;
			else {
				c += nRot;
				c = ANGLE_240(c);
				c = ANGLE_120(c);
				*ptr = (BYTE)c;
			}
		}
		dx -= dd2; dy -= dd1;
	}
}

int check_block(int Th,int nNumTh,BLOCKVECT* pFBlock,BLOCKVECT* pSBlock)
{
	int i, retVal = 0, num = 0, nDiv = 0, numF = 0, numS = 0, minNum;
	int size = pFBlock->nCol*pFBlock->nRow;
	BYTE c, *ptr, *ptrs;

	ptr = pFBlock->Data; ptrs = pSBlock->Data;

	for (i = size-1; i >= 0; i--, ptr++, ptrs++) {
		if (*ptr != 0xFF) numF++;
		if (*ptrs != 0xFF) numS++;
	}

	ptr = pFBlock->Data; ptrs = pSBlock->Data;
	
	for (i = size-1; i >= 0; i--, ptr++, ptrs++) {
		if (*ptr == 0xFF) continue;
		if (*ptrs == 0xFF) continue;
		c = abs(*ptr - *ptrs);
		num++;
		if (c > 60) c = 120 - c;
		if (c < 5){ retVal += 60; continue; }
		if (c > Th) continue; 
		retVal += 60 - c; 
	}
	
	if (num == 0) return (0);
	minNum = numS;
	if (numF < minNum) minNum = numF;
	if (nNumTh*num < minNum) return (0);

	nDiv = 60*num;
	retVal = 100*retVal/nDiv;
	return retVal;
}

int check_block2(int nNumTh,BLOCKVECT* pFBlock,BLOCKVECT* pSBlock)
{
	int i, num = 0,numF = 0, numS = 0, gVal, minNum;
	int size = pFBlock->nCol*pFBlock->nRow;
	BYTE *ptr, *ptrs;

	ptr = pFBlock->Data; ptrs = pSBlock->Data;

	for (i = size-1; i >= 0; i--, ptr++, ptrs++) {
		if (*ptr != 0xFF) numF++;
		if (*ptrs != 0xFF) numS++;
		if (*ptr != 0xFF && *ptrs != 0xFF) num++;
	}

	if (num == 0) return (0);
	minNum = numS;
	if (numF < minNum) minNum = numF;
	if (nNumTh*num < minNum) return (0);
	gVal = 2*num*100/(numF+numS);

	return gVal;
}

int coarse_matching(LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	FPVECTEX pVect;
	int i, j, cx = 0, cy = 0, dx = 0, dy = 0, rot = 0, nNum1 = 0, nNum2 = 0, nTotalNum = 0, nFCoreNum, nSCoreNum, nTh = 75;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);

	if (nFCoreNum == 0 || nSCoreNum == 0) return (0);
	if (pFile->Mp.nNumber == 0 || pSearch->Mp.nNumber == 0) return (-1);

	for (i = 0; i < nFCoreNum; i++)	{
		cx = FileCore[i].x; cy = FileCore[i].y;
		for (j = 0; j < nSCoreNum; j++)	{
			dx = SearchCore[j].x - cx;
			dy = SearchCore[j].y - cy;
			rot = SearchCore[j].dir - FileCore[i].dir;
			rot = ANGLE_0_240(rot);
			pVect = *pFile;
			transform_points(&pVect,cx,cy,rot,dx,dy);
			nTotalNum = get_min_points_number(&pVect.Mp,&pSearch->Mp);
			get_matched_points_number(&pVect.Mp,&pSearch->Mp,&nNum1,&nNum2);
			if (nNum1 > 6 && nNum1*100 > 80*nTotalNum && nTotalNum >= 10) return (1);
			if (nNum2 > 13 && nNum2*100 > nTh*nTotalNum) return (1);
		}
	}

	if (nFCoreNum == nSCoreNum && nFCoreNum == 1) {
		if (nNum1 > 6 && nNum2 >= 13 && nNum2*100 > 45*nTotalNum) {
			transform_block(rot,dx,dy,cx,cy,&pVect.Block);
			nNum1 = check_block(30,5,&pVect.Block,&pSearch->Block);
			if (nNum1 > 90) return (2);
		}
	}
	return (0);
}

int get_distance(LPFPVECTEX pVect,int CoreNumber,int MinMaxFlag)
{
	int x, y, dx, dy, nn, len1, len2;
	COREITEMEX FileCore[2];

	mch_sub_func_01(&(pVect->Core),FileCore,NULL,NULL);
	
	if (CoreNumber < 0 || CoreNumber > 1) return (0);
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

	if (MinMaxFlag == 0) { 
		nn = len1;
		if (nn > len2) nn = len2;
	}
	else {
		nn = len1;
		if (nn < len2) nn = len2;
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

	if (pFile->nType < pSearch->nType) {
		t1 = pFile->nType; t2 = pSearch->nType;
	}
	else {
		t1 = pSearch->nType; t2 = pFile->nType;
	}
	if (t1 == t2) return TRUE;
	pScore1 = pFile->Mp.quality;
	pScore2 = pSearch->Mp.quality;
	if (pScore1 == 0 || pScore2 == 0) return FALSE;

	i1 = i2 = pScore1;
	if (i1 > pScore2) i1 = pScore2;
	if (i2 < pScore2) i2 = pScore2;
	if (i1 <= 35 && i2 < 40) return TRUE;

	if (t1 < 3) {
		if (t2 == 3) return TRUE;
		if (t2 < 3)	{
			i1 = i2 = 0;
			if (nFCoreNum == 2)	{
				dx = FileCore[0].x - FileCore[1].x;
				dy = FileCore[0].y - FileCore[1].y;
				i1 = op_func_02(dx*dx+dy*dy);
			}
			if (nSCoreNum == 2) {
				dx = SearchCore[0].x - SearchCore[1].x;
				dy = SearchCore[0].y - SearchCore[1].y;
				i2 = op_func_02(dx*dx+dy*dy);
			}
			if (abs(i1-i2) <= 55) return TRUE;
			if (pFile->nType == 1 && nSCoreNum == 1) return TRUE;
			if (pSearch->nType == 1 && nFCoreNum == 1) return TRUE;
			return FALSE ;
		}
		if (t1 == 2 && (t2 == 4 || t2 == 5)) return TRUE;
		if (t2 > 3 && t2 < 8) {
			if (t1 == 0 && (t2 == 4 || t2 == 5))
			{
				if (pFile->nType == 0) {
					dx = pSearch->MainLine.points_x[0] - pSearch->MainLine.points_x[1];
					dy = pSearch->MainLine.points_y[0] - pSearch->MainLine.points_y[1];
				}
				else {
					dx = pFile->MainLine.points_x[0] - pFile->MainLine.points_x[1];
					dy = pFile->MainLine.points_y[0] - pFile->MainLine.points_y[1];
				}
				if (dx*dx+dy*dy <= 45*45) return TRUE;
			}
			if (t1 == 1) {
				if (pFile->nType == 1) i2 = get_distance(pSearch,0,0);
				else i2 = get_distance(pFile,0,0);
				if (i2 > 70) return FALSE;
				return TRUE;
			}
			return FALSE;
		}
		if (t1 != 2 && t2 == 8) {
			if (pScore2 > 35 && nFCoreNum == 2 && nSCoreNum == 1) {
				dx = FileCore[0].x - FileCore[1].x;
				dy = FileCore[0].y - FileCore[1].y;
				i1 = op_func_02(dx*dx+dy*dy);
				i2 = get_distance(pSearch,0,0);
				if (t1 == 1) {
					if (i1 < i2-40) return FALSE;
					return TRUE;
				}
				if (i1 < i2-36) return FALSE;
				if (i1 < i2-5) {
					dx = pSearch->MainLine.points_x[0] - pSearch->MainLine.points_x[1];
					dy = pSearch->MainLine.points_y[0] - pSearch->MainLine.points_y[1];
					if (dx*dx+dy*dy > 45*45) return FALSE;
				}
			}
			if (pScore1 > 35 && nFCoreNum == 1 && nSCoreNum == 2) {
				dx = SearchCore[0].x - SearchCore[1].x;
				dy = SearchCore[0].y - SearchCore[1].y;
				i1 = op_func_02(dx*dx+dy*dy);
				i2 = get_distance(pFile,0,0);
				if (t1 == 1) {
					if (i1 < i2-40) return FALSE;
					return TRUE;
				}
				if (i1 < i2-36) return FALSE;
				if (i1 < i2-5) {
					dx = pFile->MainLine.points_x[0] - pFile->MainLine.points_x[1];
					dy = pFile->MainLine.points_y[0] - pFile->MainLine.points_y[1];
					if (dx*dx+dy*dy > 45*45) return FALSE;
				}
			}
		}
		return TRUE;
	}
	if (t1 == 3) {
		if (t2 > 8 || pScore1 <= 35 || pScore2 <= 35) return TRUE;
		if (nFCoreNum == 1 && nSCoreNum == 1) {
			i1 = get_distance(pFile,0,0);
			i2 = get_distance(pSearch,0,0);
			if (pSearch->nType == 3) {
				n1 = i1; i1 = i2; i2 = n1;
			}
			if (i1 > i2-30) return TRUE;
		}
	}
	if (t1 == 4) {
		if (t2 == 5 || t2 == 6) return FALSE;
	}
	if (t1 == 5 && t2 == 6) return FALSE;
	if (t1 == 6 && t2 == 7) return FALSE;
	if (t1 == 4 || t1 == 5)	{
		if (t2 == 7) {
			if (pFile->nType == 7 && nFDeltaNum == 1) {
				i = op_func_01(FileDelta[0].x,FileDelta[0].y,FileCore[0].x,FileCore[0].y);
				i1 = i - FileCore[0].dir;
				if (i1 > 0)	{ 
					if (pSearch->nType == 4) return FALSE;
				}
				else if (i1 < 0) {
					if (pSearch->nType == 5) return FALSE;
				}
			}
			else if (pSearch->nType == 7 && nSDeltaNum == 1) {
				i = op_func_01(SearchDelta[0].x,SearchDelta[0].y,SearchCore[0].x,SearchCore[0].y);
				i2 = i - SearchCore[0].dir;
				if (i2 > 0)	{ 
					if (pFile->nType == 4) return FALSE;
				}
				else if ( i2 < 0 ) {
					if (pFile->nType == 5) return FALSE;
				}
			}
		}
	}
	if (t1 == 7 && t2 == 7 && nFDeltaNum == 1 && nSDeltaNum == 1) {
		i = op_func_01(FileDelta[0].x,FileDelta[0].y,FileCore[0].x,FileCore[0].y);
		i1 = i - FileCore[0].dir;
		i = op_func_01(SearchDelta[0].x,SearchDelta[0].y,SearchCore[0].x,SearchCore[0].y);
		i2 = i - SearchCore[0].dir;
		if (i1 > 0 && i2 < 0) return FALSE;
		if (i1 < 0 && i2 > 0) return FALSE;
	}
	return TRUE;
}

void get_tag_item(LPMPVECTEX pVect,BAR* pBar)
{
	int dir = op_func_01( pVect->item[pBar->nID2].x,
						  pVect->item[pBar->nID2].y,
						  pVect->item[pBar->nID1].x,
						  pVect->item[pBar->nID1].y);
	pBar->nSlope = dir;
	pBar->nSlope = ANGLE_120(pBar->nSlope);

	pBar->nDiff1 = dir - pVect->item[pBar->nID1].dir;
	pBar->nDiff1 = ANGLE_0(pBar->nDiff1);

	dir += 120;
	dir = ANGLE_240(dir);

	pBar->nDiff2 = dir - pVect->item[pBar->nID2].dir;
	pBar->nDiff2 = ANGLE_0(pBar->nDiff2);
}

void get_search_tag(LPFPVECTEX pSearch,BARVECT* pSearchBar,
				  int *nMaxSearchBarLen, int* SDiffField,
				  int *SArrangBarPtr,int nMinLen,int nMaxLen)
{
	MPVECTEX *pS = &pSearch->Mp;
	BAR *pBarItem = &pSearchBar->item[0];
	
	int nRealMaxNum = 300,  nMaxNumTh = MAX_SEARCH_BAR_NUM;
	int pSum[MAX_BAR_LENGTH];
	int i, j, len, sum, id, nNum = 0, nBarNum = 0;
	int	nMax = nMaxLen*nMaxLen, nMin = nMinLen*nMinLen;
	int tmp;
	BOOL flag = FALSE;

	BARTEMP pTemp[MAX_SEARCH_BAR_NUM];
	int pTmpPtr[MAX_SEARCH_BAR_NUM];
	
	pSearchBar->nNumber = 0;
	*nMaxSearchBarLen = 0;

	if (nMinLen > nMaxLen)  return;
	if (nMaxLen > MAX_BAR_LENGTH)  return;
	
	if (pS->quality <= 45)  nMaxNumTh = 350;

	for (i = 0; i < pS->nNumber-1; i++)	{
		if (pS->item[i].score < 15) continue;
		for (j = i+1; j < pS->nNumber; j++)	{
			if (pS->item[j].score < 15) continue;
			len = (pS->item[i].x-pS->item[j].x)*(pS->item[i].x-pS->item[j].x) + (pS->item[i].y-pS->item[j].y)*(pS->item[i].y-pS->item[j].y);
			if (len < nMin || len >= nMax) continue;
			pTemp[nNum].nLen = op_func_02(len);
			pTemp[nNum].nID1 = i;
			pTemp[nNum++].nID2 = j;

			if (nNum >= nMaxNumTh)  break;
		}
		if (j < pS->nNumber) break;
	}
	
	if (2*nNum >= nRealMaxNum) {
		memset( pSum, 0, sizeof(int)*MAX_BAR_LENGTH );
		for (i = 0; i < nNum; i++)  pSum[pTemp[i].nLen]++;
		for (i = 1; i < MAX_BAR_LENGTH; i++)  pSum[i] += pSum[i-1];
		for (i = 1; i < MAX_BAR_LENGTH; i++)  pSum[i]--;
		for (i = 0; i < nNum; i++) {
			len = pTemp[i].nLen;
			sum = pSum[len];
			pSum[len] = sum - 1;
			pTmpPtr[sum] = i;
		}
		if (pS->quality > 35)  nNum = nRealMaxNum / 2;
		else {
			if (nNum > nRealMaxNum)  nNum = nRealMaxNum;
		}
	}
	else {
		for (i = 0; i < nNum; i++) pTmpPtr[i] = i;
	}
	
	memset( SDiffField, 0, sizeof(int)*240 );
	nBarNum = 0;
	for (i = 0; i < nNum; i++) {
		len = pTemp[pTmpPtr[i]].nLen;
		pBarItem[nBarNum].nLen = len;
		if (*nMaxSearchBarLen < len)  *nMaxSearchBarLen = len + 1;
		pBarItem[nBarNum].nID1 = pTemp[pTmpPtr[i]].nID1;
		pBarItem[nBarNum].nID2 = pTemp[pTmpPtr[i]].nID2;
		get_tag_item( pS, &pBarItem[nBarNum] );

		if (abs(pBarItem[nBarNum].nDiff1 - pBarItem[nBarNum].nDiff2) > ANGLE) {
			if (pBarItem[nBarNum].nDiff1 - pBarItem[nBarNum].nDiff2 > ANGLE) {
				pBarItem[nBarNum].nID1 = pTemp[pTmpPtr[i]].nID2;
				pBarItem[nBarNum].nID2 = pTemp[pTmpPtr[i]].nID1;
				tmp = pBarItem[nBarNum].nDiff1;
				pBarItem[nBarNum].nDiff1 = pBarItem[nBarNum].nDiff2;
				pBarItem[nBarNum].nDiff2 = tmp;	
			}
			flag = TRUE;
		}

		id = pBarItem[nBarNum].nDiff1;
		SArrangBarPtr[id*MAX_SAMEDIFF_NUM + SDiffField[id]] = nBarNum;

		if (++SDiffField[id] == MAX_SAMEDIFF_NUM) { SDiffField[id]--; }
		if (++nBarNum >= MAX_SEARCH_BAR_NUM)  break;
		if (flag == TRUE) { flag = FALSE;  continue; }

		pBarItem[nBarNum].nLen = pTemp[pTmpPtr[i]].nLen;
		pBarItem[nBarNum].nID1 = pTemp[pTmpPtr[i]].nID2;
		pBarItem[nBarNum].nID2 = pTemp[pTmpPtr[i]].nID1;
		pBarItem[nBarNum].nDiff1 = pBarItem[nBarNum-1].nDiff2;
		pBarItem[nBarNum].nDiff2 = pBarItem[nBarNum-1].nDiff1;
		pBarItem[nBarNum].nSlope = pBarItem[nBarNum-1].nSlope;

		id = pBarItem[nBarNum].nDiff1;
		SArrangBarPtr[id*MAX_SAMEDIFF_NUM + SDiffField[id]] = nBarNum;

		if (++SDiffField[id] == MAX_SAMEDIFF_NUM) SDiffField[id]--;  
		if (++nBarNum >= MAX_SEARCH_BAR_NUM) break;
	}

	pSearchBar->nNumber = nBarNum;
}

void get_file_tag(LPFPVECTEX pFile,BARVECT* pFileBar,int* FDiffField,int *FArrangBarPtr,
				int* nFileCX,int* nFileCY,int nMinLen,int nMaxLen)
{
	MPVECTEX *pF = &pFile->Mp;
	BAR *pBarItem = &pFileBar->item[0];

	int nMaxNum = 600;
	int nMin = (nMinLen - LENGTH)*(nMinLen - LENGTH), nMax = (nMaxLen+LENGTH)*(nMaxLen+LENGTH);
	int nXMax = 0, nYMax = 0, nXMin = 10000, nYMin = 10000;
	int i, j, len, tmp, id, nBarNum = 0;

	memset( FDiffField, 0, sizeof(int)*240 );

	for (i = 0; i < pF->nNumber-1; i++)	{
		if (pF->item[i].score < 15) continue;
		for (j = i+1; j < pF->nNumber; j++)	{
			if (pF->item[j].score < 15) continue;

			len = (pF->item[i].x-pF->item[j].x)*(pF->item[i].x-pF->item[j].x) + (pF->item[i].y-pF->item[j].y)*(pF->item[i].y-pF->item[j].y);
			if (len <= nMin || len >= nMax) continue;
			pBarItem[nBarNum].nLen = op_func_02(len);

			pBarItem[nBarNum].nID1 = i;
			pBarItem[nBarNum].nID2 = j;

			get_tag_item(pF, &pBarItem[nBarNum]);

			if (pBarItem[nBarNum].nDiff1 > pBarItem[nBarNum].nDiff2) {
				tmp = pBarItem[nBarNum].nID1;
				pBarItem[nBarNum].nID1 = pBarItem[nBarNum].nID2;
				pBarItem[nBarNum].nID2 = tmp;

				tmp = pBarItem[nBarNum].nDiff1;
				pBarItem[nBarNum].nDiff1 = pBarItem[nBarNum].nDiff2;
				pBarItem[nBarNum].nDiff2 = tmp;			
			}

			id = pBarItem[nBarNum].nDiff1;
			FArrangBarPtr[id*MAX_SAMEDIFF_NUM + FDiffField[id]] = nBarNum;

			if (++FDiffField[id] == MAX_SAMEDIFF_NUM) FDiffField[id]--; 
			if (++nBarNum >= nMaxNum)  break;
		}
		if (j < pF->nNumber)  break;
	}
	pFileBar->nNumber = nBarNum;

	for (i = 0; i < pF->nNumber; i++) {
		if (nXMin > pF->item[i].x) nXMin = pF->item[i].x;
		if (nXMax < pF->item[i].x) nXMax = pF->item[i].x;
		if (nYMin > pF->item[i].y) nYMin = pF->item[i].y;
		if (nYMax < pF->item[i].y) nYMax = pF->item[i].y;
	}
	*nFileCX = (nXMin+nXMax) / 2; *nFileCY = (nYMin+nYMax) / 2;
}

int get_score_sub(LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	int i, j, x, y, dir, dx, dy, diff, val, minval, score = 0;

	for (i = 0; i < pVect1->nNumber; i++) {
		x = pVect1->item[i].x; y = pVect1->item[i].y;
		dir = pVect1->item[i].dir;
		minval = 10000;
		for (j = 0; j < pVect2->nNumber; j++) {
			dx = abs(pVect2->item[j].x - x);
			if (dx > LENGTH+3) continue;
			dy = abs(pVect2->item[j].y - y);
			if (dy > LENGTH+3) continue;
			diff = abs(pVect2->item[j].dir - dir);
			diff = IANGLE_120(diff);
			if (diff > ANGLE) continue;
			val = dx + dy + diff;
			if (minval > val) minval = val;
			if (minval < 20) break;
		}
		if (minval < 35) score += (35 - minval);
	}
	val = (pVect1->nNumber + pVect2->nNumber) / 2;
	if (val == 0) return (0);
	score = (score * 100 ) / val;
	return (score);
}

int get_point_score(LPFPVECTEX pFile,LPFPVECTEX pSearch) 
{
	int i, j, k, cx, cy, dx, dy, rot;
	int retScore = 0, score, maxscore;
	int nNumF = 0, nMaxF = 0, nMaxS = 0, coreflag = 0;
	int fcx = 0, fcy = 0, scx = 0, scy = 0, flen = 0, slen, fdir = 0, fdir0 = 0, fdir1 = 0, sdir = 0, sdir0, sdir1, difdir, difdir1;
	int scoreF, scoreS, nTh = 18, AngTh = 30;
	MINUTIAEX tmpF[7];
	MPVECTEX tmpMp = pFile->Mp;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	if ( pFile->Mp.nNumber < MIN_MINUTIA_NUMBER || pSearch->Mp.nNumber < MIN_MINUTIA_NUMBER ) return (0);

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);

	for (i = 0; i < pFile->Mp.nNumber; i++)	{
		if (nMaxF < pFile->Mp.item[i].curv) nMaxF = pFile->Mp.item[i].curv;
	}
	for (i = 0; i < pSearch->Mp.nNumber; i++) {
		if ( nMaxS < pSearch->Mp.item[i].curv ) nMaxS = pSearch->Mp.item[i].curv;
	}
	if (nMaxF > nMaxS) nMaxF = nMaxS;
	scoreF = pFile->Mp.quality; scoreS = pSearch->Mp.quality;
	if (scoreF < 35 || scoreS < 35) { nTh = 30; AngTh = 40; }

	for (i = 0, j = 0; i < pFile->Mp.nNumber; i++) {
		if (pFile->Mp.item[i].curv >= nMaxF+8) continue;
		for (k = 0; k < j; k++)	{
			dx = pFile->Mp.item[i].x - tmpF[k].x;
			dy = pFile->Mp.item[i].y - tmpF[k].y;
			if (dx*dx+dy*dy <= 30*30) break;
		}
		if (k < j) continue;
		tmpF[j++] = pFile->Mp.item[i];
		if (j >= 7) break;
	}
	nNumF = j;
	pFile->Mp = tmpMp;

	if (nFCoreNum > 0) {
		fcx = fcy = 0;
		for (i = 0; i < nFCoreNum; i++)	{
			fcx += FileCore[i].x; fcy += FileCore[i].y;
		}
		fcx /= nFCoreNum; fcy /= nFCoreNum;
		if (nFCoreNum == 2)
			fdir = op_func_01(FileCore[0].x,FileCore[0].y,FileCore[1].x,FileCore[1].y);
	}
	if (nSCoreNum > 0) {
		scx = scy = 0;
		for (i = 0; i < nSCoreNum; i++)	{
			scx += SearchCore[i].x; scy += SearchCore[i].y;
		}
		scx /= nSCoreNum; scy /= nSCoreNum;
		if (nSCoreNum == 2)
			sdir = op_func_01(SearchCore[0].x,SearchCore[0].y,SearchCore[1].x,SearchCore[1].y);
	}
	if (nFCoreNum == nSCoreNum && nFCoreNum > 0) coreflag = 1;

	for (i = 0; i<nNumF; i++) {
		cx = tmpF[i].x; cy = tmpF[i].y;
		if (coreflag == 1) flen = op_func_02((cx-fcx)*(cx-fcx)+(cy-fcy)*(cy-fcy));
		if (nFCoreNum == 1)	{
			fdir0 = tmpF[i].dir - FileCore[0].dir;
			fdir0 = ANGLE_0(fdir0);
		}
		if (nFCoreNum == 2)	{
			fdir0 = tmpF[i].dir - fdir;
			fdir0 = ANGLE_0(fdir0);
			fdir1 = fdir0 + 120;
			fdir1 = ANGLE_240(fdir1);
		}
		maxscore = 0;
		for (j = 0; j < pSearch->Mp.nNumber; j++) {
			if (abs(tmpF[i].curv-pSearch->Mp.item[j].curv) >= 6) continue;
			rot = pSearch->Mp.item[j].dir - tmpF[i].dir;
			rot = ANGLE_0(rot);
			if (nSCoreNum == 1) {
				sdir0 = pSearch->Mp.item[j].dir - SearchCore[0].dir;
				sdir0 = ANGLE_0(sdir0);
				if (nFCoreNum == 1) {
					difdir = abs(fdir0 - sdir0);
					difdir = IANGLE_120(difdir);
					if (difdir > AngTh) continue;
				}
				if (nFCoreNum == 2) {
					difdir = abs(fdir0 - sdir0);
					difdir = IANGLE_120(difdir);
					difdir1 = abs(fdir1 - sdir0);
					difdir1 = IANGLE_120(difdir1);
					if (difdir > AngTh && difdir1 > AngTh) continue;
				}
			}
			if (nSCoreNum == 2)	{
				sdir0 = pSearch->Mp.item[j].dir - sdir;
				sdir0 = ANGLE_0(sdir0);
				sdir1 = sdir0 + 120;
				sdir1 = ANGLE_240(sdir1);
				if (nFCoreNum == 1)	{
					difdir = abs(fdir0 - sdir0);
					difdir = IANGLE_120(difdir);
					difdir1 = abs(fdir0 - sdir1);
					difdir1 = IANGLE_120(difdir1);
					if (difdir > AngTh && difdir1 > AngTh) continue;
				}
				if (nFCoreNum == 2)	{
					difdir = abs(fdir0 - sdir0);
					difdir = IANGLE_120(difdir);
					difdir1 = abs(fdir0 - sdir1);
					difdir1 = IANGLE_120(difdir1);
					if (difdir > AngTh && difdir1 > AngTh) continue;
					difdir = abs(fdir1 - sdir0);
					difdir = IANGLE_120(difdir);
					difdir1 = abs(fdir1 - sdir1);
					difdir1 = IANGLE_120(difdir1);
					if (difdir > AngTh && difdir1 > AngTh) continue;
				}
			}
			if (coreflag == 1) {
				slen = op_func_02((pSearch->Mp.item[j].x-scx)*(pSearch->Mp.item[j].x-scx) +
						(pSearch->Mp.item[j].y-scy)*(pSearch->Mp.item[j].y-scy));
				if (abs(flen - slen) > nTh) continue;
			}

			dx = pSearch->Mp.item[j].x - cx; 
			dy = pSearch->Mp.item[j].y - cy;
			transform_mp(&pFile->Mp,cx,cy,rot,dx,dy);
			score = get_score_sub(&pFile->Mp,&pSearch->Mp);
			if (maxscore < score) maxscore = score;
			if (maxscore > 1700) return (retScore);
			pFile->Mp = tmpMp;
		}
		if (retScore < maxscore) retScore = maxscore;
	}
	return (retScore);
}

int get_point_score2(LPFPVECTEX pFile,LPFPVECTEX pSearch,int *nBlkScore) 
{
	int i, j, k, cx, cy, dx, dy, rot;
	int retScore = 0, score, maxscore, maxrot = 0, maxdx = 0, maxdy = 0;
	int nNumF = 0, nMaxF = 0, nMaxS = 0, scoreF, scoreS;
	int mRot = 0, mDx = 0, mDy = 0, mCx = 0, mCy = 0;
	MINUTIAEX tmpF[7];
	MPVECTEX tmpMp = pFile->Mp;
	BLOCKVECT tmpBlk;

	*nBlkScore = 0;
	if (pFile->Mp.nNumber < MIN_MINUTIA_NUMBER || pSearch->Mp.nNumber < MIN_MINUTIA_NUMBER) return (0);

	for (i = 0; i < pFile->Mp.nNumber; i++)	{
		if (nMaxF < pFile->Mp.item[i].curv) nMaxF = pFile->Mp.item[i].curv;
	}
	for (i = 0; i < pSearch->Mp.nNumber; i++) {
		if (nMaxS < pSearch->Mp.item[i].curv) nMaxS = pSearch->Mp.item[i].curv;
	}
	if (nMaxF > nMaxS) nMaxF = nMaxS;
	scoreF = pFile->Mp.quality; scoreS = pSearch->Mp.quality;
	if (scoreF < 25 || scoreS < 25) return (0); 
	for (i = 0,j = 0; i < pFile->Mp.nNumber; i++) {
		if (pFile->Mp.item[i].score < 30) continue;
		if (pFile->Mp.item[i].curv >= nMaxF+8) continue;
		for (k = 0; k < j; k++)	{
			dx = pFile->Mp.item[i].x - tmpF[k].x;
			dy = pFile->Mp.item[i].y - tmpF[k].y;
			if (dx*dx+dy*dy <= 30*30) break;
		}
		if (k < j) continue;
		tmpF[j++] = pFile->Mp.item[i];
		if (j >= 7) break;
	}
	nNumF = j;
	pFile->Mp = tmpMp;

	for (i = 0; i < nNumF; i++)	{
		cx = tmpF[i].x; cy = tmpF[i].y;
		maxscore = 0;
		for (j = 0; j < pSearch->Mp.nNumber; j++) {
			if (abs(tmpF[i].curv-pSearch->Mp.item[j].curv) >= 6) continue;
			if (pSearch->Mp.item[j].score < 30) continue;
			rot = pSearch->Mp.item[j].dir - tmpF[i].dir;
			rot = ANGLE_0(rot);
			dx = pSearch->Mp.item[j].x - cx; 
			dy = pSearch->Mp.item[j].y - cy;
			transform_mp(&pFile->Mp,cx,cy,rot,dx,dy);
			score = get_score_sub(&pFile->Mp,&pSearch->Mp);
			if (maxscore < score) {
				maxscore = score; maxrot = rot; maxdx = dx; maxdy = dy;
			}
			pFile->Mp = tmpMp;
		}
		if (retScore < maxscore) {
			retScore = maxscore;
			mRot = maxrot; mDx = maxdx; mDy = maxdy; mCx = cx; mCy = cy; 
		}
	}
	if (retScore == 0) return (retScore);
	tmpBlk = pFile->Block;
	transform_block(mRot,mDx,mDy,mCx,mCy,&tmpBlk);
	*nBlkScore = check_block(30,4,&tmpBlk,&pSearch->Block);
	return (retScore);
}

BOOL arrange_points_sub(int cx,int cy,int nRot,int nXDiff,int nYDiff,
						 LPMPVECTEX pFile,BLOCKVECT* tmpFBlock,LPFPVECTEX pSearch,
						 BARVECT *pSBar,int *nMaxSBarLen,
						 int *SBarPtr,int *SDiffField)
{
	MPVECTEX tmpFMp = *pFile;
	int i, j, k, x, y, bx, by, nNum, nRow, nCol, nDiv;
	int nNewId[MAX_MINUTIA_NUMBER];
	int nList[MAX_BLOCK_NUMBER];
  
	transform_mp(&tmpFMp,cx,cy,nRot,nXDiff,nYDiff);
	nNum = 0;
	for (i = 0; i < tmpFBlock->nCol*tmpFBlock->nRow; i++) {
		if (tmpFBlock->Data[i] == 0xFF) continue;
		if (pSearch->Block.Data[i] == 0xFF) continue;
		nList[nNum++] = i;
	}
	if (nNum == 0) return FALSE;
	for (i = 0; i < MAX_MINUTIA_NUMBER; i++) nNewId[i] = -1;
	nDiv = tmpFBlock->nCol;
	for (i = 0, k = 0; i < tmpFMp.nNumber; i++)	{
		x = tmpFMp.item[i].x; y = tmpFMp.item[i].y;
		for (j = 0; j < nNum; j++) {
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if (abs(x-bx) > LENGTH+8 || abs(y-by) > LENGTH+8) continue;
			nNewId[i] = k;
			pFile->item[k++] = pFile->item[i];
			break;
		}
	}
	pFile->nNumber = k;
	for (i = 0; i<MAX_MINUTIA_NUMBER; i++) nNewId[i] = -1;
	nDiv = pSearch->Block.nCol;
	for (i = 0, k = 0; i < pSearch->Mp.nNumber; i++) {
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		for (j = 0; j < nNum; j++) {
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if (abs(x-bx) > LENGTH+8 || abs(y-by) > LENGTH+8) continue;
			nNewId[i] = k;
			pSearch->Mp.item[k++] = pSearch->Mp.item[i];
			break;
		}
	}
	pSearch->Mp.nNumber = k;
	for (i = 0,j = 0; i < pSBar->nNumber; i++) {
		if (nNewId[pSBar->item[i].nID1] == -1) continue;
		if (nNewId[pSBar->item[i].nID2] == -1) continue;
		pSBar->item[j] = pSBar->item[i];
		pSBar->item[j].nID1 = nNewId[pSBar->item[i].nID1];
		pSBar->item[j++].nID2 = nNewId[pSBar->item[i].nID2];
	}
	if (j == 0) return FALSE;
	pSBar->nNumber = j;

	memset(SDiffField, 0, sizeof(int)*240);
	*nMaxSBarLen = 0;
	for (i = 0; i < pSBar->nNumber; i++) {
		if (*nMaxSBarLen < pSBar->item[i].nLen)
			*nMaxSBarLen = pSBar->item[i].nLen + 1;
		k = pSBar->item[i].nDiff1;
		SBarPtr[k*MAX_SAMEDIFF_NUM+SDiffField[k]] = i;
		if (++SDiffField[k] == MAX_SAMEDIFF_NUM) SDiffField[k]--; 
	}
	return TRUE;
}

int arrange_points(LPFPVECTEX pFile,LPFPVECTEX pSearch,BARVECT *pSearchBar,int *nMaxSBarLen,
				   int *SBarPtr,int *SDiffField,BOOL nCommonFlag)
{
	MPVECTEX tmpMP;
	BLOCKVECT tmpBlk;
	int i, k1, k2, nBlkScore, nBlkScore1;
	int fcx, fcy, scx, scy, dx, dy, fdir, sdir, rot;
	int maxrot = 0, maxdx = 0, maxdy = 0, cx = 0, cy = 0, score;
	COREITEMEX FileCore[2], SearchCore[2];
	int nFCoreNum, nSCoreNum;
	BOOL nTypeFlag = TRUE;
	FPVECTEX tmpFile = *pFile;

	if (pFile->Mp.nNumber == 0 || pSearch->Mp.nNumber == 0) return(-1);

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);
	if ( nFCoreNum == 0 || nSCoreNum == 0 ) return (0);

	if (nFCoreNum == nSCoreNum)	{
		fcx = fcy = scx = scy = 0;
		for (i = 0; i < nFCoreNum; i++)	{
			fcx += FileCore[i].x; fcy += FileCore[i].y;
			scx += SearchCore[i].x; scy += SearchCore[i].y;
		}
		fcx /= nFCoreNum; fcy /= nFCoreNum;
		scx /= nFCoreNum; scy /= nFCoreNum;
		dx = scx - fcx; dy = scy - fcy;
		if (nFCoreNum == 1)	{
			rot = SearchCore[0].dir - FileCore[0].dir;
			rot = ANGLE_0(rot);
		}
		else {
			fdir = op_func_01(FileCore[0].x,FileCore[0].y,FileCore[1].x,FileCore[1].y);
			sdir = op_func_01(SearchCore[0].x,SearchCore[0].y,SearchCore[1].x,SearchCore[1].y);
			rot = sdir - fdir;
			rot = ANGLE_0(rot);
			tmpMP = pFile->Mp;
			transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
			k1 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
			maxrot = rot;
			rot = 120 + rot;
			rot = ANGLE_240(rot);
			tmpMP = pFile->Mp;
			transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
			k2 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
			if (k1 > k2) rot = maxrot;
		}
		cx = fcx; cy = fcy;
	}
	else {
		if (nFCoreNum == 1)	{
			fcx = FileCore[0].x; fcy = FileCore[0].y;
			fdir = FileCore[0].dir;
			k2 = 0;
			for (i = 0; i < 2; i++)	{
				dx = SearchCore[i].x - fcx; dy = SearchCore[i].y - fcy;
				rot = SearchCore[i].dir - fdir;
				rot = ANGLE_0(rot);
				tmpMP = pFile->Mp;
				transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
				k1 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
				if (k2 < k1) {
					k2 = k1; maxrot = rot; maxdx = dx; maxdy = dy;
				}
			}
			cx = fcx; cy = fcy;
		}
		else {
			scx = SearchCore[0].x; scy = SearchCore[0].y;
			sdir = SearchCore[0].dir;
			k2 = 0;
			for (i = 0; i < 2; i++)	{
				fcx = FileCore[i].x; fcy = FileCore[i].y;
				dx = scx - fcx; dy = scy - fcy;
				rot = sdir - FileCore[i].dir;
				rot = ANGLE_0(rot);
				tmpMP = pFile->Mp;
				transform_mp(&tmpMP,fcx,fcy,rot,dx,dy);
				k1 = get_matched_mp_num(LENGTH,7,&tmpMP,&pSearch->Mp);
				if (k2 < k1) {
					k2 = k1; maxrot = rot; maxdx = dx; maxdy = dy;
					cx = fcx; cy = fcy;
				}
			}
		}
		if (k2 == 0) return 0;
		rot = maxrot; dx = maxdx; dy = maxdy;
	}

	if ( nCommonFlag) {
		tmpBlk = pFile->Block;
		transform_block(rot,dx,dy,cx,cy,&tmpBlk);
		nBlkScore = check_block(30,5,&tmpBlk,&pSearch->Block);
		if (nBlkScore > 80)	{
			if (pFile->nType == 8 || pSearch->nType == 8 || (pFile->nType == pSearch->nType && pFile->nType < 2)){
				score = get_point_score2(&tmpFile,pSearch,&nBlkScore1);
				if (nBlkScore1 >= nBlkScore && ((nBlkScore1 >= 92 && score > 700) || 
					(pFile->nType == 1 && pSearch->nType == 1 && score > 1000))) nTypeFlag = FALSE;
			}
			if (nTypeFlag)
				arrange_points_sub(cx,cy,rot,dx,dy,&pFile->Mp,&tmpBlk,pSearch,pSearchBar,nMaxSBarLen,SBarPtr,SDiffField);
		}
	}
	return(0);
}

int rotate_points(int cx,int cy,int* pAngle,BARVECT* pBar,LPFPVECTEX pVect)
{
	int i, j, sum, nMax, nMaxId, sumN1, sumN2, temp[300];
	int rotAngle, rot, nCos, nSin, x, y, angle;

	for (i = 0; i < 240; i++) { 
		sum = 0;
		for (j = i-4; j <= i+4; j++) {
			rot = j;
			rot = ANGLE_0_240(rot);
			sum += pAngle[rot];
		}
		temp[i] = sum;
	}
	memcpy( pAngle, temp, sizeof(int)*240 );
	nMax = 0; nMaxId = 0;
	for (i = 0; i < 240; i++) {
		if (pAngle[i] > nMax) {
			nMax = pAngle[i]; nMaxId = i;
		}
	}


    for (i = 0; i < 10; i++) temp[i] = pAngle[230+i];
	for (i = 0; i < 240; i++) temp[i+10] = pAngle[i];
	for (i = 0; i < 10; i++) temp[i+250] = pAngle[i];

	nMax /= 2;
	sumN1 = sumN2 = 0;
	for (i = nMaxId; i < nMaxId+20; i++) {
		if (temp[i] <= nMax) continue;
		if (temp[i] <= 20) continue;
		sumN1 += (temp[i] - nMax) * i;
		sumN2 += (temp[i] - nMax);
	}
	if (sumN2 == 0) rotAngle = 0;
	else
		rotAngle = (100*sumN1/sumN2 + 50) / 100;
	
	rotAngle -= 10;
	rotAngle = ANGLE_0_240(rotAngle);

	rot = 240 - rotAngle;
	rot = ANGLE_240(rot);
	nCos = _table_03[rot]; nSin = _table_04[rot];
	for (i = 0; i < pVect->Mp.nNumber; i++)	{
		x = (pVect->Mp.item[i].x-cx)*nCos + (pVect->Mp.item[i].y-cy)*nSin;
		x = ROUND(x);
		x = x >> 14;
		y = (pVect->Mp.item[i].y-cy)*nCos - (pVect->Mp.item[i].x-cx)*nSin;
		y = ROUND(y);
		y = y >> 14;

		pVect->Mp.item[i].x = x + cx;
		pVect->Mp.item[i].y = y + cy;

		angle = pVect->Mp.item[i].dir + rotAngle;
		angle = ANGLE_0_240(angle);
		
		pVect->Mp.item[i].dir = angle;
	}
	for (i = 0; i < pBar->nNumber; i++)	{
		angle = pBar->item[i].nSlope + rotAngle;
		angle = ANGLE_0_240(angle);
		angle = ANGLE_120(angle);
		pBar->item[i].nSlope = angle;
	}

	for (i = 0; i < pVect->Core.nNumber; i++) {
		x = (pVect->Core.item[i].x-cx)*nCos + (pVect->Core.item[i].y-cy)*nSin;
		x = ROUND(x);
		x = x >> 14;
		y = (pVect->Core.item[i].y-cy)*nCos - (pVect->Core.item[i].x-cx)*nSin;
		y = ROUND(y);
		y = y >> 14;

		pVect->Core.item[i].x = x + cx;
		pVect->Core.item[i].y = y + cy;

		angle = pVect->Core.item[i].dir + rotAngle;
		angle = ANGLE_0_240(angle);
		pVect->Core.item[i].dir = angle;
	}
	
	return(rotAngle);
}

void get_shift_param(int nTH,int nScore,BAR* pFBar,BAR* pSBar, 
				   int* XField,int* YField,LPMPVECTEX pFile,LPMPVECTEX pSearch)
{
	int nSid1, nSid2, nFid1, nFid2, dx1, dx2, dy1, dy2, dx, dy;
	if (nScore == 0) return;

	nSid1 = pSBar->nID1; nSid2 = pSBar->nID2;
	nFid1 = pFBar->nID1; nFid2 = pFBar->nID2;
	dx1 = pSearch->item[nSid1].x - pFile->item[nFid1].x;
	dx2 = pSearch->item[nSid2].x - pFile->item[nFid2].x;
	dy1 = pSearch->item[nSid1].y - pFile->item[nFid1].y;
	dy2 = pSearch->item[nSid2].y - pFile->item[nFid2].y;

	dx = abs(dx1 - dx2);
	if (dx >= nTH) return;
	dy = abs(dy1 - dy2);
	if (dy >= nTH) return;
	if (abs(dx2) >= 640) return;
	if (abs(dy2) >= 640) return;
	if (abs(dx1) >= 640) return;
	if (abs(dy1) >= 640) return;

	dx = dx1 + dx2;
	if (dx < 0) dx++;
	dx /= 2;
	XField[dx+320] += nScore;
	dy = dy1 + dy2;
	if (dy < 0) dy++;
	dy /= 2;
	YField[dy+320] += nScore;
}

void shift_points(int* nXoffset,int* nYoffset,LPFPVECTEX pVect,int* XField,int* YField)
{
	int pTmp[640];
	int i, j, starti, endi, sum, nMax, nMaxId, nSum1, nSum2;
	memset( pTmp, 0, sizeof(int)*640 );

	for (i = 5; i < 640-5; i++)	{
		sum = 0;
		for (j = i-5; j < i+5; j++) sum += XField[j];
		pTmp[i] = sum;
	}
	memcpy ( XField, pTmp, sizeof(int)*640 );
	nMax = nMaxId = 0;
	for (i = 0; i < 640; i++) {
		if (XField[i] <= nMax) continue;
		nMax = XField[i]; nMaxId = i;
	}
	nSum1 = nSum2 = 0;
	starti = ( nMaxId-10 < 0 ) ? 0 : nMaxId-10;
	endi = ( nMaxId+10 > 640-1 ) ? 640-1 : nMaxId+10;
	nMax = (nMax*2) / 3;
	for (i = starti; i < endi; i++)	{
		if (XField[i] <= nMax) continue;
		nSum1 += XField[i] * i;
		nSum2 += XField[i];
	}
	*nXoffset = ( nSum2 == 0 ) ? 0 : (nSum1/nSum2) - 320;
	
	for (i = 5; i < 640-5; i++)	{
		sum = 0;
		for (j = i-5; j < i+5; j++) sum += YField[j];
		pTmp[i] = sum;
	}
	memcpy ( YField, pTmp, sizeof(int)*640 );
	nMax = nMaxId = 0;
	for (i = 0; i < 640; i++) {
		if (YField[i] <= nMax) continue;
		nMax = YField[i]; nMaxId = i;
	}
	nSum1 = nSum2 = 0;
	starti = ( nMaxId-10 < 0 ) ? 0 : nMaxId-10;
	endi = ( nMaxId+10 > 640-1 ) ? 640-1 : nMaxId+10;
	nMax = (nMax*2) / 3;
	for (i = starti; i < endi; i++)	{
		if (YField[i] <= nMax) continue;
		nSum1 += YField[i] * i;
		nSum2 += YField[i];
	}
	*nYoffset = ( nSum2 == 0 ) ? 0 : (nSum1/nSum2) - 320;

	for (i = 0; i < pVect->Mp.nNumber; i++)	{
		pVect->Mp.item[i].x += *nXoffset; pVect->Mp.item[i].y += *nYoffset;
	}
	for (i = 0; i < pVect->Core.nNumber; i++) {
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
	for (i = 0; i < pFile->Block.nCol*pFile->Block.nRow; i++) {
		if (pFile->Block.Data[i] == 0xFF) continue;
		if (pSearch->Block.Data[i] == 0xFF) continue;
		nList[nNum++] = i;
	}
	if (nNum == 0) return FALSE;
	for (i = 0; i < MAX_MINUTIA_NUMBER; i++) nNewId[i] = -1;
	nDiv = pFile->Block.nCol;
	for (i = 0, k = 0; i < pFile->Mp.nNumber; i++) {
		x = pFile->Mp.item[i].x; y = pFile->Mp.item[i].y;
		for (j = 0; j < nNum; j++) {
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if (abs(x-bx) > 8 || abs(y-by) > 8) continue;
			nNewId[i] = k;
			pFile->Mp.item[k++] = pFile->Mp.item[i];
			break;
		}
	}
	if (k == 0) return FALSE;
	pFile->Mp.nNumber = k;

	for (i = 0; i < MAX_BAR_NUM; i++) nNewBarId[i] = -1;
	for (i = 0, j = 0; i < pFBar->nNumber; i++)	{
		if (nNewId[pFBar->item[i].nID1] == -1) continue;
		if (nNewId[pFBar->item[i].nID2] == -1) continue;
		nNewBarId[i] = j;
		pFBar->item[j] = pFBar->item[i];
		pFBar->item[j].nID1 = nNewId[pFBar->item[i].nID1];
		pFBar->item[j++].nID2 = nNewId[pFBar->item[i].nID2];
	}
	if (j == 0) return FALSE;
	pFBar->nNumber = j;

	nEmptyPairNum = 0;
	for (i = 0; i < nPairNum; i++) {
		j = nNewBarId[PairList[i].fid];
		if (j == -1) nEmptyPair[nEmptyPairNum++] = i;
		else PairList[i].fid = j;
		if (nEmptyPairNum >= MAX_BAR_NUM) break;
	}
	
	for (i = 0, k = 0; i < *nLastNum; i++) {
		for (j = 0; j < nEmptyPairNum; j++) {
			if (LastList[i] == nEmptyPair[j]) break;
		}
		if (j >= nEmptyPairNum) nNewId[k++] = LastList[i];
	}
	for (i = 0; i < k; i++) LastList[i] = nNewId[i];
	*nLastNum = k;
	for (i = 0; i < MAX_MINUTIA_NUMBER; i++) nNewId[i] = -1;
	nDiv = pSearch->Block.nCol;
	for (i = 0, k = 0; i < pSearch->Mp.nNumber; i++) {
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		for (j = 0; j < nNum; j++) {
			nRow = nList[j] / nDiv; nCol = nList[j] % nDiv;
			bx = nCol*BLOCK_SIZE + BLOCK_SIZE/2;
			by = nRow*BLOCK_SIZE + BLOCK_SIZE/2;
			if (abs(x-bx) > 8 || abs(y-by) > 8) continue;
			nNewId[i] = k;
			pSearch->Mp.item[k++] = pSearch->Mp.item[i];
			break;
		}
	}
	if (k == 0) return FALSE;
	pSearch->Mp.nNumber = k;
	
	for (i = 0; i < MAX_BAR_NUM; i++) nNewBarId[i] = -1;
	for (i = 0,j = 0; i < pSBar->nNumber; i++) {
		if (nNewId[pSBar->item[i].nID1] == -1) continue;
		if (nNewId[pSBar->item[i].nID2] == -1) continue;
		nNewBarId[i] = j;
		pSBar->item[j] = pSBar->item[i];
		pSBar->item[j].nID1 = nNewId[pSBar->item[i].nID1];
		pSBar->item[j++].nID2 = nNewId[pSBar->item[i].nID2];
	}
	if (j == 0) return FALSE;
	pSBar->nNumber = j;

	nEmptyPairNum = 0;
	for (i = 0; i < nPairNum; i++) {
		j = nNewBarId[PairList[i].sid];
		if (j == -1) nEmptyPair[nEmptyPairNum++] = i;
		else PairList[i].sid = j;
		if (nEmptyPairNum >= MAX_BAR_NUM) break;
	}
	
	for (i = 0, k = 0; i < *nLastNum; i++) {
		for (j = 0; j < nEmptyPairNum; j++)	{
			if (LastList[i] == nEmptyPair[j]) break;
		}
		if (j >= nEmptyPairNum) nNewId[k++] = LastList[i];
	}
	for (i = 0; i < k; i++) LastList[i] = nNewId[i];
	*nLastNum = k;
	return TRUE;
}

BOOL check_limit(int nTH,BAR* pFBar,BAR* pSBar,LPFPVECTEX pF,LPFPVECTEX pS,int cx,int cy)
{
	int nSid = pSBar->nID1, nFid = pFBar->nID1;
	int dx, dy, th, len, nRidge, nVal;

	nRidge = pF->nRidgeDensity;
	if (nRidge < pS->nRidgeDensity) nRidge = pS->nRidgeDensity;
	
	if (nRidge > 200) nVal = 9;
	else  nVal = 10;

	dx = pF->Mp.item[nFid].x - cx; dy = pF->Mp.item[nFid].y - cy;
	len = op_func_02(dx*dx+dy*dy);
	th = ( len <= 150 ) ? nVal+len/50 : nVal+3;

	if (abs(pS->Mp.item[nSid].x - pF->Mp.item[nFid].x) >= th) return FALSE;
	if (abs(pS->Mp.item[nSid].y - pF->Mp.item[nFid].y) >= th) return FALSE;

	nSid = pSBar->nID2; nFid = pFBar->nID2;
	dx = pF->Mp.item[nFid].x - cx; dy = pF->Mp.item[nFid].y - cy;
	len = op_func_02(dx*dx+dy*dy);
	th = (len <= 150) ? nVal+len/50 : nVal+3;

	if (abs(pS->Mp.item[nSid].x - pF->Mp.item[nFid].x) >= th) return FALSE;
	if (abs(pS->Mp.item[nSid].y - pF->Mp.item[nFid].y) >= th) return FALSE;

	return(TRUE);
}

void get_paired_template(LPMPVECTEX pVect,int nPairNum,short* pPairID,FPTEMPLATE* Template)
{
	MPVECTEX tmpMp;
	MINUTIAEX tmp;
	int i, j, k, n, minid, minval, val, radTh;
	int cx, cy, cdir;
	int num = 0;

	memset(Template, 0, sizeof(FPTEMPLATE));

	for (i = 0; i < nPairNum; i++) {
		if (pVect->item[pPairID[i]].score < 20) continue;
		cx = pVect->item[pPairID[i]].x; cy = pVect->item[pPairID[i]].y; cdir = pVect->item[pPairID[i]].dir;
		radTh = 100; n = 0;
		while (1) {
			if (radTh > 200) break;
			for (k = 0; k < pVect->nNumber; k++) {
				if (pPairID[i] == k) continue;
				val = (pVect->item[k].x - cx)*(pVect->item[k].x - cx) + (pVect->item[k].y - cy)*(pVect->item[k].y - cy);
				if (val >= radTh*radTh) continue;
				tmpMp.item[n++] = pVect->item[k];
			}
			tmpMp.nNumber = n; 
			if (n < MAX_NEIGH_NUM){ radTh += 20; continue; }
			break;
		}
		if (tmpMp.nNumber < MAX_NEIGH_NUM) continue;

		for (k = 0; k < tmpMp.nNumber-1; k++) {
			minval = (cx-tmpMp.item[k].x)*(cx-tmpMp.item[k].x) + (cy-tmpMp.item[k].y)*(cy-tmpMp.item[k].y);
			minid = k;
			for (j = k+1; j < tmpMp.nNumber; j++) {
				val = (cx-tmpMp.item[j].x)*(cx-tmpMp.item[j].x) + (cy-tmpMp.item[j].y)*(cy-tmpMp.item[j].y);
				if (val >= minval) continue;
				minid = j; minval = val;
			}
			if (minid == k) continue;
			tmp = tmpMp.item[k]; tmpMp.item[k] = tmpMp.item[minid]; tmpMp.item[minid] = tmp;
		}

		Template->chunk[num].cdir = (BYTE)cdir;

		n = 0;
		for (k = 0; k < tmpMp.nNumber; k++)	{
			val = (cx-tmpMp.item[k].x)*(cx-tmpMp.item[k].x) + (cy-tmpMp.item[k].y)*(cy-tmpMp.item[k].y);
			Template->chunk[num].radius[n] = op_func_02(val);
			val = op_func_01(tmpMp.item[k].x,tmpMp.item[k].y,cx,cy);
			val = val - cdir;
			val = ANGLE_0(val);
			Template->chunk[num].angle[n] = (BYTE)val;
			val = tmpMp.item[k].dir - cdir;
			val = ANGLE_0(val);
			Template->chunk[num].delta[n] = (BYTE)val;
			n++;
			if (n >= MAX_NEIGH_NUM) break;
		}
		num++;
		if (num >= MAX_CHUNK_NUM) break;
	}
	Template->number = num;
}

int match_paired_chunk(DATACHUNK* pChunk1,DATACHUNK* pChunk2)
{
	int i, j, num = 0, val;
	int r1, r2, a1, a2, d1, d2;

	for (i = 0; i < MAX_NEIGH_NUM; i++)	{
		r1 = pChunk1->radius[i]; a1 = pChunk1->angle[i]; d1 = pChunk1->delta[i];
		for (j = 0; j < MAX_NEIGH_NUM; j++)	{
			r2 = pChunk2->radius[j]; a2 = pChunk2->angle[j]; d2 = pChunk2->delta[j];
			if (abs(r1-r2) > 10) continue;
			val = abs(a1 - a2);
			val = IANGLE_120(val);
			if (val >= ANGLE) continue;
			val = abs(d1 - d2);
			val = IANGLE_120(val);
			if (val < ANGLE) { num++; break; }
		}
		if (num >= 5) return (1);
	}
	return (0);
}

int match_template(FPTEMPLATE* Template1, FPTEMPLATE* Template2)
{
	int i, j, ret, num = 0;
	
	for (i = 0; i < Template1->number; i++)	{
		for (j = 0; j < Template2->number; j++)	{
			ret = match_paired_chunk( &(Template1->chunk[i]), &(Template2->chunk[j]) );
			if (ret == 1){ num++; break; }
		}
	}
	return (num);
}

int GetMatchedTemplateNum(LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT* pPair)
{
	int num;
	short *nFileID = pPair->nFileID, *nSearchID = pPair->nSearchID;
	FPTEMPLATE template1, template2;

	get_paired_template(pFile,pPair->nNumber,nFileID,&template1);
	get_paired_template(pSearch,pPair->nNumber,nSearchID,&template2);
	if (template1.number == 0 || template2.number == 0) return (-1);

	num = match_template(&template1,&template2);

	return num;
}

BOOL check_core(LPCOREVECTEX pFile,LPCOREVECTEX pSearch,int nLenTh,int nAngTh)
{
	int len, flen, slen, diff, fx, fy, fdir, sx, sy, sdir;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(pFile,FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(pSearch,SearchCore,NULL,NULL);

	if (nFCoreNum==0 || nSCoreNum==0) return FALSE;

	if (nFCoreNum != nSCoreNum) return FALSE;

	if (nFCoreNum == 1)	{
		fx = FileCore[0].x; fy = FileCore[0].y;
		sx = SearchCore[0].x; sy = SearchCore[0].y;
		len = op_func_02((fx-sx)*(fx-sx)+(fy-sy)*(fy-sy));
		fdir = FileCore[0].dir; sdir = SearchCore[0].dir;
	}
	else {
		flen = op_func_02((FileCore[0].x-FileCore[1].x)*(FileCore[0].x-FileCore[1].x)
					+ (FileCore[0].y-FileCore[1].y)*(FileCore[0].y-FileCore[1].y));
		slen = op_func_02((SearchCore[0].x-SearchCore[1].x)*(SearchCore[0].x-SearchCore[1].x)
					+ (SearchCore[0].y-SearchCore[1].y)*(SearchCore[0].y-SearchCore[1].y));
		len = abs(flen-slen);
		fdir = op_func_01(FileCore[0].x,FileCore[0].y,
						  FileCore[1].x,FileCore[1].y);
		fdir = ANGLE_120(fdir);
		sdir = op_func_01(SearchCore[0].x,SearchCore[0].y,
						  SearchCore[1].x,SearchCore[1].y);
		sdir = ANGLE_120(sdir);
	}
	diff = abs(fdir-sdir);
	diff = IANGLE_120(diff);
	if (len < nLenTh && diff < nAngTh) return TRUE;
	return FALSE;
}

int adjust_score(int inScore,int globalScore,int nBlockScore,int nCommonScore,int dScore,int nPairNum,int nRidge,
						int nBarNum,int minNum,BOOL nCoreCheck,int nTempNum,BOOL qFlag,BOOL cFlag)
{
	int outScore = inScore;

		if (nBarNum < 100) {
			if (nCommonScore <= 40 && nBlockScore <= 95 && nPairNum <= 5 && !nCoreCheck) 
				outScore = (inScore * 8) / 10;
			else if (nCommonScore < 60  && nTempNum < 5 && nBlockScore <= 93 && nPairNum <= 8 && !nCoreCheck) {
				if (globalScore < 650)
					outScore = (inScore * 6) / 10;
				else
					outScore = (inScore * 9) / 10;
			}
			else if (nCommonScore < 72 && nPairNum <= 5 && nPairNum*2 <= minNum && !cFlag)
				outScore = inScore;
			else if (nPairNum*100 < minNum*50 && nPairNum < 8 && !qFlag)
				outScore = (inScore * 9) / 10;
			else if ( nPairNum < 5 || nBlockScore <= 86 || nPairNum*100<minNum*42 || (( nRidge<232 || nPairNum*2 <= minNum) && nBlockScore<92)) 
				outScore = inScore;
			else if ( nPairNum <= 6 && ((nBarNum<50 && nRidge<=236) || (nPairNum*100<=minNum*50 && (nRidge <= 233 || nCommonScore < 67))) )
				outScore = (inScore * 14) / 10;
			else if ( nPairNum*100 < minNum*43 && nRidge < 235 && nBlockScore <= 92 )
				outScore = (inScore * 15) / 10;
			else if (nTempNum <= 2 && nPairNum < 12) {
				if (nTempNum == 0 && nPairNum*100 < minNum*50)
					outScore = inScore;
				else if (nTempNum == 0 && (nBlockScore < 95 && (!nCoreCheck || nPairNum > 7 || nRidge <= 215 )))
					outScore = (inScore * 12) / 10;
				else if (nCommonScore < 60)
					outScore = (inScore * 13) / 10;
				else 
					outScore = (inScore * 16) / 10; 
			}
			else if (nCommonScore < 60 && nBlockScore < 88 && nPairNum <= 8 && !nCoreCheck)
				outScore = (inScore * 15) / 10;
			else
				outScore = inScore * 2;
		}
		else {
			if (nPairNum*100 < minNum*40 && nPairNum <= 8 && !qFlag)
				outScore = inScore;	
			else if (nCommonScore < 68 && nBlockScore <= 92 && nPairNum <= 9 && !nCoreCheck) 
				outScore = (inScore * 100) / nBarNum;
			else if (globalScore < 900 && nPairNum < 10 && nPairNum*100 < minNum*34 && !cFlag)
				outScore = (inScore * 120) / nBarNum;
			else if (nBlockScore <= 77 || (nPairNum < 16 && nPairNum*100 < minNum*34 && nTempNum < 5)) 
				outScore = (inScore * 130) / nBarNum;
			else if (nPairNum < 12 && globalScore < 850 && nTempNum == 0 && nPairNum*100 < minNum*50)
				outScore = (inScore * 130) / nBarNum;
			else if ( nPairNum < 12 && globalScore < 1190 && ( nBlockScore <= 88 || nRidge < 225 ))	{
				if (nCommonScore < 67 && nTempNum <= 1)
					outScore = (inScore * 110) / nBarNum;
				else if (nCommonScore < 80 && nTempNum <= 3)
					outScore = (inScore * 130) / nBarNum;
				else
					outScore = (inScore * 140) / nBarNum;
			}
			else if (nTempNum <= 2 && nPairNum <= 12 && nBarNum < 200) {
				if (nTempNum == 0 && nPairNum <= 8 && nPairNum*100 < minNum*45)
					outScore = (inScore * 110) / nBarNum;
				else if (nBarNum < 135)
					outScore = (inScore * 130) / nBarNum;
				else if (nTempNum >= 1 && nPairNum >= 12 && nBlockScore >= 94)
					outScore = (inScore * 190) / nBarNum;
				else if (nCommonScore < 70 && nPairNum*100 < minNum*50)
					outScore = (inScore * 140) / nBarNum;
				else 
					outScore = (inScore * 150) / nBarNum;
			}
			else if (nBarNum < 200)	{
				if (nCommonScore < 73 && nRidge < 233 && !nCoreCheck)
					outScore = (inScore * 150) / nBarNum;
				else if (nPairNum*100 < minNum*45 && nRidge <= 240)
					outScore = (inScore * 160) / nBarNum;
				else if (nBarNum <= 130 || nCommonScore <= 72)
					outScore = (inScore * 160) / nBarNum;
				else  
					outScore = (inScore * 200) / nBarNum;
			}
			else {
				if (nPairNum < 12 && nPairNum*100 < minNum*50 && nCommonScore < 75)
					outScore = (inScore * 160) / nBarNum;
				else if (nTempNum <= 1 && nPairNum <= 10 && nPairNum*100 < minNum*42)
					outScore = (inScore * 170) / nBarNum;
				else if (nTempNum <= 2 && nPairNum <= 13 && nPairNum*100 < minNum*45)
					outScore = (inScore * 165) / nBarNum;
				else if (nTempNum <= 3 && nPairNum <= 13 && nPairNum*100 < minNum*55)
					outScore = (inScore * 170) / nBarNum;
				else
					outScore = (inScore * 200) / nBarNum;
			}
		}

	outScore = (outScore*929 + 1137 + 500) / 1000;

	outScore = (outScore * nRidge) / 255;

		if ( nTempNum < 5 && globalScore < 1185 && nRidge+nBlockScore < 335 && (nBlockScore < 95 || nCommonScore < 70) && (nPairNum < 10 || 			    (outScore < nPairNum*13 && nPairNum < 15)) )
		{
			if (nBlockScore > 91 && nTempNum >= 1 && nPairNum > 5 && nBarNum < 50)
				outScore = (outScore * nRidge) / 210;
			else
				outScore = (outScore * nRidge) / 255;
		}

	return outScore;
}

BOOL check_point_kind(LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int i, fid, sid, nNum = 0;

	if (pPair->nNumber == 0) return FALSE;
	for (i = 0; i < pPair->nNumber; i++) {
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		if (pFile->item[fid].kind == pSearch->item[sid].kind) nNum++;
	}
	if (nNum == pPair->nNumber) return TRUE;
	return FALSE;
}

int dec_func_01(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair,int nTempNum)
{
	int nRetScore = score, nFileScore = 0, nSearchScore = 0;
	int maxlen, diff1, len, diflen, diff, difdir, dx, dy;
	BOOL nScoreFlag = TRUE;
	int nRot = pPair->nRot;
	int nFCoreNum, nSCoreNum, nFDeltaNum = 0, nSDeltaNum = 0;
	COREITEMEX FileCore[2], SearchCore[2];
	COREITEMEX FileDelta[2], SearchDelta[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,FileDelta,&nFDeltaNum);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,SearchDelta,&nSDeltaNum);

	nFileScore = pFile->Mp.quality;
	nSearchScore = pSearch->Mp.quality;
	if (nFileScore < 30 || nSearchScore < 30) nScoreFlag = FALSE;

	if (nFCoreNum == 1 && nSCoreNum == 1) {
		dx = FileCore[0].x - SearchCore[0].x;
		dy = FileCore[0].y - SearchCore[0].y;
		len = op_func_02(dx*dx + dy*dy);
		diff = abs(FileCore[0].dir - SearchCore[0].dir);
		diff = IANGLE_120(diff);
		if (len < LENGTH  && diff >= ANGLE*4) {
			if (nScoreFlag) nRetScore = nRetScore/2;
			else nRetScore = nRetScore*2/3;
		}
		else if (len >= LENGTH*3 && diff < ANGLE*2)	{
			if (nScoreFlag) {
				nRetScore = nRetScore/2;
			}
			else {
				if (nTempNum >= 2) nRetScore = nRetScore*80/100;
				else nRetScore = nRetScore*68/100;
			}
		}
		else if (2*len >= LENGTH*5 && diff < ANGLE && nScoreFlag) nRetScore = nRetScore*2/3;
		else if (4*len >= LENGTH*7 && diff < ANGLE && nScoreFlag) nRetScore = nRetScore*85/100; 
		else if (len > LENGTH && diff > 3*ANGLE) nRetScore = nRetScore*80/100; 
		if (nRot >= 90 && nRot < 150) {
			if (len >= LENGTH*3 && diff >= ANGLE*3 && nScoreFlag)
				nRetScore = nRetScore/2;
		}
		diff = abs(120-diff);
		if (len >= LENGTH*7 && diff >= ANGLE*3)	{
			if (nScoreFlag) nRetScore = nRetScore/2;
			else nRetScore = nRetScore*2/3;
		}
		else {
			if (len >= 120 && nScoreFlag) {
				if ((pFile->nType == 8 && (pSearch->nType == 4 || pSearch->nType == 5)) || 
					 (pSearch->nType == 8 && (pFile->nType == 4 || pFile->nType == 5))) nRetScore = nRetScore*3/5;
				else nRetScore = nRetScore/3;
			}
			else if (len > 100) nRetScore = nRetScore*2/3;
		}
	}
	if (nFCoreNum == 2 && nSCoreNum == 2) {
		dx = FileCore[0].x - FileCore[1].x;
		dy = FileCore[0].y - FileCore[1].y;
		len = op_func_02(dx*dx + dy*dy);
		dx = SearchCore[0].x - SearchCore[1].x;
		dy = SearchCore[0].y - SearchCore[1].y;
		maxlen = op_func_02(dx*dx + dy*dy);
		len = abs(len-maxlen);
		diff = op_func_01(FileCore[0].x,FileCore[0].y,
						  FileCore[1].x,FileCore[1].y);
		diff = ANGLE_120(diff);
		diff1 = op_func_01(SearchCore[0].x,SearchCore[0].y,
						   SearchCore[1].x,SearchCore[1].y);
		diff1 = ANGLE_120(diff1);
		diff = abs(diff-diff1);
		diff = IANGLE_60(diff);
		diff = op_func_02(len*len + diff*diff);
		if (diff >= 50) {
			if (nFileScore < 30 && nSearchScore < 30)
				nRetScore = nRetScore*2/3;
			else
				nRetScore = nRetScore/2;
		}
	}
	if (nFCoreNum != nSCoreNum && nFCoreNum > 0 && nSCoreNum > 0) {
		if (nFCoreNum == 1)	{
			dx = FileCore[0].x - SearchCore[0].x;
			dy = FileCore[0].y - SearchCore[0].y;
			len = dx*dx + dy*dy;
			dx = FileCore[0].x - SearchCore[1].x;
			dy = FileCore[0].y - SearchCore[1].y;
			maxlen = dx*dx + dy*dy;
			if (maxlen > len) maxlen = len;
		}
		else {
			dx = FileCore[0].x - SearchCore[0].x;
			dy = FileCore[0].y - SearchCore[0].y;
			len = dx*dx + dy*dy;
			dx = FileCore[1].x - SearchCore[0].x;
			dy = FileCore[1].y - SearchCore[0].y;
			maxlen = dx*dx + dy*dy;
			if (maxlen > len) maxlen = len;
		}
		if (maxlen >= LENGTH*LENGTH*100) nRetScore = 0;
		if (maxlen >= LENGTH*LENGTH*64) nRetScore = nRetScore/2;
	}

	if (nFDeltaNum == 0 || nSDeltaNum == 0) return nRetScore;
	if (nFDeltaNum == 2 && nSDeltaNum == 2)	{
		dx = FileDelta[0].x-FileDelta[1].x;
		dy = FileDelta[0].y-FileDelta[1].y;
		len = op_func_02(dx*dx+dy*dy);
		dx = SearchDelta[0].x-SearchDelta[1].x;
		dy = SearchDelta[0].y-SearchDelta[1].y;
		maxlen = op_func_02(dx*dx+dy*dy);
		diflen = abs(len-maxlen);
		diff = op_func_01(FileDelta[0].x,FileDelta[0].y,
						  FileDelta[1].x,FileDelta[1].y);
		diff = ANGLE_120(diff);
		diff1 = op_func_01(SearchDelta[0].x,SearchDelta[0].y,
						   SearchDelta[1].x,SearchDelta[1].y);
		diff1 = ANGLE_120(diff1);
		difdir = abs(diff-diff1);
		difdir = IANGLE_60(difdir);
		diff = op_func_02(diflen*diflen + difdir*difdir);
		if (diff >= 200 || (nScoreFlag && diff >= 100) ) nRetScore = 0;
		if (diff >= 60 || difdir > 4*ANGLE)
			nRetScore = nRetScore - nRetScore*diff/200;
		if (diflen < LENGTH && difdir < ANGLE) nRetScore = nRetScore*6/5;
	}

	if (nFCoreNum == 1 && nFDeltaNum == 1 && nSCoreNum == 1 && nSDeltaNum == 1)	{
		if (pFile->nType == pSearch->nType && ( pFile->nType == 4 || pFile->nType == 5 || pFile->nType == 7)) {
			dx = FileCore[0].x-FileDelta[0].x;
			dy = FileCore[0].y-FileDelta[0].y;
			len = op_func_02(dx*dx+dy*dy);
			dx = SearchCore[0].x-SearchDelta[0].x;
			dy = SearchCore[0].y-SearchDelta[0].y;
			maxlen = op_func_02(dx*dx+dy*dy);
			diflen = abs(len-maxlen);
			diff = op_func_01(FileCore[0].x,FileCore[0].y,
							  FileDelta[0].x,FileDelta[0].y);
			diff = ANGLE_120(diff);
			diff1 = op_func_01(SearchCore[0].x,SearchCore[0].y,
							   SearchDelta[0].x,SearchDelta[0].y);
			diff1 = ANGLE_120(diff1);
			difdir = abs(diff-diff1);
			difdir = IANGLE_60(difdir);
			if (nScoreFlag && (diflen >= 10*LENGTH)) nRetScore = 0;
			if (diflen > 4*LENGTH || difdir > 3*ANGLE) {
				if (nScoreFlag) nRetScore = nRetScore/3;
				else nRetScore = nRetScore*2/3;
			}
			if (diflen <= 10 && difdir <= 7 && score > 50) nRetScore = nRetScore*6/5;
		}
	}
	return nRetScore;
}

int dec_func_02(int score,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int nRetScore = score, nDiffNum = 0;
	int i, j, k, fid, sid, dx, dy, nFNum, nSNum;
	int nLenTh = 35*35;

	for (i = 0; i<pPair->nNumber; i++) {
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		if (pFile->item[fid].score >= 45 && pSearch->item[sid].score >= 30)	{
			nFNum = 0;
			for (j = 0; j < pFile->nNumber; j++) {
				if (pFile->item[j].score < 45) continue;
				for (k = 0; k < pPair->nNumber; k++) {
					if (j == pPair->nFileID[k]) break;
				}
				if (k < pPair->nNumber) continue;
				dx = pFile->item[fid].x - pFile->item[j].x;
				dy = pFile->item[fid].y - pFile->item[j].y;
				if (dx*dx+dy*dy < nLenTh) nFNum++;
			}
			nSNum = 0;
			for (j = 0; j < pSearch->nNumber; j++) {
				for (k = 0; k < pPair->nNumber; k++) {
					if (j == pPair->nSearchID[k]) break;
				}
				if (k < pPair->nNumber) continue;
				dx = pSearch->item[sid].x - pSearch->item[j].x;
				dy = pSearch->item[sid].y - pSearch->item[j].y;
				if (dx*dx+dy*dy < nLenTh) nSNum++;
			}
			if (nSNum == 0 && nFNum >= 3) nDiffNum++;
		}
		else if (pSearch->item[sid].score >= 45 && pFile->item[fid].score >= 30) {
			nFNum = 0;
			for (j = 0; j < pFile->nNumber; j++) {
				for (k = 0; k < pPair->nNumber; k++) {
					if (j == pPair->nFileID[k]) break;
				}
				if (k < pPair->nNumber) continue;
				dx = pFile->item[fid].x - pFile->item[j].x;
				dy = pFile->item[fid].y - pFile->item[j].y;
				if (dx*dx+dy*dy < nLenTh) nFNum++;
			}
			nSNum = 0;
			for (j = 0; j < pSearch->nNumber; j++) {
				if (pSearch->item[j].score < 45) continue;
				for (k = 0; k < pPair->nNumber; k++) {
					if (j == pPair->nSearchID[k]) break;
				}
				if (k < pPair->nNumber) continue;
				dx = pSearch->item[sid].x - pSearch->item[j].x;
				dy = pSearch->item[sid].y - pSearch->item[j].y;
				if (dx*dx+dy*dy < nLenTh) nSNum++;
			}
			if (nFNum == 0 && nSNum >= 3) nDiffNum++;
		}
	}
	if (nDiffNum >= 2) nRetScore = nRetScore*3/4;
	return nRetScore;
}

int dec_func_03(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
					  int nRidge,int nBlockScore,int nTempNum)
{
	int retScore = score, th = (LENGTH+5)*(LENGTH+5);
	int i, j, len, fx, fy, sx, sy, nFlag = 0, nRTh = 220;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);
	if ( nFCoreNum==0 || nSCoreNum==0 ) return retScore;

	for (i = 0; i < nFCoreNum; i++)	{
		fx = FileCore[i].x; fy = FileCore[i].y;
		for (j = 0; j < nSCoreNum; j++) {
			sx = SearchCore[j].x; sy = SearchCore[j].y;
			len = (fx-sx)*(fx-sx) + (fy-sy)*(fy-sy);
			if (len < th){ nFlag = 1; i = 3; break; }
		}
	}
	if (nFlag == 0) return retScore;

		if (nRidge <= nRTh && nBlockScore <= 90) retScore = retScore/2;
		else if (nRidge <= nRTh+10 && nBlockScore < 80) retScore = 2*retScore/3;
		else if (nRidge <= 222 && nBlockScore <= 96) retScore = 2*retScore/3;
		else if (nRidge <= 227 && nBlockScore < 90) retScore = 8*retScore/10;
		else if (nRidge <= 235 && nBlockScore < 95 && nTempNum <= 1) retScore = 8*retScore/10;

	return retScore;
}

BOOL check_exist(int x,int y,int dir,int nID,int nLenTh,int nAngTh,
			 LPMPVECTEX pVect,PAIRVECT *pPair,
			 BOOL nPairFlag,BOOL nSimple,BOOL nForS)
{
	int i, j, id, len, diff, th = nLenTh*nLenTh;

	for (i = 0; i < pVect->nNumber; i++) {
		if (i == nID) continue;
		if (nPairFlag) {
			for (j = 0; j < pPair->nNumber; j++) {
				id = (nForS==FALSE)?pPair->nFileID[j]:pPair->nSearchID[j];
				if (i == id) break;
			}
			if (j < pPair->nNumber) continue;
		}
		len = (x-pVect->item[i].x)*(x-pVect->item[i].x) + (y-pVect->item[i].y)*(y-pVect->item[i].y);
		diff = abs(dir - pVect->item[i].dir);
		diff = IANGLE_120(diff);
		if (len < th) {
			if (!nSimple) {
				if (diff < nAngTh) return TRUE;
				continue;
			}
			return TRUE;
		}
	}
	return FALSE;
}

int dec_func_04(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int retScore = score;
	MPVECTEX tmpF, tmpS;
	POLYGON polyF, polyS;
	short *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber, nFScore = 0, nSScore = 0;
	int i, j, x, y, dir, nFNum = 0, nSNum = 0;

	if (nPairNum < 3) return retScore;
	tmpF.nNumber = tmpS.nNumber = nPairNum;
	for (i = 0; i < nPairNum; i++) {
		tmpF.item[i] = pFile->Mp.item[nFileID[i]];
		nFScore += pFile->Mp.item[nFileID[i]].score;
		tmpS.item[i] = pSearch->Mp.item[nSearchID[i]];
		nSScore += pSearch->Mp.item[nSearchID[i]].score;
	}
	nFScore /= nPairNum; nSScore /= nPairNum;
	if (nFScore > nSScore) nFScore = nSScore;
	if (nFScore < 50) return retScore;
	if (FALSE == get_polygon_points(&tmpF,&polyF)) return retScore;
	if (FALSE == get_polygon_points(&tmpS,&polyS)) return retScore;

	for (i = 0; i < pFile->Mp.nNumber; i++)	{
		if (pFile->Mp.item[i].score < 40) continue;
		for (j = 0; j < nPairNum; j++){ if (i==nFileID[j]) break; }
		if (j < nPairNum) continue;
		x = pFile->Mp.item[i].x; y = pFile->Mp.item[i].y;
		dir = pFile->Mp.item[i].dir;
		if (FALSE == check_in_polygon(x,y,&polyF,0)) continue;
		if (FALSE == check_exist(x,y,dir,-1,20,15,&pSearch->Mp,pPair,TRUE,FALSE,TRUE)) nFNum++;
	}
	for (i = 0; i < pSearch->Mp.nNumber; i++) {
		if (pSearch->Mp.item[i].score < 40) continue;
		for (j = 0; j < nPairNum; j++) {
			if (i==nSearchID[j]) break; 
		}
		if (j < nPairNum) continue;
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		dir = pSearch->Mp.item[i].dir;
		if (FALSE == check_in_polygon(x,y,&polyS,0)) continue;
		if (FALSE == check_exist(x,y,dir,-1,20,15,&pFile->Mp,pPair,TRUE,FALSE,FALSE)) nSNum++;
	}
	nFNum = nFNum + nSNum;
	if (nFNum >= 5) retScore = retScore/2;
	else if (nFNum >= 3) retScore = retScore - nFNum*5;
	return retScore;
}

BOOL check_overlap(LPCOREVECTEX pFile,LPCOREVECTEX pSearch)
{
	int i, j, len, diff, fx, fy, fdir, sx, sy, sdir;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(pFile,FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(pSearch,SearchCore,NULL,NULL);

	if (nFCoreNum==0 || nSCoreNum==0) return FALSE;
	for (i = 0; i < nFCoreNum; i++)	{
		fx = FileCore[i].x; fy = FileCore[i].y;
		fdir = FileCore[i].dir;
		for (j = 0; j < nSCoreNum; j++)	{
			sx = SearchCore[j].x; sy = SearchCore[j].y;
			sdir = SearchCore[j].dir;
			len = (fx-sx)*(fx-sx) + (fy-sy)*(fy-sy);
			diff = abs(fdir-sdir);
			diff = IANGLE_120(diff);
			if (len < 16*16 && diff < 7) return TRUE;
		}
	}
	return FALSE;
}

int dec_func_05(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int retScore = score;
	short *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber;
	BOOL nCoreFlag = TRUE, nFalse = FALSE;
	int sid, fid, fx, fy, sx, sy, flen, slen, sid1, fid1, sid2, fid2;
	int fdir1, fdir2, fdir, sdir1, sdir2, sdir, diff;
	int i, j, k;

	if (nPairNum < 3 || nPairNum > 10) return retScore;
	if (pFile->Core.nNumber==0 || pSearch->Core.nNumber == 0)
		nCoreFlag = FALSE;
	else {
		if (FALSE == check_overlap(&pFile->Core,&pSearch->Core))	nCoreFlag = FALSE;
	}
	for (i = 0; i < nPairNum; i++) {
		sid = nSearchID[i]; fid = nFileID[i];
		fx = pFile->Mp.item[fid].x; fy = pFile->Mp.item[fid].y;
		sx = pSearch->Mp.item[sid].x; sy = pSearch->Mp.item[sid].y;
		if (pFile->Mp.item[fid].score < 30) continue;
		if (pSearch->Mp.item[sid].score < 30) continue;
		for (j = 0; j < nPairNum; j++) {
			if (i == j) continue;
			sid1 = nSearchID[j]; fid1 = nFileID[j];
			if (pFile->Mp.item[fid1].score < 20) continue;
			if (pSearch->Mp.item[sid1].score < 20) continue;
			flen = (fx-pFile->Mp.item[fid1].x)*(fx-pFile->Mp.item[fid1].x) 
				 + (fy-pFile->Mp.item[fid1].y)*(fy-pFile->Mp.item[fid1].y);
			slen = (sx-pSearch->Mp.item[sid1].x)*(sx-pSearch->Mp.item[sid1].x) 
				 + (sy-pSearch->Mp.item[sid1].y)*(sy-pSearch->Mp.item[sid1].y);
			if (flen >= 9*LENGTH*LENGTH || slen >= 9*LENGTH*LENGTH) continue;
			fdir1 = op_func_01(fx,fy,pFile->Mp.item[fid1].x,pFile->Mp.item[fid1].y);
			sdir1 = op_func_01(sx,sy,pSearch->Mp.item[sid1].x,pSearch->Mp.item[sid1].y);
			for (k = 0; k < nPairNum; k++) {
				if (k == j || k == i) continue;
				sid2 = nSearchID[k]; fid2 = nFileID[k];
				if (pFile->Mp.item[fid2].score < 20) continue;
				if (pSearch->Mp.item[sid2].score < 20) continue;
				flen = (fx-pFile->Mp.item[fid2].x)*(fx-pFile->Mp.item[fid2].x) 
					 + (fy-pFile->Mp.item[fid2].y)*(fy-pFile->Mp.item[fid2].y);
				slen = (sx-pSearch->Mp.item[sid2].x)*(sx-pSearch->Mp.item[sid2].x) 
					 + (sy-pSearch->Mp.item[sid2].y)*(sy-pSearch->Mp.item[sid2].y);
				if (flen >= 9*LENGTH*LENGTH || slen >= 9*LENGTH*LENGTH) continue;
				fdir2 = op_func_01(fx,fy,pFile->Mp.item[fid2].x,pFile->Mp.item[fid2].y);
				sdir2 = op_func_01(sx,sy,pSearch->Mp.item[sid2].x,pSearch->Mp.item[sid2].y);
				fdir = abs(fdir1 - fdir2);
				fdir = IANGLE_120(fdir);
				sdir = abs(sdir1 - sdir2);
				sdir = IANGLE_120(sdir);
				diff = abs(fdir - sdir);
				diff = IANGLE_120(diff);
				if (diff >= ANGLE*3) {
					j = 100; i = 100; nFalse = TRUE; break;
				}
			}
		}
	}
	if (nFalse == TRUE)	{
		if (!nCoreFlag) retScore = retScore*8/10;
		else retScore = retScore*7/10;
	}
	return retScore;
}

int dec_func_06(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
						PAIRVECT *pPair,int nMatchScoreTh,int nRidge,
						int nBlockScore,int score_old)
{
	int retScore = score;
	MPVECTEX tmpF, tmpS;
	POLYGON polyF, polyS;
	short *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber, nFScore, nSScore;
	int i, j, x, y, dir, nFNum = 0, nSNum = 0;
	BOOL nCoreFlag, nScoreFlag = TRUE;

	if (pFile->Mp.nNumber==0 || pSearch->Mp.nNumber==0) return 0;
	if (nPairNum == 0) return 0;
	if (nPairNum > 11 || score >= nMatchScoreTh*3 || score_old >= nMatchScoreTh*4) return retScore;
	
		if (nRidge >= 245 && nBlockScore >= 95) return retScore;
		if (nPairNum > 9 && nRidge >= 248 && nBlockScore >= 92) return retScore;

	nCoreFlag = check_core(&pFile->Core,&pSearch->Core,16,7);
	nFScore = pFile->Mp.quality; nSScore = pSearch->Mp.quality;
	if (nFScore < 35 || nSScore < 35) nScoreFlag = FALSE;

	tmpF.nNumber = tmpS.nNumber = nPairNum;
	for (i = 0; i < nPairNum; i++) {
		tmpF.item[i] = pFile->Mp.item[nFileID[i]];
		tmpS.item[i] = pSearch->Mp.item[nSearchID[i]];
	}
	if (get_polygon_points(&tmpF,&polyF) == FALSE) return retScore;
	if (get_polygon_points(&tmpS,&polyS) == FALSE) return retScore;

	for (i = 0; i < pFile->Mp.nNumber; i++)	{
		if (pFile->Mp.item[i].score < 30 ) continue;
		for (j = 0; j < nPairNum; j++) {
			if (i==nFileID[j]) break; 
		}
		if (j < nPairNum) continue;
		x = pFile->Mp.item[i].x; y = pFile->Mp.item[i].y;
		dir = pFile->Mp.item[i].dir;
		if (FALSE == check_in_polygon(x,y,&polyF,0)) continue;
		if (FALSE == check_exist(x,y,dir,-1,20,20,&pSearch->Mp,NULL,FALSE,FALSE,TRUE)) nFNum++;
	}
	for (i = 0; i < pSearch->Mp.nNumber; i++) {
		if (pSearch->Mp.item[i].score < 30) continue;
		for (j = 0; j < nPairNum; j++) { if (i==nSearchID[j]) break; }
		if (j < nPairNum) continue;
		x = pSearch->Mp.item[i].x; y = pSearch->Mp.item[i].y;
		dir = pSearch->Mp.item[i].dir;
		if (FALSE == check_in_polygon(x,y,&polyS,0)) continue;
		if (FALSE == check_exist(x,y,dir,-1,20,20,&pFile->Mp,NULL,FALSE,FALSE,FALSE)) nSNum++;
	}
	nFNum = nFNum + nSNum;
	if (nFNum >= 5)	{
		if (nCoreFlag && nScoreFlag)
			retScore = retScore/2;
		else
			retScore = retScore*7/10;
	}
	else if (nFNum >= 3) {
		if (nScoreFlag) retScore -= nFNum*5;
		else retScore -= nFNum*3;
	}
	return retScore;
}

int dec_func_07(int score,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int nRetScore = score;
	int i, fid, sid, nNum = 0, maxscore, minscore;

	if (pPair->nNumber == 0) return 0;
	for (i = 0; i < pPair->nNumber; i++) {
		fid = pPair->nFileID[i]; sid = pPair->nSearchID[i];
		maxscore = pFile->item[fid].score;
		if (maxscore < pSearch->item[sid].score) maxscore = pSearch->item[sid].score;
		minscore = pFile->item[fid].score;
		if (minscore > pSearch->item[sid].score) minscore = pSearch->item[sid].score;
		if (maxscore < 40 && minscore < 25)  continue;
		if (pFile->item[fid].kind != pSearch->item[sid].kind) nNum++;
	}
	if (pPair->nNumber >= 5) {
		if (nNum*10 >= pPair->nNumber*7) nRetScore = nRetScore/2;
		else if (nNum*3 >= pPair->nNumber*2) nRetScore = nRetScore*6/8;
		else if (nNum*2 >= pPair->nNumber) nRetScore = nRetScore*65/80;
	}
	return nRetScore;
}

int dec_func_08(int score,int nPairNeigh, LPMPVECTEX pFile,LPMPVECTEX pSearch,
						PAIRVECT *pPair,int nBlockScore,int nRidge,int globalScore,BOOL nCoreFlag,int nCommonScore,int localScore)
{
	int nRetScore = score, nLenTh = (LENGTH+5)*(LENGTH+5);
	int i, j, k, dx, dy, diff, num0, num1, dScore = 3;
	int nNum = pPair->nNumber, num, nScoreFTh = 30, nScoreSTh = 30;
	int nMinNum = pFile->nNumber;
	BOOL qflag = TRUE;

	if (abs(pFile->quality-pSearch->quality) >= 12) qflag = FALSE;
	
	if (nMinNum > pSearch->nNumber) nMinNum = pSearch->nNumber;

	if (nNum > 18 && nRidge >= 235 && nPairNeigh > 0) return nRetScore;
	if (nNum*100 >= nMinNum*44 && nBlockScore >= 98 && nNum > 8 && globalScore > 1400) return nRetScore;

	if (nNum >= 10 && nNum*100 > nMinNum*30 && nCommonScore > 78) {
		num = 0;
		for (i = 0; i < pFile->nNumber; i++) {
			if (pFile->item[i].score < 35) continue;
			for (j = 0; j < nNum; j++) {
				if (i == pPair->nFileID[j]) break;
			}
			if (j < nNum) continue;
			if (TRUE == check_exist(pFile->item[i].x,pFile->item[i].y,pFile->item[i].dir,
						-1,20,ANGLE,pSearch,pPair,FALSE,FALSE,TRUE)) num++;
		}
		num = nNum + num;
		if (100*num > nMinNum*60 && nNum >= 13 && nBlockScore > 91 && nRidge > 237 && globalScore > 850 && nPairNeigh > 0 && qflag) return nRetScore;
		if (100*num > nMinNum*80 && nNum >= 11 && nBlockScore > 90 && nRidge > 215 && nPairNeigh > 0) return nRetScore;
		if (nCoreFlag && globalScore > 900 && nNum >= 13 && nBlockScore > 92 && nRidge > 250 && nPairNeigh > 0) return nRetScore;
	}
	num = 0;
	if (nNum*100 >= nMinNum*44)	{
		if (nBlockScore >= 96 && nNum > 6 && nRidge > 245) dScore = 2;
		if (nBlockScore > 92 && nNum >= 10 && nRidge >= 250) dScore = 2;
	}
	if (nNum <= 6) {
		if (nNum <= 5 && nBlockScore <= 92) dScore = 4;
		if (nBlockScore < 90) dScore = 4;
		if (nBlockScore < 90 && nNum*100 < nMinNum*38) dScore = 5;
	}
	if (globalScore > 1350 && nNum >= 12 && nNum*100 >= nMinNum*33) dScore = 2;
	if (globalScore < 350 && nNum <= 7) dScore = 4; 
	if (nPairNeigh == 0) dScore = 5;
	if (globalScore > 1150 && nNum >= 15) dScore = 4;
	if (nPairNeigh == 1 || nCoreFlag) dScore = 4;
	if ((nPairNeigh == 3 && nBlockScore > 93) || (nPairNeigh == 5 && nNum*100 > nMinNum*34)) dScore = 2;
	if (nPairNeigh == 4 && nBlockScore >= 92 && nNum*100 > nMinNum*70) dScore = 3;
	if (nPairNeigh == 2 && nNum*100 > nMinNum*50 && globalScore > 750 && (nNum < 8 || nCoreFlag)) dScore = 2;
	if (nPairNeigh >= 2 && globalScore > 900 && nCoreFlag && nBlockScore >= 88 && nRidge > 250) dScore--;
	if (nPairNeigh == 2 && globalScore < 400 && nBlockScore < 88) dScore++;

	if (nPairNeigh == 0 && globalScore > 1100 && nNum >= 9 && nBlockScore > 90) dScore = 4;
	if (nPairNeigh == 0 && nNum >= 18 && nNum*100 > nMinNum*45) dScore = 4;
	if (nPairNeigh <= 1 && ((nNum > 12 && globalScore < 900) || (globalScore < 650 && nRidge < 235 && nBlockScore < 95))) dScore++;
	else if (nPairNeigh == 0 && globalScore < 550 && nNum >= 7) dScore++;
	if (nPairNeigh <= 1 && nNum < 10 && globalScore < 700 && (nBlockScore <= 85 || (nBlockScore <= 92 && nNum*100 < nMinNum*48))) dScore++;
	if (nCoreFlag && globalScore > 1160 && nBlockScore > 93) dScore = 3;  
	
	if (nCommonScore < 70 && !nCoreFlag) dScore++;
	if (nCommonScore <= 75 && nBlockScore < 90 && nNum*100 < nMinNum*45 && !nCoreFlag) dScore++;

	num1 = 0;

	if (pFile->quality < 20) nScoreFTh = 20;
	else if (pFile->quality <= 26) nScoreFTh = 26;
	if (pSearch->quality < 20) nScoreSTh = 20;
	else if (pSearch->quality <= 26) nScoreSTh = 26;

	if (!qflag)	{
		if (nScoreFTh < 30) nScoreFTh = pFile->quality;
		if (nScoreSTh < 30) nScoreSTh = pSearch->quality;
		if (dScore >= 5) dScore++;
		else dScore = 5;
	}

	if (nPairNeigh == 0 && nNum >= 8 && nNum*100 < nMinNum*45) {
		if (pFile->quality >= 27) nScoreFTh = 28;
		if (pSearch->quality >= 27) nScoreSTh = 28;
	}

	for (i = 0; i < pFile->nNumber; i++) {
		if (pFile->item[i].score < nScoreFTh) continue;
		num++;
		for (j = 0; j < pPair->nNumber; j++) {
			if (i == pPair->nFileID[j]) break;
		}
		if (j < pPair->nNumber) continue;
		num0 = 0;
		for (j = 0; j < pSearch->nNumber; j++) {
			if (pSearch->item[j].score < nScoreSTh) continue;
			for (k = 0; k < pPair->nNumber; k++) {
				if (j == pPair->nSearchID[k]) break;
			}
			if (k < pPair->nNumber) continue;
			dx = pFile->item[i].x - pSearch->item[j].x;
			dy = pFile->item[i].y - pSearch->item[j].y;
			dx = dx*dx + dy*dy;
			if (dx >= nLenTh) continue;
			num0++;
	 		diff = abs(pFile->item[i].dir-pSearch->item[j].dir);
			diff = IANGLE_120(diff);
			if (diff >= 8*ANGLE){ nRetScore -= 4; continue; }
			if (diff >= 4*ANGLE){ nRetScore -= 2; continue; }
			if (diff >= 2*ANGLE){ nRetScore -= 1; continue; }
		}
		if (num0 == 0) num1++;
	}
	if (num1 > 0) {
		if ((nNum <= 5 && nBlockScore <= 70) || (nNum <= 4 && nBlockScore < 85)) nRetScore /= (num1*dScore);
		else {
			if (nNum*100 > nMinNum*30 && nBlockScore > 96 && num1 < nNum && globalScore > 400) nRetScore -= (dScore-1)*num1;	
			else nRetScore -= dScore*num1;
		}
	}
	if (nRetScore < 0) nRetScore = 0;
	return nRetScore;
}

int dec_func_09(int score,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair,int dScore)
{
	POLYGON polyF, polyS;
	int retScore = score, nTh = 40*40, nPairNum = pPair->nNumber;
	short *nFileID = &(pPair->nFileID[0]);
	short *nSearchID = &(pPair->nSearchID[0]);
	int i, j, k, h, fid1, fid2, sid1, sid2, flen, slen, len;
	int fx1, fy1, fx2, fy2, sx1, sy1, sx2, sy2, fx, fy, sx, sy;
	int x, y, dir, fnum = 0, snum = 0, th = 9;
	int flist[MAX_MINUTIA_NUMBER], slist[MAX_MINUTIA_NUMBER];


	if (FALSE == get_polygon_points(pSearch,&polyS)) return retScore;
	if (FALSE == get_polygon_points(pFile,&polyF)) return retScore;

	for (i = 0; i < nPairNum-1; i++) {
		fid1 = nFileID[i]; sid1 = nSearchID[i];
		if (pFile->item[fid1].score < 45) continue;
		if (pSearch->item[sid1].score < 45) continue;
		fx1 = pFile->item[fid1].x; fy1 = pFile->item[fid1].y;
		sx1 = pSearch->item[sid1].x; sy1 = pSearch->item[sid1].y;
		for (j = i+1; j < nPairNum; j++) {
			fid2 = nFileID[j]; sid2 = nSearchID[j];
			if (pFile->item[fid2].score < 45) continue;
			if (pSearch->item[sid2].score < 45) continue;
			fx2 = pFile->item[fid2].x; fy2 = pFile->item[fid2].y;
			sx2 = pSearch->item[sid2].x; sy2 = pSearch->item[sid2].y;
			flen = (fx1-fx2)*(fx1-fx2) + (fy1-fy2)*(fy1-fy2);
			slen = (sx1-sx2)*(sx1-sx2) + (sy1-sy2)*(sy1-sy2);
			fx = (fx1 + fx2)/2; fy = (fy1 + fy2)/2;
			sx = (sx1 + sx2)/2; sy = (sy1 + sy2)/2;
			if (flen < slen) flen = slen;
			if (flen >= nTh) continue;

			for (k = 0; k < pFile->nNumber; k++) {
				if (k == fid1 || k == fid2) continue;
				if (pFile->item[k].score < 45) continue;
				for (h = 0; h < fnum; h++) {
					if (flist[h] == k) break;
				}
				if (h < fnum) continue;
				for (h = 0; h < nPairNum; h++) {
					if (k == nFileID[h]) break;
				}
				if (h < nPairNum) continue;
				len = (fx-pFile->item[k].x)*(fx-pFile->item[k].x)
					+ (fy-pFile->item[k].y)*(fy-pFile->item[k].y);
				if (len > nTh) continue;
				x = sx - (fx - pFile->item[k].x);
				y = sy - (fy - pFile->item[k].y);
				dir = pFile->item[k].dir;
				if (FALSE == check_in_polygon(x,y,&polyS,-1)) continue;
				if (FALSE == check_exist(x,y,dir,-1,LENGTH+10,25,pSearch,pPair,TRUE,FALSE,TRUE)) {
					retScore = retScore*th/10;
					flist[fnum++] = k;
				}
			}
			for (k = 0; k < pSearch->nNumber; k++) {
				if (k == sid1 || k == sid2) continue;
				if (pSearch->item[k].score < 45) continue;
				for (h = 0; h < snum; h++) {
					if (slist[h] == k) break;
				}
				if (h < snum) continue;
				for (h = 0; h < nPairNum; h++) {
					if (k == nSearchID[h]) break;
				}
				if (h < nPairNum) continue;
				len = (sx-pSearch->item[k].x)*(sx-pSearch->item[k].x)
					+ (sy-pSearch->item[k].y)*(sy-pSearch->item[k].y);
				if (len > nTh) continue;
				x = fx - (sx - pSearch->item[k].x);
				y = fy - (sy - pSearch->item[k].y);
				dir = pSearch->item[k].dir;
				if (FALSE == check_in_polygon(x,y,&polyF,0)) continue;
				if (FALSE == check_exist(x,y,dir,-1,LENGTH+10,25,pFile,pPair,TRUE,FALSE,FALSE))	{
					retScore = retScore*th/10;
					slist[snum++] = k;
				}
			}
		}
	}

	return retScore;
}

int dec_func_10(int score,int nTempNum,LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT* pPair,
					   int nBlockScore,int nRidge,int globalScore,BOOL nCoreFlag,int nCommonScore,int dscore)
{
	int retScore = score;
	int nNum = pPair->nNumber, nMinNum = pFile->nNumber;
	
	if (nMinNum > pSearch->nNumber) nMinNum = pSearch->nNumber;

	if (nNum == 0 ) return (0);
	if (nRidge >= 235 && ((nNum >= 15 && nTempNum >= 2 && nNum*100 >= nMinNum*50) || nNum > 18)) return (retScore);
	if (nRidge >= 235 && nNum >= 14 && nTempNum >= 3 && globalScore > 850 && nNum*100 >= nMinNum*50) return (retScore);

	if (nTempNum <= 6) {
		if (nTempNum == 5 && nNum >= 8 && globalScore > 1035) return retScore;
		if (nTempNum == 4 && nRidge > 240 && nCommonScore > 80 && nBlockScore > 81 && (nNum > 10 || (nNum*100 > nMinNum*35 && nNum > 4))) 
			return retScore;
		if (nTempNum == 4 && nRidge >= 252 && nCommonScore > 82 && nNum > 17 && nNum*100 > nMinNum*43) return retScore;
		if (nTempNum == 2 && nNum >= 10 && nBlockScore > 91 && globalScore >= 1100 && nCoreFlag) return retScore;
		if (nTempNum == 3 && nCommonScore > 75 && nNum >= 7 && nBlockScore >= 85 && globalScore > 1000 && nRidge > 250 && nNum*100 > nMinNum*40)
			return retScore;
		if (nTempNum == 3 && nCommonScore > 90 && nNum >= 15 && nBlockScore >= 83 && globalScore > 1000 && nRidge > 252 && nNum*100 > nMinNum*44) 				return retScore;
		if (nTempNum == 0 && nNum > 6 && nBlockScore < 83) retScore -= 2*nNum;
		retScore = dec_func_08(retScore,nTempNum,pFile,pSearch,pPair,nBlockScore,nRidge,globalScore,nCoreFlag,nCommonScore,dscore);
	}
	if (nTempNum >= 5 && nNum*100 > nMinNum*34) return (retScore);
	if (nTempNum >= 2 && pPair->nNumber >= 9 && globalScore > 930 && nRidge > 252 && nCoreFlag) return (retScore);
	if (nTempNum >= 3 && pPair->nNumber >= 7 && globalScore > 850 && nRidge > 250) return (retScore);
	retScore = dec_func_09(retScore,pFile,pSearch,pPair,dscore);
	if (retScore < 0) retScore = 0;

	return retScore;
}

BOOL check_paired_mp(LPMPVECTEX pFile,LPMPVECTEX pSearch,PAIRVECT *pPair)
{
	int nPairNum = pPair->nNumber;
	short *nFileID = &(pPair->nFileID[0]);
	short *nSearchID = &(pPair->nSearchID[0]);
	int i, j, len;
	int fx1 = 0, fy1 = 0;
	int fx2 = 0, fy2 = 0;
	int num = 0, nMinNum = pFile->nNumber;

	if (pFile->quality < 35 || pSearch->quality < 35) return TRUE;
	if (nMinNum > pSearch->nNumber) nMinNum = pSearch->nNumber;
	if (nPairNum*100 > nMinNum*50 ) return TRUE;
	
	for (i = 0; i < nPairNum; i++) {
		fx1 += pFile->item[nFileID[i]].x; fy1 += pFile->item[nFileID[i]].y;
	}
	fx1 /= nPairNum; fy1 /= nPairNum;

	if (pFile->nNumber > pSearch->nNumber) {
		for (i = 0; i < pFile->nNumber; i++) {
			if ( pFile->item[i].score < 30 ) continue;
			for (j = 0; j < nPairNum; j++) {
				if (i == nFileID[j]) break;
			}
			if (j < nPairNum) continue;
			fx2 += pFile->item[i].x; fy2 += pFile->item[i].y;
			num++;
		}
		if (num == 0) return TRUE;
		fx2 /= num; fy2 /= num;
	}
	else {
		for (i = 0; i < pSearch->nNumber; i++) {
			if ( pSearch->item[i].score < 30 ) continue;
			for (j = 0; j < nPairNum; j++) {
				if (i == nSearchID[j]) break;
			}
			if (j < nPairNum) continue;
			fx2 += pSearch->item[i].x; fy2 += pSearch->item[i].y;
			num++;
		}
		if (num == 0) return TRUE;
		fx2 /= num; fy2 /= num;
	}
	len = (fx1-fx2)*(fx1-fx2) + (fy1-fy2)*(fy1-fy2);
	len = op_func_02(len);

	if (len > 6*LENGTH) return FALSE;

	return TRUE;
}

void get_neighbor(int cx,int cy,LPMPVECTEX pVect,short *pPairID,int nPairNum,
				 BOOL nPairFlag,int nLenTh,BOOL nScoreFlag,int nScoreTh,
				 BOOL nNumFlag,int nNumTh,LPMPVECTEX pNewVect)
{
	int i, x, y, len, th = nLenTh*nLenTh;
	int j, nNum = 0, list[MAX_MINUTIA_NUMBER], lenlist[MAX_MINUTIA_NUMBER], nMinIdx, nMinValue, tmp;
	int nn = 0;

	pNewVect->nNumber = 0;
	for (i = 0; i < pVect->nNumber; i++) {
		if (nScoreFlag)	{
			if (pVect->item[i].score < nScoreTh) continue;
		}
		x = pVect->item[i].x; y = pVect->item[i].y;
		if (x == cx && y == cy) continue;
		if (nPairFlag) {
			for (j = 0; j < nPairNum; j++) {
				if (i == pPairID[j]) break;
			}
			if (j < nPairNum) continue;
		}
		len = (x-cx)*(x-cx) + (y-cy)*(y-cy);
		if (len >= th) continue;
		list[nNum] = i; lenlist[nNum++] = len;
		pNewVect->item[nn++] = pVect->item[i];
	}
	if (nNumFlag) {
		if (nNum > nNumTh) {
			for (i = 0; i < nNum-1; i++) {
				nMinIdx = i;
				nMinValue = lenlist[i];
				for (j = i+1; j < nNum; j++) {
					if (lenlist[j] >= nMinValue) continue;
					nMinIdx = j; nMinValue = lenlist[j];
				}
				if (nMinIdx == i) continue;
				tmp = list[i]; list[i] = list[nMinIdx]; list[nMinIdx] = tmp;
				tmp = lenlist[i]; lenlist[i] = lenlist[nMinIdx]; lenlist[nMinIdx] =tmp;
			}
			nn = 0;
			for (i = 0; i < nNumTh; i++) {
				pNewVect->item[nn++] = pVect->item[list[i]];
			}
		}
	}
	pNewVect->nNumber = nn;
}

BOOL check_neighbor(int nFid,int nSid,LPMPVECTEX tmpF,LPMPVECTEX tmpS,
				   LPFPVECTEX pFile,LPFPVECTEX pSearch)
{
	POLYGON polyF, polyS;
	int i, x, y, dir, nFNum = 0, nSNum = 0;

	if (FALSE == get_polygon_points(&pFile->Mp,&polyF)) return TRUE;
	if (FALSE == get_polygon_points(&pSearch->Mp,&polyS)) return TRUE;

	for (i = 0; i < tmpF->nNumber; i++)	{
		x = tmpF->item[i].x; y = tmpF->item[i].y;
		dir = tmpF->item[i].dir;
		if (FALSE == check_in_polygon(x,y,&polyS,0)) continue;
		if (FALSE == check_exist(x,y,dir,nSid,20,20,&pSearch->Mp,NULL,FALSE,FALSE,FALSE)) nFNum++;			
	}
	for (i = 0; i < tmpS->nNumber; i++)	{
		x = tmpS->item[i].x; y = tmpS->item[i].y;
		dir = tmpS->item[i].dir;
		if (FALSE == check_in_polygon(x,y,&polyF,0)) continue;
		if (FALSE == check_exist(x,y,dir,nFid,20,20,&pFile->Mp,NULL,FALSE,FALSE,FALSE))	nSNum++;
	}
	if (nSNum > 0 && nSNum == tmpS->nNumber) return FALSE;
	if (nFNum > 0 && nFNum == tmpF->nNumber) return FALSE;
	return TRUE;
}

int dec_func_11(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair,int nRidge, int mScore, int dScore)
{
	int retScore = score;
	short *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber;
	BOOL nCoreFlag, nScoreFlag = TRUE;
	MPVECTEX tmpF, tmpS;
	int scoreF, scoreS, nTotalNum, nFalseNum;
	int i, fid, sid, fx, fy, sx, sy;

	if (nPairNum > 13) return retScore;
	if (pFile->Mp.nNumber==0 || pSearch->Mp.nNumber==0) return 0;
	nCoreFlag = check_overlap(&pFile->Core,&pSearch->Core);
	scoreF = pFile->Mp.quality; scoreS = pSearch->Mp.quality;
	if (scoreF < 35 || scoreS < 35) nScoreFlag = FALSE;

	nTotalNum = 0; nFalseNum = 0;
	for (i = 0; i < nPairNum; i++) {
		fid = nFileID[i]; sid = nSearchID[i];
		if (pFile->Mp.item[fid].score < 30) continue;
		if (pSearch->Mp.item[sid].score < 30) continue;
		fx = pFile->Mp.item[fid].x; fy = pFile->Mp.item[fid].y;
		sx = pSearch->Mp.item[sid].x; sy = pSearch->Mp.item[sid].y;
		nTotalNum++;
		get_neighbor(fx,fy,&pFile->Mp,nFileID,nPairNum,TRUE,50,TRUE,30,FALSE,0,&tmpF);
		get_neighbor(sx,sy,&pSearch->Mp,nSearchID,nPairNum,TRUE,50,TRUE,30,FALSE,0,&tmpS);
		if (FALSE == check_neighbor(fid,sid,&tmpF,&tmpS,pFile,pSearch)) nFalseNum++;
	}

	if (nTotalNum > 0) {
		if (nFalseNum >= 5) {
				if (nCoreFlag && nScoreFlag)
					retScore = retScore/2;
				else
					retScore = retScore*33/50;
		}
		else if (nFalseNum >= 4) {
			if (nScoreFlag) retScore -= nFalseNum*6;
			else retScore -= nFalseNum*5;
		}
		else if (nFalseNum >= 3) {
			if (nScoreFlag || nFalseNum*4 >= nTotalNum*3) retScore -= nFalseNum*5;
			else retScore -= 13*nFalseNum/3;
		}
		else if (nFalseNum >= 2) {
			retScore -= 5;
		}
		else if (nFalseNum >= 1) 
			retScore -= 3;
	}
	else {
		retScore = retScore*80/100;
	}

	return retScore;
}

int dec_func_12(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair)
{
	int retScore = score;
	int fx = 0, fy = 0, fdir, sx = 0, sy = 0, sdir, len, diff, nFlag = 0;
	int fx1, fy1, fdir1, sx1, sy1, sdir1;
	short *nFileID = &(pPair->nFileID[0]), *nSearchID = &(pPair->nSearchID[0]);
	int nPairNum = pPair->nNumber, num;
	POLYGON polyF, polyS;
	int i, j, k, th = (LENGTH+5)*(LENGTH+5);
	int nFCoreNum, nSCoreNum;
	MPVECTEX *pFMp, *pSMp;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&pFile->Core,FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&pSearch->Core,SearchCore,NULL,NULL);

	pFMp = &(pFile->Mp); pSMp = &(pSearch->Mp);

	if (nFCoreNum==0 || nSCoreNum==0) return retScore;
	for (i = 0; i < nFCoreNum; i++)	{
		fx = FileCore[i].x; fy = FileCore[i].y;
		fdir = FileCore[i].dir;
		for (j = 0; j < nSCoreNum; j++) {
			sx = SearchCore[j].x; sy = SearchCore[j].y;
			sdir = SearchCore[j].dir;
			diff = abs(fdir - sdir);
			diff = IANGLE_120(diff);
			if (diff > ANGLE) continue;
			len = (fx-sx)*(fx-sx) + (fy-sy)*(fy-sy);
			if (len < LENGTH*LENGTH){ nFlag = 1; i = 3; break; }
		}
	}
	if (nFlag == 0) return retScore;
	if (FALSE == get_polygon_points(pFMp,&polyF)) return retScore;
	if (FALSE == get_polygon_points(pSMp,&polyS)) return retScore;
	for (i = 0; i < pFMp->nNumber; i++)	{
		if (pFMp->item[i].score < 20) continue;
		for (j = 0; j < nPairNum; j++) {
			if (i == nFileID[j]) break; 
		}
		if (j < nPairNum) continue;
		fx1 = pFMp->item[i].x; fy1 = pFMp->item[i].y;
		fdir1 = pFMp->item[i].dir;
		len = (fx-fx1)*(fx-fx1) + (fy-fy1)*(fy-fy1);
		if (len < LENGTH*LENGTH) continue;
		if (len >= 50*50) continue;
		if (check_in_polygon(fx1,fy1,&polyS,0) == FALSE) continue;
		num = 0;
		for (j = 0; j < pSMp->nNumber; j++)	{
			for (k = 0; k < nPairNum; k++) {
				if (j == nSearchID[k]) break; 
			}
			if (k < nPairNum) continue;
			sx1 = pSMp->item[j].x; sy1 = pSMp->item[j].y;
			sdir1 = pSMp->item[j].dir;
			len = (sx-sx1)*(sx-sx1) + (sy-sy1)*(sy-sy1);
			if (len < LENGTH*LENGTH) continue;
			if (len >= 50*50) continue;
			len = (fx1-sx1)*(fx1-sx1) + (fy1-sy1)*(fy1-sy1);
			if (len >= th) continue;
			num++;
	 		diff = abs(fdir1 - sdir1);
			diff = IANGLE_120(diff);
			if (diff >= 4*ANGLE) { retScore -= 5; continue; }
			if (2*diff >= 5*ANGLE) { retScore -= 2; continue; }
		}
		if (num == 0) {
			if (pFMp->item[i].score < 40) continue;
			retScore -= 7;
		}
	}
	if (retScore < 0) retScore = 0;
	return retScore;
}

int dec_func_13(int score,LPFPVECTEX pFile,LPFPVECTEX pSearch,
					   PAIRVECT *pPair,int nMatchScoreTh,int nRidge,int nBlockScore,int dScore)
{
	int retScore = score;
	int i, dx, dy, len, nFalseNum = 0;

		if (pPair->nNumber >= 10) return retScore;
		if (nRidge >= 250 && nBlockScore >= 90) return retScore;

	for (i = 0; i < pPair->nNumber; i++) {
		if (pFile->Mp.item[pPair->nFileID[i]].score < 30) continue;
		if (pSearch->Mp.item[pPair->nSearchID[i]].score < 30) continue;
		dx = pFile->Mp.item[pPair->nFileID[i]].x-pSearch->Mp.item[pPair->nSearchID[i]].x;
		dy = pFile->Mp.item[pPair->nFileID[i]].y-pSearch->Mp.item[pPair->nSearchID[i]].y;
		len = op_func_02(dx*dx+dy*dy);
		if (len >= 8) nFalseNum++;
	}

	if (nFalseNum >= 5)
		retScore = retScore/2;
	else if (nFalseNum >= 3) {
			if (nFalseNum >= 4) retScore -= nFalseNum*5;
			else retScore -= nFalseNum*4;
	}
	else if (nFalseNum >= 2) {
		retScore -= 3;
	}
	else if (nFalseNum >= 1) retScore -= 2;

		if (nFalseNum >= 1 && nRidge <= 228 && nBlockScore <= 86) retScore = retScore*2/3;

	return retScore;
}

int point_matching(LPFPVECTEX pFile,LPFPVECTEX pSearch,PAIRVECT *pPair,BOOL nCollectFlag,BOOL nCommonFlag,int nScoreTh,int *nGlobalScore) 
{
	BARVECT FileBar, SearchBar;
	PAIRBAR PairList[1000], PairTemp[200];
	int AngleField[240], SDiffField[240], FDiffField[240];
	int XField[640], YField[640], LastList[500];
	int SArrangBarPtr[240*MAX_SAMEDIFF_NUM], FArrangBarPtr[240*MAX_SAMEDIFF_NUM];
	BAR FBar, SBar;
	int nPairNum, nLastNum, nMaxSearchBarLen;
	int nSid1, nSid2, nFid1, nFid2;
	int rdiff1, rdiff2, curv1, curv2, dslope, nCheck;
	int nFileCX = 0, nFileCY = 0, nRot, nXoffset, nYoffset, nRidge;
	int score_sum, score, score_max, score_old, nBlockScore = 0, nCommonScore = 0;
	int i, j, h, k, nUpper, nLower, id, diff1, diff2, lendiff;
	int nNum, snum = 0, fnum = 0, tmpNum, nSid, nFid, id1 = 0, id2;
	int newCX, newCY, nMatchScoreTh = nScoreTh;
	BOOL nRepeat, hflag, sflag, fflag, flag, nflag = FALSE;
	BOOL nMpScoreFlag = TRUE, qflag = TRUE;
	int  nTh = ((ANGLE+LENGTH+3)*600)/1000;
	int x, y, fdir, sdir;
	FPVECTEX saveFile = *pFile, saveSearch = *pSearch;
	int nLimitTh = 7*nMatchScoreTh/2;
	BYTE dScoreList[MAX_MINUTIA_NUMBER*MAX_MINUTIA_NUMBER];
	int dscore1 = 0, dscore2 = 0, dscore = 0;

	memset(dScoreList,0xFF,MAX_MINUTIA_NUMBER*MAX_MINUTIA_NUMBER);

	nRidge = 255 - abs(pSearch->nRidgeDensity - pFile->nRidgeDensity);

	if (pFile->Mp.quality < 35 || pSearch->Mp.quality < 35) nMpScoreFlag = FALSE;

	get_search_tag(pSearch,&SearchBar,&nMaxSearchBarLen,SDiffField,SArrangBarPtr,20,200);
	if (SearchBar.nNumber <= 0) return 0;

	arrange_points(pFile,pSearch,&SearchBar,&nMaxSearchBarLen,SArrangBarPtr,SDiffField,nCommonFlag);

	get_file_tag(pFile,&FileBar,FDiffField,FArrangBarPtr,&nFileCX,&nFileCY,20,nMaxSearchBarLen);
	if (FileBar.nNumber <= 0) return 0;

	score_sum = 0; nPairNum = 0;
	memset(AngleField, 0, sizeof(int)*240);
	for (i = 0; i < 240; i++) {
		nUpper = i + ANGLE; nLower = i - ANGLE;
		for (j = SDiffField[i]-1; j >= 0; j--) {
			score = 0; score_max = 0;
			id1 = SArrangBarPtr[i*MAX_SAMEDIFF_NUM+j];
			SBar = SearchBar.item[id1];
			nSid1 = SBar.nID1; nSid2 = SBar.nID2;
			for (h = nLower; h < nUpper; h++) {
				id = h - 240;
				if (id < -240) id += 480;
				else {
					if (id < 0) id += 240;
				}
				for (k = FDiffField[id]-1; k >= 0; k--)	{
					id2 = FArrangBarPtr[id*MAX_SAMEDIFF_NUM+k];
					FBar = FileBar.item[id2];
					nFid1 = FBar.nID1; nFid2 = FBar.nID2;
					lendiff = abs(FBar.nLen - SBar.nLen);
					if (lendiff >= LENGTH) continue;
					diff1 = abs(FBar.nDiff1 - SBar.nDiff1);
					diff1 = IANGLE_120(diff1);
					if (diff1 >= ANGLE) continue;
					diff2 = abs(FBar.nDiff2 - SBar.nDiff2);
					diff2 = IANGLE_120(diff2);
					if (diff2 >= ANGLE) continue;
					
					rdiff1 = pSearch->Mp.item[nSid1].dir - pFile->Mp.item[nFid1].dir;
					rdiff1 = ANGLE_0(rdiff1);
					rdiff2 = pSearch->Mp.item[nSid2].dir - pFile->Mp.item[nFid2].dir;
					rdiff2 = ANGLE_0(rdiff2);
					
					score = nTh - op_func_02(diff1*diff1+lendiff*lendiff+diff2*diff2);
					if (score < 0) continue;
					curv1 = 30 - abs(pSearch->Mp.item[nSid1].curv - pFile->Mp.item[nFid1].curv);
					if (curv1 < 0) continue;
					curv2 = 30 - abs(pSearch->Mp.item[nSid2].curv - pFile->Mp.item[nFid2].curv);
					if (curv2 < 0) continue;
					score = (score*curv1*curv2)/900;
					if (pFile->Mp.quality > 32 || pSearch->Mp.quality > 32) {
						if (pFile->Mp.item[nFid1].kind != pSearch->Mp.item[nSid1].kind) score = score*9/10;
						if (pFile->Mp.item[nFid2].kind != pSearch->Mp.item[nSid2].kind) score = score*9/10;
					}					
					
					if (pFile->Mp.item[nFid1].score < 25 || pSearch->Mp.item[nSid1].score < 25) score--;
					if (pFile->Mp.item[nFid2].score < 25 || pSearch->Mp.item[nSid2].score < 25) score--;

					if (nMpScoreFlag )	{
						flag = TRUE;
						x = (pFile->Mp.item[nFid1].x+pFile->Mp.item[nFid2].x)/2;
						y = (pFile->Mp.item[nFid1].y+pFile->Mp.item[nFid2].y)/2;
						x /= BLOCK_SIZE; y /= BLOCK_SIZE;
						fdir = pFile->Block.Data[y*MAX_BLOCK_COL+x];
						if (fdir == 0xFF) flag = FALSE; 
						fdir = abs(fdir - FBar.nSlope);
						fdir = IANGLE_60(fdir);
						x = (pSearch->Mp.item[nSid1].x+pSearch->Mp.item[nSid2].x)/2;
						y = (pSearch->Mp.item[nSid1].y+pSearch->Mp.item[nSid2].y)/2;
						x /= BLOCK_SIZE; y /= BLOCK_SIZE;
						sdir = pSearch->Block.Data[y*MAX_BLOCK_COL+x];
						if (sdir == 0xFF) flag = FALSE; 
						sdir = abs(sdir - SBar.nSlope);
						sdir = IANGLE_60(sdir);
						diff1 = abs(fdir-sdir);
						if (flag) {
							if (diff1 > 2*ANGLE) continue;
							else if (diff1 > ANGLE) score = score*8/10;
						}
					}

					if (score > 0) {
						AngleField[rdiff1] += score;
						AngleField[rdiff2] += score;
						PairList[nPairNum].score = score;
						PairList[nPairNum].fid = id2;
						PairList[nPairNum].sid = id1;
						if (++nPairNum >= 1000)	{
							i = h = 1024; j = k = -1; break; 
						}
					}
					if (score > score_max) score_max = score;
				}
			}
			score_sum += score_max;
		}
	}

	if (nCommonFlag && 2*score_sum < nPairNum*5 && score_sum < nMatchScoreTh*3) return 0;
	if (SearchBar.nNumber > 100) 
		score_sum = (score_sum * 200) / SearchBar.nNumber;
	else 
		score_sum = score_sum * 2;
	score_sum = (score_sum*929 + 1137 + 500) / 1000;

	if (score_sum < nMatchScoreTh) return 0;

	nRot = rotate_points(nFileCX,nFileCY,AngleField,&FileBar,pFile);

	memset( XField, 0, sizeof(int)*640 );
	memset( YField, 0, sizeof(int)*640 );
	nLastNum = 0; score_sum = score_max = 0;
	id = PairList[0].sid;
	for (i = 0; i < nPairNum; i++) {
		score = PairList[i].score;
		FBar = FileBar.item[PairList[i].fid];
		SBar = SearchBar.item[PairList[i].sid];
		dslope = abs(FBar.nSlope - SBar.nSlope);
		dslope = IANGLE_60(dslope);
		if (dslope >= ANGLE) continue;
		get_shift_param(LENGTH,score,&FBar,&SBar,XField,YField,&pFile->Mp,&pSearch->Mp);
		LastList[nLastNum++] = i;
		if (nLastNum == 500) break; 
		if (id != PairList[i].sid) {
			score_sum += score_max; score_max = 0;
			id = PairList[i].sid;
		}
		if (score_max < PairList[i].score) score_max = PairList[i].score;
	}
	score_sum += score_max;
	if (SearchBar.nNumber > 100) 
		score_sum = (score_sum * 200) / SearchBar.nNumber;
	else 
		score_sum = score_sum * 2;
	score_sum = (score_sum*929 + 1137 + 500) / 1000;

	if (score_sum < nMatchScoreTh) return 0;

	shift_points(&nXoffset,&nYoffset,pFile,XField,YField);
	transform_block(nRot,nXoffset,nYoffset,nFileCX,nFileCY,&pFile->Block);
	nBlockScore = check_block(30,5,&pFile->Block,&pSearch->Block);
	if (nBlockScore < 75 && nPairNum < 110) return 0;
	nCommonScore = check_block2(4,&pFile->Block,&pSearch->Block);
	if (nCommonScore < 30) return 0;
	if (nBlockScore < 80 && nCommonScore < 50) return 0;

	if (FALSE == re_arrange_point(PairList,LastList,nPairNum,&nLastNum,pFile,pSearch,&FileBar,&SearchBar))
		return 0;

	newCX = 0; newCY = 0;
	for (i = 0; i < pFile->Mp.nNumber; i++)	{
		newCX += pFile->Mp.item[i].x; newCY += pFile->Mp.item[i].y;
	}
	newCX /= pFile->Mp.nNumber; newCY /= pFile->Mp.nNumber;

	score_sum = 0;
	for (i = 0,j = 0; i < nLastNum; i++) {
		FBar = FileBar.item[PairList[LastList[i]].fid];
		SBar = SearchBar.item[PairList[LastList[i]].sid];
		nCheck = check_limit(LENGTH,&FBar,&SBar,pFile,pSearch,newCX,newCY);
		if (!nCheck) continue;
		LastList[j++] = LastList[i];
	}
	if(j <= 0) return 0;
	nLastNum = j;
	for (i = 0; i < nLastNum; i++) {
		if (LastList[i] == -1) continue;
		score_max = PairList[LastList[i]].score;
		PairTemp[0] = PairList[LastList[i]];
		nNum = 1;
		do {
			nRepeat = FALSE;
			for (j = 0; j < nLastNum; j++) {
				if (i == j) continue;
				if (LastList[j] == -1) continue;
				nSid = PairList[LastList[j]].sid;
				nFid = PairList[LastList[j]].fid;
				for (k = 0; k < nNum; k++) {
					if (PairTemp[k].sid == nSid) break;
					if (PairTemp[k].fid == nFid) break;
				}
				if (k >= nNum) continue;
				nRepeat = TRUE;
				PairTemp[nNum++] = PairList[LastList[j]];
				if (nNum >= 200-1) { nRepeat = FALSE; break; }
				if (score_max < PairList[LastList[j]].score)
					score_max = PairList[LastList[j]].score;
				LastList[j] = -1;
			}
		} while (nRepeat == TRUE);
		if (nNum <= 3) { score_sum += score_max; continue; }
		else {
			snum = fnum = 0;
			for (j = 0; j < nNum; j++) {
				sflag = fflag = TRUE;
				nSid = PairTemp[j].sid; nFid = PairTemp[j].fid;
				for (k = 0; k < j; k++)	{
					if (PairTemp[k].sid == nSid) sflag = FALSE;
					if (PairTemp[k].fid == nFid) fflag = FALSE;
				}
				if (sflag == TRUE) snum++;
				if (fflag == TRUE) fnum++;
			}
			if (snum > fnum){ flag = TRUE; snum = fnum; }
			else flag = FALSE;
			for (j = 0; j < snum; j++) {
				hflag = TRUE;
				for (k = 0; k < nNum; k++) {
					if (PairTemp[k].score == -1) continue;
					if (hflag == TRUE) {
						if (flag == TRUE) {
							id1 = PairTemp[k].fid;
							score_max = PairTemp[k].score; hflag = FALSE;
						}
						else {
							id1 = PairTemp[k].sid;
							score_max = PairTemp[k].score; hflag = FALSE;
						}
					}
					else {
						if (flag == TRUE) id2 = PairTemp[k].fid;
						else id2 = PairTemp[k].sid;
						if (id1 != id2) continue;
						if (score_max < PairTemp[k].score) 
							score_max = PairTemp[k].score;
					}
					PairTemp[k].score = -1;
				}
				score_sum += score_max;
			}
		}
	}

	if (pPair != NULL) {
		id = 0;
		for (i = 0; i < nLastNum; i++) {
			if (LastList[i] == -1) continue;
			FBar = FileBar.item[PairList[LastList[i]].fid];
			nFid = FBar.nID1;
			SBar = SearchBar.item[PairList[LastList[i]].sid];
			nSid = SBar.nID1;
			for (j = 0; j < id; j++) {
				if (nFid == pPair->nFileID[j]) break;
				if (nSid == pPair->nSearchID[j]) break;
			}
			if (j >= id) {
				pPair->nFileID[id] = nFid; 
				pPair->nSearchID[id++] = nSid;
			}
			nFid = FBar.nID2; nSid = SBar.nID2;
			for (j = 0; j < id; j++) {
				if (nFid == pPair->nFileID[j]) break;
				if (nSid == pPair->nSearchID[j]) break;
			}
			if (j >= id) {
				pPair->nFileID[id] = nFid; 
				pPair->nSearchID[id++] = nSid;
			}
		}
		pPair->nNumber = id; pPair->nRot = nRot; 
		pPair->nXOffset = nXoffset; pPair->nYOffset = nYoffset;
		pPair->nXC = nFileCX; pPair->nYC = nFileCY;
	}

	if (score_sum < nMatchScoreTh/2) return 0;
	nNum = pFile->Mp.nNumber;
	if (nNum > pSearch->Mp.nNumber) nNum = pSearch->Mp.nNumber;

	tmpNum = GetMatchedTemplateNum(&pFile->Mp,&pSearch->Mp,pPair);
	if (tmpNum < 0) return 0;

	if (*nGlobalScore > 0) score = *nGlobalScore;
	else {
		score = get_point_score(&saveFile,&saveSearch);
		*nGlobalScore = score;
		if (score > 1700) return (nMatchScoreTh*2);
	}

	score_old = score_sum;
	sflag = check_core(&pFile->Core,&pSearch->Core,10,ANGLE);

	if (abs(pFile->Mp.quality - pSearch->Mp.quality) >= 12) qflag = FALSE;
	fflag = check_paired_mp(&pFile->Mp,&pSearch->Mp,pPair);

	score_sum = adjust_score(score_sum,score,nBlockScore,nCommonScore,dscore,
				pPair->nNumber,nRidge,SearchBar.nNumber,nNum,sflag,tmpNum,qflag,fflag);

		dslope = 0;
		if (pFile->nFrequency != 0 && pSearch->nFrequency != 0) {
			fnum = snum = pFile->nFrequency;
			if (pSearch->nFrequency > fnum) fnum = pSearch->nFrequency;
			if (pSearch->nFrequency < snum) snum = pSearch->nFrequency;
			dslope = fnum - snum;
		}

		if (pPair->nNumber > 16 && sflag && score_old > 240) score_sum = score_old;

		if (dslope > 21 && score < 1200) {
			score_sum = score_sum*snum/fnum;
			if ((pFile->nType == 7 && (pSearch->nType==4 || pSearch->nType==5))
				|| (pSearch->nType == 7 && (pFile->nType==4 || pFile->nType==5))) {
				score_sum = score_sum*snum/fnum;
			}
		}

		if (score_sum > nMatchScoreTh && nBlockScore >= 83)	{
			if (score > 1060 && pPair->nNumber*100 > nNum*50 && nRidge > 220) {
				if (pPair->nNumber >= 12 && tmpNum > 0 && nBlockScore >= 88) score_sum = score_sum*2;
			}
			if (pPair->nNumber > 13 && pPair->nNumber*100 >= 56*nNum && nRidge > 241 && tmpNum > 0)
				return (score_sum);		
		}
		if (tmpNum > 5 && pPair->nNumber >= 12 && pPair->nNumber*100 > nNum*44 && nRidge > 250) score_sum *= 2;

		if (pFile->Mp.quality < 32 && pSearch->Mp.quality < 32)	{
		}
		else {
				score_sum = score_sum*nBlockScore*nBlockScore/10000;
		}

	if (score_sum < nMatchScoreTh) return 0;

	
		if (nCommonFlag) {
			if (tmpNum == 0 && score_sum < nMatchScoreTh*12/10 && nRidge+nBlockScore < 335 && nBlockScore < 95) {
				if (score < 1100 || pPair->nNumber < 9 || pPair->nNumber*100 < 40*nNum)
					return 0;
			}
		}

	flag = TRUE;
	fflag = hflag = FALSE;
	
		if (pFile->nType == 8 || pSearch->nType == 8) flag = TRUE;
		else if (pFile->nType != pSearch->nType) flag = FALSE;
		fflag = ((TRUE == check_core(&pFile->Core,&pSearch->Core,16,7)) && nRidge+nBlockScore < 330 && score < 1250);
		hflag = (score < 1350 || pPair->nNumber <= 10 || pPair->nNumber*100 < nNum*38);
		nflag = (score > 1000 && pPair->nNumber > 10 && nBlockScore > 92 && sflag && flag);
		flag = (score<1170 && tmpNum<=3 && (pPair->nNumber<12 || nCommonScore<80)) || (score<930 && tmpNum<=2) || (score<460 && nCommonScore<75);


	if (score_sum < nMatchScoreTh*4) {
		score_sum = dec_func_01(score_sum,pFile,pSearch,pPair,tmpNum);
		if ((tmpNum < 5 && (fflag || hflag))) {
			score_sum = dec_func_03(score_sum,pFile,pSearch,nRidge,nBlockScore,tmpNum);
		}
	}

	if (pFile->Mp.quality < 35 && pSearch->Mp.quality < 35)	{
		nLimitTh = 5*nMatchScoreTh/2;
	}

	if (score_sum <= nLimitTh) {
		score_sum = dec_func_07(score_sum,&pFile->Mp,&pSearch->Mp,pPair);
		if ((hflag || fflag) && score < 1500)
			score_sum = dec_func_10(score_sum,tmpNum,&pFile->Mp,&pSearch->Mp,pPair,nBlockScore,nRidge,score,sflag,nCommonScore,dscore);
		if ((!nflag && hflag && flag))
			score_sum = dec_func_11(score_sum,pFile,pSearch,pPair,nRidge,score,dscore);
		if ((nMpScoreFlag && score < 1100 && nBlockScore < 96) )
			score_sum = dec_func_12(score_sum,pFile,pSearch,pPair);
		if ((tmpNum < 4 && pPair->nNumber*100 < nNum*50))
			score_sum = dec_func_13(score_sum,pFile,pSearch,pPair,nMatchScoreTh,nRidge,nBlockScore,dscore);
	}

	return score_sum;
}

int matching_main(LPFPFEATURE pFeatureVect1,LPFPFEATURE pFeatureVect2,int securitylevel)
{
	FPVECTEX Vect1, Vect2, TmpVect1, TmpVect2;
	PAIRVECT pPair;
	int score, nCoarse = 0, nGlobalScore = 0, nTh = MATCH_TH_MEDIUM;

	if ( securitylevel == HIGH_LEVEL )	nTh = MATCH_TH_HIGH;
	if ( securitylevel == LOW_LEVEL )	nTh = MATCH_TH_LOW;

	mch_sub_func_02(pFeatureVect1,&Vect1);
	TmpVect1 = Vect1;
	mch_sub_func_02(pFeatureVect2,&Vect2);
	TmpVect2 = Vect2;

	if (mch_sub_func_03(&Vect1) == FALSE) return (-1);
	if (mch_sub_func_03(&Vect2) == FALSE) return (-1);

	nCoarse = coarse_matching(&Vect1,&Vect2);

	if (nCoarse == 1) {
		return (1);
	}

	if (type_matching(&Vect1,&Vect2) == FALSE) {
		return (-1);
	}

	score = point_matching(&Vect1,&Vect2,&pPair,FALSE,TRUE,nTh,&nGlobalScore);
	if (score < nTh) {
		Vect1 = TmpVect1; Vect2 = TmpVect2;
		score = point_matching(&Vect2,&Vect1,&pPair,FALSE,TRUE,nTh,&nGlobalScore);
		if (score < nTh) {
			if (nCoarse == 2 && score > nTh/2) return (1);
			return (-1);
		}
	}
	return (1);
}

void sch_sub_func_03(LPMPVECTEX pVect,int cx,int cy,int nAngle,int nDiffX,int nDiffY)
{
	int i, x, y, angle, rot, nCos, nSin, dx, dy;
	int nX = cx + nDiffX, nY = cy + nDiffY;

	rot = 240 - nAngle;
	rot = ANGLE_240(rot);

	nCos = _table_03[rot]; nSin=_table_04[rot];

	for (i = 0; i < pVect->nNumber; i++) {
		dx = pVect->item[i].x-cx; dy = pVect->item[i].y-cy;
		x = dx*nCos + dy*nSin;
		x = x >> 14;
		y = dy*nCos - dx*nSin;
		y = y >> 14;
		pVect->item[i].x = x + nX;
		pVect->item[i].y = y + nY;

		angle = pVect->item[i].dir + nAngle;
		angle = ANGLE_0_240(angle);
		pVect->item[i].dir = angle;
	}
}

int sch_sub_func_04(LPMPVECTEX pVect1,LPMPVECTEX pVect2)
{
	int i, j, x, y, dir, dx, dy, diff, val, minval, score = 0;

	for (i = 0; i < pVect1->nNumber; i++) {
		x = pVect1->item[i].x; y = pVect1->item[i].y;
		dir = pVect1->item[i].dir;
		minval = 10000;
		for (j = 0; j < pVect2->nNumber; j++) {
			dx = abs(pVect2->item[j].x - x);
			if (dx > LENGTH+3) continue;
			dy = abs(pVect2->item[j].y - y);
			if (dy > LENGTH+3) continue;
			diff = abs(pVect2->item[j].dir - dir);
			diff = IANGLE_120(diff);
			if (diff > ANGLE) continue;
			val = dx + dy + diff;
			if (minval > val) minval = val;
			if (minval < 20) break;
		}
		if (minval < 35) score += (35 - minval);
	}
	val = (pVect1->nNumber + pVect2->nNumber) / 2;
	if (val == 0) return (0);
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
	
	for (i = 0; i < nFCoreNum; i++)	{
		cx0 = FileCore[i].x; cy0 = FileCore[i].y;
		dir0 = FileCore[i].dir;
		radTh = 100; nn = 0;
		while (1) {
			if (radTh > 200) break;
			for (k = 0; k < pFPEx0->Mp.nNumber; k++) {
				dx = pFPEx0->Mp.item[k].x - cx0;
				dy = pFPEx0->Mp.item[k].y - cy0;
				if (dx*dx+dy*dy < radTh*radTh) {
					tmpMPEx0.item[nn] = pFPEx0->Mp.item[k];
					nn++;
				}
			}
			if (nn < nNumTh){ radTh += 20; continue; }
			tmpMPEx0.nNumber = nn; break;
		}
		for (j = 0; j < nSCoreNum; j++)	{
			cx1 = SearchCore[j].x; cy1 = SearchCore[j].y;
			dir1 = SearchCore[j].dir;

			rot = dir0 - dir1;
			rot = ANGLE_0(rot);
			radTh = 100; nn = 0;
			while (1) {
				if (radTh > 200) break;
				for (k = 0; k < pFPEx1->Mp.nNumber; k++) {
					dx = pFPEx1->Mp.item[k].x - cx1;
					dy = pFPEx1->Mp.item[k].y - cy1;
					if (dx*dx+dy*dy < radTh*radTh) {
						tmpMPEx1.item[nn] = pFPEx1->Mp.item[k];
						nn++;
					}
				}
				if (nn < nNumTh) { radTh += 20; continue; }
				tmpMPEx1.nNumber = nn; break;
			}
			dx = cx0 - cx1; dy = cy0 - cy1;
			sch_sub_func_03(&tmpMPEx1,cx1,cy1,rot,dx,dy);
			nn = sch_sub_func_04(&tmpMPEx0,&tmpMPEx1);
			if (score < nn) score = nn;
		}
	}
	return(score);
}

int sch_sub_func_01(LPFPVECTEX pFile,LPFPVECTEX pSearch) 
{
	int i, j, k, cx, cy, dx, dy, rot;
	int retScore = 0, score, maxscore, nNumF = 0, nMaxF = 0, nMaxS = 0;
	MINUTIAEX tmpF[5];
	MPVECTEX tmpVect = pFile->Mp;
	int nFCoreNum, nSCoreNum;
	COREITEMEX FileCore[2], SearchCore[2];

	nFCoreNum = mch_sub_func_01(&(pFile->Core),FileCore,NULL,NULL);
	nSCoreNum = mch_sub_func_01(&(pSearch->Core),SearchCore,NULL,NULL);

	if (pFile->Mp.nNumber < MIN_MINUTIA_NUMBER || pSearch->Mp.nNumber < MIN_MINUTIA_NUMBER) return (0);

	if (nFCoreNum != 0 && nSCoreNum != 0) {
		retScore = sch_sub_func_05(pSearch,pFile);
		return (retScore);
	}

	for (i = 0; i < pFile->Mp.nNumber; i++)	{
		if (nMaxF < pFile->Mp.item[i].curv) nMaxF = pFile->Mp.item[i].curv;
	}
	for (i = 0; i < pSearch->Mp.nNumber; i++) {
		if (nMaxS < pSearch->Mp.item[i].curv) nMaxS = pSearch->Mp.item[i].curv;
	}
	if (nMaxF > nMaxS) nMaxF = nMaxS;

	for (i = 0, j = 0; i <pFile->Mp.nNumber; i++) {
		if (pFile->Mp.item[i].curv >= nMaxF+8) continue;
		for (k = 0; k < j; k++)	{
			dx = pFile->Mp.item[i].x - tmpF[k].x;
			dy = pFile->Mp.item[i].y - tmpF[k].y;
			if (dx*dx+dy*dy <= 40*40) break;
		}
		if (k < j) continue;
		tmpF[j++] = pFile->Mp.item[i];
		if (j >= 5) break;
	}
	nNumF = j;
	pFile->Mp = tmpVect;
	for (i = 0; i < nNumF; i++)	{
		cx = tmpF[i].x; cy = tmpF[i].y;
		maxscore = 0;
		for (j = 0; j < pSearch->Mp.nNumber; j++) {
			if (abs(tmpF[i].curv-pSearch->Mp.item[j].curv) >= 6) continue;
			rot = pSearch->Mp.item[j].dir - tmpF[i].dir;
			rot = ANGLE_0(rot);
			dx = pSearch->Mp.item[j].x - cx; dy = pSearch->Mp.item[j].y - cy;
			sch_sub_func_03(&pFile->Mp,cx,cy,rot,dx,dy);
			score = sch_sub_func_04(&pFile->Mp,&pSearch->Mp);
			if (maxscore < score) maxscore = score;
			pFile->Mp = tmpVect;
		}
		if (retScore < maxscore) retScore = maxscore;
	}
	return (retScore);
}

void sch_sub_func_02(int *pScore,int nSize,short *pIndex)
{
	int i, j, tmp, nMin;

	for (i = 0; i < nSize; i++) pIndex[i] = i;
	
	nMin = nSize - 1;
	if (nMin > 10) nMin = 10;

	for (i = 0; i < nMin; i++) {
		for (j = i+1; j < nSize; j++) {
			if (pScore[pIndex[i]] < pScore[pIndex[j]]) {
				tmp = pIndex[i]; pIndex[i] = pIndex[j]; pIndex[j] = tmp;
			}
		}
	}
}


BOOL sch_sub_func(LPFPFEATURE pFeatureVect,LPFPFEATURE pDataBaseVects,int nDataBaseSize,short* pIndex) 
{
	int i, *pScore;
	FPVECTEX FileFPEx, SearchFPEx;
	
	pScore = (int*)malloc(sizeof(int)*nDataBaseSize);
	if (pScore == NULL) return FALSE;

	mch_sub_func_02(pFeatureVect,&SearchFPEx);

	for (i = 0; i < nDataBaseSize; i++)	{
		mch_sub_func_02(&(pDataBaseVects[i]),&FileFPEx);
		pScore[i] = sch_sub_func_01(&FileFPEx,&SearchFPEx);
	}

	sch_sub_func_02(pScore,nDataBaseSize,pIndex);
	free(pScore);
	return TRUE;
}

  
// exported function //

/*
 *	fingerprint matching function.
 *  parameter ;
 *		pFeature1 : pointer to first fingerprint feature template buffer
 *		pFeature2 : pointer to second fingerprint feature template buffer
 *		securitylevel : value of security level (default - MEDIUM_LEVEL)
 *  return value ;
 *		if success, return 1.
 *		else , return other value
 */
int finger_match(BYTE* pFeature1,BYTE* pFeature2,int securitylevel)
{
	int res;
	LPFPFEATURE pVect1 = (LPFPFEATURE)(pFeature1);
	LPFPFEATURE pVect2 = (LPFPFEATURE)(pFeature2);
	int nTh = MATCH_TH_MEDIUM;

	if ( securitylevel == HIGH_LEVEL )	nTh = MATCH_TH_HIGH;
	if ( securitylevel == LOW_LEVEL )	nTh = MATCH_TH_LOW;


	if ( pFeature1 == NULL || pFeature2 == NULL ) return (0);
	
	res = matching_main(pVect1,pVect2,securitylevel);
	/* O retorno da função matching_main é 1 ou -1!!!
	if (res < 0) res = 0;
	if (res > 1000) res = 1000;

	if ( res >= nTh ) return ERR_OK;
	return ERR_MATCH_FAILED;
	*/
	return (res);
}


/*
 *	fingerprint searching function.
 *  parameter ;
 *		pFeature : pointer to inputed fingerprint feature template buffer
 *		pDBFeature : pointer to buffer of fingerprint features(size is 512 * nDBSize bytes)
 *		nDBSize : number of fingerprint feature in database
 *		securitylevel : value of security level (default - MEDIUM_LEVEL)
 *  return value ;
 *		if success, return the index value of matched fingerprint feature (0 ~ nDBSize).
 *		if failed, return error code ( < 0 ).
 */
int finger_search(BYTE* pFeature,BYTE* pDBFeature,int nDBSize,int securitylevel)
{
	int i, idx, res, nMin;
	short *pIndexArray;
	LPFPFEATURE pVect = (LPFPFEATURE)(pFeature);
	LPFPFEATURE pDBVect = (LPFPFEATURE)(pDBFeature);
	int nTh = MATCH_TH_MEDIUM;

	if ( securitylevel == HIGH_LEVEL )	nTh = MATCH_TH_HIGH;
	if ( securitylevel == LOW_LEVEL )	nTh = MATCH_TH_LOW;


	if ( nDBSize < 1 ) return ERR_GENERAL_ERROR;
	pIndexArray = (short*)malloc(sizeof(int)*nDBSize);
	if ( pIndexArray == NULL ) return ERR_CAN_NOT_ALLOC_MEMORY;

	if ( nDBSize == 1 ){ pIndexArray[0] = 0; }
	else{
		if ( ERR_OK != sch_sub_func(pVect,pDBVect,nDBSize,pIndexArray) ){
			free(pIndexArray); return ERR_GENERAL_ERROR;
		}
	}

	nMin = nDBSize;
	/* Percorre o banco inteiro de templates
	if (nMin > 10) nMin = 10;
	*/
	
	for (i = 0; i < nMin; i++) {
		idx = pIndexArray[i];
		if ( idx < 0 || idx >= nDBSize ) continue;
		res = matching_main(&pDBVect[idx],pVect,securitylevel);
		if ( res == 1 ) {
			free(pIndexArray);
			return (idx);
		}
	}		

	free(pIndexArray);
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
	int ret;

	LPFPFEATURE pVect = (LPFPFEATURE)(pFeature);

	if ( pImage == NULL || pFeature == NULL ) return (ERR_GENERAL_ERROR);
	if ( nWidth < 0 || nWidth > MAX_IMG_WIDTH || nHeight < 0 || nHeight > MAX_IMG_HEIGHT ) return ERR_INVALID_IMAGESIZE;

	ret = ext_main(pImage,nWidth,nHeight,pVect);
	
	return ret;
}
