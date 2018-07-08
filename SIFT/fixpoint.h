#ifndef fixpoint
#define fixpoint
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"
#include "GUSSIANdelta.h"
#include "deletepoint.h"
using namespace cv;
using namespace std;


void FixPoint(CvMat* spoint,int *num,delta dogdelta)
{
	int xoct=dogdelta.longoct;
	int xinte=dogdelta.longinte;
	int i,j;
	int m,n;
	int xrow;
	int xcol;
	int x,y;
	int h=1;
	double dxx,dxy,dyy,dx,dy;
	double extrmx,dex;
	double extrmy,dey;
	double sdey;
	double D_XY;
	double tdey;
	double tr,det;
	double inr=10;      //====曲率因子
	double fininr,fink;
	CvMat *del=cvCreateMat(xoct,10000,CV_64FC1);
	int dele[10];
	int numk=0;
	//=======计算一阶导与二阶导=======
	for (i=0;i<xoct;i++)
	{
		m=2*i;
		CvMat *fixing=cvCreateMat(dogdelta.oct[i].inter[1].Gussianimg->rows,dogdelta.oct[i].inter[1].Gussianimg->cols,CV_64FC1);
		fixing=cvCloneMat(dogdelta.oct[i].inter[1].Gussianimg);
		xrow=fixing->rows;
  		xcol=fixing->cols;
		dele[i]=0;
		
		for (n=0;n<*(num+i);n++)
		{
			x=*(spoint->data.db+m*spoint->step/8+n); 
			y=*(spoint->data.db+(m+1)*spoint->step/8+n);

			dx=(*(fixing->data.db+x*fixing->step/8+y+1)-*(fixing->data.db+x*fixing->step/8+y-1))/(double)(2.0*h);
			dy=(*(fixing->data.db+(x+1)*fixing->step/8+y)-*(fixing->data.db+(x-1)*fixing->step/8+y))/(double)(2.0*h);
			dxx=(*(fixing->data.db+x*fixing->step/8+y+1)+*(fixing->data.db+x*fixing->step/8+y-1)-*(fixing->data.db+x*fixing->step/8+y)*2)/(double)(pow(h,2.0));
			dyy=(*(fixing->data.db+(x+1)*fixing->step/8+y)+*(fixing->data.db+(x-1)*fixing->step/8+y)-*(fixing->data.db+x*fixing->step/8+y)*2)/(double)(pow(h,2.0));
  			dxy=((*(fixing->data.db+(x+1)*fixing->step/8+y+1)+*(fixing->data.db+(x-1)*fixing->step/8+y-1))-(*(fixing->data.db+(x-1)*fixing->step/8+y+1)+*(fixing->data.db+(x+1)*fixing->step/8+y-1)))/(double)(4.0*pow(h,2.0));

			//==========消除曲率响应=============
			tr=dxx+dyy;
			det=dxx*dyy-dxy*dxy;
			fininr=(inr+1)*(inr+1)/inr;
			fink=tr*tr/det;
			if (fink>=fininr)
			{
				*(del->data.db+i*del->step/8+dele[i])=n;
				dele[i]++;
			}
			//避免做检测时有特征点邻域溢出====
			int l=40;
			if (x>xrow-l*1.4)
			{
				*(del->data.db+i*del->step/8+dele[i])=n;
				dele[i]++;
			}

			if (x<l*1.4)
			{
				*(del->data.db+i*del->step/8+dele[i])=n;
				dele[i]++;
			}

			if (y>xcol-l*1.4)
			{
				*(del->data.db+i*del->step/8+dele[i])=n;
				dele[i]++;
			}

			if (y<l*1.4)
			{
				*(del->data.db+i*del->step/8+dele[i])=n;
				dele[i]++;
			}
			//===========================

			if (dxy!=0) //&& dx!=0 && dy!=0)
			{
				dex=(dy-dx*dyy/dxy)/(dxx*dyy/dxy-dxy);
				dey=(-dx-dex*dxx)/(dxy);
				sdey=(-dy-dex*dxy)/(dyy);
				tdey=(dy-dx*dxy/dxx)/(dxy*dxy/dxx-dyy);
				D_XY=dex*dx+dey*dy+0.5*dex*dex*dxx+0.5*dey*dey*dyy+dex*dey*dxy;


				if (abs(D_XY)<0.03)
				{
					*(del->data.db+i*del->step/8+dele[i])=n;
					dele[i]++;
				}
	
				*(spoint->data.db+m*spoint->step/8+n)=*(spoint->data.db+m*spoint->step/8+n)+dex;
				*(spoint->data.db+(m+1)*spoint->step/8+n)=*(spoint->data.db+(m+1)*spoint->step/8+n)+dey;
				



			}
			else
			{
	/*			*(del->data.db+i*del->step/8+dele[i])=n;
				dele[i]++;
*/
			}
		}
		

	}


	DelPoint(spoint,del,num,dele);
	

	//for (int o=0;o<xoct;o++)
	//{
	//	MatOUT(spoint,2,*(num+o),1);
	//}




}



#endif