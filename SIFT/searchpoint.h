#ifndef searchpoint
#define searchpoint
#include "GUSSIANdelta.h"
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
using namespace cv;
using namespace std;
//num为每个oct检测出来的数目 在外部 *num需要malloc 在栈上开好空间
CvMat *SearchP(delta dogdelta,int *num)
{
	int xoct=dogdelta.longoct;
	int xinte=dogdelta.longinte;
	int i,j,k;
	int p=0;
	int h=0;
	CvMat *search[30];
	int searchx[30];
	double r=0;
	for (k=0;k<xoct;k++)
	{
		searchx[k]=0;//每个oct检测出来的点的个数
		search[k]=cvCreateMat(2,100000,CV_64FC1);//用来存放点
		CvMat *spc;
		CvMat *spb;
		CvMat *spa;
		spc=dogdelta.oct[k].inter[2].Gussianimg;
		spb=dogdelta.oct[k].inter[1].Gussianimg;
		spa=dogdelta.oct[k].inter[0].Gussianimg;
  		for (i=1;i<spb->rows-1;i++)
		{
			for (j=1;j<spb->cols-1;j++)
			{
				double a[3][3];
				double b[3][3];
				double c[3][3];

				b[0][0]=*(spb->data.db+(i-1)*spb->step/8+j-1);
				b[0][1]=*(spb->data.db+(i-1)*spb->step/8+j);
				b[0][2]=*(spb->data.db+(i-1)*spb->step/8+j+1);
				b[1][0]=*(spb->data.db+(i)*spb->step/8+j-1);
				b[1][1]=*(spb->data.db+(i)*spb->step/8+j);
				b[1][2]=*(spb->data.db+(i)*spb->step/8+j+1);
				b[2][0]=*(spb->data.db+(i+1)*spb->step/8+j-1);
				b[2][1]=*(spb->data.db+(i+1)*spb->step/8+j);
				b[2][2]=*(spb->data.db+(i+1)*spb->step/8+j+1);

				a[0][0]=*(spa->data.db+(i-1)*spa->step/8+j-1);
				a[0][1]=*(spa->data.db+(i-1)*spa->step/8+j);
				a[0][2]=*(spa->data.db+(i-1)*spa->step/8+j+1);
				a[1][0]=*(spa->data.db+(i)*spa->step/8+j-1);
				a[1][1]=*(spa->data.db+(i)*spa->step/8+j);
				a[1][2]=*(spa->data.db+(i)*spa->step/8+j+1);
				a[2][0]=*(spa->data.db+(i+1)*spa->step/8+j-1);
				a[2][1]=*(spa->data.db+(i+1)*spa->step/8+j);
				a[2][2]=*(spa->data.db+(i+1)*spa->step/8+j+1);

				c[0][0]=*(spc->data.db+(i-1)*spc->step/8+j-1);
				c[0][1]=*(spc->data.db+(i-1)*spc->step/8+j);
				c[0][2]=*(spc->data.db+(i-1)*spc->step/8+j+1);
				c[1][0]=*(spc->data.db+(i)*spc->step/8+j-1);
				c[1][1]=*(spc->data.db+(i)*spc->step/8+j);
				c[1][2]=*(spc->data.db+(i)*spc->step/8+j+1);
				c[2][0]=*(spc->data.db+(i+1)*spc->step/8+j-1);
				c[2][1]=*(spc->data.db+(i+1)*spc->step/8+j);
				c[2][2]=*(spc->data.db+(i+1)*spc->step/8+j+1);

 				p=0;
				h=0;
				int all=0;
				int allx;
				for (int m=0;m<3;m++)
				{
					for (int n=0;n<3;n++)
					{
						//all+=a[m][n]+b[m][n]+c[m][n];
						//all+=abs(a[m][n]-b[1][1])+abs(b[m][n]-b[1][1])+abs(c[m][n]-b[1][1]);
						if (b[1][1]<a[m][n]-r)
						{
							p++;
						}
						if (b[1][1]<b[m][n]-r)
						{
							p++;
						}
						if (b[1][1]<c[m][n]-r)
						{
							p++;
						}

						if (b[1][1]>a[m][n]+r)
						{
							h++;
						}
						if (b[1][1]>b[m][n]+r)
						{
							h++;
						}
						if (b[1][1]>c[m][n]+r)
						{
							h++;
						}
					}
				}
				allx=all-b[1][1];
				
				//if (allx>5)
				//{
					if (p==26||h==26)
					{
  						*(search[k]->data.db+0*search[k]->step/8+searchx[k])=i;
						*(search[k]->data.db+1*search[k]->step/8+searchx[k])=j;
						searchx[k]++;
				//	}

					}
				
			}
		}



	}
	//for (int l=0;l<xoct;l++)
	//{
	//	MatOUT(search[l],2,searchx[l],1);
	//}

	CvMat *allsearch=cvCreateMat(2*xoct,100000,CV_64FC1);
	for(int q=0;q<2*xoct;q=q+2)
	{
		for(int w=0;w<searchx[(int)(q/2)];w++)
		{

  			*(allsearch->data.db+q*allsearch->step/8+w)=*(search[(int)(q/2)]->data.db+0*search[(int)(q/2)]->step/8+w);
			*(allsearch->data.db+(q+1)*allsearch->step/8+w)=*(search[(int)(q/2)]->data.db+1*search[(int)(q/2)]->step/8+w);
		}
		*(num+(int)(q/2))=searchx[(int)(q/2)];
	}

	return allsearch;
	



}











#endif