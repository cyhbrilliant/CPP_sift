#ifndef fixmatch
#define fixmatch
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"
#include "GUSSIANdelta.h"
//#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace cv;
using namespace std;



//==========================================================双线性插值，返回像素间的灰度值
#define ImMat(ROW,COL) ((float *)(im->data.fl + im->step/sizeof(float) *(ROW)))[(COL)]

float getPixelBI(CvMat * im, float col, float row) 
{


	int irow, icol;
	float rfrac, cfrac;
	float row1 = 0, row2 = 0;
	int width=im->cols;
	int height=im->rows;


	irow = (int) row;
	icol = (int) col;

	if (irow < 0 || irow >= height|| icol < 0 || icol >= width)
		return 0;
	if (row > height - 1)
		row = height - 1;
	if (col > width - 1)
		col = width - 1;
	rfrac = 1.0 - (row - (float) irow);
	cfrac = 1.0 - (col - (float) icol);
	if (cfrac < 1) 
	{
		row1 = cfrac * ImMat(irow,icol) + (1.0 - cfrac) * ImMat(irow,icol+1);
	} 
	else 
	{
		row1 = ImMat(irow,icol);
	}
	if (rfrac < 1) 
	{
		if (cfrac < 1) 
		{
			row2 = cfrac * ImMat(irow+1,icol) + (1.0 - cfrac) * ImMat(irow+1,icol+1);
		} else 
		{
			row2 = ImMat(irow+1,icol);
		}
	}
	return rfrac * row1 + (1.0 - rfrac) * row2;
}



class histo
{
public:
	CvMat *eight_histo;
};
class point_Feature
{
public:
	histo blocks[4][4];
	CvMat *Feature_Vector;
};
class img_Feature
{
public:
	int Feature_num;
	point_Feature p_Feature[10000];
};
class Feature
{
public:
	int oct;
	img_Feature oct_imgFeature[10];
};






//===============//
//传参数时 x y要反着传参数//
//===============//







//====================================================================================

//初始化特征类============================================================================
void  init_Feature(Feature *F,int xoct,int *num)
{
	int i,j,k,l,m;


	F->oct=xoct;
	for (i=0;i<xoct;i++)
	{
		F->oct_imgFeature[i].Feature_num=*(num+i);
		for(j=0;j<*(num+i);j++)
		{
			F->oct_imgFeature[i].p_Feature[j].Feature_Vector=cvCreateMat(1,128,CV_64FC1);

			for (k=0;k<4;k++)
			{
				for (l=0;l<4;l++)
				{
					F->oct_imgFeature[i].p_Feature[j].blocks[k][l].eight_histo=cvCreateMat(1,8,CV_64FC1);
					for (m=0;m<8;m++)
					{
						*(F->oct_imgFeature[i].p_Feature[j].blocks[k][l].eight_histo->data.db+m)=0;
					}
				}
			}
		}

	}
}

//=====================================================================================


//计算模值================================================================================
double Mod(CvMat *NowMat,int lx,int ly)// l
{
	//直接在fixpoint里面去除这样的特征点
	//if (NowMat->rows-x<2||x<2||NowMat->cols-y<2||y<2)
	//{
	//	return 0;
	//}

	double xp1_y=*(NowMat->data.db+(lx+1)*NowMat->step/8+ly);
	double xs1_y=*(NowMat->data.db+(lx-1)*NowMat->step/8+ly);
	double x_yp1=*(NowMat->data.db+(lx)*NowMat->step/8+ly+1);
	double x_ys1=*(NowMat->data.db+(lx)*NowMat->step/8+ly-1);

	double pointmod=sqrt(pow(xp1_y-xs1_y,2)+pow(x_yp1-x_ys1,2));

	return pointmod;
}

//======================================================================================

double gussian_value(double x,double y,double sigma)
{
	double gv=pow(2.71,-(x*x+y*y)/(2*sigma*sigma));
	return gv;
}


//计算梯度方向即旋转角度 返回角度值 360制式
double Grad_sita(CvMat *NowMat,int Ix,int Iy)//角度按照图像坐标系写 减90即为正常坐标
{

	double xp1_y;
	double xs1_y;
	double x_yp1;
	double x_ys1;

	//if (NowMat->rows-x<2||x<2||NowMat->cols-y<2||y<2)
	//{
	//	return 0;
	//}

	xp1_y=*(NowMat->data.db+(Ix+1)*NowMat->step/8+Iy);
	xs1_y=*(NowMat->data.db+(Ix-1)*NowMat->step/8+Iy);
	x_yp1=*(NowMat->data.db+(Ix)*NowMat->step/8+Iy+1);
	x_ys1=*(NowMat->data.db+(Ix)*NowMat->step/8+Iy-1);

	double Grad=(double)(x_yp1-x_ys1)/(double)(xp1_y-xs1_y);

	if (((x_yp1-x_ys1)==0)&&((xp1_y-xs1_y)==0))
	{
		return 0;
	}
	else if((x_yp1-x_ys1)==0)
	{
		if ((xp1_y-xs1_y)>0)
		{
			return 180;
		} 
		else
		{
			return 0;
		}

	}
	else if((xp1_y-xs1_y)==0)
	{
		if ((x_yp1-x_ys1)>0)
		{
			return 270;
		} 
		else
		{
			return 90;
		}

	}
	else
	{
		if (Grad>0)
		{
			if ((xp1_y-xs1_y)>0)
			{
				double sita=atan(Grad)*180/3.14159+180;
				return sita;
			}
			else
			{
				double sita=atan(Grad)*180/3.14159;

				if (sita<0)
				{
					return sita+360.0;
				} 
				else
				{
					return sita;
				}
			}
		}

		if (Grad<0)
		{
			if ((xp1_y-xs1_y)>0)
			{
				double sita=atan(Grad)*180/3.14159+180;
				return sita;
			} 
			else
			{
				double sita=atan(Grad)*180/3.14159;
				if (sita<0)
				{
					return sita+360.0;
				} 
				else
				{
					return sita;
				}
				
			}
		}

	}

}







//创建128维方形特征描述==========================================================================
Feature *histo_mat(delta Delta,CvMat *spoint,int *num,char *name)
{
	//循环变量===
	int i,j,l,k,p;
	//========

	int xoct=Delta.longoct;//delta的层数

	CvMat *main_direction=cvCreateMat(xoct,100000,CV_64FC1);//=========================存放每个特征点的主方向 序号对应为spoint的序号


	//创建圆形邻域===
	int neibor[2][300];//======================储存每个邻域点的x y 调用时需要用中心点加上这些值
	int neibor_num=0;//======================储存邻域内一共有多少点
	for (i=0;i<17;i++)//======================8*2+1个像素点为长和宽的正方形里面寻找半径为8的圆
	{
		for (j=0;j<17;j++)
		{
			int nx=-8+i;
			int ny=-8+j;
			int dx=abs(nx);
			int dy=abs(ny);
			double modxy=sqrt(dx*dx+dy*dy);
			if (modxy<=8)
			{
				if (dx==0||dy==0)
				{
				}
				else
				{
					neibor[0][neibor_num]=nx;
					neibor[1][neibor_num]=ny;
					neibor_num++;

				}
				
     		}

		}
	}
	//===========
	cout<<"点个数"<<neibor_num<<endl;


	//对圆形邻域内的点进行分块 分成36块===
	int neiborbin[36][100];//======================储存每个扇形柱内点的序号，对应的是neibor[2][300]，即要寻找坐标，用序号去在neibor二维数组用对应的序号查找
	int neiborbin_num[36];//====================储存每个扇形柱内点的数量
	for (i=0;i<36;i++)
	{
		neiborbin_num[i]=0;
	}

	for (i=0;i<neibor_num;i++)
	{
		int nx=neibor[0][i];
		int ny=neibor[1][i];
		double k=(double)(ny)/(double)(nx);
		int sitabin;

		if ((ny==0)||(nx==0))//====================中心点不计入
		{
			/*sitabin=0;
			neiborbin[sitabin][neiborbin_num[sitabin]]=i;*/
		}
		else if(ny==0)
		{
			if (nx>0)
			{
				sitabin=0;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			} 
			else
			{
				sitabin=180/10;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			}

		}
		else if(nx==0)
		{
			if (ny>0)
			{
				sitabin=90/10;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			} 
			else
			{
				sitabin=270/10;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			}

		}
		else
		{
			if (k>0)
			{
				if (nx>0)
				{
					double sita=atan(k)*180/3.14159;
					if (sita<0)
					{
						sita=sita+360;
					}
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				}
				else
				{
					double sita=atan(k)*180/3.14159+180;
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				}
			}

			if (k<0)
			{
				if (nx>0)
				{
					double sita=atan(k)*180/3.14159;
					if (sita<0)
					{
						sita=sita+360;
					}
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				} 
				else
				{
					double sita=atan(k)*180/3.14159+180;
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				}
			}

		}


	}

	for(i=0;i<36;i++)
	{
		int x=neiborbin_num[i];
		cout<<x<<" ";
	}

	//==================================================

	for (i=0;i<xoct;i++)
	{

		CvMat *NowMat=Delta.oct[i].inter[0].Gussianimg;//===========================需要处理的图像矩阵
		int point_num=*(num+i);
		for(j=0;j<point_num;j++)
		{
			double doux=*(spoint->data.db+(2*i)*spoint->step/8+j);
			double douy=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
			int x=cvRound(doux);//======================================特征点的坐标x   整形坐标
			int y=cvRound(douy);//======================================特征点的坐标y   整形坐标


			//=======================================计算主方向与80%的辅助方向添加

			double direction[36];//==========================存放角度与模值
			for (l=0;l<36;l++)//============================初始化direction
			{
				direction[l]=0;
			}


			for (l=0;l<neibor_num;l++)
			{
				int dx=neibor[0][l];//==========================邻域中相对主坐标的相对坐标
				int dy=neibor[1][l];

				int nx=x+dx;//=============================围绕主坐标的邻域坐标
				int ny=y+dy;

				double mod;//=============================邻域坐标的模值
				double sita;//==============================邻域坐标的角度梯度

				mod=Mod(NowMat,nx,ny);//gussian_value(dx,dy,0.5);
				sita=Grad_sita(NowMat,nx,ny);

				int sitabin;//==============================存放每个角度在direction里面的位置
				sitabin=(int)sita/10;
				
				direction[sitabin]+=mod;
			}

			int topbin=0;//================================存放最大值的角度柱
			double topmod=direction[0];//=============================存放最大值的总模值
			int secbin=0;//================================存放辅助值的角度柱
			double secmod=direction[0];//=============================存放辅助值的总模值

			for (l=0;l<36;l++)//=============================找主方向
			{
				if (direction[l]>topmod)
				{
					topmod=direction[l];
					topbin=l;
				}

			}

			//for (l=0;l<36;l++)//==============================找辅助方向
			//{
			//	if (l==topbin)
			//	{
			//	}
			//	else if (direction[l]>secmod)
			//	{
			//		secmod=direction[l];
			//		secbin=l;
			//	}

			//}
		
			//if ((double)(secmod)>=(double)(0.8*topmod))
			//{
			//	*(spoint->data.db+(2*i)*spoint->step/8+*(num+i))=doux;
			//	*(spoint->data.db+(2*i+1)*spoint->step/8+*(num+i))=douy;
			//	*(main_direction->data.db+i*main_direction->step/8+*(num+i))=secbin*10;
			//	*(num+i)=*(num+i)+1;
			//	
			//}


			*(main_direction->data.db+i*main_direction->step/8+j)=topbin*10;


			//==============================================================主方向梯度计算完成

		}
	}















#endif