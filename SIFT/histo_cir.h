#ifndef histo_cir
#define histo_cir

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
#include "his.h"
using namespace cv;
using namespace std;

//建立特征点数据对象================================================================
////
////class histo
////{
////public:
////	CvMat *eight_histo;
////};
////class point_Feature
////{
////public:
////	histo blocks[9];
////	CvMat *Feature_Vector;
////};
////class img_Feature
////{
////public:
////	int Feature_num;
////	point_Feature p_Feature[10000];
////};
////class Feature
////{
////public:
////	int oct;
////	img_Feature oct_imgFeature[10];
////};
//
//
////====================================================================================
//
////初始化特征类============================================================================
//void  init_Feature(Feature *F,int xoct,int *num)
//{
//	int i,j,k,m;
//
//
//	F->oct=xoct;
//	for (i=0;i<xoct;i++)
//	{
//		F->oct_imgFeature[i].Feature_num=*(num+i);
//		for(j=0;j<*(num+i);j++)
//		{
//			F->oct_imgFeature[i].p_Feature[j].Feature_Vector=cvCreateMat(1,72,CV_64FC1);
//
//			for (k=0;k<9;k++)
//			{
//				F->oct_imgFeature[i].p_Feature[j].blocks[k].eight_histo=cvCreateMat(1,8,CV_64FC1);
//				for (m=0;m<8;m++)
//				{
//					*(F->oct_imgFeature[i].p_Feature[j].blocks[k].eight_histo->data.db+m)=0;
//				}
//				
//			}
//		}
//
//	}
//}
//#define ImMat(ROW,COL) ((float *)(im->data.fl + im->step/sizeof(float) *(ROW)))[(COL)]
//float getPixelBI(CvMat * im, float col, float row) 
//{
//	int irow, icol;
//	float rfrac, cfrac;
//	float row1 = 0, row2 = 0;
//	int width=im->cols;
//	int height=im->rows;
//
//
//	irow = (int) row;
//	icol = (int) col;
//
//	if (irow < 0 || irow >= height|| icol < 0 || icol >= width)
//		return 0;
//	if (row > height - 1)
//		row = height - 1;
//	if (col > width - 1)
//		col = width - 1;
//	rfrac = 1.0 - (row - (float) irow);
//	cfrac = 1.0 - (col - (float) icol);
//	if (cfrac < 1) 
//	{
//		row1 = cfrac * ImMat(irow,icol) + (1.0 - cfrac) * ImMat(irow,icol+1);
//	} 
//	else 
//	{
//		row1 = ImMat(irow,icol);
//	}
//	if (rfrac < 1) 
//	{
//		if (cfrac < 1) 
//		{
//			row2 = cfrac * ImMat(irow+1,icol) + (1.0 - cfrac) * ImMat(irow+1,icol+1);
//		} else 
//		{
//			row2 = ImMat(irow+1,icol);
//		}
//	}
//	return rfrac * row1 + (1.0 - rfrac) * row2;
//}




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



//初始化类Feature============
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






//
////创建72维扇形特征描述==========================================================================
//Feature *histocir(delta Delta,CvMat *spoint,int *num,char *name)
//{
//	//循环变量===
//	int i,j,l,k,p;
//	//========
//
//	int xoct=Delta.longoct;//delta的层数
//
//	CvMat *main_direction=cvCreateMat(xoct,100000,CV_64FC1);//=========================存放每个特征点的主方向 序号对应为spoint的序号
//
//
//	//创建圆形邻域===
//	int neibor[2][300];//======================储存每个邻域点的x y 调用时需要用中心点加上这些值
//	int neibor_num=0;//======================储存邻域内一共有多少点
//	for (i=0;i<17;i++)//======================8*2+1个像素点为长和宽的正方形里面寻找半径为8的圆
//	{
//		for (j=0;j<17;j++)
//		{
//			int nx=-8+i;
//			int ny=-8+j;
//			int dx=abs(nx);
//			int dy=abs(ny);
//			double modxy=sqrt(dx*dx+dy*dy);
//			if (modxy<=8)
//			{
//				if (dx==0||dy==0)
//				{
//				}
//				else
//				{
//					neibor[0][neibor_num]=nx;
//					neibor[1][neibor_num]=ny;
//					neibor_num++;
//
//				}
//				
//     		}
//
//		}
//	}
//	//===========
//	cout<<"点个数"<<neibor_num<<endl;
//
//
//	//对圆形邻域内的点进行分块 分成36块===
//	int neiborbin[36][100];//======================储存每个扇形柱内点的序号，对应的是neibor[2][300]，即要寻找坐标，用序号去在neibor二维数组用对应的序号查找
//	int neiborbin_num[36];//====================储存每个扇形柱内点的数量
//	for (i=0;i<36;i++)
//	{
//		neiborbin_num[i]=0;
//	}
//
//	for (i=0;i<neibor_num;i++)
//	{
//		int nx=neibor[0][i];
//		int ny=neibor[1][i];
//		double k=(double)(ny)/(double)(nx);
//		int sitabin;
//
//		if ((ny==0)||(nx==0))//====================中心点不计入
//		{
//			/*sitabin=0;
//			neiborbin[sitabin][neiborbin_num[sitabin]]=i;*/
//		}
//		else if(ny==0)
//		{
//			if (nx>0)
//			{
//				sitabin=0;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			} 
//			else
//			{
//				sitabin=180/10;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			}
//
//		}
//		else if(nx==0)
//		{
//			if (ny>0)
//			{
//				sitabin=90/10;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			} 
//			else
//			{
//				sitabin=270/10;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			}
//
//		}
//		else
//		{
//			if (k>0)
//			{
//				if (nx>0)
//				{
//					double sita=atan(k)*180/3.14159;
//					if (sita<0)
//					{
//						sita=sita+360;
//					}
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				}
//				else
//				{
//					double sita=atan(k)*180/3.14159+180;
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				}
//			}
//
//			if (k<0)
//			{
//				if (nx>0)
//				{
//					double sita=atan(k)*180/3.14159;
//					if (sita<0)
//					{
//						sita=sita+360;
//					}
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				} 
//				else
//				{
//					double sita=atan(k)*180/3.14159+180;
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				}
//			}
//
//		}
//
//
//	}
//
//	for(i=0;i<36;i++)
//	{
//		int x=neiborbin_num[i];
//		cout<<x<<" ";
//	}
//
//	//==================================================
//
//	for (i=0;i<xoct;i++)
//	{
//
//		CvMat *NowMat=Delta.oct[i].inter[0].Gussianimg;//===========================需要处理的图像矩阵
//		int point_num=*(num+i);
//		for(j=0;j<point_num;j++)
//		{
//			double doux=*(spoint->data.db+(2*i)*spoint->step/8+j);
//			double douy=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
//			int x=cvRound(doux);//======================================特征点的坐标x   整形坐标
//			int y=cvRound(douy);//======================================特征点的坐标y   整形坐标
//
//
//			//=======================================计算主方向与80%的辅助方向添加
//
//			double direction[36];//==========================存放角度与模值
//			for (l=0;l<36;l++)//============================初始化direction
//			{
//				direction[l]=0;
//			}
//
//
//			for (l=0;l<neibor_num;l++)
//			{
//				int dx=neibor[0][l];//==========================邻域中相对主坐标的相对坐标
//				int dy=neibor[1][l];
//
//				int nx=x+dx;//=============================围绕主坐标的邻域坐标
//				int ny=y+dy;
//
//				double mod;//=============================邻域坐标的模值
//				double sita;//==============================邻域坐标的角度梯度
//
//				mod=Mod(NowMat,nx,ny);//gussian_value(dx,dy,0.5);
//				sita=Grad_sita(NowMat,nx,ny);
//
//				int sitabin;//==============================存放每个角度在direction里面的位置
//				sitabin=(int)sita/10;
//				
//				direction[sitabin]+=mod;
//			}
//
//			int topbin=0;//================================存放最大值的角度柱
//			double topmod=direction[0];//=============================存放最大值的总模值
//			int secbin=0;//================================存放辅助值的角度柱
//			double secmod=direction[0];//=============================存放辅助值的总模值
//
//			for (l=0;l<36;l++)//=============================找主方向
//			{
//				if (direction[l]>topmod)
//				{
//					topmod=direction[l];
//					topbin=l;
//				}
//
//			}
//
//			//for (l=0;l<36;l++)//==============================找辅助方向
//			//{
//			//	if (l==topbin)
//			//	{
//			//	}
//			//	else if (direction[l]>secmod)
//			//	{
//			//		secmod=direction[l];
//			//		secbin=l;
//			//	}
//
//			//}
//		
//			//if ((double)(secmod)>=(double)(0.8*topmod))
//			//{
//			//	*(spoint->data.db+(2*i)*spoint->step/8+*(num+i))=doux;
//			//	*(spoint->data.db+(2*i+1)*spoint->step/8+*(num+i))=douy;
//			//	*(main_direction->data.db+i*main_direction->step/8+*(num+i))=secbin*10;
//			//	*(num+i)=*(num+i)+1;
//			//	
//			//}
//
//
//			*(main_direction->data.db+i*main_direction->step/8+j)=topbin*10;
//
//
//			//==============================================================主方向梯度计算完成
//
//		}
//	}
//	
//	//==============================================================计算所有邻域点的72维直方图
//	//==============================================================x轴对齐梯度方向
//	//==============================================================所计算出邻域点的梯度角度也要和相应的主方向进行相对化
//
//	Feature *F;
//	F=(Feature*)malloc(sizeof(Feature));
//	init_Feature(F,xoct,num);//===============================================初始化特征值的类
//
//	for (i=0;i<xoct;i++)
//	{
//
//		CvMat *NowMat=Delta.oct[i].inter[2].Gussianimg;//===========================需要处理的图像矩阵
//		int point_num=*(num+i);
//		for(j=0;j<point_num;j++)
//		{
//			double doux=*(spoint->data.db+(2*i)*spoint->step/8+j);
//			double douy=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
//			int x=cvRound(doux);//======================================特征点的坐标x   整形坐标
//			int y=cvRound(douy);//======================================特征点的坐标y   整形坐标
//
//			double direction=*(main_direction->data.db+i*main_direction->step/8+j);//===============================提取特征点主方向
//
//			int neibor_point[9][100];//=====================================存放按主方向为x轴的划区内邻域点的序号
//			int neibor_pointnum[9];//======================================存放每个区域内邻域点的数量
//
//			for (p=0;p<9;p++)//========================================初始化neibor_pointnum
//			{
//				neibor_pointnum[p]=0;
//			}
//
//			for (p=0;p<9;p++)//=========================================将从主方向划分九个区块 将邻域内点的序号放入9个区块中
//			{
//				for (l=p*4;l<p*4+4;l++)
//				{
//					for (k=0;k<neiborbin_num[((int)direction/10+l)%36];k++)//==================有可能加过 需要取模运算
//					{
//						//int tew=((int)direction/10+l)%36;
//						neibor_point[p][neibor_pointnum[p]]=neiborbin[((int)direction/10+l)%36][k];
//						neibor_pointnum[p]++;
//					}
//
//				}
//			}
//			
//			//for (p=0;p<9;p++)
//			//{
//			//	cout<<neibor_pointnum[p]<<endl;
//			//}
//
//
//			for (p=0;p<9;p++)
//			{
//				for (l=0;l<neibor_pointnum[p];l++)
//				{
//					int number=neibor_point[p][l];
//					int dx=neibor[0][number];
//					int dy=neibor[1][number];
//
//					int exact_x=x+dx;//==============================图像中邻域的真实点
//					int exact_y=y+dy;
//
//					double exact_sita=Grad_sita(NowMat,exact_x,exact_y);//=========================================绝对的方向 也就是图像中直接像素点的方向
//					double exact_mod=Mod(NowMat,exact_x,exact_y);//===========================================绝对模值
//					double maind=*(main_direction->data.db+i*main_direction->step/8+j);
//
//					double ralative_sita=(int)(exact_sita-*(main_direction->data.db+i*main_direction->step/8+j)+360.0)%360;//=========相对方向 也就是图像中的邻域点的梯度方向-特征点主方向 保证了旋转不变性
//					int eight_direction=ralative_sita/45;//=====================================================相对方向在八方图中的表示
//					double addmod=exact_mod*pow(2.71,-(dx*dx+dy*dy)/(2.0*0.5*0.5));
//	/*				if (addmod>0.2)
//					{
//						addmod=0.2;
//					}*/
//			/*		if (addmod<0.001)
//					{
//						addmod=0;
//					}*/
//
//					*(F->oct_imgFeature[i].p_Feature[j].blocks[p].eight_histo->data.db+eight_direction)+=addmod;
//					
//
//				}
//			}
//
//
//			for (p=0;p<9;p++)
//			{
//				for (l=0;l<8;l++)
//				{
//					*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+8*p+l)=*(F->oct_imgFeature[i].p_Feature[j].blocks[p].eight_histo->data.db+l);//将9*8放到72维向量
//				}
//			}
//
//			//======================================================向量归一化
//			double add2=0;
//			for (p=0;p<72;p++)
//			{
//				add2+=(*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+p))*(*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+p));
//			}
//			add2=sqrt(add2);
//			for (p=0;p<72;p++)
//			{
//				double te=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+p);
//				*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+p)=te/add2;
//			}
//
//
//
//		}
//	}
//
//	char file[100];
//	sprintf_s(file,"G:\\opencv\\1\\test10\\SIFT\\SIFT\\%s.txt",name);
//	
//	FILE *fl;
//	fl=fopen(file,"w");
//	
////	of.open("G:\\opencv\\1\\test10\\SIFT\\SIFT\\img1.txt",ios::app|ios::out);
//
//
//	int point_num=*(num+0);
//	for(j=0;j<point_num;j++)
//	{
//		fprintf(fl,"特征点%d     %f\n",j,*(main_direction->data.db+j));
////		of<<"特征点"<<j<<endl;
//		CvMat *e=F->oct_imgFeature[0].p_Feature[j].Feature_Vector;
//		for (p=0;p<9;p++)
//		{
//			for (l=0;l<8;l++)
//			{
//				//of<<*(e->data.db+8*p+l);
//				fprintf(fl,"%10f   ",*(e->data.db+8*p+l));
//			}
//			//of<<endl;
//			fprintf(fl,"\n");
//		}
//		//of<<endl;
//		fprintf(fl,"\n");
//	}
////	of.close();
//
//
//	Feature *F1;
//
//
//
//
//
//	return F;
//
//
//}


Feature *histmat(delta Delta,CvMat *spoint,int *num,char *name)
{

	int xoct=Delta.longoct;//delta的层数
	Feature *F;
	F=(Feature*)malloc(sizeof(Feature));
	init_Feature(F,xoct,num);

	//循环变量===
	int i,j,l,k,p;
	//========

	

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
	
	//int neicir[2][300];
	//int neicir_num=0;
	//int R=8;
	//int bigmat_R=R*sqrt(2)+1;
	//for (i=0;i<2*bigmat_R+1;i++)//======================8*2+1个像素点为长和宽的正方形里面寻找半径为8*sqrt(2)的圆
	//{
	//	for (j=0;j<2*bigmat_R+1;j++)
	//	{
	//		int nx=-bigmat_R+i;
	//		int ny=-bigmat_R+j;
	//		int dx=abs(nx);
	//		int dy=abs(ny);
	//		double modxy=sqrt(dx*dx+dy*dy);
	//		if (modxy<=bigmat_R)
	//		{
	//			if (dx==0||dy==0)
	//			{
	//			}
	//			else
	//			{
	//				neicir[0][neicir_num]=nx;
	//				neicir[1][neicir_num]=ny;
	//				neicir_num++;
	//			}
	//		}
	//	}
	//}
	//cout<<neicir_num<<endl;
	//for (i=0;i<xoct;i++)
	//{
	//	CvMat *NowMat=Delta.oct[i].inter[0].Gussianimg;//===========================需要处理的图像矩阵
	//	int point_num=*(num+i);
	//	for(j=0;j<point_num;j++)
	//	{
	//	}
	//}

	for (i=0;i<xoct;i++)
	{

		CvMat *NowMat=Delta.oct[i].inter[0].Gussianimg;//===========================需要处理的图像矩阵
		IplImage *ax=cvCreateImage(cvGetSize(NowMat),IPL_DEPTH_8U,1);
		cvConvert(NowMat,ax);
		int point_num=*(num+i);
		for(j=0;j<point_num;j++)
		{

			/*Mat Matnow=cvarrToMat(ax,true);*/
			float ft[128];
			double x=*(spoint->data.db+2*i*spoint->step/8+j);
			double y=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
			Point2f p(y,x);
			float direct=*(main_direction->data.db+i*main_direction->step/8+j);
			
			//cvNamedWindow("a");
			//cvShowImage("a",ax);
			//cvWaitKey(0);
			/*namedWindow("as");
			imshow("as",Matnow);
			waitKey(0);*/
			calcSIFTDescriptor(NowMat,p,direct,ft);

			for (k=0;k<128;k++)
			{
				CvMat *z=F->oct_imgFeature[i].p_Feature[j].Feature_Vector;
				*(z->data.db+k)=ft[k];
			}
			


			}
	}

	//char file[100];
	//sprintf_s(file,"G:\\opencv\\1\\test10\\SIFT\\SIFT\\%s.txt",name);

	//FILE *fl;
	//fl=fopen(file,"w");

	//	int point_num=*(num+0);
	//	for(j=0;j<point_num;j++)
	//	{
	//		fprintf(fl,"特征点%d     %f\n",j,*(main_direction->data.db+j));
	////		of<<"特征点"<<j<<endl;
	//		CvMat *e=F->oct_imgFeature[0].p_Feature[j].Feature_Vector;
	//		for (p=0;p<9;p++)
	//		{
	//			for (l=0;l<8;l++)
	//			{
	//				//of<<*(e->data.db+8*p+l);
	//				fprintf(fl,"%10f   ",*(e->data.db+8*p+l));
	//			}
	//			//of<<endl;
	//			fprintf(fl,"\n");
	//		}
	//		//of<<endl;
	//		fprintf(fl,"\n");
	//	}
	////	of.close();
	////

	//cout<<endl;








	return F;

}







//
//
//
//CvMat *hi(delta Delta,CvMat *spoint,int *num)
//{
//	//循环变量===
//	int i,j,l,k,p;
//	//========
//
//	int xoct=Delta.longoct;//delta的层数
//
//	CvMat *main_direction=cvCreateMat(xoct,100000,CV_64FC1);//=========================存放每个特征点的主方向 序号对应为spoint的序号
//
//
//	//创建圆形邻域===
//	int neibor[2][300];//======================储存每个邻域点的x y 调用时需要用中心点加上这些值
//	int neibor_num=0;//======================储存邻域内一共有多少点
//	for (i=0;i<17;i++)//======================8*2+1个像素点为长和宽的正方形里面寻找半径为8的圆
//	{
//		for (j=0;j<17;j++)
//		{
//			int nx=-8+i;
//			int ny=-8+j;
//			int dx=abs(nx);
//			int dy=abs(ny);
//			double modxy=sqrt(dx*dx+dy*dy);
//			if (modxy<=8)
//			{
//				neibor[0][neibor_num]=nx;
//				neibor[1][neibor_num]=ny;
//				neibor_num++;
//     		}
//
//		}
//	}
//	//===========
//
//
//
//	//对圆形邻域内的点进行分块 分成36块===
//	int neiborbin[36][100];//======================储存每个扇形柱内点的序号，对应的是neibor[2][300]，即要寻找坐标，用序号去在neibor二维数组用对应的序号查找
//	int neiborbin_num[36];//====================储存每个扇形柱内点的数量
//	for (i=0;i<36;i++)
//	{
//		neiborbin_num[i]=0;
//	}
//
//	for (i=0;i<neibor_num;i++)
//	{
//		int nx=neibor[0][i];
//		int ny=neibor[1][i];
//		double k=(double)(ny)/(double)(nx);
//		int sitabin;
//
//		if ((ny==0)&&(nx==0))//====================中心点不计入
//		{
//			/*sitabin=0;
//			neiborbin[sitabin][neiborbin_num[sitabin]]=i;*/
//		}
//		else if(ny==0)
//		{
//			if (nx>0)
//			{
//				sitabin=0;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			} 
//			else
//			{
//				sitabin=180/10;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			}
//
//		}
//		else if(nx==0)
//		{
//			if (ny>0)
//			{
//				sitabin=90/10;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			} 
//			else
//			{
//				sitabin=270/10;
//				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//				neiborbin_num[sitabin]++;
//			}
//
//		}
//		else
//		{
//			if (k>0)
//			{
//				if (nx>0)
//				{
//					double sita=atan(k)*180/3.14159;
//					if (sita<0)
//					{
//						sita=sita+360;
//					}
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				}
//				else
//				{
//					double sita=atan(k)*180/3.14159+180;
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				}
//			}
//
//			if (k<0)
//			{
//				if (nx>0)
//				{
//					double sita=atan(k)*180/3.14159;
//					if (sita<0)
//					{
//						sita=sita+360;
//					}
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				} 
//				else
//				{
//					double sita=atan(k)*180/3.14159+180;
//					sitabin=(int)(sita/10);
//					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
//					neiborbin_num[sitabin]++;
//				}
//			}
//
//		}
//
//
//	}
//
//	for(i=0;i<36;i++)
//	{
//		int x=neiborbin_num[i];
//		cout<<x<<" ";
//	}
//
//	//==================================================
//
//	for (i=0;i<xoct;i++)
//	{
//
//		CvMat *NowMat=Delta.oct[i].inter[0].Gussianimg;//===========================需要处理的图像矩阵
//		int point_num=*(num+i);
//		for(j=0;j<point_num;j++)
//		{
//			double doux=*(spoint->data.db+(2*i)*spoint->step/8+j);
//			double douy=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
//			int x=cvRound(doux);//======================================特征点的坐标x   整形坐标
//			int y=cvRound(douy);//======================================特征点的坐标y   整形坐标
//
//
//			//=======================================计算主方向与80%的辅助方向添加
//
//			double direction[36];//==========================存放角度与模值
//			for (l=0;l<36;l++)//============================初始化direction
//			{
//				direction[l]=0;
//			}
//
//
//			for (l=0;l<neibor_num;l++)
//			{
//				int dx=neibor[0][l];//==========================邻域中相对主坐标的相对坐标
//				int dy=neibor[1][l];
//
//				int nx=x+dx;//=============================围绕主坐标的邻域坐标
//				int ny=y+dy;
//
//				double mod;//=============================邻域坐标的模值
//				double sita;//==============================邻域坐标的角度梯度
//
//				mod=Mod(NowMat,nx,ny);//gussian_value(dx,dy,0.5);
//				sita=Grad_sita(NowMat,nx,ny);
//
//				int sitabin;//==============================存放每个角度在direction里面的位置
//				sitabin=(int)sita/10;
//				
//				direction[sitabin]+=mod;
//			}
//
//			int topbin=0;//================================存放最大值的角度柱
//			double topmod=direction[0];//=============================存放最大值的总模值
//			int secbin=0;//================================存放辅助值的角度柱
//			double secmod=direction[0];//=============================存放辅助值的总模值
//
//			for (l=0;l<36;l++)//=============================找主方向
//			{
//				if (direction[l]>topmod)
//				{
//					topmod=direction[l];
//					topbin=l;
//				}
//
//			}
//
//			//for (l=0;l<36;l++)//==============================找辅助方向
//			//{
//			//	if (l==topbin)
//			//	{
//			//	}
//			//	else if (direction[l]>secmod)
//			//	{
//			//		secmod=direction[l];
//			//		secbin=l;
//			//	}
//
//			//}
//		
//			//if ((double)(secmod)>=(double)(0.8*topmod))
//			//{
//			//	*(spoint->data.db+(2*i)*spoint->step/8+*(num+i))=doux;
//			//	*(spoint->data.db+(2*i+1)*spoint->step/8+*(num+i))=douy;
//			//	*(main_direction->data.db+i*main_direction->step/8+*(num+i))=secbin*10;
//			//	*(num+i)=*(num+i)+1;
//			//	
//			//}
//
//
//			*(main_direction->data.db+i*main_direction->step/8+j)=topbin*10;
//
//
//			//==============================================================主方向梯度计算完成
//
//		}
//	}
//
//
//
//
//
//
//
//	return main_direction;
//}
//	
#endif