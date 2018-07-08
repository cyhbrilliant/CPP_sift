#ifndef histogram
#define histogram
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"
#include "GUSSIANdelta.h"
#include <stdlib.h>
#include <math.h>
using namespace cv;
using namespace std;

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


//计算模值=====================
double Mod(CvMat *NowMat,double x,double y)
{
	int Ix=cvRound(x);
	int Iy=cvRound(y);    //四舍五入到整数点
	if (NowMat->rows-x<2||x<2||NowMat->cols-y<2||y<2)
	{
		return 0;
	}

	double xp1_y=*(NowMat->data.db+(Ix+1)*NowMat->step/8+Iy);
	double xs1_y=*(NowMat->data.db+(Ix-1)*NowMat->step/8+Iy);
	double x_yp1=*(NowMat->data.db+(Ix)*NowMat->step/8+Iy+1);
	double x_ys1=*(NowMat->data.db+(Ix)*NowMat->step/8+Iy-1);

	double pointmod=sqrt(pow(xp1_y-xs1_y,2)+pow(x_yp1-x_ys1,2));

	return pointmod;
}

double xp1_y;
double xs1_y;
double x_yp1;
double x_ys1;


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

//用X,Y,SIGMA三个参数计算高斯权重==============
double Gaussian_value(int i,int j,double sigma)
{
	double gaussianv=pow(2.72,-(i*i+j*j)/(2*sigma*sigma));
	return gaussianv;
}



//建立直方图=====================
CvMat* BuildHistogram(delta Delta,CvMat *spoint,int *num,int R_BasicDrection)
{
	int xoct=Delta.longoct;
	int i,j,k,l;
	double Top_num=0;
	int top_k=0;
	double *sigma=(double *)malloc(sizeof(double)*xoct);
	for(i=0;i<xoct;i++)
	{
		*(sigma+i)=pow(2,i)*1.6;
	}

	
	int NUM_R_BasicDrection=R_BasicDrection*R_BasicDrection;


	CvMat *Point_GradMod=cvCreateMat(spoint->rows,spoint->cols*NUM_R_BasicDrection,CV_64FC1);  //存放每个特征点以及邻域点的基本特征方向矩阵信息
	CvMat *Spoint_GradMod=cvCreateMat(spoint->rows,spoint->cols,CV_64FC1); //若有辅助方向 新建的numGM也要加上数量也要相应加上数量
	double thirtysix_bin[36];  //存放36柱的直方图

	//bfR_XY未进行基础方向旋转的坐标
	CvMat* bfR_XY=cvCreateMat(2,10000,CV_64FC1);
	int num_abfR_XY=0;

	//afR_XY进行基础方向旋转后的坐标
	CvMat* afR_XY=cvCreateMat(2,10000,CV_64FC1);


	for (i=0;i<xoct;i++)
	{
		for (j=0;j<*(num+i);j++)
		{

			/*R_BasicDrection=3*(int)*(sigma+i)*1.5;
			NUM_R_BasicDrection=R_BasicDrection*R_BasicDrection;*/

			Top_num=0;
			top_k=0;
			for(k=0;k<36;k++)
			{
				thirtysix_bin[k]=0;
			}

			double x=*(spoint->data.db+2*i*spoint->step/8+j);
			double y=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
			int Intx=cvRound(x);
			int Inty=cvRound(y); 

			CvMat *NowMat=Delta.oct[i].inter[2].Gussianimg;
			

			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k-R_BasicDrection/2,Inty+l-R_BasicDrection/2);
					//cout<<*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)<<"  "; 
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k-R_BasicDrection/2,Inty+l-R_BasicDrection/2)*Gaussian_value(k-R_BasicDrection/2,l-R_BasicDrection/2,*(sigma+i)*1.5);
				}
			}




			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+1*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k-R_BasicDrection/2,Inty+1+l);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+1*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k-R_BasicDrection/2,Inty+1+l)*Gaussian_value(k-R_BasicDrection/2,1+l,*(sigma+i)*1.5);
				}
			}


			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+2*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k+1,Inty+l-R_BasicDrection/2);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+2*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k+1,Inty+l-R_BasicDrection/2)*Gaussian_value(k+1,l-R_BasicDrection/2,*(sigma+i)*1.5);
				}
			}



			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+3*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k+1,Inty+l+1);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+3*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k+1,Inty+l+1)*Gaussian_value(k+1,l+1,*(sigma+i)*1.5);
				}
			}

			//建立初始特征点方向直方图
			for (k=0;k<NUM_R_BasicDrection;k++)
			{
				double now_sita=*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+k);
				double now_mod=*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+k);

				if (now_sita<0)
				{
					thirtysix_bin[(int)(now_sita/10)+36]+=now_mod;
				} 
				else if(now_sita>0)
				{
					thirtysix_bin[(int)(now_sita/10)]+=now_mod;
				}
				
			}


			//对直方图进行高斯平滑
			thirtysix_bin[0]=(thirtysix_bin[0]*6+thirtysix_bin[1]*4+thirtysix_bin[2])/11;
			thirtysix_bin[1]=(thirtysix_bin[0]*4+thirtysix_bin[1]*6+thirtysix_bin[2]*4+thirtysix_bin[3])/15;
			for (k=2;k<34;k++)
			{
				thirtysix_bin[k]=(thirtysix_bin[k-2]+thirtysix_bin[k-1]*4+thirtysix_bin[k]*6+thirtysix_bin[k+1]*4+thirtysix_bin[k+2])/16;
			}
			thirtysix_bin[34]=(thirtysix_bin[35]*4+thirtysix_bin[34]*6+thirtysix_bin[33]*4+thirtysix_bin[32])/15;
			thirtysix_bin[35]=(thirtysix_bin[35]*6+thirtysix_bin[34]*4+thirtysix_bin[33])/11;


			


			Top_num=thirtysix_bin[0];
			for (k=0;k<36;k++)
			{
				if (thirtysix_bin[k]>Top_num)
				{
					Top_num=thirtysix_bin[k];
					top_k=k;
				}
			}

			*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j)=top_k*10;
			*(Spoint_GradMod->data.db+(2*i+1)*Spoint_GradMod->step/8+j)=Top_num;

			}

		}
			

	return Spoint_GradMod;

	
}


//128维特征向量=========================
Feature* calFeature(delta Delta,CvMat *spoint,int *num,int R_BasicDrection,int R_FeatureGRAD)
{
	int xoct=Delta.longoct;
	int i,j,k,l;
	double Top_num=0;
	int top_k=0;
	double *sigma=(double *)malloc(sizeof(double)*xoct);
	for(i=0;i<xoct;i++)
	{
		*(sigma+i)=pow(2,i)*1.6;
	}
	//================================================================================================================
	//初步测试 使用8*8探测基本特征点方向 16*16探测梯度

	
	int NUM_R_BasicDrection=R_BasicDrection*R_BasicDrection;
	int NUM_R_FeatureGRAD=R_FeatureGRAD*R_FeatureGRAD;

	CvMat *Point_GradMod=cvCreateMat(spoint->rows,spoint->cols*NUM_R_BasicDrection,CV_64FC1);  //存放每个特征点以及邻域点的基本特征方向矩阵信息
	CvMat *Spoint_GradMod=cvCreateMat(spoint->rows,spoint->cols,CV_64FC1); //若有辅助方向 新建的numGM也要加上数量也要相应加上数量
	double thirtysix_bin[36];  //存放36柱的直方图

	//bfR_XY未进行基础方向旋转的坐标
	CvMat* bfR_XY=cvCreateMat(2,1000,CV_64FC1);
	int num_abfR_XY=0;

	//afR_XY进行基础方向旋转后的坐标
	CvMat* afR_XY=cvCreateMat(2,1000,CV_64FC1);


	//计算领域所需要的半径点坐标
	double R=(1.4*R_FeatureGRAD+1)/2;
	cout<<"检测半径:"<<R<<endl;

	for (k=-R;k<=R;k++)
	{
		for (l=-R;l<=R;l++)
		{
			if (sqrt(pow(k,2)+pow(l,2))<R)
			{
				if (k!=0&&l!=0)
				{
					*(bfR_XY->data.db+0*bfR_XY->step/8+num_abfR_XY)=k;
					*(bfR_XY->data.db+1*bfR_XY->step/8+num_abfR_XY)=l;
					//cout<<"x:"<<k<<" "<<"y:"<<l<<endl;
					num_abfR_XY++; 
				}

			}

		}
	}

	cout<<"领域个数"<<num_abfR_XY<<endl;



	//储存每个特征点的特征向量 也就是最终的SIFT特征向量
	//每一层有两行 第一行为角度 第二行为幅值
	//CvMat *Feature=cvCreateMat(spoint->rows,128*spoint->cols,CV_64FC1);



	//================初始化======================================
	Feature *F;
	F=(Feature*)malloc(sizeof(Feature));
	init_Feature(F,xoct,num);
	//============================================================



	for (i=0;i<xoct;i++)
	{
		for (j=0;j<*(num+i);j++)
		{

			/*R_BasicDrection=3*(int)*(sigma+i)*1.5;
			NUM_R_BasicDrection=R_BasicDrection*R_BasicDrection;*/

			Top_num=0;
			top_k=0;
			for(k=0;k<36;k++)
			{
				thirtysix_bin[k]=0;
			}

			double x=*(spoint->data.db+2*i*spoint->step/8+j);
			double y=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
			int Intx=cvRound(x);
			int Inty=cvRound(y); 

			CvMat *NowMat=Delta.oct[i].inter[1].Gussianimg;
			

			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k-R_BasicDrection/2,Inty+l-R_BasicDrection/2);
					//cout<<*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)<<"  "; 
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k-R_BasicDrection/2,Inty+l-R_BasicDrection/2)*Gaussian_value(k-R_BasicDrection/2,l-R_BasicDrection/2,*(sigma+i)*1.5);
				}
			}




			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+1*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k-R_BasicDrection/2,Inty+1+l);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+1*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k-R_BasicDrection/2,Inty+1+l)*Gaussian_value(k-R_BasicDrection/2,1+l,*(sigma+i)*1.5);
				}
			}


			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+2*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k+1,Inty+l-R_BasicDrection/2);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+2*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k+1,Inty+l-R_BasicDrection/2)*Gaussian_value(k+1,l-R_BasicDrection/2,*(sigma+i)*1.5);
				}
			}



			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+3*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k+1,Inty+l+1);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+3*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k+1,Inty+l+1)*Gaussian_value(k+1,l+1,*(sigma+i)*1.5);
				}
			}

			//建立初始特征点方向直方图
			for (k=0;k<NUM_R_BasicDrection;k++)
			{
				double now_sita=*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+k);
				double now_mod=*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+k);

				if (now_sita<0)
				{
					thirtysix_bin[(int)(now_sita/10)+36]+=now_mod;
				} 
				else if(now_sita>0)
				{
					thirtysix_bin[(int)(now_sita/10)]+=now_mod;
				}
				
			}


			//对直方图进行高斯平滑
			thirtysix_bin[0]=(thirtysix_bin[0]*6+thirtysix_bin[1]*4+thirtysix_bin[2])/11;
			thirtysix_bin[1]=(thirtysix_bin[0]*4+thirtysix_bin[1]*6+thirtysix_bin[2]*4+thirtysix_bin[3])/15;
			for (k=2;k<34;k++)
			{
				thirtysix_bin[k]=(thirtysix_bin[k-2]+thirtysix_bin[k-1]*4+thirtysix_bin[k]*6+thirtysix_bin[k+1]*4+thirtysix_bin[k+2])/16;
			}
			thirtysix_bin[34]=(thirtysix_bin[35]*4+thirtysix_bin[34]*6+thirtysix_bin[33]*4+thirtysix_bin[32])/15;
			thirtysix_bin[35]=(thirtysix_bin[35]*6+thirtysix_bin[34]*4+thirtysix_bin[33])/11;


			


			Top_num=thirtysix_bin[0];
			for (k=0;k<36;k++)
			{
				if (thirtysix_bin[k]>Top_num)
				{
					Top_num=thirtysix_bin[k];
					top_k=k;
				}
			}

			*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j)=top_k*10;
			*(Spoint_GradMod->data.db+(2*i+1)*Spoint_GradMod->step/8+j)=Top_num;

			//printf("角度%10lf 幅值%10lf",*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j),*(Spoint_GradMod->data.db+(2*i+1)*Spoint_GradMod->step/8+j));
			//printf("角度%10lf\n",*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j));

			

			//对每个种子点领域像素进行坐标方向旋转
			CvMat *grad_point=cvCreateMat(2,num_abfR_XY,CV_64FC1);//sita mod
			CvMat *point_block=cvCreateMat(2,num_abfR_XY,CV_64FC1);//x y
			int inblock[1000];//在领域内的点的序号
			int num_inblock=0;

			int q,w;
			for (k=0;k<num_abfR_XY;k++)
			{
				double point_sita=(*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j))*3.14/180;
				double bfx=*(bfR_XY->data.db+0*bfR_XY->step/8+k);
				double bfy=*(bfR_XY->data.db+1*bfR_XY->step/8+k);
				*(afR_XY->data.db+0*afR_XY->step/8+k)=bfx*cos(point_sita)-bfy*sin(point_sita);
				*(afR_XY->data.db+1*afR_XY->step/8+k)=bfx*sin(point_sita)+bfy*cos(point_sita);


				int afX_block,afY_block;

				*(grad_point->data.db+0*grad_point->step/8+k)=Grad_sita(NowMat,Intx+bfx,Inty+bfy);
				*(grad_point->data.db+1*grad_point->step/8+k)=Mod(NowMat,Intx+bfx,Inty+bfy);


				if (*(afR_XY->data.db+0*afR_XY->step/8+k)>0&&*(afR_XY->data.db+1*afR_XY->step/8+k)<0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)((*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1);
					*(point_block->data.db+1*point_block->step/8+k)=(int)((*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1);
				}
				if (*(afR_XY->data.db+0*afR_XY->step/8+k)>0&&*(afR_XY->data.db+1*afR_XY->step/8+k)>0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)((*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1);
					*(point_block->data.db+1*point_block->step/8+k)=(int)((*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1);
				}
				if (*(afR_XY->data.db+0*afR_XY->step/8+k)<0&&*(afR_XY->data.db+1*afR_XY->step/8+k)<0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)((*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1);
					*(point_block->data.db+1*point_block->step/8+k)=(int)((*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1);
				}
				if (*(afR_XY->data.db+0*afR_XY->step/8+k)<0&&*(afR_XY->data.db+1*afR_XY->step/8+k)>0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)((*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1);
					*(point_block->data.db+1*point_block->step/8+k)=(int)((*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1);
				}

			/*	if (abs(*(point_block->data.db+0*point_block->step/8+k))<=2&&abs(*(point_block->data.db+1*point_block->step/8+k))<=2)
				{
					inblock[num_inblock]=k;
					num_inblock++;
				}*/
				if (abs(*(afR_XY->data.db+0*afR_XY->step/8+k))<=(double)(R_FeatureGRAD+1.0)/2.0&&abs(*(afR_XY->data.db+1*afR_XY->step/8+k))<=(double)(R_FeatureGRAD+1.0)/2.0)
				{
					inblock[num_inblock]=k;
					num_inblock++;
				}

			}
			cout<<"sdsdsdsdsd"<<num_inblock<<endl;


			for (l=0;l<num_inblock;l++)
			{
				int now_n=inblock[l];
				double pix_x=*(afR_XY->data.db+0*afR_XY->step/8+now_n);
				double pix_y=*(afR_XY->data.db+1*afR_XY->step/8+now_n);
				double sita=*(grad_point->data.db+0*grad_point->step/8+now_n);
				double MID_sitabin=sita/45;
				double mod=*(grad_point->data.db+1*grad_point->step/8+now_n);
				int UP_sitabin=(int)(sita/45)+1;
				int DOWN_sitabin=(int)(sita/45);

				int blockx=*(point_block->data.db+0*point_block->step/8+now_n);
				int blocky=*(point_block->data.db+1*point_block->step/8+now_n);

				//cout<<"blockx:"<<blockx<< " blocky"<<blocky<<endl;
				//====

				//方向赋值与计算
				//========================================
				if (blockx<0&&blocky<0)
				{
					if (abs(blockx)>2&&abs(blocky)<=2)
					{
						blockx=-2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)<=2&&abs(blocky)>2)
					{
						blocky=-2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)>2&&abs(blocky)>2)
					{
						blockx=-2;
						blocky=-2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);

					}

					if (abs(blockx)<=2&&abs(blocky)<=2)
					{
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
				}
				//========================================


				//========================================
				if (blockx<0&&blocky>0)
				{
					if (abs(blockx)>2&&abs(blocky)<=2)
					{
						blockx=-2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)<=2&&abs(blocky)>2)
					{
						blocky=2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)>2&&abs(blocky)>2)
					{
						blockx=-2;
						blocky=2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);

					}

					if (abs(blockx)<=2&&abs(blocky)<=2)
					{
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx+1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+2][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
				}
				//================================================




				//================================================
				if (blockx>0&&blocky<0)
				{
					if (abs(blockx)>2&&abs(blocky)<=2)
					{
						blockx=2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)<=2&&abs(blocky)>2)
					{
						blocky=-2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)>2&&abs(blocky)>2)
					{
						blockx=2;
						blocky=-2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}

					if (abs(blockx)<=2&&abs(blocky)<=2)
					{
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky+1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+2].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
				}
				//=======================================================


				//=======================================================

				if (blockx>0&&blocky>0)
				{
					if (abs(blockx)>2&&abs(blocky)<=2)
					{
						blockx=2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)<=2&&abs(blocky)>2)
					{
						blocky=2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
					if (abs(blockx)>2&&abs(blocky)>2)
					{
						blockx=2;
						blocky=2;
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);

					}

					if (abs(blockx)<=2&&abs(blocky)<=2)
					{
						//===测试代码
					/*	double _dx=1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4);
						double _dy=1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4);
						double _dupsita=1.0-(double)UP_sitabin+MID_sitabin;
						double _ddownsita=1.0-(double)DOWN_sitabin+MID_sitabin;
						double upxmod=abs(mod*_dx*_dy*_dupsita);
						double downmod=abs(mod*_dx*_dy*_ddownsita);*/
						double up_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)UP_sitabin+MID_sitabin));
						double down_mod=abs(mod*(1.0-abs(pix_x-(blockx-1.0))/(R_FeatureGRAD/4))*(1.0-abs(pix_y-(blocky-1.0))/(R_FeatureGRAD/4))*(1.0-(double)DOWN_sitabin+MID_sitabin));
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+UP_sitabin)+=up_mod*Gaussian_value(pix_x,pix_y,0.5);
						*(F->oct_imgFeature[i].p_Feature[j].blocks[blockx+1][blocky+1].eight_histo->data.db+DOWN_sitabin)+=down_mod*Gaussian_value(pix_x,pix_y,0.5);
					}
				}

				//=======================================================
			}
			
			for (k=0;k<4;k++)
			{
				for (l=0;l<4;l++)
				{
					for (int q=0;q<8;q++)
					{
						double nmnm=*(F->oct_imgFeature[i].p_Feature[j].blocks[k][l].eight_histo->data.db+q);
						*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k*32+l*8+q)=*(F->oct_imgFeature[i].p_Feature[j].blocks[k][l].eight_histo->data.db+q);
						
					}
				}	
			}


			//归一化128维算子
			double add_mod=0;
			for (k=0;k<128;k++)
			{
				add_mod=add_mod+*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k);
			}
			for (k=0;k<128;k++)
			{
				//*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)/add_mod;
				*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)/sqrt(add_mod);
				double tsd=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k);
			}
			double psd=0;
			for (k=0;k<128;k++)
			{
				psd+=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k);
			}
			
		}
	}
	return F;

	
}


//128维特征向量=============================fix===
Feature* calFeature_fix(delta Delta,CvMat *spoint,int *num,int R_BasicDrection,int R_FeatureGRAD)
{
	int xoct=Delta.longoct;
	int i,j,k,l;
	double Top_num=0;
	int top_k=0;
	double *sigma=(double *)malloc(sizeof(double)*xoct);
	for(i=0;i<xoct;i++)
	{
		*(sigma+i)=pow(2,i)*1.6;
	}
	//================================================================================================================
	//初步测试 使用8*8探测基本特征点方向 24*24探测梯度

	
	int NUM_R_BasicDrection=R_BasicDrection*R_BasicDrection;
	int NUM_R_FeatureGRAD=R_FeatureGRAD*R_FeatureGRAD;

	CvMat *Point_GradMod=cvCreateMat(spoint->rows,spoint->cols*NUM_R_BasicDrection,CV_64FC1);  //存放每个特征点以及邻域点的基本特征方向矩阵信息
	CvMat *Spoint_GradMod=cvCreateMat(spoint->rows,spoint->cols,CV_64FC1); //若有辅助方向 新建的numGM也要加上数量也要相应加上数量
	double thirtysix_bin[36];  //存放36柱的直方图

	//bfR_XY未进行基础方向旋转的坐标
	CvMat* bfR_XY=cvCreateMat(2,10000,CV_64FC1);//坐标中心为特征点坐标
	int num_abfR_XY=0;//邻域内点的个数

	//afR_XY进行基础方向旋转后的坐标
	CvMat* afR_XY=cvCreateMat(2,10000,CV_64FC1);


	//计算领域所需要的半径点坐标
	double R=(1.4*R_FeatureGRAD+1)/2;
	cout<<"检测半径:"<<R<<endl;

	for (k=-R;k<=R;k++)
	{
		for (l=-R;l<=R;l++)
		{
			if (sqrt(pow(k,2)+pow(l,2))<R)
			{
				if (k!=0&&l!=0)
				{
					*(bfR_XY->data.db+0*bfR_XY->step/8+num_abfR_XY)=k;
					*(bfR_XY->data.db+1*bfR_XY->step/8+num_abfR_XY)=l;
					//cout<<"x:"<<k<<" "<<"y:"<<l<<endl;
					num_abfR_XY++; 
				}

			}

		}
	}

	cout<<"领域个数"<<num_abfR_XY<<endl;



	//储存每个特征点的特征向量 也就是最终的SIFT特征向量
	//每一层有两行 第一行为角度 第二行为幅值
	//CvMat *Feature=cvCreateMat(spoint->rows,128*spoint->cols,CV_64FC1);



	//================初始化======================================

	Feature *F;
	F=(Feature*)malloc(sizeof(Feature));
	init_Feature_fix(F,xoct,num);
	//============================================================



	for (i=0;i<xoct;i++)
	{
		for (j=0;j<*(num+i);j++)
		{

			/*R_BasicDrection=3*(int)*(sigma+i)*1.5;
			NUM_R_BasicDrection=R_BasicDrection*R_BasicDrection;*/

			Top_num=0;
			top_k=0;
			for(k=0;k<36;k++)
			{
				thirtysix_bin[k]=0;
			}

			double x=*(spoint->data.db+2*i*spoint->step/8+j);
			double y=*(spoint->data.db+(2*i+1)*spoint->step/8+j);
			int Intx=cvRound(x);
			int Inty=cvRound(y); 

			CvMat *NowMat=Delta.oct[i].inter[0].Gussianimg;
			

			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k-R_BasicDrection/2,Inty+l-R_BasicDrection/2);
					//cout<<*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)<<"  "; 
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k-R_BasicDrection/2,Inty+l-R_BasicDrection/2)*Gaussian_value(k-R_BasicDrection/2,l-R_BasicDrection/2,*(sigma+i)*1.5);
				}
			}




			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+1*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k-R_BasicDrection/2,Inty+1+l);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+1*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k-R_BasicDrection/2,Inty+1+l)*Gaussian_value(k-R_BasicDrection/2,1+l,*(sigma+i)*1.5);
				}
			}


			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+2*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k+1,Inty+l-R_BasicDrection/2);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+2*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k+1,Inty+l-R_BasicDrection/2)*Gaussian_value(k+1,l-R_BasicDrection/2,*(sigma+i)*1.5);
				}
			}



			for (k=0;k<R_BasicDrection/2;k++)
			{
				for (l=0;l<R_BasicDrection/2;l++)
				{
					*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+3*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Grad_sita(NowMat,Intx+k+1,Inty+l+1);
					*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+3*NUM_R_BasicDrection/4+R_BasicDrection/2*k+l)=Mod(NowMat,Intx+k+1,Inty+l+1)*Gaussian_value(k+1,l+1,*(sigma+i)*1.5);
				}
			}

			//建立初始特征点方向直方图
			for (k=0;k<NUM_R_BasicDrection;k++)
			{
				double now_sita=*(Point_GradMod->data.db+2*i*Point_GradMod->step/8+j*NUM_R_BasicDrection+k);
				double now_mod=*(Point_GradMod->data.db+(2*i+1)*Point_GradMod->step/8+j*NUM_R_BasicDrection+k);

				if (now_sita<0)
				{
					thirtysix_bin[(int)(now_sita/10)+36]+=now_mod;
				} 
				else if(now_sita>0)
				{
					thirtysix_bin[(int)(now_sita/10)]+=now_mod;
				}
				
			}


			//对直方图进行高斯平滑
			thirtysix_bin[0]=(thirtysix_bin[0]*6+thirtysix_bin[1]*4+thirtysix_bin[2])/11;
			thirtysix_bin[1]=(thirtysix_bin[0]*4+thirtysix_bin[1]*6+thirtysix_bin[2]*4+thirtysix_bin[3])/15;
			for (k=2;k<34;k++)
			{
				thirtysix_bin[k]=(thirtysix_bin[k-2]+thirtysix_bin[k-1]*4+thirtysix_bin[k]*6+thirtysix_bin[k+1]*4+thirtysix_bin[k+2])/16;
			}
			thirtysix_bin[34]=(thirtysix_bin[35]*4+thirtysix_bin[34]*6+thirtysix_bin[33]*4+thirtysix_bin[32])/15;
			thirtysix_bin[35]=(thirtysix_bin[35]*6+thirtysix_bin[34]*4+thirtysix_bin[33])/11;


			


			Top_num=thirtysix_bin[0];
			for (k=0;k<36;k++)
			{
				if (thirtysix_bin[k]>Top_num)
				{
					Top_num=thirtysix_bin[k];
					top_k=k;
				}
			}

			*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j)=top_k*10;
			*(Spoint_GradMod->data.db+(2*i+1)*Spoint_GradMod->step/8+j)=Top_num;

			//printf("角度%10lf 幅值%10lf",*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j),*(Spoint_GradMod->data.db+(2*i+1)*Spoint_GradMod->step/8+j));
			//printf("角度%10lf\n",*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j));

//旋转坐标探测完成===========================================================================================================================

			//对每个种子点领域像素进行坐标方向旋转
			CvMat *grad_point=cvCreateMat(2,10000,CV_64FC1);//每个特征点邻域坐标的sita mod       
			CvMat *point_block=cvCreateMat(2,10000,CV_64FC1);//每个领域坐标在旋转后坐标的方块坐标x y
			int inblock[10000];//在领域内的点的序号
			int num_inblock=0;

			int q,w;  
			for (k=0;k<num_abfR_XY;k++)
			{
				double point_sita=(*(Spoint_GradMod->data.db+2*i*Spoint_GradMod->step/8+j))*3.14/180;
				double bfx=*(bfR_XY->data.db+0*bfR_XY->step/8+k);
				double bfy=*(bfR_XY->data.db+1*bfR_XY->step/8+k);
				*(afR_XY->data.db+0*afR_XY->step/8+k)=bfx*cos(point_sita)-bfy*sin(point_sita);
				*(afR_XY->data.db+1*afR_XY->step/8+k)=bfx*sin(point_sita)+bfy*cos(point_sita);

				//check OK===================================================================

				int afX_block,afY_block;//像素所在旋转后方块的方块坐标

				*(grad_point->data.db+0*grad_point->step/8+k)=Grad_sita(NowMat,Intx+bfx,Inty+bfy);
				*(grad_point->data.db+1*grad_point->step/8+k)=Mod(NowMat,Intx+bfx,Inty+bfy);

				//check OK============grad_point 第一行为领域序号的角度(360) 第二行为领域序号的幅值


				if (*(afR_XY->data.db+0*afR_XY->step/8+k)>0&&*(afR_XY->data.db+1*afR_XY->step/8+k)<0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)(*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1;
					*(point_block->data.db+1*point_block->step/8+k)=(int)(*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1;
				}
				if (*(afR_XY->data.db+0*afR_XY->step/8+k)>0&&*(afR_XY->data.db+1*afR_XY->step/8+k)>0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)(*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1;
					*(point_block->data.db+1*point_block->step/8+k)=(int)(*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1;
				}
				if (*(afR_XY->data.db+0*afR_XY->step/8+k)<0&&*(afR_XY->data.db+1*afR_XY->step/8+k)<0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)(*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1;
					*(point_block->data.db+1*point_block->step/8+k)=(int)(*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1;
				}
				if (*(afR_XY->data.db+0*afR_XY->step/8+k)<0&&*(afR_XY->data.db+1*afR_XY->step/8+k)>0)
				{
					*(point_block->data.db+0*point_block->step/8+k)=(int)(*(afR_XY->data.db+0*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)-1;
					*(point_block->data.db+1*point_block->step/8+k)=(int)(*(afR_XY->data.db+1*afR_XY->step/8+k))/(int)(R_FeatureGRAD/4)+1;
				}

			/*	if (abs(*(point_block->data.db+0*point_block->step/8+k))<=2&&abs(*(point_block->data.db+1*point_block->step/8+k))<=2)
				{
					inblock[num_inblock]=k;
					num_inblock++;
				}*/
				if (abs(*(afR_XY->data.db+0*afR_XY->step/8+k))<=(double)(R_FeatureGRAD+1.0)/2.0&&abs(*(afR_XY->data.db+1*afR_XY->step/8+k))<=(double)(R_FeatureGRAD+1.0)/2.0)
				{
					inblock[num_inblock]=k;
					num_inblock++;
				}

			}
			//cout<<"sdsdsdsdsd"<<num_inblock<<endl;


			for (l=0;l<num_inblock;l++)
			{
				int now_n=inblock[l];
				double pix_x=*(afR_XY->data.db+0*afR_XY->step/8+now_n);
				double pix_y=*(afR_XY->data.db+1*afR_XY->step/8+now_n);
				double sita=*(grad_point->data.db+0*grad_point->step/8+now_n);
				double MID_sitabin=sita/45;
				double mod=*(grad_point->data.db+1*grad_point->step/8+now_n);
				int UP_sitabin=(int)(sita/45)+1;
				if (UP_sitabin==8)
				{
					UP_sitabin=0;
				}
				int DOWN_sitabin=(int)(sita/45);

				int blockx=*(point_block->data.db+0*point_block->step/8+now_n);
				int blocky=*(point_block->data.db+1*point_block->step/8+now_n);


				double lu;
				double ld;
				double ru;
				double rd;
				double ruld_add=0;

				double lu_d;
				double ld_d;
				double ru_d;
				double rd_d;

				int lu_X;
				int lu_Y;
				int ld_X;
				int ld_Y;
				int ru_X;
				int ru_Y;
				int rd_X;
				int rd_Y;

				ru_X=lu_X=blockx;
				ld_Y=lu_Y=blocky;
				rd_X=ld_X=blockx-blockx/abs(blockx);
				rd_Y=ru_Y=blocky-blocky/abs(blocky);

				lu_d=sqrt(pow(pix_x-lu_X,2)+pow(pix_y-lu_Y,2));
				ld_d=sqrt(pow(pix_x-ld_X,2)+pow(pix_y-ld_Y,2));
				ru_d=sqrt(pow(pix_x-ru_X,2)+pow(pix_y-ru_Y,2));
				rd_d=sqrt(pow(pix_x-rd_X,2)+pow(pix_y-rd_Y,2));

				ruld_add=lu_d+ld_d+ru_d+rd_d;

				lu=lu_d/ruld_add;
				ld=lu_d/ruld_add;
				ru=ru_d/ruld_add;
				rd=rd_d/ruld_add;

				double downs=MID_sitabin-(double)(DOWN_sitabin);
				double ups=1.0-downs;
				

				//lu加成===============================================================
				if (lu_X+2>=0&&lu_X+2<5&&lu_Y+2>=0&&lu_Y+2<5)
				{
					*(F->oct_imgFeature[i].p_Feature[j].blocks[lu_X+2][lu_Y+2].eight_histo->data.db+UP_sitabin)=mod*lu*ups*Gaussian_value(pix_x,pix_y,0.5);
					*(F->oct_imgFeature[i].p_Feature[j].blocks[lu_X+2][lu_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*lu*downs*Gaussian_value(pix_x,pix_y,0.5);
				}
				
				
				//ld加成===============================================================
				if (ld_X+2>=0&&ld_X+2<5&&ld_Y+2>=0&&ld_Y+2<5)
				{
					*(F->oct_imgFeature[i].p_Feature[j].blocks[ld_X+2][ld_Y+2].eight_histo->data.db+UP_sitabin)=mod*ld*ups*Gaussian_value(pix_x,pix_y,0.5);
					*(F->oct_imgFeature[i].p_Feature[j].blocks[ld_X+2][ld_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*ld*downs*Gaussian_value(pix_x,pix_y,0.5);
				}
				

				//ru加成===============================================================
				if (ru_X+2>=0&&ru_X+2<5&&ru_Y+2>=0&&ru_Y+2<5)
				{
					*(F->oct_imgFeature[i].p_Feature[j].blocks[ru_X+2][ru_Y+2].eight_histo->data.db+UP_sitabin)=mod*ru*ups*Gaussian_value(pix_x,pix_y,0.5);
					*(F->oct_imgFeature[i].p_Feature[j].blocks[ru_X+2][ru_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*ru*downs*Gaussian_value(pix_x,pix_y,0.5);
				}
				

				//rd加成===============================================================
				if (rd_X+2>=0&&rd_X+2<5&&rd_Y+2>=0&&rd_Y+2<5)
				{
					*(F->oct_imgFeature[i].p_Feature[j].blocks[rd_X+2][rd_Y+2].eight_histo->data.db+UP_sitabin)=mod*rd*ups*Gaussian_value(pix_x,pix_y,0.5);
					*(F->oct_imgFeature[i].p_Feature[j].blocks[rd_X+2][rd_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*rd*downs*Gaussian_value(pix_x,pix_y,0.5);
				}
				//=======================================================

			////	lu加成===============================================================
			//	if (lu_X+2>=0&&lu_X+2<5&&lu_Y+2>=0&&lu_Y+2<5)
			//	{
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[lu_X+2][lu_Y+2].eight_histo->data.db+UP_sitabin)=mod*lu*ups;
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[lu_X+2][lu_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*lu*downs;
			//	}


			////	ld加成===============================================================
			//	if (ld_X+2>=0&&ld_X+2<5&&ld_Y+2>=0&&ld_Y+2<5)
			//	{
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[ld_X+2][ld_Y+2].eight_histo->data.db+UP_sitabin)=mod*ld*ups;
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[ld_X+2][ld_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*ld*downs;
			//	}


			////	ru加成===============================================================
			//	if (ru_X+2>=0&&ru_X+2<5&&ru_Y+2>=0&&ru_Y+2<5)
			//	{
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[ru_X+2][ru_Y+2].eight_histo->data.db+UP_sitabin)=mod*ru*ups;
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[ru_X+2][ru_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*ru*downs;
			//	}


			////	rd加成===============================================================
			//	if (rd_X+2>=0&&rd_X+2<5&&rd_Y+2>=0&&rd_Y+2<5)
			//	{
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[rd_X+2][rd_Y+2].eight_histo->data.db+UP_sitabin)=mod*rd*ups;
			//		*(F->oct_imgFeature[i].p_Feature[j].blocks[rd_X+2][rd_Y+2].eight_histo->data.db+DOWN_sitabin)=mod*rd*downs;
			//	}
			}
			
			for (k=0;k<5;k++)
			{
				for (l=0;l<5;l++)
				{
					for (int q=0;q<8;q++)
					{
						double nmnm=*(F->oct_imgFeature[i].p_Feature[j].blocks[k][l].eight_histo->data.db+q);
						*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k*40+l*8+q)=*(F->oct_imgFeature[i].p_Feature[j].blocks[k][l].eight_histo->data.db+q);
						
					}
				}	
			}


			//归一化200维算子
			double add_mod=0;
			for (k=0;k<200;k++)
			{
				add_mod=add_mod+pow(*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k),2);
			}
			for (k=0;k<200;k++)
			{
				//*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)/add_mod;
				*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k)/sqrt(add_mod);
				double tsd=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k);
			}
			double psd=0;
			for (k=0;k<200;k++)
			{
				psd+=*(F->oct_imgFeature[i].p_Feature[j].Feature_Vector->data.db+k);
			}
			
		}
	}
	return F;

	
}


#endif

