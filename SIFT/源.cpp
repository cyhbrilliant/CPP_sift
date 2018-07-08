#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <iostream>
#include <math.h>
#include "MatOUT.h"
#include "GUSSIANmodel.h"
#include "GDgetmodel.h"
#include "GUSSIANdelta.h"
#include "DeltaOUT.h"
#include "searchpoint.h"
#include "SearchOUT.h"
#include "fixpoint.h"
//#include "histogram;.h"
#include "histo_cir.h"
#include "SIFT.h"
#include "O_distance.h"
#include "mul_img.h"

using namespace cv;
using namespace std;
 


int main()
{
	int oct_F1=0;
	int oct_F2=0;

	int R_BasicDrection=8;
	int R_FeatureGRAD_img1_1=32;
//	int R_FeatureGRAD_img1_2=40;

	int R_FeatureGRAD_img2_1=32;
//	int R_FeatureGRAD_img2_2=40;

	int nameimg1=14;
	int nameimg2=15;
	
	double _sita=15;
	double _ratio=0.6;
	
	char img1_name[50];
	char img2_name[50];


	sprintf_s(img1_name,"test%d.jpg",nameimg1);
	sprintf_s(img2_name,"test%d.jpg",nameimg2);

	IplImage *img1= cvLoadImage(img1_name,0); 
	IplImage *img2=cvLoadImage(img2_name,0);
	//IplImage *img2=cvCreateImage(cvGetSize(img2x),IPL_DEPTH_8U,1);

	//CvMat *qw=cvCreateMat(img2x->height,img2->width,CV_64FC1);
	//cvConvert(img2x,qw);
	//for (int i=0;i<qw->rows;i++)
	//{
	//	for (int j=0;j<qw->cols;j++)
	//	{
	//		double t=*(qw->data.db+i*qw->step/8+j);
	//		*(qw->data.db+i*qw->step/8+j)=(uchar)(t+20);
	//	}
	//}
	//cvConvert(qw,img2);
	//

	//IplImage *img1=cvCreateImage(cvSize(img2->width/2,img2->height/2),IPL_DEPTH_8U,1);
	//cvPyrDown(img1x,img1);

	CvMat *position_img1=cvCreateMat(2,50000,CV_64FC1);
	CvMat *position_img2=cvCreateMat(2,50000,CV_64FC1);

	Feature *F1,*F2;
	F1=sift(img1,img1_name,position_img1,oct_F1,R_BasicDrection,R_FeatureGRAD_img1_1);
	F2=sift(img2,img2_name,position_img2,oct_F2,R_BasicDrection,R_FeatureGRAD_img2_1);

	cout<<"Feature finish"<<endl;

	CvMat *dist;
	dist=O_dis_fix(F1,F2,oct_F1,oct_F2);



	CvMat *point_pos1=cvCreateMat(2,50000,CV_64FC1);//第一行img1的特征点序号 第二行img2的特征点序号
	CvMat *point_pos2=cvCreateMat(2,50000,CV_64FC1);
	CvMat *point_pos=cvCreateMat(2,50000,CV_64FC1);




	int numx=0;
	int numx1=0;
	int numx2=0;


	double odis1=0;
	double odis2=0;
	int i,j,k,l;
	double comp1=0;
	double comp2=0;
	int reco1=0;
	int reco2=0;



	//========================================================================================
	//for (i=0;i<F1[0]->oct_imgFeature[oct_F1].Feature_num;i++)
	//{
	//	comp1=0;
	//	comp2=0;
	//	for (j=0;j<F2[0]->oct_imgFeature[oct_F2].Feature_num;j++)
	//	{
	//		odis1=*(dist1->data.db+i*dist1->step/8+j);
	//		if (j==0)
	//		{
	//			comp1=odis1;
	//			reco1=0;
	//		}
	//		else if(odis1<comp1)
	//		{
	//			comp1=odis1;
	//			reco1=j;

	//		}
	//	}


	//	*(point_pos1->data.db+0*point_pos1->step/8+numx1)=i;
	//	*(point_pos1->data.db+1*point_pos1->step/8+numx1)=reco1;
	//	numx1++;
	//	
	//	
	//}

	//for (i=0;i<F1[1]->oct_imgFeature[oct_F1].Feature_num;i++)
	//{
	//	comp1=0;
	//	comp2=0;
	//	for (j=0;j<F2[1]->oct_imgFeature[oct_F2].Feature_num;j++)
	//	{
	//		odis1=*(dist2->data.db+i*dist2->step/8+j);
	//		if (j==0)
	//		{
	//			comp1=odis1;
	//			reco1=0;
	//		}
	//		else if(odis1<comp1)
	//		{
	//			comp1=odis1;
	//			reco1=j;

	//		}
	//	}


	//	*(point_pos2->data.db+0*point_pos2->step/8+numx2)=i;
	//	*(point_pos2->data.db+1*point_pos2->step/8+numx2)=reco1;
	//	numx2++;


	//}

	//cout<<numx1<<endl<<numx2<<endl;
	//for (i=0;i<numx1;i++)
	//	for(j=0;j<numx2;j++)
	//		{

	//			if(*(point_pos2->data.db+0*point_pos2->step/8+j)==*(point_pos1->data.db+0*point_pos1->step/8+i))
	//			{
	//				if (*(point_pos2->data.db+1*point_pos2->step/8+j)==*(point_pos1->data.db+1*point_pos1->step/8+i))
	//				{
	//					//平行约束======================================================
	//				int img1_pos=*(point_pos1->data.db+0*point_pos1->step/8+i);
	//				int img2_pos=*(point_pos1->data.db+1*point_pos1->step/8+i);
	//				double x_img1=*(position_img1->data.db+0*position_img1->step/8+img1_pos);
	//				double y_img1=*(position_img1->data.db+1*position_img1->step/8+img1_pos);
	//				double x_img2=*(position_img2->data.db+0*position_img2->step/8+img2_pos);
	//				double y_img2=*(position_img2->data.db+1*position_img2->step/8+img2_pos);
	//				double d_x=x_img1-x_img2;
	//				double d_y=y_img1-y_img2;
	//				double sita=atan(d_x/d_y)*180.0/3.14;
	//				if (abs(sita)<_sita)
	//				{
	//					*(point_pos->data.db+0*point_pos->step/8+numx)=*(point_pos1->data.db+0*point_pos1->step/8+i);
	//					*(point_pos->data.db+1*point_pos->step/8+numx)=*(point_pos1->data.db+1*point_pos1->step/8+i);
	//					numx++;
	//				}

	//					*(point_pos->data.db+0*point_pos->step/8+numx)=*(point_pos1->data.db+0*point_pos1->step/8+i);
	//					*(point_pos->data.db+1*point_pos->step/8+numx)=*(point_pos1->data.db+1*point_pos1->step/8+i);
	//					numx++;
	//				}
	//			}

	//		}

	//		
	//cout<<numx<<endl;
	//		=====================================================================================================
			//if(*(point_pos2->data.db+0*point_pos2->step/8+j)==*(point_pos1->data.db+0*point_pos1->step/8+i))
			//{
			//	if (*(point_pos2->data.db+1*point_pos2->step/8+j)==*(point_pos1->data.db+1*point_pos1->step/8+i))
			//	{
			//		//平行约束======================================================
			//		int img1_pos=*(point_pos1->data.db+0*point_pos1->step/8+i);
			//		int img2_pos=*(point_pos1->data.db+1*point_pos1->step/8+i);
			//		double x_img1=*(position_img1->data.db+0*position_img1->step/8+img1_pos);
			//		double y_img1=*(position_img1->data.db+1*position_img1->step/8+img1_pos);
			//		double x_img2=*(position_img2->data.db+0*position_img2->step/8+img2_pos);
			//		double y_img2=*(position_img2->data.db+1*position_img2->step/8+img2_pos);
			//		double d_x=x_img1-x_img2;
			//		double d_y=y_img1-y_img2;
			//		double sita=atan(d_x/d_y)*180.0/3.14;
			//		if (abs(sita)<_sita)
			//		{
			//			*(point_pos->data.db+0*point_pos->step/8+numx)=*(point_pos1->data.db+0*point_pos1->step/8+i);
			//			*(point_pos->data.db+1*point_pos->step/8+numx)=*(point_pos1->data.db+1*point_pos1->step/8+i);
			//			numx++;
			//		}
			//		//*(point_pos->data.db+0*point_pos->step/8+numx)=*(point_pos1->data.db+0*point_pos1->step/8+i);
			//		//*(point_pos->data.db+1*point_pos->step/8+numx)=*(point_pos1->data.db+1*point_pos1->step/8+i);
			//		//numx++;
			//	}
			//}



	//========================================================
	for (i=0;i<F1->oct_imgFeature[oct_F1].Feature_num;i++)
	{
		comp1=0;
		comp2=0;
		for (j=0;j<F2->oct_imgFeature[oct_F2].Feature_num;j++)
		{
			odis1=*(dist->data.db+i*dist->step/8+j);
			if (j==0)
			{
				comp1=odis1;
				reco1=0;
			}
			else if(odis1<comp1)
			{
				comp1=odis1;
				reco1=j;

			}
		}

		for (j=0;j<F2->oct_imgFeature[oct_F2].Feature_num;j++)
		{
			odis2=*(dist->data.db+i*dist->step/8+j);
			if (j==0)
			{
				comp2=odis2;
				reco2=0;
			}
			else if(j==reco1)
			{

			}
			else if(odis2<comp2)
			{
				comp2=odis2;
				reco2=j;

			}
		}
		
		//if (comp<0.05)
		//{
		//	*(point_pos1->data.db+0*point_pos1->step/8+numx1)=i;
		//	*(point_pos1->data.db+1*point_pos1->step/8+numx1)=reco;
		//	numx1++;
		//}

		double ratio=comp1/comp2;
		if (ratio<_ratio)
		{
			/*if (comp1<)
			{
			}*/
			*(point_pos->data.db+0*point_pos->step/8+numx)=i;
			*(point_pos->data.db+1*point_pos->step/8+numx)=reco1;
			numx++;
		}
		//*(point_pos->data.db+0*point_pos->step/8+numx)=i;
		//*(point_pos->data.db+1*point_pos->step/8+numx)=reco1;
		//numx++;
		
		
	}

	////for (i=0;i<F2->oct_imgFeature[oct_F2].Feature_num;i++)
	////{
	////	comp=0;
	////	for (j=0;j<F1->oct_imgFeature[oct_F1].Feature_num;j++)
	////	{
	////		odis=*(dist->data.db+i*dist->step/8+j);
	////		if (j==0)
	////		{
	////			comp=odis;
	////			reco=0;
	////		}
	////		else if(odis<comp)
	////		{
	////			comp=odis;
	////			reco=j;

	////		}
	////		/*if (odis<0.05)
	////		{
	////			*(point_pos->data.db+0*point_pos->step/8+numx)=i;
	////			*(point_pos->data.db+1*point_pos->step/8+numx)=j;
	////			numx++;
	////			
	////		}*/
	////	}
	/////*	if (comp<0.05)
	////	{
	////		*(point_pos2->data.db+0*point_pos2->step/8+numx2)=i;
	////		*(point_pos2->data.db+1*point_pos2->step/8+numx2)=reco;
	////		numx2++;
	////	}*/
	////	*(point_pos2->data.db+0*point_pos2->step/8+numx2)=i;
	////	*(point_pos2->data.db+1*point_pos2->step/8+numx2)=reco;
	////	numx2++;
	////}
//	for (i=0;i<F2->oct_imgFeature[oct_F2].Feature_num;i++)
//	{
//		comp1=0;
//		comp2=0;
//		for (j=0;j<F1->oct_imgFeature[oct_F1].Feature_num;j++)
//		{
//			odis1=*(dist->data.db+j*dist->step/8+i);
//			if (j==0)
//			{
//				comp1=odis1;
//				reco1=0;
//			}
//			else if(odis1<comp1)
//			{
//				comp1=odis1;
//				reco1=j;
//
//			}
//		}
//
//		for (j=0;j<F1->oct_imgFeature[oct_F1].Feature_num;j++)
//		{
//			odis2=*(dist->data.db+j*dist->step/8+i);
//			if (j==0)
//			{
//				comp2=odis2;
//				reco2=0;
//			}
//			else if(j==reco1)
//			{
//
//			}
//			else if(odis2<comp2)
//			{
//				comp2=odis2;
//				reco2=j;
//
//			}
//		}
//		//
//		//if (comp<0.05)
//		//{
//		//	*(point_pos1->data.db+0*point_pos1->step/8+numx1)=i;
//		//	*(point_pos1->data.db+1*point_pos1->step/8+numx1)=reco;
//		//	numx1++;
//		//}
//
//		ratio=comp1/comp2;
//		if (ratio<_ratio)
//		{
//			/*if (comp1<)
//			{
//			}*/
//			*(point_pos2->data.db+0*point_pos2->step/8+numx2)=i;
//			*(point_pos2->data.db+1*point_pos2->step/8+numx2)=reco1;
//			numx2++;
//		}
//		//*(point_pos->data.db+0*point_pos->step/8+numx)=i;
//		//*(point_pos->data.db+1*point_pos->step/8+numx)=reco1;
//		//numx++;
//		
//		
//	}
////================================================================
//	for (i=0;i<numx1;i++)
//	{
//		for(j=0;j<numx2;j++)
//		{
//			if(*(point_pos2->data.db+0*point_pos2->step/8+j)==*(point_pos1->data.db+0*point_pos1->step/8+i))
//			{
//				//if (*(point_pos2->data.db+1*point_pos2->step/8+j)==*(point_pos1->data.db+1*point_pos1->step/8+i))
//				//{
//				//	//平行约束======================================================
//				//	int img1_pos=*(point_pos1->data.db+0*point_pos1->step/8+i);
//				//	int img2_pos=*(point_pos1->data.db+1*point_pos1->step/8+i);
//				//	double x_img1=*(position_img1->data.db+0*position_img1->step/8+img1_pos);
//				//	double y_img1=*(position_img1->data.db+1*position_img1->step/8+img1_pos);
//				//	double x_img2=*(position_img2->data.db+0*position_img2->step/8+img2_pos);
//				//	double y_img2=*(position_img2->data.db+1*position_img2->step/8+img2_pos);
//				//	double d_x=x_img1-x_img2;
//				//	double d_y=y_img1-y_img2;
//				//	double sita=atan(d_x/d_y)*180.0/3.14;
//				//	if (abs(sita)<_sita)
//				//	{
//				//		*(point_pos->data.db+0*point_pos->step/8+numx)=*(point_pos1->data.db+0*point_pos1->step/8+i);
//				//		*(point_pos->data.db+1*point_pos->step/8+numx)=*(point_pos1->data.db+1*point_pos1->step/8+i);
//				//		numx++;
//				//	}
//				*(point_pos->data.db+0*point_pos->step/8+numx)=*(point_pos1->data.db+0*point_pos1->step/8+i);
//				*(point_pos->data.db+1*point_pos->step/8+numx)=*(point_pos1->data.db+1*point_pos1->step/8+i);
//				numx++;
//				}
//			}
//
//		}
	//}

	cout<<"sift finish"<<endl;

	char file1[100];
	sprintf_s(file1,"G:\\opencv\\1\\test10\\SIFT\\SIFT\\x%s.txt",img1_name);

	FILE *fl1;
	fl1=fopen(file1,"w");


	char file2[100];
	sprintf_s(file2,"G:\\opencv\\1\\test10\\SIFT\\SIFT\\x%s.txt",img2_name);

	FILE *fl2;
	fl2=fopen(file2,"w");

	for (i=0;i<numx;i++)
	{
		int img1_pos=*(point_pos->data.db+0*point_pos->step/8+i);
		int img2_pos=*(point_pos->data.db+1*point_pos->step/8+i);
		fprintf(fl1,"特征点有：  %d   ",img1_pos);
		fprintf(fl2,"特征点有：  %d   ",img2_pos);

	}
	fclose(fl1);
	fclose(fl2);

	
	//连线=============
	double x_img1=0;
	double y_img1=0;
	double x_img2=0;
	double y_img2=0;
	int x1,y1,x2,y2;
	IplImage *out;
	out=mulimg(img1_name,img2_name,img1->height,img1->width,img2->height,img2->width);

	for (i=0;i<numx;i++)
	{
		int img1_pos=*(point_pos->data.db+0*point_pos->step/8+i);
		int img2_pos=*(point_pos->data.db+1*point_pos->step/8+i);
		x_img1=*(position_img1->data.db+0*position_img1->step/8+img1_pos);
		y_img1=*(position_img1->data.db+1*position_img1->step/8+img1_pos);
		x_img2=*(position_img2->data.db+0*position_img2->step/8+img2_pos);
		y_img2=*(position_img2->data.db+1*position_img2->step/8+img2_pos);
		x1=cvRound(x_img1);
		y1=cvRound(y_img1);
		x2=cvRound(x_img2);
		y2=cvRound(y_img2);
		//cout<<"1:"<<x1<<"_"<<y1<<" 2:"<<x2<<"_"<<y2<<endl;
		cvLine(out,cvPoint(y1,x1),cvPoint(y2+img1->width,x2),CV_RGB(0,0,255),2);

		


	}
	cvNamedWindow("x",0);
	cvShowImage("x",out);
	cvWaitKey(0);






	cout<<numx<<endl;


	system("pause");



	return 0;
}


//
////==============搜索点时row为x col为y  opencv默认col为x row为y
//
//
//
////读取图像
//IplImage *bfimg= cvLoadImage("test9.jpg",0); 
//
//
//
////图像转换成矩阵
//CvMat* bf=cvCreateMat(bfimg->height,bfimg->width,CV_64FC1);
//cvConvert(bfimg,bf);
//
//
//
////创建高斯金字塔
//delta Delta;
//Delta=CreateDelta(bf,6);      //未升采样
////Delta=upCreateDelta(bf,6); //升采样
//
//
//
//
////显示高斯金字塔
////   Delta_OUT(Delta);
//
//
//
//
////高斯金字塔差分为DOG金字塔
//delta dogdelta;
//dogdelta=DOG(Delta);
//
//
//
////显示DOG金字塔
////	Delta_OUT(dogdelta);
//
//
//
////创建num为每个尺度空间中特征点个数的栈数组
//int *num=(int*)malloc(sizeof(int)*dogdelta.longoct);
//
//
//
////创建allsearch为每个尺度空间中的特征点
//CvMat *allsearch=cvCreateMat(2*dogdelta.longoct,100000,CV_64FC1);
//
//
//
////SearchP函数为找金字塔找点函数
//allsearch=SearchP(dogdelta,num);
//
//
//
////点的修正与精确
//FixPoint(allsearch,num,dogdelta);
//
//
//
////创建矩阵存放每个特征点的初始旋转参数 sita mod
//CvMat *spoint_gradsita=cvCreateMat(allsearch->rows,allsearch->cols,CV_64FC1);
//
//
//
////特征点初始旋转参数
//spoint_gradsita=BuildHistogram(Delta,allsearch,num);
//
//
//
//
//
////128维特征向量
//Feature *F;
//F=calFeature(Delta,allsearch,num);
//
//
//
//
//
////输出点在图上的位置
//search_OUT(dogdelta,allsearch,num,Delta,spoint_gradsita);
//
//
//
//
//
//
//
//
//free(num);

//cout<<int(-1.5)<<endl;

////=======高斯模糊=============
////GussianModel为计算高斯高斯函数的模板 但未保证和值为一 也就是图像亮度会因为权值而变化
////Gussian2one为高斯归一化函数 让模板和值为一 保证亮度不变性
//CvMat *gsm=cvCreateMat(5,5,CV_64FC1);
//gsm=GussianModel(5,0.6);
////cvNormalize(gsm,gsm);
//gsm=Gussian2one(gsm,5,5);
//cout<<"高斯模板"<<endl;
//MatOUT(gsm,5,5,1);
//cout<<endl;
////===========================



////======检测高斯算子==============

//IplImage *frame_gray,*frame_gussian;
//frame_gray=cvLoadImage("test2.jpg",0);
//frame_gussian=cvCreateImage(cvSize(frame_gray->width,frame_gray->height),frame_gray->depth,frame_gray->nChannels);
////int i,j;
////double g00,g01,g02,g03,g04,g10,g11,g12,g13,g14,g20,g21,g22,g23,g24,g30,g31,g32,g33,g34,g40,g41,g42,g43,g44;
////g00=*(gsm->data.db+0*gsm->step/8+0);
////g01=*(gsm->data.db+0*gsm->step/8+1);
////g02=*(gsm->data.db+0*gsm->step/8+2);
////g03=*(gsm->data.db+0*gsm->step/8+3);
////g04=*(gsm->data.db+0*gsm->step/8+4);

////g10=*(gsm->data.db+1*gsm->step/8+0);
////g11=*(gsm->data.db+1*gsm->step/8+1);
////g12=*(gsm->data.db+1*gsm->step/8+2);
////g13=*(gsm->data.db+1*gsm->step/8+3);
////g14=*(gsm->data.db+1*gsm->step/8+4);

////g20=*(gsm->data.db+2*gsm->step/8+0);
////g21=*(gsm->data.db+2*gsm->step/8+1);
////g22=*(gsm->data.db+2*gsm->step/8+2);
////g23=*(gsm->data.db+2*gsm->step/8+3);
////g24=*(gsm->data.db+2*gsm->step/8+4);

////g30=*(gsm->data.db+3*gsm->step/8+0);
////g31=*(gsm->data.db+3*gsm->step/8+1);
////g32=*(gsm->data.db+3*gsm->step/8+2);
////g33=*(gsm->data.db+3*gsm->step/8+3);
////g34=*(gsm->data.db+3*gsm->step/8+4);

////g40=*(gsm->data.db+4*gsm->step/8+0);
////g41=*(gsm->data.db+4*gsm->step/8+1);
////g42=*(gsm->data.db+4*gsm->step/8+2);
////g43=*(gsm->data.db+4*gsm->step/8+3);
////g44=*(gsm->data.db+4*gsm->step/8+4);

////cout<<g00<<endl<<g01<<endl<<g02<<endl<<g10<<endl<<g11<<endl<<g12<<endl<<g20<<endl<<g21<<endl<<g22<<endl;
////unsigned char a00,a01,a02,a03,a04,a10,a11,a12,a13,a14,a20,a21,a22,a23,a24,a30,a31,a32,a33,a34,a40,a41,a42,a43,a44;
////for (j=2;j<frame_gray->height-2;j++)
////{
////	for(i=2;i<frame_gray->width-2;i++)
////	{
////		a00=cvGet2D(frame_gray,j-2,i-2).val[0];
////		a01=cvGet2D(frame_gray,j-2,i-1).val[0];
////		a02=cvGet2D(frame_gray,j-2,i).val[0];
////		a03=cvGet2D(frame_gray,j-2,i+1).val[0];
////		a04=cvGet2D(frame_gray,j-2,i+2).val[0];
////		
////		a10=cvGet2D(frame_gray,j-1,i-2).val[0];
////		a11=cvGet2D(frame_gray,j-1,i-1).val[0];
////		a12=cvGet2D(frame_gray,j-1,i).val[0];
////		a13=cvGet2D(frame_gray,j-1,i+1).val[0];
////		a14=cvGet2D(frame_gray,j-1,i+2).val[0];

////		a20=cvGet2D(frame_gray,j,i-2).val[0];
////		a21=cvGet2D(frame_gray,j,i-1).val[0];
////		a22=cvGet2D(frame_gray,j,i).val[0];
////		a23=cvGet2D(frame_gray,j,i+1).val[0];
////		a24=cvGet2D(frame_gray,j,i+2).val[0];

////		a30=cvGet2D(frame_gray,j+1,i-2).val[0];
////		a31=cvGet2D(frame_gray,j+1,i-1).val[0];
////		a32=cvGet2D(frame_gray,j+1,i).val[0];
////		a33=cvGet2D(frame_gray,j+1,i+1).val[0];
////		a34=cvGet2D(frame_gray,j+1,i+2).val[0];

////		a40=cvGet2D(frame_gray,j+2,i-2).val[0];
////		a41=cvGet2D(frame_gray,j+2,i-1).val[0];
////		a42=cvGet2D(frame_gray,j+2,i).val[0];
////		a43=cvGet2D(frame_gray,j+2,i+1).val[0];
////		a44=cvGet2D(frame_gray,j+2,i+2).val[0];

////		double ux=(a00*g00+a01*g01+a02*g02+a03*g03+a04*g04)+(a10*g10+a11*g11+a12*g12+a13*g13+a14*g14)+(a20*g20+a21*g21+a22*g22+a23*g23+a24*g24)+(a30*g30+a31*g31+a32*g32+a33*g33+a34*g34)+(a40*g40+a41*g41+a42*g42+a43*g43+a44*g44);
////		CvScalar cencor;
////		cencor.val[0]=ux;
////		cvSet2D(frame_gussian,j,i,cencor);

////	}

////}
//CvMat *bfgu=cvCreateMat(frame_gray->height,frame_gray->width,CV_64FC1);
//cvConvert(frame_gray,bfgu);
//bfgu=gussianmat(bfgu,gsm);
//cvConvert(bfgu,frame_gussian);

//cvNamedWindow("gray",0);
//cvNamedWindow("gussian",0);
//cvShowImage("gray",frame_gray);
//cvShowImage("gussian",frame_gussian);

//cvWaitKey(0);
////==================================================
////===========图像金字塔================================
////====隔点采样=======
////行列变为原来一半 清晰度减为四分之一
//IplImage *bfimg= cvLoadImage("test8.jpg",0);
//
//
//cvNamedWindow("GDgetmodelBF",0);
//cvShowImage("GDgetmodelBF",bfimg);
//
//CvMat* bf=cvCreateMat(bfimg->height,bfimg->width,CV_64FC1);
//cvConvert(bfimg,bf);
////bf=GDGM(GDGM(bf));
//
//CvMat* af=cvCreateMat(bfimg->height/2,bfimg->width/2,CV_64FC1);
//af=GDGM(bf);
//IplImage *afimg=cvCreateImage(cvGetSize(af),IPL_DEPTH_8U,1);
//cvConvert(af,afimg);
//cvNamedWindow("GDgetmodel",CV_WINDOW_AUTOSIZE);
//cvShowImage("GDgetmodel",afimg);
//cvWaitKey(0);




	//====================
	/*IplImage *imgx=cvCreateImage(cvGetSize(dogdelta.oct[0].inter[2].Gussianimg),IPL_DEPTH_8U,1);
	IplImage *t;
	CvMat *k=dogdelta.oct[0].inter[2].Gussianimg;
	t=cvCreateImage(cvGetSize(k),IPL_DEPTH_8U,1); 
	cvConvert(k,t);
	cvCanny(t,imgx,100,100,3);
	cvNamedWindow("ccc",0);
	cvShowImage("ccc",imgx);
	cvWaitKey(0);*/
