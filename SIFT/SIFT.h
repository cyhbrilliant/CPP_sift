#ifndef SIFT
#define SIFT


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
#include "O_distance.h"

using namespace cv;
using namespace std;

Feature *sift(IplImage *bfimg,char name[],CvMat *position,int oct,int R_BasicDrection,int R_FeatureGRAD1)
{




	//==============搜索点时row为x col为y  opencv默认col为x row为y

	



	//图像转换成矩阵
	CvMat* bf=cvCreateMat(bfimg->height,bfimg->width,CV_64FC1);
	cvConvert(bfimg,bf);



	//创建高斯金字塔
	delta Delta;
	Delta=CreateDelta(bf,6);      //未升采样
	//Delta=upCreateDelta(bf,6); //升采样




	//显示高斯金字塔
	  // Delta_OUT(Delta);




	//高斯金字塔差分为DOG金字塔
	delta dogdelta;
	dogdelta=DOG(Delta);



	//显示DOG金字塔
		Delta_OUT(dogdelta);
//checkfix OK===========================================================


	//创建num为每个尺度空间中特征点个数的栈数组
	int *num=(int*)malloc(sizeof(int)*dogdelta.longoct);



	//创建allsearch为每个尺度空间中的特征点
	CvMat *allsearch=cvCreateMat(2*dogdelta.longoct,100000,CV_64FC1);



	//SearchP函数为找金字塔找点函数
	allsearch=SearchP(dogdelta,num);



	//点的修正与精确
	FixPoint(allsearch,num,dogdelta);

//	search_O(dogdelta,allsearch,num,Delta);

	//====================================================================================

	//创建矩阵存放每个特征点的初始旋转参数 sita mod
	CvMat *spoint_gradsita=cvCreateMat(allsearch->rows/2,allsearch->cols,CV_64FC1);



	////特征点初始旋转参数
//	spoint_gradsita=hi(Delta,allsearch,num);

	//===========================================================================================
	//返回主函数特征点坐标
	int i,j;
	for (i=0;i<*(num+oct);i++)
	{
		*(position->data.db+0*position->step/8+i)=*(allsearch->data.db+2*oct*allsearch->step/8+i);
		*(position->data.db+1*position->step/8+i)=*(allsearch->data.db+(2*oct+1)*allsearch->step/8+i);
	}


	//search_OUT(dogdelta,allsearch,num,Delta,spoint_gradsita,name);
	//128维特征向量
	Feature *F1;
	F1=histmat(Delta,allsearch,num,name);
//	F1=calFeature_fix(Delta,allsearch,num,R_BasicDrection,R_FeatureGRAD1);
//	F1=calFeature_fix(Delta,allsearch,num,R_BasicDrection,R_FeatureGRAD1);

	//输出点在图上的位置




	free(num);



	return F1;
}


#endif