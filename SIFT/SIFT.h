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




	//==============������ʱrowΪx colΪy  opencvĬ��colΪx rowΪy

	



	//ͼ��ת���ɾ���
	CvMat* bf=cvCreateMat(bfimg->height,bfimg->width,CV_64FC1);
	cvConvert(bfimg,bf);



	//������˹������
	delta Delta;
	Delta=CreateDelta(bf,6);      //δ������
	//Delta=upCreateDelta(bf,6); //������




	//��ʾ��˹������
	  // Delta_OUT(Delta);




	//��˹���������ΪDOG������
	delta dogdelta;
	dogdelta=DOG(Delta);



	//��ʾDOG������
		Delta_OUT(dogdelta);
//checkfix OK===========================================================


	//����numΪÿ���߶ȿռ��������������ջ����
	int *num=(int*)malloc(sizeof(int)*dogdelta.longoct);



	//����allsearchΪÿ���߶ȿռ��е�������
	CvMat *allsearch=cvCreateMat(2*dogdelta.longoct,100000,CV_64FC1);



	//SearchP����Ϊ�ҽ������ҵ㺯��
	allsearch=SearchP(dogdelta,num);



	//��������뾫ȷ
	FixPoint(allsearch,num,dogdelta);

//	search_O(dogdelta,allsearch,num,Delta);

	//====================================================================================

	//����������ÿ��������ĳ�ʼ��ת���� sita mod
	CvMat *spoint_gradsita=cvCreateMat(allsearch->rows/2,allsearch->cols,CV_64FC1);



	////�������ʼ��ת����
//	spoint_gradsita=hi(Delta,allsearch,num);

	//===========================================================================================
	//��������������������
	int i,j;
	for (i=0;i<*(num+oct);i++)
	{
		*(position->data.db+0*position->step/8+i)=*(allsearch->data.db+2*oct*allsearch->step/8+i);
		*(position->data.db+1*position->step/8+i)=*(allsearch->data.db+(2*oct+1)*allsearch->step/8+i);
	}


	//search_OUT(dogdelta,allsearch,num,Delta,spoint_gradsita,name);
	//128ά��������
	Feature *F1;
	F1=histmat(Delta,allsearch,num,name);
//	F1=calFeature_fix(Delta,allsearch,num,R_BasicDrection,R_FeatureGRAD1);
//	F1=calFeature_fix(Delta,allsearch,num,R_BasicDrection,R_FeatureGRAD1);

	//�������ͼ�ϵ�λ��




	free(num);



	return F1;
}


#endif