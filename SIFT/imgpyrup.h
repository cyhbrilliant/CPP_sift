#ifndef imgpyrup
#define imgpyrup

#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <iostream>
#include "MatOUT.h"

using namespace cv;
using namespace std;

//========双线性插值升采样======
CvMat *Img_pyrup(CvMat *bfmat)//双线性插值
{
	int xrow=bfmat->rows;
	int xcol=bfmat->cols;
	CvMat *afmat=cvCreateMat(2*xrow,2*xcol,CV_64FC1);
	int i,j;
	for (i=0;i<xrow*2;i=i+2)//像素插入
	{
		for (j=0;j<xcol*2;j=j+2)
		{
			*(afmat->data.db+i*afmat->step/8+j)=*(bfmat->data.db+(i/2)*bfmat->step/8+(j/2));

		}
	}

	//====源像素斜四点中心像素计算===
	for(i=1;i<xrow*2-1;i=i+2)
	{
		for (j=1;j<xcol*2-1;j=j+2)
		{
			double allcen=*(afmat->data.db+(i-1)*afmat->step/8+j-1)+*(afmat->data.db+(i-1)*afmat->step/8+j+1)+*(afmat->data.db+(i+1)*afmat->step/8+j-1)+*(afmat->data.db+(i+1)*afmat->step/8+j+1);
			*(afmat->data.db+i*afmat->step/8+j)=(double)(allcen/4);
		}
	}

	//====混合像素正四点中心像素计算===
	for (i=1;i<xrow*2-2;i=i+2)
	{
		for (j=2;j<xcol*2-3;j=j+2)
		{
			double fallcen=*(afmat->data.db+(i-1)*afmat->step/8+j)+*(afmat->data.db+i*afmat->step/8+j-1)+*(afmat->data.db+i*afmat->step/8+j+1)+*(afmat->data.db+(i+1)*afmat->step/8+j);
			*(afmat->data.db+i*afmat->step/8+j)=(double)(fallcen/4);
		}
	}
	for (i=2;i<xrow*2-3;i=i+2)
	{
		for (j=1;j<xcol*2-2;j=j+2)
		{
			double fallcen_=*(afmat->data.db+(i-1)*afmat->step/8+j)+*(afmat->data.db+i*afmat->step/8+j-1)+*(afmat->data.db+i*afmat->step/8+j+1)+*(afmat->data.db+(i+1)*afmat->step/8+j);
			*(afmat->data.db+i*afmat->step/8+j)=(double)(fallcen_/4);
		}
	}
	//====边混合像素三点像素计算===
	for (i=1;i<xrow*2-2;i=i+2)
	{
		double lcen=*(afmat->data.db+(i-1)*afmat->step/8)+*(afmat->data.db+(i+1)*afmat->step/8)+*(afmat->data.db+i*afmat->step/8+1);
		*(afmat->data.db+i*afmat->step/8)=(double)(lcen/3);

		double rcen=*(afmat->data.db+(i-1)*afmat->step/8+xcol*2-2)+*(afmat->data.db+(i+1)*afmat->step/8+xcol*2-2)+*(afmat->data.db+i*afmat->step/8+xcol*2-3);
		*(afmat->data.db+i*afmat->step/8+xcol*2-2)=(double)(rcen/3);
	}

	for (j=1;j<xcol*2-2;j=j+2)
	{
		double ucen=*(afmat->data.db+j-1)+*(afmat->data.db+j+1)+*(afmat->data.db+1*afmat->step/8+j);
		*(afmat->data.db+j)=(double)(ucen/3);

		double dcen=*(afmat->data.db+(2*xrow-2)*afmat->step/8+j-1)+*(afmat->data.db+(2*xrow-2)*afmat->step/8+j+1)+*(afmat->data.db+(2*xrow-3)*afmat->step/8+j);
		*(afmat->data.db+j)=(double)(dcen/3);

	}
	//====两底像素赋值===
	for (j=0;j<2*xcol-1;j++)
	{
		*(afmat->data.db+(xrow*2-1)*afmat->step/8+j)=*(afmat->data.db+(xrow*2-4)*afmat->step/8+j);

	}
	for (i=0;i<xrow*2;i++)
	{
		*(afmat->data.db+i*afmat->step/8+2*xcol-1)=*(afmat->data.db+i*afmat->step/8+2*xcol-4);
	}

	return afmat;
}

/*
	*(afmat->data.db+(xrow*2-1)*afmat->step/8+0)=(*(afmat->data.db+(xrow*2-2)*afmat->step/8+0)+*(afmat->data.db+(xrow*2-1)*afmat->step/8+1))/2;
	*(afmat->data.db+(xrow*2-1)*afmat->step/8+0)*/

















#endif