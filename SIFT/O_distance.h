#ifndef O_distance
#define O_distance

#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <math.h>
#include <iostream>
#include "MatOUT.h"
#include "SIFT.h"

#include "histo_cir.h"
using namespace cv;
using namespace std;

//F1为匹配模板 F2为被匹配图像
CvMat *O_dis(Feature *F1,Feature *F2,int F1_oct,int F2_oct)
{
	int i,j,k,l;
	int F1_num=F1->oct_imgFeature[F1_oct].Feature_num;
	int F2_num=F2->oct_imgFeature[F2_oct].Feature_num;

	//O_dist存放计算好的欧式距离 行为F1的每个特征点 列为每个F1特征点对于与F2特征点的欧式距离
	CvMat *O_dist=cvCreateMat(F1_num,F2_num,CV_64FC1);
	double O=0;
	double temp=0;

	for (i=0;i<F1_num;i++)
	{
		for (j=0;j<F2_num;j++)
		{
			temp=0;
			CvMat *p1=F1->oct_imgFeature[F1_oct].p_Feature[i].Feature_Vector;
			CvMat *p2=F2->oct_imgFeature[F2_oct].p_Feature[j].Feature_Vector;

			for (k=0;k<128;k++)
			{
				double lp1=*(p1->data.db+k);
				double lp2=*(p2->data.db+k);
				if((lp1-lp2)>0.01)
				{
					temp+=pow(lp1-lp2,2);
				}
				
			}

			O=sqrt(temp);
			*(O_dist->data.db+i*O_dist->step/8+j)=O;

		}
	}
	return O_dist;


}





CvMat *O_dis_fix(Feature *F1,Feature *F2,int F1_oct,int F2_oct)
{
	int i,j,k,l;
	int F1_num=F1->oct_imgFeature[F1_oct].Feature_num;
	int F2_num=F2->oct_imgFeature[F2_oct].Feature_num;

	//O_dist存放计算好的欧式距离 行为F1的每个特征点 列为每个F1特征点对于与F2特征点的欧式距离
	CvMat *O_dist=cvCreateMat(F1_num,F2_num,CV_64FC1);
	double O=0;
	double temp=0;

	for (i=0;i<F1_num;i++)
	{
		for (j=0;j<F2_num;j++)
		{
			temp=0;
			CvMat *p1=F1->oct_imgFeature[F1_oct].p_Feature[i].Feature_Vector;
			CvMat *p2=F2->oct_imgFeature[F2_oct].p_Feature[j].Feature_Vector;

			for (k=0;k<128;k++)
			{
				double lp1=*(p1->data.db+k);
				double lp2=*(p2->data.db+k);
				temp+=pow(lp1-lp2,2.0);

			}

			O=sqrt(temp);
			//O=temp;
			*(O_dist->data.db+i*O_dist->step/8+j)=O;

		}
	}
	return O_dist;


}





















#endif