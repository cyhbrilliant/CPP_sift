#ifndef GUSSIANmodel
#define GUSSIANmodel

#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <iostream>
#include "MatOUT.h"

using namespace cv;
using namespace std;

CvMat *GussianModel(int row,double stdless)
{
	double x,y;
	double gs;
	int arow=row/2;
	CvMat *gsm=cvCreateMat(row,row,CV_64FC1);
	for (int i=0;i<arow+1;i++)
	{
		for(int j=arow;j<row;j++)
		{
			/*x=j-2.0;
			y=2-i;*/
			gs=(1/(2*3.14*stdless*stdless))*(pow(2.72,-((j-arow)*(j-arow)+(i-arow)*(i-arow))/(2*stdless*stdless)));
			cout<<gs<<" "<<endl;
			*(gsm->data.db+i*gsm->step/8+j)=gs;

		}
	}
	for (int i=0;i<arow+1;i++)
	{
		for (int j=0;j<arow;j++)
		{
			*(gsm->data.db+i*gsm->step/8+j)=*(gsm->data.db+i*gsm->step/8+(4-j));
		}
	}
	for (int i=arow+1;i<row;i++)
	{
		for (int j=0;j<arow;j++)
		{
			*(gsm->data.db+i*gsm->step/8+j)=*(gsm->data.db+j*gsm->step/8+i);
		}
	}
	for (int i=arow;i<row;i++)
	{
		for (int j=arow;j<row;j++)
		{
			*(gsm->data.db+i*gsm->step/8+j)=*(gsm->data.db+(4-i)*gsm->step/8+j);
		}
	}
	cout<<"未归一化之前高斯模板----------------------------------------"<<endl;
	MatOUT(gsm,row,row,1);
	cout<<endl;
	return gsm;

}


CvMat *Gussian2one(CvMat *gsm,int row,int col)
{
	

	double ad=0;
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<col;j++)
		{
			ad+=(*(gsm->data.db+i*gsm->step/8+j));
		}
	}
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<col;j++)
		{
			double k=(*(gsm->data.db+i*gsm->step/8+j));
			(*(gsm->data.db+i*gsm->step/8+j))=k/ad;
		}
	}
	cout<<"归一化之后的高斯模板----------------------------------------"<<endl;
	MatOUT(gsm,row,row,1);
	cout<<endl;
	return gsm;


}

CvMat *gussianmat(CvMat *bfgussian,CvMat *gsm)//只支持模板为5*5的高斯模板
{
	int i,j;
	//thirtysix_bin[0]=(thirtysix_bin[0]*6+thirtysix_bin[1]*4+thirtysix_bin[2])/11;
	//thirtysix_bin[1]=(thirtysix_bin[0]*4+thirtysix_bin[1]*6+thirtysix_bin[2]*4+thirtysix_bin[3])/15;
	//for (k=2;k<34;k++)
	//{
	//	thirtysix_bin[k]=(thirtysix_bin[k-2]+thirtysix_bin[k-1]*4+thirtysix_bin[k]*6+thirtysix_bin[k+1]*4+thirtysix_bin[k+2])/16;
	//}
	//thirtysix_bin[34]=(thirtysix_bin[35]*4+thirtysix_bin[34]*6+thirtysix_bin[33]*4+thirtysix_bin[32])/15;
	//thirtysix_bin[35]=(thirtysix_bin[35]*6+thirtysix_bin[34]*4+thirtysix_bin[33])/11;
	double g00,g01,g02,g03,g04,g10,g11,g12,g13,g14,g20,g21,g22,g23,g24,g30,g31,g32,g33,g34,g40,g41,g42,g43,g44;
	g00=*(gsm->data.db+0*gsm->step/8+0);
	g01=*(gsm->data.db+0*gsm->step/8+1);
	g02=*(gsm->data.db+0*gsm->step/8+2);
	g03=*(gsm->data.db+0*gsm->step/8+3);
	g04=*(gsm->data.db+0*gsm->step/8+4);

	g10=*(gsm->data.db+1*gsm->step/8+0);
	g11=*(gsm->data.db+1*gsm->step/8+1);
	g12=*(gsm->data.db+1*gsm->step/8+2);
	g13=*(gsm->data.db+1*gsm->step/8+3);
	g14=*(gsm->data.db+1*gsm->step/8+4);

	g20=*(gsm->data.db+2*gsm->step/8+0);
	g21=*(gsm->data.db+2*gsm->step/8+1);
	g22=*(gsm->data.db+2*gsm->step/8+2);
	g23=*(gsm->data.db+2*gsm->step/8+3);
	g24=*(gsm->data.db+2*gsm->step/8+4);

	g30=*(gsm->data.db+3*gsm->step/8+0);
	g31=*(gsm->data.db+3*gsm->step/8+1);
	g32=*(gsm->data.db+3*gsm->step/8+2);
	g33=*(gsm->data.db+3*gsm->step/8+3);
	g34=*(gsm->data.db+3*gsm->step/8+4);

	g40=*(gsm->data.db+4*gsm->step/8+0);
	g41=*(gsm->data.db+4*gsm->step/8+1);
	g42=*(gsm->data.db+4*gsm->step/8+2);
	g43=*(gsm->data.db+4*gsm->step/8+3);
	g44=*(gsm->data.db+4*gsm->step/8+4);

	IplImage *frame_gray=cvCreateImage(cvGetSize(bfgussian),IPL_DEPTH_8U,1);
	cvConvert(bfgussian,frame_gray);
	IplImage *frame_gussian=cvCreateImage(cvGetSize(bfgussian),IPL_DEPTH_8U,1);
	unsigned char a00,a01,a02,a03,a04,a10,a11,a12,a13,a14,a20,a21,a22,a23,a24,a30,a31,a32,a33,a34,a40,a41,a42,a43,a44;
	//double  a00,a01,a02,a03,a04,a10,a11,a12,a13,a14,a20,a21,a22,a23,a24,a30,a31,a32,a33,a34,a40,a41,a42,a43,a44;
	for (j=2;j<frame_gray->height-2;j++)
	{
		for(i=2;i<frame_gray->width-2;i++)
		{
			a00=cvGet2D(frame_gray,j-2,i-2).val[0];
			a01=cvGet2D(frame_gray,j-2,i-1).val[0];
			a02=cvGet2D(frame_gray,j-2,i).val[0];
			a03=cvGet2D(frame_gray,j-2,i+1).val[0];
			a04=cvGet2D(frame_gray,j-2,i+2).val[0];

			a10=cvGet2D(frame_gray,j-1,i-2).val[0];
			a11=cvGet2D(frame_gray,j-1,i-1).val[0];
			a12=cvGet2D(frame_gray,j-1,i).val[0];
			a13=cvGet2D(frame_gray,j-1,i+1).val[0];
			a14=cvGet2D(frame_gray,j-1,i+2).val[0];

			a20=cvGet2D(frame_gray,j,i-2).val[0];
			a21=cvGet2D(frame_gray,j,i-1).val[0];
			a22=cvGet2D(frame_gray,j,i).val[0];
			a23=cvGet2D(frame_gray,j,i+1).val[0];
			a24=cvGet2D(frame_gray,j,i+2).val[0];

			a30=cvGet2D(frame_gray,j+1,i-2).val[0];
			a31=cvGet2D(frame_gray,j+1,i-1).val[0];
			a32=cvGet2D(frame_gray,j+1,i).val[0];
			a33=cvGet2D(frame_gray,j+1,i+1).val[0];
			a34=cvGet2D(frame_gray,j+1,i+2).val[0];

			a40=cvGet2D(frame_gray,j+2,i-2).val[0];
			a41=cvGet2D(frame_gray,j+2,i-1).val[0];
			a42=cvGet2D(frame_gray,j+2,i).val[0];
			a43=cvGet2D(frame_gray,j+2,i+1).val[0];
			a44=cvGet2D(frame_gray,j+2,i+2).val[0];

			double ux=(a00*g00+a01*g01+a02*g02+a03*g03+a04*g04)+(a10*g10+a11*g11+a12*g12+a13*g13+a14*g14)+(a20*g20+a21*g21+a22*g22+a23*g23+a24*g24)+(a30*g30+a31*g31+a32*g32+a33*g33+a34*g34)+(a40*g40+a41*g41+a42*g42+a43*g43+a44*g44);
			CvScalar cencor;
			cencor.val[0]=ux;
			cvSet2D(frame_gussian,j,i,cencor);
		}
	}
	//===========消除边缘影响============
	for (j=0;j<2;j++)
	{
		for (i=0;i<frame_gray->width;i++)
		{
			unsigned char a;
			a=cvGet2D(frame_gray,2,i).val[0];
			CvScalar cen;
			cen.val[0]=(double)a;
			cvSet2D(frame_gussian,j,i,cen);
		}
	}
	for (j=frame_gray->height-3;j<frame_gray->height;j++)
	{
		for (i=0;i<frame_gray->width;i++)
		{
			unsigned char a;
			a=cvGet2D(frame_gray,frame_gray->height-4,i).val[0];
			CvScalar cen;
			cen.val[0]=(double)a;
			cvSet2D(frame_gussian,j,i,cen);
		}
	}
	for (j=0;j<frame_gray->height;j++)
	{
		for (i=0;i<2;i++)
		{
			unsigned char a;
			a=cvGet2D(frame_gray,j,2).val[0];
			CvScalar cen;
			cen.val[0]=(double)a;
			cvSet2D(frame_gussian,j,i,cen);
		}
	}

	for (j=0;j<frame_gray->height;j++)
	{
		for (i=frame_gray->width-3;i<frame_gray->width;i++)
		{
			unsigned char a;
			a=cvGet2D(frame_gray,j,frame_gray->width-4).val[0];
			CvScalar cen;
			cen.val[0]=(double)a;
			cvSet2D(frame_gussian,j,i,cen);
		}
	}
	//===============================




	CvMat *afgussian=cvCreateMat(bfgussian->rows,bfgussian->cols,CV_64FC1);
	cvConvert(frame_gussian,afgussian);
	return afgussian;

}







#endif