#ifndef mul_img
#define  mul_img


#include "cv.h"
#include "highgui.h"
#include <math.h>
#include <iostream>

int max(int x,int y)
{
	if (x>y)
	{
		return x;
	} 
	else
	{
		return y;
	}
}


IplImage* mulimg(char n1[],char n2[],int row1,int col1,int row2,int col2)
{
	IplImage *img1,*img2,*dst1,*dst2,*dst_big; //img1 img2 原图 dst1、dst2放缩后的图 dst_big 大图
	CvRect rect1=cvRect(0,0,col1,row1); //两个ROI区域
	CvRect rect2=cvRect(col1,0,col1+col2,row2);
	img1=cvLoadImage(n1);
	img2=cvLoadImage(n2);
	dst1=cvCreateImage(cvSize(col1,row1),img1->depth,3);
	dst2=cvCreateImage(cvSize(col2,row2),img2->depth,3);
	dst_big=cvCreateImage(cvSize(col1+col2,max(row1,row2)),img2->depth,3);
	cvResize(img1,dst1); //放缩
	cvResize(img2,dst2);
	cvSetImageROI(dst_big,rect1); //设置ROI
	cvCopy(dst1,dst_big);
	cvSetImageROI(dst_big,rect2);
	cvCopy(dst2,dst_big);
	cvResetImageROI(dst_big);//释放ROI

	return dst_big;
}

#endif