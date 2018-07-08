#ifndef cvmat_mat
#define cvmat_mat
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"
using namespace cv;
using namespace std;

Mat CvmatToMat(CvMat *cvmat)
{
	int i,j;
	
	int xrow=cvmat->rows;
	int xcol=cvmat->cols;
	Mat mat(cvGetSize(cvmat),CV_64FC1);
	int matchannal=mat.channels();
	for (i=0;i<xrow;i++)
	{
		for (j=0;j<xcol;j++)
		{
			
 			*(mat.data+i*mat.step/8+j)=*(cvmat->data.db+i*cvmat->step/8+j);
		}
	}
	return mat;
}


CvMat *MatToCvmat(Mat mat)
{
	int xrow=mat.rows;
	int xcol=mat.cols;
	CvMat *cvmat=cvCreateMat(mat.rows,mat.cols,CV_64FC1);
	int matchannal=mat.channels();
	int i,j;
	for (i=0;i<xrow;i++)
	{
		for (j=0;j<xcol;j++)
		{

			*(cvmat->data.db+i*cvmat->step/8+j)=*(mat.data+i*mat.step/8+j);
		}
	}

	return cvmat;

}





#endif