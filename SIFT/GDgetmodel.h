#ifndef GDgetmodel
#define GDgetmodel
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"




using namespace cv;
using namespace std;


CvMat *GDGM(CvMat *bf)
{
	CvMat *af=cvCreateMat(bf->rows/2,bf->cols/2,CV_64FC1);
	int m,n;
	for (int i=0;i<bf->rows-1;i=i+2)
	{
		for (int j=0;j<bf->cols-1;j=j+2)
		{
			m=i/2;
			n=j/2;
			*(af->data.db+m*af->step/8+n)=*(bf->data.db+i*bf->step/8+j);
		}
	}
	return af;


}









#endif