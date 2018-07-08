
#ifndef DeltaOUT
#define DeltaOUT
#include "GUSSIANdelta.h"
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
using namespace cv;
using namespace std;


void Delta_OUT(delta Delta)
{

	char namewin[20];
	for (int i=0;i<Delta.longoct;i++)
	{
		for(int j=0;j<Delta.longinte;j++)
		{
			sprintf_s(namewin,"%d%dµÄ´°¿Ú",i,j);
			cvNamedWindow(namewin,CV_WINDOW_AUTOSIZE);
			IplImage *t;
			CvMat *k=Delta.oct[i].inter[j].Gussianimg;
			t=cvCreateImage(cvGetSize(k),IPL_DEPTH_8U,1); 
			cvConvert(k,t);
 			cvShowImage(namewin,t);
			cvWaitKey(0);
		}
	}

}




#endif