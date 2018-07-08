#ifndef SearchOUT
#define SearchOUT
#include "GUSSIANdelta.h"
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
using namespace cv;
using namespace std;

//void search_OUT(delta dogdelta,CvMat *allsearch,int *num,delta Delta)
//{
//	int xoct=dogdelta.longoct;
//	int xinte=dogdelta.longinte;
//	int p=0;
//	for (int i=0;i<xoct;i++)
//	{
//		p=0;
//		IplImage *searchshow;
//		searchshow=cvCreateImage(cvGetSize(Delta.oct[i].inter[0].Gussianimg),IPL_DEPTH_8U,1);
//		cvConvert(Delta.oct[i].inter[1].Gussianimg,searchshow);
//		Mat showse=cvarrToMat(searchshow,true);
//		cout<<"长"<<searchshow->width<<endl<<"宽"<<searchshow->height<<endl;
// 		for (int j=0;j<*(num+i)*2;j++)
//		{
//			int x=cvRound(*(allsearch->data.db+2*i*allsearch->step/8+j));
//			int y=cvRound(*(allsearch->data.db+(2*i+1)*allsearch->step/8+j));
//			
//
//
//
//			//cvCircle(searchshow,cvPoint(x,y),2,CV_RGB(0,0,0));
//  			circle(showse,cvPoint(y,x),2,CV_RGB(0,0,0),2,8,0);
//			p++;
//
//
//
//
//
//		}
//		char namewin[20];
// 		sprintf_s(namewin,"%d的窗口",i);
//		/*cvNamedWindow(namewin,0);
//		cvShowImage(namewin,searchshow);*/
//		namedWindow(namewin,0);
//		imshow(namewin,showse);
//		cout<<"p"<<p<<endl<<"num"<<*(num+i)<<endl;
//		cvWaitKey(0);
//	}
//
//}




void search_OUT(delta dogdelta,CvMat *allsearch,int *num,delta Delta,CvMat *spoint_gradsita,char name[])
{
	int xoct=dogdelta.longoct;
	int xinte=dogdelta.longinte;
	int p=0;
	for (int i=0;i<xoct;i++)
	{
		p=0;
		IplImage *searchshow;
		searchshow=cvCreateImage(cvGetSize(Delta.oct[i].inter[0].Gussianimg),IPL_DEPTH_8U,1);
		cvConvert(Delta.oct[i].inter[1].Gussianimg,searchshow);
		Mat showse=cvarrToMat(searchshow,true);
		cout<<"长"<<searchshow->width<<endl<<"宽"<<searchshow->height<<endl;
 		for (int j=0;j<*(num+i)*2;j++)
		{
			int x=cvRound(*(allsearch->data.db+2*i*allsearch->step/8+j));
			int y=cvRound(*(allsearch->data.db+(2*i+1)*allsearch->step/8+j));
			

			//cvCircle(searchshow,cvPoint(x,y),2,CV_RGB(0,0,0));
  			circle(showse,cvPoint(y,x),2,CV_RGB(0,0,0),2,8,0);
			p++;


			double sita=*(spoint_gradsita->data.db+i*spoint_gradsita->step/8+j);
			//double mod=*(spoint_gradsita->data.db+(2*i+1)*spoint_gradsita->step/8+j);
			int mod=25;
			int smod=sin(sita*3.14/180)*mod;
			int cmod=cos(sita*3.14/180)*mod;
			
			line(showse,cvPoint(y,x),cvPoint(y+smod,x+cmod),CV_RGB(0,0,0));  //============


		}
		char namewin[20];
 		sprintf_s(namewin,"%d的窗口",i);
		strcat(name,namewin);
		/*cvNamedWindow(namewin,0);
		cvShowImage(namewin,searchshow);*/
		namedWindow(name,0);
		imshow(name,showse);
		cout<<"p"<<p<<endl<<"num"<<*(num+i)<<endl;
		cvWaitKey(0);
	}

}



void search_O(delta dogdelta,CvMat *allsearch,int *num,delta Delta)
{
	int xoct=dogdelta.longoct;
	int xinte=dogdelta.longinte;
	int p=0;
	for (int i=0;i<xoct;i++)
	{
		p=0;
		IplImage *searchshow;
		searchshow=cvCreateImage(cvGetSize(Delta.oct[i].inter[0].Gussianimg),IPL_DEPTH_8U,1);
		cvConvert(Delta.oct[i].inter[1].Gussianimg,searchshow);
		Mat showse=cvarrToMat(searchshow,true);
		cout<<"长"<<searchshow->width<<endl<<"宽"<<searchshow->height<<endl;
 		for (int j=0;j<*(num+i)*2;j++)
		{
			int x=cvRound(*(allsearch->data.db+2*i*allsearch->step/8+j));
			int y=cvRound(*(allsearch->data.db+(2*i+1)*allsearch->step/8+j));
			

			//cvCircle(searchshow,cvPoint(x,y),2,CV_RGB(0,0,0));
  			circle(showse,cvPoint(y,x),2,CV_RGB(0,0,0),2,8,0);
			p++;



		}
		char namewin[20];
 		sprintf_s(namewin,"%d的窗口",i);
	
		/*cvNamedWindow(namewin,0);
		cvShowImage(namewin,searchshow);*/
		namedWindow(namewin,0);
		imshow(namewin,showse);
		cout<<"p"<<p<<endl<<"num"<<*(num+i)<<endl;
		cvWaitKey(0);
	}

}


#endif