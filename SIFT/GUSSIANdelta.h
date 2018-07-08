
#ifndef GUSSIANdeltal
#define GUSSIANdeltal
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <iostream>
#include <math.h>
#include "GDgetmodel.h"
#include "GUSSIANmodel.h"
#include "cvmat_mat.h"
#include "imgpyrup.h"

using namespace cv;
using namespace std;

class interval
{
public:
	CvMat *Gussianimg;
};
class octave
{
public:
	interval inter[10];
	
};

class delta
{
public:
	octave oct[20];
	int longoct;
	int longinte;
};

//

void sub_imagex(CvMat *x,CvMat *y,CvMat *z)
{
	IplImage *image_first=cvCreateImage(cvGetSize(x),IPL_DEPTH_8U,1);
	IplImage *image_second=cvCreateImage(cvGetSize(x),IPL_DEPTH_8U,1);
	IplImage *ksc=cvCreateImage(cvGetSize(x),IPL_DEPTH_8U,1);
	cvConvert(x,image_first);
	cvConvert(y,image_second);
	cvConvert(z,ksc);
	for(int i = 0; i<image_first->height;i++)
		for(int j=0;j<image_first->width;j++)
		{
			int temp = (int)((uchar*)(image_first->imageData + image_first->width*i))[j];
			int temp_sec = (int)((uchar*)(image_second->imageData + image_second->width*i))[j];
			((uchar*)(ksc->imageData + ksc->width*i))[j] = temp - temp_sec;
		}

}


//==============未进行升采样====================
delta CreateDelta(CvMat *DeltaMat,int Ginterval)
{
	delta Delta;
	int longm;
	if (DeltaMat->width<DeltaMat->height)
	{
		longm=DeltaMat->width;
	} 
	else
	{
		longm=DeltaMat->height;
	}
	int k;
	k=log10(longm)/log10(2)/2-1;

	for (int i=0;i<k;i++)
	{
		for(int j=0;j<Ginterval;j++)
		{
			Delta.oct[i].inter[j].Gussianimg=cvCreateMat(DeltaMat->height/(pow(2,i)),DeltaMat->width/(pow(2,i)),CV_64FC1);
		}
	}

//	CvMat *temp=cvCreateMat(DeltaMat->height,DeltaMat->width,CV_64FC1);
	int x=0,y=0;
	int m,n;
//	double v=2.0;
	double v=pow(2.0,1.0/((double)Ginterval-3.0));
//	double v=pow(2,1.0/(Ginterval));
	for (m=0;m<k;m++)
	{
		
		if (m==0)
		{
  			//temp=cvCloneMat(DeltaMat);

			for(n=0;n<Ginterval;n++)
			{

				if (n==0)
				{
					Delta.oct[m].inter[0].Gussianimg=cvCloneMat(DeltaMat);
				}
				else
				{
					
					CvMat *gsm=cvCreateMat(5,5,CV_64FC1);
					double sigmapre=1.6*pow(v,n-1);
					double sigmanow=sigmapre*v;
				//	double sigma=1.6*pow(v,n);
			//		double sigma=sqrt(1.6*pow(v,n)*1.6*pow(v,n)*(v*v-1));
					double sigma=sqrt(sigmanow*sigmanow-sigmapre*sigmapre);
				//	cout<<sigma<<endl;
			//		double sigma=1.6*pow(v,n)*(v-1);
			//		gsm=GussianModel(5,sigma);
			//		gsm=Gussian2one(gsm,5,5);
			//		Delta.oct[m].inter[n].Gussianimg=gussianmat(Delta.oct[m].inter[n-1].Gussianimg,gsm);
					IplImage *a=cvCreateImage(cvGetSize(Delta.oct[m].inter[n-1].Gussianimg),IPL_DEPTH_8U,1);
					IplImage *b=cvCreateImage(cvGetSize(Delta.oct[m].inter[n].Gussianimg),IPL_DEPTH_8U,1);
					cvConvert(Delta.oct[m].inter[n-1].Gussianimg,a);
					cvSmooth(a,b,CV_GAUSSIAN,5,5,sqrt(sigma),sqrt(sigma));
					cvConvert(b,Delta.oct[m].inter[n].Gussianimg);
				//cvSmooth(Delta.oct[m].inter[n-1].Gussianimg,Delta.oct[m].inter[n].Gussianimg,CV_GAUSSIAN,5,5,(sigma),(sigma));
				
					/*Mat src=cvarrToMat(Delta.oct[m].inter[n-1].Gussianimg,true);
					Mat dst=cvarrToMat(Delta.oct[m].inter[n].Gussianimg,false);
					GaussianBlur(src,dst,CvSize(5,5),sqrt(sigma),0);
					CvMat tempcopy=dst;
					cvCopy(&tempcopy,Delta.oct[m].inter[n].Gussianimg);*/
					cvReleaseImage(&a);
					cvReleaseImage(&b);
					}
			}


		
		}
		else
		{
			for(n=0;n<Ginterval;n++)
			{

				if (n==0)
				{
					cout<<"行 "<<Delta.oct[m].inter[0].Gussianimg->rows<<"   列 "<<Delta.oct[m].inter[0].Gussianimg->cols<<endl;
					cout<<"行 "<<Delta.oct[x].inter[y].Gussianimg->rows<<"   列 "<<Delta.oct[x].inter[y].Gussianimg->cols<<endl;
					cout<<"行 "<<GDGM(Delta.oct[x].inter[y].Gussianimg)->rows<<"   列 "<<GDGM(Delta.oct[x].inter[y].Gussianimg)->cols<<endl;
					Delta.oct[m].inter[0].Gussianimg=cvCloneMat(GDGM(Delta.oct[x].inter[y].Gussianimg));
				}
				else
				{
					CvMat *gsm=cvCreateMat(5,5,CV_64FC1);
					double sigmapre=1.6*pow(v,n-1);
					double sigmanow=sigmapre*v;
				//	double sigma=1.6*pow(v,n);
					//double sigma=sqrt(1.6*pow(v,n)*1.6*pow(v,n)*(v*v-1));
					double sigma=sqrt(sigmanow*sigmanow-sigmapre*sigmapre);
				//	cout<<sigma<<endl;
				//	double sigma=1.6*pow(v,n)*(v-1);
				//	gsm=GussianModel(5,sigma);
				//	gsm=Gussian2one(gsm,5,5);
				//	Delta.oct[m].inter[n].Gussianimg=gussianmat(Delta.oct[m].inter[n-1].Gussianimg,gsm);
					IplImage *a=cvCreateImage(cvGetSize(Delta.oct[m].inter[n-1].Gussianimg),IPL_DEPTH_8U,1);
					IplImage *b=cvCreateImage(cvGetSize(Delta.oct[m].inter[n].Gussianimg),IPL_DEPTH_8U,1);
					cvConvert(Delta.oct[m].inter[n-1].Gussianimg,a);
					cvSmooth(a,b,CV_GAUSSIAN,5,5,sqrt(sigma),sqrt(sigma));
					cvConvert(b,Delta.oct[m].inter[n].Gussianimg);

					cvReleaseImage(&a);
					cvReleaseImage(&b);
					//cvSmooth(Delta.oct[m].inter[n-1].Gussianimg,Delta.oct[m].inter[n].Gussianimg,CV_GAUSSIAN,5,5,(sigma),(sigma));
				}

			}
		}
		
		x=m;
		y=Ginterval-3;
//		temp=GDGM(Delta.oct[m].inter[Ginterval-3].Gussianimg);
	

	}

	Delta.longinte=Ginterval;
	Delta.longoct=k;

	return Delta;



}





//===================进行升采样=====================
delta upCreateDelta(CvMat *DeltaMat,int Ginterval)
{
	delta Delta;
	int longm;
	if (DeltaMat->width<DeltaMat->height)
	{
		longm=DeltaMat->width;
	} 
	else
	{
		longm=DeltaMat->height;
	}
	int k;
	k=log10(longm)/log10(2)/2;//===

	for (int i=1;i<k;i++)
	{
		for(int j=0;j<Ginterval;j++)
		{
			Delta.oct[i].inter[j].Gussianimg=cvCreateMat(DeltaMat->height/(pow(2,i-1)),DeltaMat->width/(pow(2,i-1)),CV_64FC1);
		}
	}
	for(int j=0;j<Ginterval;j++)
	{
		Delta.oct[0].inter[j].Gussianimg=cvCreateMat(DeltaMat->height*2,DeltaMat->width*2,CV_64FC1);
	}
	
	int x=0,y=0;
	int m,n;

//	double v=2;
	double v=pow(2,1.0/(Ginterval-3.0));
//	double v=pow(2,1.0/(Ginterval));

	for (m=0;m<k;m++)
	{
		
		if (m==0)//===
		{
  			//temp=cvCloneMat(DeltaMat);

			for(n=0;n<Ginterval;n++)
			{

				if (n==0)
				{
					CvMat *pymat=cvCreateMat(2*DeltaMat->rows,2*DeltaMat->cols,CV_64FC1);
					 pymat=Img_pyrup(DeltaMat);
					Delta.oct[m].inter[0].Gussianimg=cvCloneMat(pymat);
				}
				else
				{
					
					CvMat *gsm=cvCreateMat(5,5,CV_64FC1);
					double sigmapre=1.6*pow(v,n-1);
					double sigmanow=sigmapre*v;
				//	double sigma=1.6*(v,n);
				//	double sigma=sqrt(1.6*pow(v,n)*1.6*pow(v,n)*(v*v-1));
					double sigma=sqrt(sigmanow*sigmanow-sigmapre*sigmapre);
				//	double sigma=1.6*pow(v,n)*(v-1);
				/*	gsm=GussianModel(5,sigma);
					gsm=Gussian2one(gsm,5,5);
					Delta.oct[m].inter[n].Gussianimg=gussianmat(Delta.oct[m].inter[n-1].Gussianimg,gsm);*/
				cvSmooth(Delta.oct[m].inter[n-1].Gussianimg,Delta.oct[m].inter[n].Gussianimg,CV_GAUSSIAN,5,5,(sigma),(sigma));
					/*Mat src=cvarrToMat(Delta.oct[m].inter[n-1].Gussianimg,true);
					Mat dst=cvarrToMat(Delta.oct[m].inter[n].Gussianimg,false);
					GaussianBlur(src,dst,CvSize(5,5),sqrt(sigma),sqrt(sigma));
					CvMat tempcopy=dst;
					cvCopy(&tempcopy,Delta.oct[m].inter[n].Gussianimg);*/
					}
			}




		}
		else if (m==1)//===
		{
  			//temp=cvCloneMat(DeltaMat);

			for(n=0;n<Ginterval;n++)
			{

				if (n==0)
				{
					Delta.oct[m].inter[0].Gussianimg=cvCloneMat(DeltaMat);
				}
				else
				{
					
					CvMat *gsm=cvCreateMat(5,5,CV_64FC1);
					double sigmapre=1.6*pow(v,n-1);
					double sigmanow=sigmapre*v;
				//	double sigma=1.6*(v,n);
				//	double sigma=sqrt(1.6*pow(v,n)*1.6*pow(v,n)*(v*v-1));
					double sigma=sqrt(sigmanow*sigmanow-sigmapre*sigmapre);
				//	double sigma=1.6*pow(v,n)*(v-1);
				/*	gsm=GussianModel(5,sigma);
					gsm=Gussian2one(gsm,5,5);*/
				//	Delta.oct[m].inter[n].Gussianimg=gussianmat(Delta.oct[m].inter[n-1].Gussianimg,gsm);
				
					cvSmooth(Delta.oct[m].inter[n-1].Gussianimg,Delta.oct[m].inter[n].Gussianimg,CV_GAUSSIAN,5,5,(sigma),(sigma));
					/*Mat src=cvarrToMat(Delta.oct[m].inter[n-1].Gussianimg,true);
					Mat dst=cvarrToMat(Delta.oct[m].inter[n].Gussianimg,false);
					GaussianBlur(src,dst,CvSize(5,5),sqrt(sigma),sqrt(sigma));
					CvMat tempcopy=dst;
					cvCopy(&tempcopy,Delta.oct[m].inter[n].Gussianimg);*/
					}
			}


		
		}
		else
		{
			for(n=0;n<Ginterval;n++)
			{

				if (n==0)
				{
					cout<<"行 "<<Delta.oct[m].inter[0].Gussianimg->rows<<"   列 "<<Delta.oct[m].inter[0].Gussianimg->cols<<endl;
					cout<<"行 "<<Delta.oct[x].inter[y].Gussianimg->rows<<"   列 "<<Delta.oct[x].inter[y].Gussianimg->cols<<endl;
					cout<<"行 "<<GDGM(Delta.oct[x].inter[y].Gussianimg)->rows<<"   列 "<<GDGM(Delta.oct[x].inter[y].Gussianimg)->cols<<endl;
					Delta.oct[m].inter[0].Gussianimg=cvCloneMat(GDGM(Delta.oct[x].inter[y].Gussianimg));
				}
				else
				{
					CvMat *gsm=cvCreateMat(5,5,CV_64FC1);
					double sigmapre=1.6*pow(v,n-1);
					double sigmanow=sigmapre*v;
				//	double sigma=1.6*(v,n);
				//	double sigma=sqrt(1.6*pow(v,n)*1.6*pow(v,n)*(v*v-1));
					double sigma=sqrt(sigmanow*sigmanow-sigmapre*sigmapre);
				//	double sigma=1.6*pow(v,n)*(v-1);
				/*	gsm=GussianModel(5,sigma);
					gsm=Gussian2one(gsm,5,5);*/
				//	Delta.oct[m].inter[n].Gussianimg=gussianmat(Delta.oct[m].inter[n-1].Gussianimg,gsm);
					cvSmooth(Delta.oct[m].inter[n-1].Gussianimg,Delta.oct[m].inter[n].Gussianimg,CV_GAUSSIAN,5,5,sqrt(sigma),sqrt(sigma));
				}

			}
		}
		
		x=m;
		y=Ginterval-3;
//		temp=GDGM(Delta.oct[m].inter[Ginterval-3].Gussianimg);
	

	}

	Delta.longinte=Ginterval;
	Delta.longoct=k;

	return Delta;



}






//void builddelta(delta Delta,int oct,int inte)
//{
//	for (int i=0;i<oct;i++)
//	{
//		for(int j=0;j<inte-1;j++)
//		{
//			Delta.oct[i].inter[j].Gussianimg=cvCreateMat(oct,inte,CV_64FC1);
//		}
//	}
//}

delta DOG(delta Delta)
{
	int xrow,xcol;
	delta DogDelta;
	for (int i=0;i<Delta.longoct;i++)
	{
		for(int j=0;j<Delta.longinte-1;j++)
		{
			xrow=Delta.oct[i].inter[j].Gussianimg->rows;
			xcol=Delta.oct[i].inter[j].Gussianimg->cols;
			DogDelta.oct[i].inter[j].Gussianimg=cvCreateMat(xrow,xcol,CV_64FC1);
		}
	}

	for (int i=0;i<Delta.longoct;i++)
	{
		for(int j=0;j<Delta.longinte-1;j++)
		{
			CvMat *bf1=Delta.oct[i].inter[j+1].Gussianimg;
			CvMat *bf2=Delta.oct[i].inter[j].Gussianimg;
			CvMat *dst=DogDelta.oct[i].inter[j].Gussianimg;

			for(int i = 0; i<bf1->height;i++)
				for(int j=0;j<bf1->width;j++)
				{
					uchar temp = *(bf1->data.db+i*bf1->step/8+j);
					uchar temp_sec =*(bf2->data.db+i*bf2->step/8+j);
					*(dst->data.db+i*dst->step/8+j)=(uchar)(temp-temp_sec);
				}

			/*cvAbsDiff(Delta.oct[i].inter[j+1].Gussianimg,Delta.oct[i].inter[j].Gussianimg,DogDelta.oct[i].inter[j].Gussianimg);*/
			//sub_imagex(Delta.oct[i].inter[j+1].Gussianimg,Delta.oct[i].inter[j].Gussianimg,DogDelta.oct[i].inter[j].Gussianimg);
			//CvMat *aft=cvCreateMat(Delta.oct[i].inter[j+1].Gussianimg->rows,Delta.oct[i].inter[j+1].Gussianimg->cols,CV_64FC1);
			/*Mat af=cvarrToMat(Delta.oct[i].inter[j+1].Gussianimg,true);
			Mat bf=cvarrToMat(Delta.oct[i].inter[j].Gussianimg,true);
			Mat sub=af-bf;
			CvMat tempcopy=sub;
			cvCopy(&tempcopy,DogDelta.oct[i].inter[j].Gussianimg);*/



		}
	}
	DogDelta.longoct=Delta.longoct;
	DogDelta.longinte=Delta.longinte-1;
	return DogDelta;


}





#endif