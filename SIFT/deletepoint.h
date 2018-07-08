#ifndef deletepoint
#define deletepoint
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"
using namespace cv;
using namespace std;

void DelPoint(CvMat *spoint,CvMat *delp,int *num,int delnum[])
{
	int i,j,k;
	int xnum;
	int p;
	for (i=0;i<spoint->rows/2;i++)
	{
		for (j=0;j<delnum[i];j++)
		{
			xnum=*(delp->data.db+i*delp->step/8+j);
			/*for(k=xnum+1;k<*(num+i);k++)
			{
				*(spoint->data.db+2*i*spoint->step/8+k-1)=*(spoint->data.db+2*i*spoint->step/8+k);
				*(spoint->data.db+(2*i+1)*spoint->step/8+k-1)=*(spoint->data.db+(2*i+1)*spoint->step/8+k);
				

			}*/


			//*(num+i)=*(num+i)-1;

			*(spoint->data.db+2*i*spoint->step/8+xnum)=-1;
			*(spoint->data.db+(2*i+1)*spoint->step/8+xnum)=-1;


		}

	}



	for (i=0;i<spoint->rows/2;i++)
	{

		for (j=0;j<*(num+i);j++)
		{
			
			if (*(spoint->data.db+2*i*spoint->step/8+j)==-1)
			{
					for(k=j+1;k<*(num+i);k++)
					{
						*(spoint->data.db+2*i*spoint->step/8+k-1)=*(spoint->data.db+2*i*spoint->step/8+k);
						*(spoint->data.db+(2*i+1)*spoint->step/8+k-1)=*(spoint->data.db+(2*i+1)*spoint->step/8+k);
				

					}
					j--;
					*(num+i)=*(num+i)-1;
			}




		}



	}






}



#endif