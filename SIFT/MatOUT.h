
#ifndef Mat_OUT
#define Mat_OUT
#include <math.h>
#include <opencv2/opencv.hpp>
#include <iostream>

 


using namespace std;
void MatOUT(CvMat *mat,int row,int col,int channal)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		for (j=0;j<col;j++)
		{
			printf("%10lf",*(mat->data.db+channal*j+i*mat->step/8));
		}
		cout<<endl;
	}
}







#endif