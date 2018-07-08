#ifndef histo_cir
#define histo_cir
#include <opencv2/opencv.hpp>
#include <highgui.hpp>
#include <cv.h>
#include <iostream>
#include "MatOUT.h"
#include "GUSSIANdelta.h"
#include "histogram;.h"
#include <stdlib.h>
#include <math.h>
using namespace cv;
using namespace std;

//�������������ݶ���================================================================

class histo
{
public:
	CvMat *eight_histo;
};
class point_Feature
{
public:
	histo blocks[9];
	CvMat *Feature_Vector;
};
class img_Feature
{
public:
	int Feature_num;
	point_Feature p_Feature[10000];
};
class Feature
{
public:
	int oct;
	img_Feature oct_imgFeature[10];
};






//===============//
//������ʱ x yҪ���Ŵ�����//
//===============//







//====================================================================================

//��ʼ��������============================================================================
void  init_Feature(Feature *F,int xoct,int *num)
{
	int i,j,k,m;


	F->oct=xoct;
	for (i=0;i<xoct;i++)
	{
		F->oct_imgFeature[i].Feature_num=*(num+i);
		for(j=0;j<*(num+i);j++)
		{
			F->oct_imgFeature[i].p_Feature[j].Feature_Vector=cvCreateMat(1,72,CV_64FC1);

			for (k=0;k<9;k++)
			{
				F->oct_imgFeature[i].p_Feature[j].blocks[k].eight_histo=cvCreateMat(1,8,CV_64FC1);
				for (m=0;m<8;m++)
				{
					*(F->oct_imgFeature[i].p_Feature[j].blocks[k].eight_histo->data.db+m)=0;
				}
				
			}
		}

	}
}

//=====================================================================================


//����ģֵ================================================================================
double Mod(CvMat *NowMat,int lx,int ly)// l
{
	//ֱ����fixpoint����ȥ��������������
	//if (NowMat->rows-x<2||x<2||NowMat->cols-y<2||y<2)
	//{
	//	return 0;
	//}

	double xp1_y=*(NowMat->data.db+(lx+1)*NowMat->step/8+ly);
	double xs1_y=*(NowMat->data.db+(lx-1)*NowMat->step/8+ly);
	double x_yp1=*(NowMat->data.db+(lx)*NowMat->step/8+ly+1);
	double x_ys1=*(NowMat->data.db+(lx)*NowMat->step/8+ly-1);

	double pointmod=sqrt(pow(xp1_y-xs1_y,2)+pow(x_yp1-x_ys1,2));

	return pointmod;
}

//======================================================================================



//�����ݶȷ�����ת�Ƕ� ���ؽǶ�ֵ 360��ʽ
double Grad_sita(CvMat *NowMat,int Ix,int Iy)//�ǶȰ���ͼ������ϵд ��90��Ϊ��������
{

	double xp1_y;
	double xs1_y;
	double x_yp1;
	double x_ys1;

	//if (NowMat->rows-x<2||x<2||NowMat->cols-y<2||y<2)
	//{
	//	return 0;
	//}

	xp1_y=*(NowMat->data.db+(Ix+1)*NowMat->step/8+Iy);
	xs1_y=*(NowMat->data.db+(Ix-1)*NowMat->step/8+Iy);
	x_yp1=*(NowMat->data.db+(Ix)*NowMat->step/8+Iy+1);
	x_ys1=*(NowMat->data.db+(Ix)*NowMat->step/8+Iy-1);

	double Grad=(x_yp1-x_ys1)/(xp1_y-xs1_y);

	if (((x_yp1-x_ys1)==0)&&((xp1_y-xs1_y)==0))
	{
		return 0;
	}
	else if((x_yp1-x_ys1)==0)
	{
		if ((xp1_y-xs1_y)>0)
		{
			return 180;
		} 
		else
		{
			return 0;
		}

	}
	else if((xp1_y-xs1_y)==0)
	{
		if ((x_yp1-x_ys1)>0)
		{
			return 270;
		} 
		else
		{
			return 90;
		}

	}
	else
	{
		if (Grad>0)
		{
			if ((xp1_y-xs1_y)>0)
			{
				double sita=atan(Grad)*180/3.14+180;
				return sita;
			}
			else
			{
				double sita=atan(Grad)*180/3.14;
				return sita;
			}
		}

		if (Grad<0)
		{
			if ((xp1_y-xs1_y)>0)
			{
				double sita=atan(Grad)*180/3.14+180;
				return sita;
			} 
			else
			{
				double sita=atan(Grad)*180/3.14;
				return sita;
			}
		}

	}

}







//����72ά������������==========================================================================
Feature *histocir(delta Delta,CvMat *spoint,int *num)
{
	//ѭ������===
	int i,j,k,l;
	//========

	int xoct=Delta.longoct;//delta�Ĳ���



	//����Բ������===
	int neibor[2][300];//======================����ÿ��������x y ����ʱ��Ҫ�����ĵ������Щֵ
	int neibor_num=0;//======================����������һ���ж��ٵ�
	for (i=0;i<17;i++)//======================8*2+1�����ص�Ϊ���Ϳ������������Ѱ�Ұ뾶Ϊ8��Բ
	{
		for (j=0;j<17;j++)
		{
			int nx=-8+i;
			int ny=-8+j;
			int dx=abs(nx);
			int dy=abs(ny);
			double modxy=sqrt(dx*dx+dy*dy);
			if (modxy<=8)
			{
				neibor[0][neibor_num]=nx;
				neibor[1][neibor_num]=ny;
				neibor_num++;
     		}

		}
	}
	//===========



	//��Բ�������ڵĵ���зֿ� �ֳ�36��===
	int neiborbin[36][100];//======================����ÿ���������ڵ����ţ���Ӧ����neibor[2][300]����ҪѰ�����꣬�����ȥ��neibor��ά�����ö�Ӧ����Ų���
	int neiborbin_num[36];//====================����ÿ���������ڵ������
	for (i=0;i<36;i++)
	{
		neiborbin_num[i]=0;
	}

	for (i=0;i<neibor_num;i++)
	{
		int nx=neibor[0][i];
		int ny=neibor[1][i];
		double k=ny/nx;
		int sitabin;

		if ((ny==0)&&(nx==0))//====================���ĵ㲻����
		{
			/*sitabin=0;
			neiborbin[sitabin][neiborbin_num[sitabin]]=i;*/
		}
		else if(ny==0)
		{
			if (nx>0)
			{
				sitabin=180/10;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			} 
			else
			{
				sitabin=0;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			}

		}
		else if(nx==0)
		{
			if (ny>0)
			{
				sitabin=270/10;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			} 
			else
			{
				sitabin=90/10;
				neiborbin[sitabin][neiborbin_num[sitabin]]=i;
				neiborbin_num[sitabin]++;
			}

		}
		else
		{
			if (k>0)
			{
				if (nx>0)
				{
					double sita=atan(k)*180/3.14+180;
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				}
				else
				{
					double sita=atan(k)*180/3.14;
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				}
			}

			if (k<0)
			{
				if (nx>0)
				{
					double sita=atan(k)*180/3.14+180;
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				} 
				else
				{
					double sita=atan(k)*180/3.14;
					sitabin=(int)(sita/10);
					neiborbin[sitabin][neiborbin_num[sitabin]]=i;
					neiborbin_num[sitabin]++;
				}
			}

		}


	}

	for(i=0;i<36;i++)
	{
		int x=neiborbin_num[i];
		cout<<x<<" ";
	}

	//==================================================




	Feature *F;
	init_Feature(F,xoct,num);




	return F;


}





#endif