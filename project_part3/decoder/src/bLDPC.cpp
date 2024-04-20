#include "bLDPC.h"
#include "function4Binary.h"
#include <iostream>
#include <algorithm>
#include <bitset>
using namespace std;
#include <stdio.h>
#include <math.h>

qcldpc::qcldpc(int Hbase_rows, int Hbase_cols, int cpm_size, int q)
{
	infoLen = (Hbase_cols - Hbase_rows)*cpm_size;
	codeLen = Hbase_cols*cpm_size;
	base_rows = Hbase_rows;
	base_cols = Hbase_cols;
	CPM_size = cpm_size;
	ROWS = base_rows*CPM_size;
	cols = base_cols*CPM_size;
	now_rows = ROWS;
	qnumber = q;//q-array
	G = new int[(cols - ROWS)*ROWS];
	rearrange = new int[cols];
	H = new int[ROWS*cols];

	for (int j = 0; j<cols; j++)
		rearrange[j] = 0;

	for (int i = 0; i<ROWS; i++)
		for (int j = 0; j<cols; j++)
			H[i*cols + j] = 0;

	for (int i = 0; i<cols - ROWS; i++)
		for (int j = 0; j<ROWS; j++)
			G[i*ROWS + j] = 0;
}
qcldpc::~qcldpc()
{
	delete[]H;
	delete[]G;
	delete[]rearrange;
	if (destroy_crossList())
	{
		cout << "succssfully delete the crossList !" << endl;
	}
}

bool qcldpc::Contained(int *vector, int ele, int length)
{
	for (int i = 0; i != length; ++i)
	{
		if (ele == vector[i])
			return true;
	}
	return false;
}
bool qcldpc::Contained(vector<int> &vector, int ele)
{
	for (int i = 0; i != vector.size(); ++i)
	{
		if (ele == vector[i])
			return true;
	}
	return false;
}


// decoding themes
int qcldpc::log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected)
{
	const int qrM = 8;
	const int qrN = 17;
	const int Usize = 256;
	int Hqr[qrM][qrN] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
		0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
		0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
		0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
		0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
		0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
		0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
		0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };
	int hm[qrN] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };

	/*int hm[23] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };
	double **U = new double*[2048];
	for (int i = 0; i != 2048; ++i)
	{
	U[i] = new double[24];
	}
	for (int i = 0; i != 2048; ++i)
	{
	if (i == 0)
	U[i][0] = 1.0;
	else
	U[i][0] = 0.0;
	}*/

	double **U = new double*[Usize];
	for (int i = 0; i != Usize; ++i)
	{
		U[i] = new double[qrN + 1];
	}
	for (int i = 0; i != Usize; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}



	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	int qr_new[qrN] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	int qr_old[qrN] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	int qr_most[qrN] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };



	int SuccConceutive = 0;
	int SuccConceutiveMost = -1;

	for (int i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
	for (int i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps->Zmn = Fn[ps->col_num];
			ps = ps->right;
		}
	}

	int judgeFlag = 0;

	ps = M.rhead[checkNodeSelected[0]];
	int idx = 0;
	while (ps != NULL)
	{
		if (pn0[ps->col_num] > pn1[ps->col_num])
			qr_old[idx] = 0;
		else
			qr_old[idx] = 1;
		ps = ps->right;    //ָ����һ�����
		idx++;
	}

	vector<int> qrVNode;
	ps = M.rhead[checkNodeSelected[0]];
	while (ps != NULL)
	{
		qrVNode.push_back(ps->col_num);
		ps = ps->right;
	}

	for (int itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		for (int i = 0; i < ROWS; i++)
		{
			if (Contained(checkNodeSelected, i))
			{
				ps = M.rhead[i];
				for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < Usize; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}
					ps = ps->right;
				}

				qs = M.rhead[i];
				for (int m = 0; m < qrN; m++)
				{
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;

						//qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
						qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}
					//cout << qs->Lt_extr << " ";
					//cout << log(qs->outp0 / qs->outp1) << " ";
					qs = qs->right;
				}

				//bool isLegal = true;
				//for (int i = 0; i != 8; ++i)
				//{
				//	int sum = 0;
				//	for (int j = 0; j != 17; ++j)
				//	{
				//		sum ^= qr[j] * Hqr[i][j];
				//	}
				//	if (sum == 1)
				//	{
				//		isLegal = false;
				//		break;
				//	}
				//}
				//if (isLegal)
				//{
				//	//cout << "Legal" << endl;
				//	qs = M.rhead[i];
				//	for (int m = 0; m < 17; m++)
				//	{
				//	    qs->Lt_extr *= 1;
				//		qs = qs->right; //ָ����һ�����
				//	}
				//}
				//else
				//{
				//	//cout << "nonLegal" << endl;
				//	qs = M.rhead[i];
				//	for (int m = 0; m < 17; m++)
				//	{
				//		qs->Lt_extr *= 1;
				//		qs = qs->right; //ָ����һ�����
				//	}
				//}
			}
			else
			{
				//cout << "Normal Check Row" << i << endl;
				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					double min = 100000.0;
					ts = M.rhead[i];
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;

						}
						ts = ts->right;
					}

					//����������,variable nodes get value from check nodes
					qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
					if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;
					qs = qs->right;
				}
			}
		}



		//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (int j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				// transment ps->Zmn as inputs posibility in the next iteration 
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;

				ps = ps->down; //ָ����һ�����  
			}
		}

		ps = M.rhead[checkNodeSelected[0]];
		int idx2 = 0;
		while (ps != NULL)
		{
			qr_new[idx2] = vhat[ps->col_num];
			ps = ps->right;    //ָ����һ�����
			idx2++;
		}

		bool isEqual = true;
		for (int j = 0; j != qrN; ++j)
		{
			if (qr_old[j] != qr_new[j])
			{
				isEqual = false;
				break;
			}
		}
		if (isEqual)
			SuccConceutive++;
		else
			SuccConceutive = 0;

		if (SuccConceutive > SuccConceutiveMost)
		{
			SuccConceutiveMost = SuccConceutive;
			for (int j = 0; j != qrN; ++j)
			{
				qr_most[j] = qr_new[j];
			}
		}

		for (int j = 0; j != qrN; ++j)
		{
			qr_old[j] = qr_new[j];
		}

		//cout << SuccConceutive << endl;

		judgeFlag = otherLDPC1.judgeZero(vhat);
		//judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
			decodingflag++;
			//if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/

	 //cout << "SuccConceutive = " << SuccConceutive << endl;
	 //cout << "judgeFlag = " << judgeFlag << endl;
	 /*
	 if (judgeFlag == 0 && SuccConceutive >= 10)
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (Contained(checkNodeSelected, i))
	 {
	 ps = M.rhead[i];
	 for (int j = 1; j < 18; j++)  //����״̬�����ȸ�����(s)����(n)
	 {
	 for (int m = 0; m < 256; m++)
	 {
	 U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
	 }
	 ps = ps->right;    //ָ����һ�����
	 }
	 qs = M.rhead[i];
	 for (int m = 0; m < 17; m++)
	 {
	 if (qs->qmn0 != qs->qmn1)
	 {
	 qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][17] / U[0][17]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
	 if (qs->outp0 < MIN) qs->outp0 = MIN;
	 qs->outp1 = 1 - qs->outp0;
	 if (qs->outp1 < MIN) qs->outp1 = MIN;

	 //qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	 qs->Lt_extr = coefficient * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	 //qs->Lt_extr = 0;

	 if (qs->Lt_extr > MAXLLR2)
	 qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr < -MAXLLR2)
	 qs->Lt_extr = -MAXLLR2;
	 }
	 //cout << qs->Lt_extr << " ";
	 //cout << log(qs->outp0 / qs->outp1) << " ";
	 qs = qs->right; //ָ����һ�����
	 }//cout<<endl;
	 //cout << endl;


	 }
	 else
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }
	 //Tmn *= min;
	 //if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 //if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 //qs->Lt_extr = Tmn;
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }

	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;

	 // transment ps->Zmn as inputs posibility in the next iteration
	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 } // for(j=0;j<cols;j++)


	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 */
	 /*
	 if (judgeFlag == 0 && SuccConceutive >= 10)
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if(!Contained(checkNodeSelected, i))
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }
	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 if (!Contained(qrVNode, j))
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }
	 }

	 for (j = 0; j != qrVNode.size(); ++j)
	 {
	 ps = M.chead[j];
	 vhat[qrVNode[j]] = qr_new[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 if (qr_new[j] == 0)
	 ps->Zmn = MAXLLR2;
	 else
	 ps->Zmn = -MAXLLR2;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }

	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < 10))
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (!Contained(checkNodeSelected, i))
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }
	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 if (!Contained(qrVNode, j))
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }
	 }

	 for (j = 0; j != qrVNode.size(); ++j)
	 {
	 ps = M.chead[j];
	 vhat[qrVNode[j]] = qr_most[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 if (qr_new[j] == 0)
	 ps->Zmn = MAXLLR2;
	 else
	 ps->Zmn = -MAXLLR2;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }

	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 if (judgeFlag == 0 && SuccConceutive == 0)
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (Contained(checkNodeSelected, i))
	 {
	 qs = M.rhead[i];
	 for (int m = 0; m < 17; m++)
	 {
	 qs->Lt_extr = 0;

	 if (qs->Lt_extr > MAXLLR2)
	 qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr < -MAXLLR2)
	 qs->Lt_extr = -MAXLLR2;

	 qs = qs->right; //ָ����һ�����
	 }//cout<<endl;
	 //cout << endl;


	 }
	 else
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }

	 //qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }

	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;

	 // transment ps->Zmn as inputs posibility in the next iteration
	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;




	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 } // for(j=0;j<cols;j++)


	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 */
	if (judgeFlag == 0 && SuccConceutive >= 10)
	{
		for (int i = 0; i < ROWS; i++)
		{
			ps = M.rhead[i];
			while (ps != NULL)
			{
				ps->qmn0 = pn0[ps->col_num];
				ps->qmn1 = pn1[ps->col_num];
				ps->Zmn = Fn[ps->col_num];
				ps = ps->right;
			}
		}
		for (int itr_num = 1; itr_num <= 50; itr_num++)
		{
			for (int i = 0; i < ROWS; i++)
			{
				if (!Contained(checkNodeSelected, i))
				{
					//cout << "Normal Check Row" << i << endl;
					//2.(SPC) check nodes transmit to variable nodes 
					qs = M.rhead[i];
					while (qs != NULL)
					{
						double Tmn = 1.0;
						double min = 100000.0;
						ts = M.rhead[i]; //��������������
						while (ts != NULL)
						{
							if ((ts->col_num) != (qs->col_num))
							{
								Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
								if (Tmn > MAXLLR2) Tmn = MAXLLR2;
								if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
							}
							ts = ts->right;
						}

						//����������,variable nodes get value from check nodes
						qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
						if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

						//cout << qs->Lt_extr << " ";
						qs = qs->right;
					}
				}
			}//end-for(0;Ldpc_Row)  



			 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
			for (int j = 0; j < cols; j++)
			{
				if (!Contained(qrVNode, j))
				{
					ps = M.chead[j];
					while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
					{
						double prod_rmn = 0.0;

						qs = M.chead[j];
						while (qs != NULL)
						{
							if ((qs->row_num) != (ps->row_num))
							{
								prod_rmn = prod_rmn + qs->Lt_extr;
							}
							qs = qs->down;
						}
						ps->Zmn = Fn[j] + prod_rmn;
						if (ps->Zmn > MAXLLR2)
							ps->Zmn = MAXLLR2;
						if (ps->Zmn < -MAXLLR2)
							ps->Zmn = -MAXLLR2;


						ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
						if (ps->qmn1 < MIN)
							ps->qmn1 = MIN;

						ps->qmn0 = 1 - ps->qmn1;
						if (ps->qmn0 < MIN)
							ps->qmn0 = MIN;

						ps->Zn = ps->Zmn + ps->Lt_extr;
						if (ps->Zn > MAXLLR2)
							ps->Zn = MAXLLR2;
						if (ps->Zn < -MAXLLR2)
							ps->Zn = -MAXLLR2;

						if (ps->Zn < 0)
							vhat[j] = 1;
						else
							vhat[j] = 0;

						ps = ps->down; //ָ����һ�����  

					}//cout<<endl;//end- while(qs!=NULL)
				}
			}

			for (int j = 0; j != qrVNode.size(); ++j)
			{
				ps = M.chead[qrVNode[j]];
				vhat[qrVNode[j]] = qr_new[j];
				while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
				{
					double prod_rmn = 0.0;

					if (qr_new[j] == 0)
						ps->Zmn = MAXLLR2;
					else
						ps->Zmn = -MAXLLR2;

					ps = ps->down; //ָ����һ�����  

				}//cout<<endl;//end- while(qs!=NULL)
			}

			judgeFlag = otherLDPC1.judgeZero(vhat);
			//judgeFlag = judgeZero(vhat);

			if (judgeFlag)
			{  //---����У����˳�,���򷵻ؼ���������
			   //decodingflag = 1;//������ȷ
				decodingflag++;
				break;
			}
		}
	}
	if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < 10))
	{
		for (int i = 0; i < ROWS; i++)//42  rows/3
		{
			ps = M.rhead[i];
			while (ps != NULL)
			{

				ps->qmn0 = pn0[ps->col_num];
				ps->qmn1 = pn1[ps->col_num];
				ps->Zmn = Fn[ps->col_num];
				ps = ps->right;
			}
		}
		for (int itr_num = 1; itr_num <= 50; itr_num++)
		{
			for (int i = 0; i < ROWS; i++)
			{
				if (!Contained(checkNodeSelected, i))
				{
					//cout << "Normal Check Row" << i << endl;
					//2.(SPC) check nodes transmit to variable nodes 
					qs = M.rhead[i];
					while (qs != NULL)
					{
						double Tmn = 1.0;
						double min = 100000.0;
						ts = M.rhead[i]; //��������������
						while (ts != NULL)
						{
							if ((ts->col_num) != (qs->col_num))
							{
								Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
								if (Tmn > MAXLLR2) Tmn = MAXLLR2;
								if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
							}
							ts = ts->right;
						}
						//����������,variable nodes get value from check nodes
						qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
						if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

						//cout << qs->Lt_extr << " ";
						qs = qs->right;
					}//end-while(qs!=NULL)	
					 //cout << endl;
				}
			}//end-for(0;Ldpc_Row)  



			 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
			for (int j = 0; j < cols; j++)//147
			{
				if (!Contained(qrVNode, j))
				{
					ps = M.chead[j];
					while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
					{
						double prod_rmn = 0.0;

						qs = M.chead[j];
						while (qs != NULL)
						{
							if ((qs->row_num) != (ps->row_num))
							{
								prod_rmn = prod_rmn + qs->Lt_extr;
							}
							qs = qs->down;
						}
						ps->Zmn = Fn[j] + prod_rmn;
						if (ps->Zmn > MAXLLR2)
							ps->Zmn = MAXLLR2;
						if (ps->Zmn < -MAXLLR2)
							ps->Zmn = -MAXLLR2;


						ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
						if (ps->qmn1 < MIN)
							ps->qmn1 = MIN;

						ps->qmn0 = 1 - ps->qmn1;
						if (ps->qmn0 < MIN)
							ps->qmn0 = MIN;

						ps->Zn = ps->Zmn + ps->Lt_extr;
						if (ps->Zn > MAXLLR2)
							ps->Zn = MAXLLR2;
						if (ps->Zn < -MAXLLR2)
							ps->Zn = -MAXLLR2;

						if (ps->Zn < 0)
							vhat[j] = 1;
						else
							vhat[j] = 0;

						ps = ps->down; //ָ����һ�����  

					}//cout<<endl;//end- while(qs!=NULL)
				}
			}

			for (int j = 0; j != qrVNode.size(); ++j)
			{
				ps = M.chead[qrVNode[j]];
				vhat[qrVNode[j]] = qr_most[j];
				while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
				{
					double prod_rmn = 0.0;

					if (qr_most[j] == 0)
						ps->Zmn = MAXLLR2;
					else
						ps->Zmn = -MAXLLR2;

					ps = ps->down; //ָ����һ�����  

				}//cout<<endl;//end- while(qs!=NULL)
			}

			judgeFlag = otherLDPC1.judgeZero(vhat);
			//judgeFlag = judgeZero(vhat);

			if (judgeFlag)
			{
				decodingflag++;
				break;
			}
		}
	}
	if (judgeFlag == 0 && SuccConceutive == 0)
	{
		for (int i = 0; i < ROWS; i++)//42  rows/3
		{
			ps = M.rhead[i];
			while (ps != NULL)
			{

				ps->qmn0 = pn0[ps->col_num];
				ps->qmn1 = pn1[ps->col_num];
				ps->Zmn = Fn[ps->col_num];
				ps = ps->right;
			}
		}
		for (int itr_num = 1; itr_num <= 50; itr_num++)
		{
			for (int i = 0; i < ROWS; i++)
			{
				if (Contained(checkNodeSelected, i))
				{
					qs = M.rhead[i];
					while (qs != NULL)
					{
						qs->Lt_extr = 0;

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;

						qs = qs->right;
					}


				}
				else
				{
					//cout << "Normal Check Row" << i << endl;
					//2.(SPC) check nodes transmit to variable nodes 
					qs = M.rhead[i];
					while (qs != NULL)
					{
						double Tmn = 1.0;
						double min = 100000.0;
						ts = M.rhead[i]; //��������������
						while (ts != NULL)
						{
							if ((ts->col_num) != (qs->col_num))
							{
								Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
								if (Tmn > MAXLLR2) Tmn = MAXLLR2;
								if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;

							}
							ts = ts->right;
						}
						////����������,variable nodes get value from check nodes
						qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
						if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

						//cout << qs->Lt_extr << " ";
						qs = qs->right;
					}//end-while(qs!=NULL)	
					 //cout << endl;
				}

			}//end-for(0;Ldpc_Row)  



			 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
			for (int j = 0; j < cols; j++)//147
			{
				ps = M.chead[j];
				while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
				{
					double prod_rmn = 0.0;

					qs = M.chead[j];
					while (qs != NULL)
					{
						if ((qs->row_num) != (ps->row_num))
						{
							prod_rmn = prod_rmn + qs->Lt_extr;
						}
						qs = qs->down;
					}
					ps->Zmn = Fn[j] + prod_rmn;
					if (ps->Zmn > MAXLLR2)
						ps->Zmn = MAXLLR2;
					if (ps->Zmn < -MAXLLR2)
						ps->Zmn = -MAXLLR2;

					/** transment ps->Zmn as inputs posibility in the next iteration **/
					ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
					if (ps->qmn1 < MIN)
						ps->qmn1 = MIN;

					ps->qmn0 = 1 - ps->qmn1;
					if (ps->qmn0 < MIN)
						ps->qmn0 = MIN;

					ps->Zn = ps->Zmn + ps->Lt_extr;
					if (ps->Zn > MAXLLR2)
						ps->Zn = MAXLLR2;
					if (ps->Zn < -MAXLLR2)
						ps->Zn = -MAXLLR2;

					if (ps->Zn < 0)
						vhat[j] = 1;
					else
						vhat[j] = 0;




					ps = ps->down; //ָ����һ�����  

				}//cout<<endl;//end- while(qs!=NULL)
			} // for(j=0;j<cols;j++)


			judgeFlag = otherLDPC1.judgeZero(vhat);
			//judgeFlag = judgeZero(vhat);

			if (judgeFlag)
			{  //---����У����˳�,���򷵻ؼ���������
			   //decodingflag = 1;//������ȷ
				decodingflag++;
				//if(decodingflag==2)
				break;
			}
		}
	}



	for (int i = 0; i != Usize; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}

int qcldpc::min_sum(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected)
{
	const int NQR = 1;
	const int qrM = 8;
	const int qrN = 17;
	const int Usize = 256;
	int Hqr[qrM][qrN] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
		0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
		0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
		0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
		0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
		0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
		0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
		0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };
	int hm[qrN] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };



	/*int hm[23] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };
	double **U = new double*[2048];
	for (int i = 0; i != 2048; ++i)
	{
	U[i] = new double[24];
	}
	for (int i = 0; i != 2048; ++i)
	{
	if (i == 0)
	U[i][0] = 1.0;
	else
	U[i][0] = 0.0;
	}*/

	double **U = new double*[Usize];
	for (int i = 0; i != Usize; ++i)
	{
		U[i] = new double[qrN + 1];
	}
	for (int i = 0; i != Usize; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}



	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	int qr_new[NQR][qrN] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	int qr_old[NQR][qrN] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	int qr_most[NQR][qrN] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };



	int SuccConceutive[NQR] = { 0 };
	int SuccConceutiveMost[NQR] = { -1 };

	for (int i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
	for (int i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps->Zmn = Fn[ps->col_num];
			ps = ps->right;
		}
	}

	int judgeFlag = 0;

	for (int i = 0; i != NQR; ++i)
	{
		ps = M.rhead[checkNodeSelected[i]];
		int idx = 0;
		while (ps != NULL)
		{
			if (pn0[ps->col_num] > pn1[ps->col_num])
				qr_old[i][idx] = 0;
			else
				qr_old[i][idx] = 1;
			ps = ps->right;    //ָ����һ�����
			idx++;
		}
	}

	int **qrVNode = new int*[NQR];
	for (int i = 0; i != NQR; ++i)
	{
		qrVNode[i] = new int[qrN];
	}



	for (int i = 0; i != NQR; ++i)
	{
		ps = M.rhead[checkNodeSelected[i]];
		int idx = 0;
		while (ps != NULL)
		{

			qrVNode[i][idx] = ps->col_num;
			idx++;
			ps = ps->right;
		}

	}



	int itr_num;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		for (int i = 0; i < ROWS; i++)
		{
			if (Contained(checkNodeSelected, i))
			{
				ps = M.rhead[i];
				for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < Usize; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}
					ps = ps->right;
				}

				qs = M.rhead[i];
				for (int m = 0; m < qrN; m++)
				{
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;

						//qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
						qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}
					//cout << qs->Lt_extr << " ";
					//cout << log(qs->outp0 / qs->outp1) << " ";
					qs = qs->right;
				}

				//bool isLegal = true;
				//for (int i = 0; i != 8; ++i)
				//{
				//	int sum = 0;
				//	for (int j = 0; j != 17; ++j)
				//	{
				//		sum ^= qr[j] * Hqr[i][j];
				//	}
				//	if (sum == 1)
				//	{
				//		isLegal = false;
				//		break;
				//	}
				//}
				//if (isLegal)
				//{
				//	//cout << "Legal" << endl;
				//	qs = M.rhead[i];
				//	for (int m = 0; m < 17; m++)
				//	{
				//	    qs->Lt_extr *= 1;
				//		qs = qs->right; //ָ����һ�����
				//	}
				//}
				//else
				//{
				//	//cout << "nonLegal" << endl;
				//	qs = M.rhead[i];
				//	for (int m = 0; m < 17; m++)
				//	{
				//		qs->Lt_extr *= 1;
				//		qs = qs->right; //ָ����һ�����
				//	}
				//}
			}
			else
			{
				//cout << "Normal Check Row" << i << endl;
				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					double min = 100000.0;
					ts = M.rhead[i];
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							/*Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;*/
							double temp = abs(ts->Zmn);
							if (temp < min)
							{
								min = temp;
							}
							if (ts->Zmn >= 0)
							{
								Tmn *= 1.0;
							}
							else
							{
								Tmn *= -1.0;
							}
						}
						ts = ts->right;
					}
					Tmn *= min;
					if (Tmn > MAXLLR2) Tmn = MAXLLR2;
					if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					qs->Lt_extr = Tmn;
					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
					if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;
					qs = qs->right;
				}
			}
		}



		//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (int j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				// transment ps->Zmn as inputs posibility in the next iteration 
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;

				ps = ps->down; //ָ����һ�����  
			}
		}

		for (int j = 0; j != NQR; ++j)
		{
			ps = M.rhead[checkNodeSelected[j]];
			int idx2 = 0;
			while (ps != NULL)
			{
				qr_new[j][idx2] = vhat[ps->col_num];
				ps = ps->right;    //ָ����һ�����
				idx2++;
			}
		}

		bool isEqual[NQR];
		for (int i = 0; i != NQR; ++i)
		{
			isEqual[i] = true;
			for (int j = 0; j != qrN; ++j)
			{
				if (qr_old[i][j] != qr_new[i][j])
				{
					isEqual[i] = false;
					break;
				}
			}
		}
		for (int i = 0; i != NQR; ++i)
		{
			if (isEqual[i])
				SuccConceutive[i]++;
			else
				SuccConceutive[i] = 0;
		}

		for (int i = 0; i != NQR; ++i)
		{
			if (SuccConceutive[i] > SuccConceutiveMost[i])
			{
				SuccConceutiveMost[i] = SuccConceutive[i];
				for (int j = 0; j != qrN; ++j)
				{
					qr_most[i][j] = qr_new[i][j];
				}
			}
		}

		for (int i = 0; i != NQR; ++i)
		{
			for (int j = 0; j != qrN; ++j)
			{
				qr_old[i][j] = qr_new[i][j];
			}
		}

		//cout << SuccConceutive << endl;

		/*cout << "after iteration:" << endl;
		for (int i = 0; i != cols; ++i)
		{
		if (vhat[i] != 0)
		{
		cout << i << " ";
		}
		}
		cout << endl;*/

		judgeFlag = otherLDPC1.judgeZero(vhat);
		//judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
			decodingflag++;
			//if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/

	 /*cout << "SuccConceutive = " << SuccConceutive << endl;
	 cout << "judgeFlag = " << judgeFlag << endl;*/
	 /*
	 if (judgeFlag == 0 && SuccConceutive >= 10)
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (Contained(checkNodeSelected, i))
	 {
	 ps = M.rhead[i];
	 for (int j = 1; j < 18; j++)  //����״̬�����ȸ�����(s)����(n)
	 {
	 for (int m = 0; m < 256; m++)
	 {
	 U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
	 }
	 ps = ps->right;    //ָ����һ�����
	 }
	 qs = M.rhead[i];
	 for (int m = 0; m < 17; m++)
	 {
	 if (qs->qmn0 != qs->qmn1)
	 {
	 qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][17] / U[0][17]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
	 if (qs->outp0 < MIN) qs->outp0 = MIN;
	 qs->outp1 = 1 - qs->outp0;
	 if (qs->outp1 < MIN) qs->outp1 = MIN;

	 //qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	 qs->Lt_extr = coefficient * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	 //qs->Lt_extr = 0;

	 if (qs->Lt_extr > MAXLLR2)
	 qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr < -MAXLLR2)
	 qs->Lt_extr = -MAXLLR2;
	 }
	 //cout << qs->Lt_extr << " ";
	 //cout << log(qs->outp0 / qs->outp1) << " ";
	 qs = qs->right; //ָ����һ�����
	 }//cout<<endl;
	 //cout << endl;


	 }
	 else
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }
	 //Tmn *= min;
	 //if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 //if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 //qs->Lt_extr = Tmn;
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }

	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;

	 // transment ps->Zmn as inputs posibility in the next iteration
	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 } // for(j=0;j<cols;j++)


	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 */
	 /*
	 if (judgeFlag == 0 && SuccConceutive >= 10)
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if(!Contained(checkNodeSelected, i))
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }
	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 if (!Contained(qrVNode, j))
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }
	 }

	 for (j = 0; j != qrVNode.size(); ++j)
	 {
	 ps = M.chead[j];
	 vhat[qrVNode[j]] = qr_new[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 if (qr_new[j] == 0)
	 ps->Zmn = MAXLLR2;
	 else
	 ps->Zmn = -MAXLLR2;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }

	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < 10))
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (!Contained(checkNodeSelected, i))
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }
	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 if (!Contained(qrVNode, j))
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }
	 }

	 for (j = 0; j != qrVNode.size(); ++j)
	 {
	 ps = M.chead[j];
	 vhat[qrVNode[j]] = qr_most[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 if (qr_new[j] == 0)
	 ps->Zmn = MAXLLR2;
	 else
	 ps->Zmn = -MAXLLR2;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }

	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 if (judgeFlag == 0 && SuccConceutive == 0)
	 {
	 coefficient = 2;
	 for (i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (Contained(checkNodeSelected, i))
	 {
	 qs = M.rhead[i];
	 for (int m = 0; m < 17; m++)
	 {
	 qs->Lt_extr = 0;

	 if (qs->Lt_extr > MAXLLR2)
	 qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr < -MAXLLR2)
	 qs->Lt_extr = -MAXLLR2;

	 qs = qs->right; //ָ����һ�����
	 }//cout<<endl;
	 //cout << endl;


	 }
	 else
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 }
	 ts = ts->right;
	 }

	 //qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }

	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (j = 0; j < cols; j++)//147
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;

	 // transment ps->Zmn as inputs posibility in the next iteration
	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;




	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 } // for(j=0;j<cols;j++)


	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 */

	 /*cout << "after first decoding:" << endl;
	 for (int i = 0; i != cols; ++i)
	 {
	 if (vhat[i] != 0)
	 {
	 cout << i << " ";
	 }
	 }
	 cout << endl;*/

	 /*
	 if (judgeFlag == 0 && SuccConceutive >= itr_num/5*4)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {
	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (int itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (!Contained(checkNodeSelected, i))
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;

	 }
	 ts = ts->right;
	 }

	 //����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }
	 }
	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (int j = 0; j < cols; j++)
	 {
	 if (!Contained(qrVNode, j))
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }
	 }

	 for (int j = 0; j != qrVNode.size(); ++j)
	 {
	 ps = M.chead[qrVNode[j]];
	 vhat[qrVNode[j]] = qr_new[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 if (qr_new[j] == 0)
	 ps->Zmn = MAXLLR2;
	 else
	 ps->Zmn = -MAXLLR2;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }



	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 break;
	 }
	 }
	 }
	 if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < itr_num/5*4))
	 {
	 for (int i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (int itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (!Contained(checkNodeSelected, i))
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;

	 }
	 ts = ts->right;
	 }

	 //����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }
	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (int j = 0; j < cols; j++)//147
	 {
	 if (!Contained(qrVNode, j))
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }
	 }

	 for (int j = 0; j != qrVNode.size(); ++j)
	 {
	 ps = M.chead[qrVNode[j]];
	 vhat[qrVNode[j]] = qr_most[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 if (qr_most[j] == 0)
	 ps->Zmn = MAXLLR2;
	 else
	 ps->Zmn = -MAXLLR2;

	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 }

	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {
	 decodingflag++;
	 break;
	 }
	 }
	 }
	 if (judgeFlag == 0 && SuccConceutive == 0)
	 {
	 for (int i = 0; i < ROWS; i++)//42  rows/3
	 {
	 ps = M.rhead[i];
	 while (ps != NULL)
	 {

	 ps->qmn0 = pn0[ps->col_num];
	 ps->qmn1 = pn1[ps->col_num];
	 ps->Zmn = Fn[ps->col_num];
	 ps = ps->right;
	 }
	 }
	 for (int itr_num = 1; itr_num <= 50; itr_num++)
	 {
	 for (int i = 0; i < ROWS; i++)
	 {
	 if (Contained(checkNodeSelected, i))
	 {

	 ps = M.rhead[i];
	 for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
	 {
	 for (int m = 0; m < Usize; m++)
	 {
	 U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
	 }
	 ps = ps->right;
	 }

	 qs = M.rhead[i];
	 for (int m = 0; m < qrN; m++)
	 {
	 if (qs->qmn0 != qs->qmn1)
	 {
	 qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
	 if (qs->outp0 < MIN) qs->outp0 = MIN;
	 qs->outp1 = 1 - qs->outp0;
	 if (qs->outp1 < MIN) qs->outp1 = MIN;

	 //qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	 qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

	 if (qs->Lt_extr > MAXLLR2)
	 qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr < -MAXLLR2)
	 qs->Lt_extr = -MAXLLR2;
	 }
	 //cout << qs->Lt_extr << " ";
	 //cout << log(qs->outp0 / qs->outp1) << " ";
	 qs = qs->right;
	 }
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 if (qs->col_num<cols - qrM)
	 qs->Lt_extr = 0;

	 qs = qs->right;
	 }


	 }
	 else
	 {
	 //cout << "Normal Check Row" << i << endl;
	 //2.(SPC) check nodes transmit to variable nodes
	 qs = M.rhead[i];
	 while (qs != NULL)
	 {
	 double Tmn = 1.0;
	 double min = 100000.0;
	 ts = M.rhead[i]; //��������������
	 while (ts != NULL)
	 {
	 if ((ts->col_num) != (qs->col_num))
	 {
	 Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
	 if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;

	 }
	 ts = ts->right;
	 }

	 ////����������,variable nodes get value from check nodes
	 qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
	 if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //cout << qs->Lt_extr << " ";
	 qs = qs->right;
	 }//end-while(qs!=NULL)
	 //cout << endl;
	 }

	 }//end-for(0;Ldpc_Row)



	 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 for (int j = 0; j < cols; j++)//147
	 {
	 ps = M.chead[j];
	 while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 {
	 double prod_rmn = 0.0;

	 qs = M.chead[j];
	 while (qs != NULL)
	 {
	 if ((qs->row_num) != (ps->row_num))
	 {
	 prod_rmn = prod_rmn + qs->Lt_extr;
	 }
	 qs = qs->down;
	 }
	 ps->Zmn = Fn[j] + prod_rmn;
	 if (ps->Zmn > MAXLLR2)
	 ps->Zmn = MAXLLR2;
	 if (ps->Zmn < -MAXLLR2)
	 ps->Zmn = -MAXLLR2;


	 ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 if (ps->qmn1 < MIN)
	 ps->qmn1 = MIN;

	 ps->qmn0 = 1 - ps->qmn1;
	 if (ps->qmn0 < MIN)
	 ps->qmn0 = MIN;

	 ps->Zn = ps->Zmn + ps->Lt_extr;
	 if (ps->Zn > MAXLLR2)
	 ps->Zn = MAXLLR2;
	 if (ps->Zn < -MAXLLR2)
	 ps->Zn = -MAXLLR2;

	 if (ps->Zn < 0)
	 vhat[j] = 1;
	 else
	 vhat[j] = 0;




	 ps = ps->down; //ָ����һ�����

	 }//cout<<endl;//end- while(qs!=NULL)
	 } // for(j=0;j<cols;j++)


	 judgeFlag = otherLDPC1.judgeZero(vhat);
	 //judgeFlag = judgeZero(vhat);

	 if (judgeFlag)
	 {  //---����У����˳�,���򷵻ؼ���������
	 //decodingflag = 1;//������ȷ
	 decodingflag++;
	 //if(decodingflag==2)
	 break;
	 }
	 }
	 }
	 */

	for (int i = 0; i != NQR; ++i)
	{
		delete[]qrVNode[i];
	}
	delete[]qrVNode;
	for (int i = 0; i != Usize; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}
int qcldpc::min_sum_v2(int *CodeWord, int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected)
{
	const int qrM = 8;
	const int qrN = 17;
	const int Usize = 1 << qrM;
	int Hqr[qrM][qrN] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
		0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
		0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
		0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
		0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
		0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
		0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
		0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };

	int hm[qrN] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };

	/*int Hqr[qrM][qrN] = {	1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1 };
	int hm[qrN] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };*/


	double **U = new double*[Usize];
	for (int i = 0; i != Usize; ++i)
	{
		U[i] = new double[qrN + 1];
	}
	for (int i = 0; i != Usize; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}



	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	int qr_new[qrN] = { 0 };
	int qr_old[qrN] = { 0 };
	int qr_most[qrN] = { 0 };
	int qr[qrN] = { 0 };



	int SuccConceutive = 0;
	int SuccConceutiveMost = -1;

	for (int i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
	for (int i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps->Zmn = Fn[ps->col_num];
			ps = ps->right;
		}
	}

	int judgeFlag = 0;

	ps = M.rhead[checkNodeSelected[0]];
	int idx = 0;
	while (ps != NULL)
	{
		if (pn0[ps->col_num] > pn1[ps->col_num])
			qr_old[idx] = 0;
		else
			qr_old[idx] = 1;
		ps = ps->right;    //ָ����һ�����
		idx++;
	}

	vector<int> qrVNode;
	ps = M.rhead[checkNodeSelected[0]];
	while (ps != NULL)
	{
		qrVNode.push_back(ps->col_num);
		ps = ps->right;
	}

	ofstream os2;
	os2.open("GLDPC_new_design\\information.txt");

	const int q = 31;
	int **QR_arr_new;
	QR_arr_new = new int*[q];
	for (int i = 0; i != q; ++i)
	{
		QR_arr_new[i] = new int[qrN];
	}

	int **QR_arr_old;
	QR_arr_old = new int*[q];
	for (int i = 0; i != q; ++i)
	{
		QR_arr_old[i] = new int[qrN];
	}
	int qrPartIdx = 0;
	for (int i = 0; i < ROWS; i++)
	{
		if (Contained(checkNodeSelected, i))
		{
			qs = M.rhead[i];
			for (int m = 0; m < qrN; m++)
			{
				if (Fn[qs->col_num] > 0)
				{
					QR_arr_old[qrPartIdx][m] = 0;
				}
				else
				{
					QR_arr_old[qrPartIdx][m] = 1;
				}
				qs = qs->right;
			}
			qrPartIdx++;
		}

	}

	int S[q] = { 0 }; //��ʾ�仯��

	vector< vector< bitset<qrN> > > qrCode(q, vector< bitset<qrN> >());
	vector< vector< int > > fre(q, vector<int>());

	int itr_num;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		os2 << itr_num << " ";
		qrPartIdx = 0;
		for (int i = 0; i < ROWS; i++)
		{

			if (Contained(checkNodeSelected, i))
			{

				bitset<qrN> qrCodeSub = 0;

				ps = M.rhead[i];

				for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < Usize; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}
					ps = ps->right;

				}
				int qr[17] = { 0 };

				qs = M.rhead[i];
				for (int m = 0; m < qrN; m++)
				{
					qs->Lt_extr = 0;
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;


						qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}
					if (qs->Lt_extr > 0)
						qr[m] = 0;
					else
						qr[m] = 1;
					if (qs->Lt_extr > 0)
					{
						QR_arr_new[qrPartIdx][m] = 0;
					}
					else
					{
						QR_arr_new[qrPartIdx][m] = 1;
					}

					if (qr[m] == 1)
						qrCodeSub ^= (1 << m);

					qs = qs->right;
				}

				int isCheckOld = 0;
				for (int m = 0; m != qrM; ++m)
				{
					int value = 0;
					for (int n = 0; n != qrN; ++n)
					{
						value ^= QR_arr_old[qrPartIdx][n] * Hqr[m][n];
					}
					if (value != 0)
					{
						isCheckOld++;
					}
				}

				int isCheckNew = 0;
				for (int m = 0; m != qrM; ++m)
				{
					int value = 0;
					for (int n = 0; n != qrN; ++n)
					{
						value ^= QR_arr_new[qrPartIdx][n] * Hqr[m][n];
					}
					if (value != 0)
					{
						isCheckNew++;
					}
				}
				if (isCheckNew == 0)
				{
					bool isThere = false;
					bool atIdx = -1;
					for (int m = 0; m != qrCode[qrPartIdx].size(); ++m)
					{
						if (qrCode[qrPartIdx][m] == qrCodeSub)
						{
							isThere = true;
							fre[qrPartIdx][m]++;
							atIdx = m;
							break;
						}
					}

					if (!isThere)
					{
						qrCode[qrPartIdx].push_back(qrCodeSub);
						fre[qrPartIdx].push_back(1);
					}
				}

				bool isEqualOld2New = true;
				for (int m = 0; m < qrN; m++)
				{
					if (QR_arr_old[qrPartIdx][m] != QR_arr_new[qrPartIdx][m])
					{
						isEqualOld2New = false;
						break;
					}
				}

				for (int m = 0; m < qrN; m++)
				{
					QR_arr_old[qrPartIdx][m] = QR_arr_new[qrPartIdx][m];
				}

				//if (isEqualOld2New)
				if (isCheckNew != 0 && isCheckOld != 0)
					S[qrPartIdx] += 1;
				else
					S[qrPartIdx] += 0;

				qrPartIdx++;

				int isOriIdx = 0;
				bool isOri = true;
				qs = M.rhead[i];
				for (int n = 0; n != qrN; ++n)
				{
					if (qr[n] != CodeWord[qs->col_num])
					{
						isOri = false;
						isOriIdx++;
					}
					qs = qs->right;
				}

				if (isOri)
					os2 << i << "->right " << " ";
				else
					os2 << i << "->wrong " << " ";

				qs = M.rhead[i];
				while (qs != NULL)
				{
					os2 << CodeWord[qs->col_num];
					qs = qs->right;
				}
				os2 << " ";



				for (int n = 0; n != qrN; ++n)
				{
					os2 << qr[n];
				}
				os2 << " ";

				if (isCheckNew == 0)
					os2 << i << "->check right " << " " << isCheckNew << " ";
				else
					os2 << i << "->check wrong " << " " << isCheckNew << " ";

				os2 << " S = " << S[qrPartIdx - 1] << " ";

				/*qs = M.rhead[i];
				for (int m = 0; m < qrN; m++)
				{
				os2 << qs->Lt_extr << " ";
				qs = qs->right;
				}

				os2 << " ***** ";*/

				double coefficient[8 + 1] = { 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2 };
				double coefficientOri[17 + 1] = { 1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15 };
				//if (isCheck == 0)
				////if(isOriIdx != 0)
				//{
				//	qs = M.rhead[i];
				//	for (int m = 0; m < qrN; m++)
				//	{
				//		qs->Lt_extr *= (coefficient[isCheck] + S[qrPartIdx-1]*0.002);
				//		//qs->Lt_extr *= coefficientOri[isOriIdx];
				//		qs = qs->right;
				//	}
				//}

				if (isCheckNew != 0)
					//if(isOriIdx != 0)
				{
					qs = M.rhead[i];
					for (int m = 0; m < qrN; m++)
					{
						qs->Lt_extr *= coefficient[isCheckNew];
						//qs->Lt_extr *= coefficientOri[isOriIdx];
						qs = qs->right;
					}
				}


				/*qs = M.rhead[i];
				for (int m = 0; m < qrN; m++)
				{
				os2 << qs->Lt_extr << " ";
				qs = qs->right;
				}*/

			}
			else
			{

				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					double min = 100000.0;
					ts = M.rhead[i];
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							/*Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;*/
							double temp = abs(ts->Zmn);
							if (temp < min)
							{
								min = temp;
							}
							if (ts->Zmn >= 0)
							{
								Tmn *= 1.0;
							}
							else
							{
								Tmn *= -1.0;
							}
						}
						ts = ts->right;
					}
					Tmn *= min;
					if (Tmn > MAXLLR2) Tmn = MAXLLR2;
					if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					qs->Lt_extr = Tmn;
					////����������,variable nodes get value from check nodes
					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
					//if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					//if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;
					qs = qs->right;
				}
			}
		}
		os2 << endl;


		//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (int j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)
			{
				double prod_rmn = 0.0;
				double temp = ps->Zmn;
				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}

				ps->Zmn = Fn[j] + prod_rmn;

				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				// transment ps->Zmn as inputs posibility in the next iteration 
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn <= 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;

				ps = ps->down; //ָ����һ�����  
			}
		}





		ps = M.rhead[checkNodeSelected[0]];
		int idx2 = 0;
		while (ps != NULL)
		{
			qr_new[idx2] = vhat[ps->col_num];
			ps = ps->right;    //ָ����һ�����
			idx2++;
		}

		bool isEqual = true;
		for (int j = 0; j != qrN; ++j)
		{
			if (qr_old[j] != qr_new[j])
			{
				isEqual = false;
				break;
			}
		}
		if (isEqual)
			SuccConceutive++;
		else
			SuccConceutive = 0;

		if (SuccConceutive > SuccConceutiveMost)
		{
			SuccConceutiveMost = SuccConceutive;
			for (int j = 0; j != qrN; ++j)
			{
				qr_most[j] = qr_new[j];
			}
		}

		for (int j = 0; j != qrN; ++j)
		{
			qr_old[j] = qr_new[j];
		}



		judgeFlag = otherLDPC1.judgeZero(vhat);


		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
		   //if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/
	 /*
	 if (!judgeFlag)
	 {
	 int qrPartIdx = 0;
	 for (int i = 0; i < ROWS; i++)
	 {

	 if (Contained(checkNodeSelected, i))
	 {
	 int theMostReliable = -1;
	 int theMostReliableIdx = -1;
	 for (int j = 0; j != fre[qrPartIdx].size(); ++j)
	 {
	 if (theMostReliable < fre[qrPartIdx][j])
	 {
	 theMostReliable = fre[qrPartIdx][j];
	 theMostReliableIdx = j;
	 }
	 }



	 if (theMostReliableIdx != -1)
	 {
	 cout << i << " -> " << qrCode[qrPartIdx][theMostReliableIdx] << endl;
	 qs = M.rhead[i];
	 bitset<qrN> one = 1;
	 for (int j = 0; j != qrN; ++j)
	 {
	 if (((one << j) & (qrCode[qrPartIdx][theMostReliableIdx])) != 0)
	 {
	 vhat[qs->col_num] = 1;
	 }
	 else
	 {
	 vhat[qs->col_num] = 0;
	 }

	 qs = qs->right;
	 }
	 }

	 qrPartIdx++;


	 }

	 }
	 }
	 */
	for (int i = 0; i != q; ++i)
	{
		delete[]QR_arr_new[i];
	}
	delete[]QR_arr_new;

	decodingflag = 1;

	/*

	vector<int> idxMin;
	vector<double> valueMin;
	const int numMin = 9;
	for (int i = 0; i != cols; ++i)
	{
	vector<double>::iterator iter1 = valueMin.begin();
	vector<int>::iterator iter2 = idxMin.begin();
	while (iter1 != valueMin.end() && abs(M.chead[i]->Zn) >= *iter1)
	{
	iter1++;
	iter2++;
	}
	valueMin.insert(iter1, abs(M.chead[i]->Zn));
	idxMin.insert(iter2, i);
	while (idxMin.size() > numMin)
	{
	idxMin.pop_back();
	valueMin.pop_back();
	}
	}
	for (int i = 0; i != idxMin.size(); ++i)
	{
	valueMin[i] = M.chead[idxMin[i]]->Zn;
	}

	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn_pre = ps->Zmn;
	ps->qmn0_pre = ps->qmn0;
	ps->qmn1_pre = ps->qmn1;
	ps = ps->right;
	}
	}



	if (judgeFlag == 0 && SuccConceutive >= itr_num / 5 * 4)
	{
	decodingflag = 2;

	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps = ps->right;
	}
	}




	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (!Contained(checkNodeSelected, i))
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i];
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	//////����������,variable nodes get value from check nodes
	//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	qs = qs->right;
	}

	}
	}



	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	if (!Contained(qrVNode, j))
	{
	ps = M.chead[j];
	while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;


	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;

	ps = ps->down;

	}
	}
	}

	for (int j = 0; j != qrVNode.size(); ++j)
	{
	ps = M.chead[qrVNode[j]];
	vhat[qrVNode[j]] = qr_new[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	if (qr_new[j] == 0)
	ps->Zmn = MAXLLR2;
	else
	ps->Zmn = -MAXLLR2;
	ps = ps->down;
	}
	}


	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	ps->Zmn = -MAXLLR2;
	else
	ps->Zmn = MAXLLR2;

	ps = ps->down;
	}




	judgeFlag = otherLDPC1.judgeZero(vhat);
	//judgeFlag = judgeZero(vhat);



	if (judgeFlag)
	{
	break;
	}

	}
	idxForMin++;
	}
	}

	if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < itr_num / 5 * 4))
	{
	decodingflag = 3;



	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps = ps->right;
	}
	}



	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (!Contained(checkNodeSelected, i))
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i]; //��������������
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	////����������,variable nodes get value from check nodes
	//qs->Lt_extr = coefficient + (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;


	qs = qs->right;
	}
	}
	}



	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	if (!Contained(qrVNode, j))
	{
	ps = M.chead[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;


	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;

	ps = ps->down; //ָ����һ�����

	}
	}
	}

	for (int j = 0; j != qrVNode.size(); ++j)
	{
	ps = M.chead[qrVNode[j]];
	vhat[qrVNode[j]] = qr_most[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;

	if (qr_most[j] == 0)
	ps->Zmn = MAXLLR2;
	else
	ps->Zmn = -MAXLLR2;

	ps = ps->down; //ָ����һ�����

	}
	}

	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	ps->Zmn = -MAXLLR2;
	else
	ps->Zmn = MAXLLR2;

	ps = ps->down; //ָ����һ�����
	}

	judgeFlag = otherLDPC1.judgeZero(vhat);
	//judgeFlag = judgeZero(vhat);



	if (judgeFlag)
	{
	break;
	}
	}
	idxForMin++;
	}
	}

	if (judgeFlag == 0 && SuccConceutive == 0)
	{
	decodingflag = 4;

	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps->qmn0 = ps->qmn0_pre;
	ps->qmn1 = ps->qmn1_pre;
	ps = ps->right;
	}
	}



	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (Contained(checkNodeSelected, i))
	{

	ps = M.rhead[i];
	for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
	{
	for (int m = 0; m < Usize; m++)
	{
	U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
	}
	ps = ps->right;
	}

	qs = M.rhead[i];
	for (int m = 0; m < qrN; m++)
	{
	if (qs->qmn0 != qs->qmn1)
	{
	qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
	if (qs->outp0 < MIN) qs->outp0 = MIN;
	qs->outp1 = 1 - qs->outp0;
	if (qs->outp1 < MIN) qs->outp1 = MIN;

	//qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

	if (qs->Lt_extr > MAXLLR2)
	qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2)
	qs->Lt_extr = -MAXLLR2;
	}

	qs = qs->right;
	}
	qs = M.rhead[i];
	while (qs != NULL)
	{
	if (qs->col_num < cols - qrM)
	qs->Lt_extr = 0;

	qs = qs->right;
	}
	}
	else
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i];
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	//////����������,variable nodes get value from check nodes
	//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	qs = qs->right;
	}
	}

	}
	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	ps = M.chead[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;

	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;




	ps = ps->down; //ָ����һ�����

	}
	}

	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	{
	ps->Zmn = -MAXLLR2;
	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;
	}
	else
	{
	ps->Zmn = MAXLLR2;
	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;
	}

	ps = ps->down;
	}

	judgeFlag = otherLDPC1.judgeZero(vhat);

	if (judgeFlag)
	{
	break;
	}
	}
	idxForMin++;
	}
	}

	*/

	//delete[]variableDis;

	for (int i = 0; i != Usize; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}
int qcldpc::min_sum_v3(int *CodeWord, int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected)
{
	/*
	const int qrM = 8;
	const int qrN = 17;
	const int Usize = 1 << qrM;
	int Hqr[qrM][qrN] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
	0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
	0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
	0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
	0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
	0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
	0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
	0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };

	int hm[qrN] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };
	*/
	const int qrROWS = 3;
	const int qrCOLS = 7;
	const int Usize = 1 << qrROWS;
	int Hqr[qrROWS][qrCOLS] = { 1,0,0,1,0,1,1,
		0,1,0,1,1,1,0,
		0,0,1,0,1,1,1 };
	int hm[qrCOLS] = { 4, 2, 1, 6, 3, 7, 5 };


	/*int Hqr[qrM][qrN] = {	1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1 };
	int hm[qrN] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };*/


	double **U = new double*[Usize];
	for (int i = 0; i != Usize; ++i)
	{
		U[i] = new double[qrCOLS + 1];
	}
	for (int i = 0; i != Usize; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}



	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	int qr_new[qrCOLS] = { 0 };
	int qr_old[qrCOLS] = { 0 };
	int qr_most[qrCOLS] = { 0 };
	int qr[qrCOLS] = { 0 };



	int SuccConceutive = 0;
	int SuccConceutiveMost = -1;

	for (int i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
	for (int i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps->Zmn = Fn[ps->col_num];
			ps = ps->right;
		}
	}

	int judgeFlag = 0;

	vector< vector<int> > QRPart;
	vector< int > QRFlag;
	vector< int > QRrow;

	int itr_num;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		QRPart.clear();
		QRFlag.clear();
		QRrow.clear();
		for (int i = 0; i < ROWS; i++)
		{

			if (Contained(checkNodeSelected, i))
			{
				ps = M.rhead[i];

				for (int j = 1; j < qrCOLS + 1; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < Usize; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}
					ps = ps->right;

				}



				qs = M.rhead[i];
				for (int m = 0; m < qrCOLS; m++)
				{
					qs->Lt_extr = 0;
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrCOLS] / U[0][qrCOLS]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;


						qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}
					qs = qs->right;
				}
			}
			else
			{

				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					double min = 100000.0;
					ts = M.rhead[i];
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							/*Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;*/
							double temp = abs(ts->Zmn);
							if (temp < min)
							{
								min = temp;
							}
							if (ts->Zmn >= 0)
							{
								Tmn *= 1.0;
							}
							else
							{
								Tmn *= -1.0;
							}
						}
						ts = ts->right;
					}
					Tmn *= min;
					if (Tmn > MAXLLR2) Tmn = MAXLLR2;
					if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					qs->Lt_extr = Tmn;
					////����������,variable nodes get value from check nodes
					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
					//if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					//if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;
					qs = qs->right;
				}
			}
		}



		//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (int j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)
			{
				double prod_rmn = 0.0;
				double temp = ps->Zmn;
				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}

				ps->Zmn = Fn[j] + prod_rmn;



				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				// transment ps->Zmn as inputs posibility in the next iteration 
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn <= 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;

				ps = ps->down; //ָ����һ�����  
			}
		}



		judgeFlag = otherLDPC1.judgeZero(vhat);

		/*
		for (int i = 0; i<ROWS; i++)
		{
		ps = M.rhead[i];
		int sum = 0;
		while (ps != NULL)
		{

		if (vhat[ps->col_num] != 0)
		{
		sum ^= 1;//add[sum][mul[ps->e][vs[ps->col_num]]];
		}
		ps = ps->right;
		}

		// ��ʱk4��h(i,:)*vhat;
		if (sum != 0)
		{
		judgeFlag = 0;
		cout << i << " ";
		}
		}
		cout << endl;

		if (itr_num % 10 == 0)
		{
		for (int idx = 0; idx != QRFlag.size(); ++idx)
		{
		if (QRFlag[idx] == 1)
		{
		cout << "&&&&&&&&&&&&&" << endl;
		cout << QRrow[idx] << endl;
		ps = M.rhead[QRrow[idx]];
		for (int m = 0; m != qrN; ++m)
		{
		qs = M.chead[ps->col_num];
		while (qs != NULL)
		{
		qs->Zmn = (QRPart[idx][m] == 0) ? MAXLLR2 : -MAXLLR2;
		qs = qs->down;
		}
		ps = ps->right;
		}
		cout << endl;
		}
		}

		for (int i = 0; i != QRPart.size(); ++i)
		{
		for (int j = 0; j != QRPart[i].size(); ++j)
		{
		cout << QRPart[i][j];
		}
		cout << endl;

		ps = M.rhead[QRrow[i]];
		for (int m = 0; m != qrN; ++m)
		{
		cout << CodeWord[ps->col_num];
		ps = ps->right;
		}
		cout << endl << endl;
		}
		}
		*/

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
		   //if(decodingflag==2)

			break;
		}
	}/*********���һ֡������********/

	vector<int> disVar(cols, 0);
	for (int i = 0; i != ROWS; ++i)
	{
		if (!Contained(checkNodeSelected, i))
		{
			int value = 0;
			ps = M.rhead[i];
			while (ps != NULL)
			{
				value ^= vhat[ps->col_num];
				disVar[ps->col_num]++;
				ps = ps->right;
			}
			/*	if (value != 0)
			cout << i << " ";*/
		}
	}
	//cout << endl;

	//for (int j = 0; j != cols; ++j)
	//{
	//	cout << disVar[j] << " ";
	//}
	//cout << endl;

	decodingflag = 1;

	/*

	vector<int> idxMin;
	vector<double> valueMin;
	const int numMin = 9;
	for (int i = 0; i != cols; ++i)
	{
	vector<double>::iterator iter1 = valueMin.begin();
	vector<int>::iterator iter2 = idxMin.begin();
	while (iter1 != valueMin.end() && abs(M.chead[i]->Zn) >= *iter1)
	{
	iter1++;
	iter2++;
	}
	valueMin.insert(iter1, abs(M.chead[i]->Zn));
	idxMin.insert(iter2, i);
	while (idxMin.size() > numMin)
	{
	idxMin.pop_back();
	valueMin.pop_back();
	}
	}
	for (int i = 0; i != idxMin.size(); ++i)
	{
	valueMin[i] = M.chead[idxMin[i]]->Zn;
	}

	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn_pre = ps->Zmn;
	ps->qmn0_pre = ps->qmn0;
	ps->qmn1_pre = ps->qmn1;
	ps = ps->right;
	}
	}



	if (judgeFlag == 0 && SuccConceutive >= itr_num / 5 * 4)
	{
	decodingflag = 2;

	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps = ps->right;
	}
	}




	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (!Contained(checkNodeSelected, i))
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i];
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	//////����������,variable nodes get value from check nodes
	//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	qs = qs->right;
	}

	}
	}



	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	if (!Contained(qrVNode, j))
	{
	ps = M.chead[j];
	while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;


	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;

	ps = ps->down;

	}
	}
	}

	for (int j = 0; j != qrVNode.size(); ++j)
	{
	ps = M.chead[qrVNode[j]];
	vhat[qrVNode[j]] = qr_new[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	if (qr_new[j] == 0)
	ps->Zmn = MAXLLR2;
	else
	ps->Zmn = -MAXLLR2;
	ps = ps->down;
	}
	}


	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	ps->Zmn = -MAXLLR2;
	else
	ps->Zmn = MAXLLR2;

	ps = ps->down;
	}




	judgeFlag = otherLDPC1.judgeZero(vhat);
	//judgeFlag = judgeZero(vhat);



	if (judgeFlag)
	{
	break;
	}

	}
	idxForMin++;
	}
	}

	if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < itr_num / 5 * 4))
	{
	decodingflag = 3;



	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps = ps->right;
	}
	}



	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (!Contained(checkNodeSelected, i))
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i]; //��������������
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	////����������,variable nodes get value from check nodes
	//qs->Lt_extr = coefficient + (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;


	qs = qs->right;
	}
	}
	}



	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	if (!Contained(qrVNode, j))
	{
	ps = M.chead[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;


	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;

	ps = ps->down; //ָ����һ�����

	}
	}
	}

	for (int j = 0; j != qrVNode.size(); ++j)
	{
	ps = M.chead[qrVNode[j]];
	vhat[qrVNode[j]] = qr_most[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;

	if (qr_most[j] == 0)
	ps->Zmn = MAXLLR2;
	else
	ps->Zmn = -MAXLLR2;

	ps = ps->down; //ָ����һ�����

	}
	}

	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	ps->Zmn = -MAXLLR2;
	else
	ps->Zmn = MAXLLR2;

	ps = ps->down; //ָ����һ�����
	}

	judgeFlag = otherLDPC1.judgeZero(vhat);
	//judgeFlag = judgeZero(vhat);



	if (judgeFlag)
	{
	break;
	}
	}
	idxForMin++;
	}
	}

	if (judgeFlag == 0 && SuccConceutive == 0)
	{
	decodingflag = 4;

	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps->qmn0 = ps->qmn0_pre;
	ps->qmn1 = ps->qmn1_pre;
	ps = ps->right;
	}
	}



	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (Contained(checkNodeSelected, i))
	{

	ps = M.rhead[i];
	for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
	{
	for (int m = 0; m < Usize; m++)
	{
	U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
	}
	ps = ps->right;
	}

	qs = M.rhead[i];
	for (int m = 0; m < qrN; m++)
	{
	if (qs->qmn0 != qs->qmn1)
	{
	qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
	if (qs->outp0 < MIN) qs->outp0 = MIN;
	qs->outp1 = 1 - qs->outp0;
	if (qs->outp1 < MIN) qs->outp1 = MIN;

	//qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

	if (qs->Lt_extr > MAXLLR2)
	qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2)
	qs->Lt_extr = -MAXLLR2;
	}

	qs = qs->right;
	}
	qs = M.rhead[i];
	while (qs != NULL)
	{
	if (qs->col_num < cols - qrM)
	qs->Lt_extr = 0;

	qs = qs->right;
	}
	}
	else
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i];
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	//////����������,variable nodes get value from check nodes
	//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	qs = qs->right;
	}
	}

	}
	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	ps = M.chead[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;

	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;




	ps = ps->down; //ָ����һ�����

	}
	}

	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	{
	ps->Zmn = -MAXLLR2;
	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;
	}
	else
	{
	ps->Zmn = MAXLLR2;
	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;
	}

	ps = ps->down;
	}

	judgeFlag = otherLDPC1.judgeZero(vhat);

	if (judgeFlag)
	{
	break;
	}
	}
	idxForMin++;
	}
	}

	*/

	//delete[]variableDis;

	for (int i = 0; i != Usize; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}
int qcldpc::min_sum_v4(int *CodeWord, int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected)
{
	const int qrM = 8;
	const int qrN = 17;
	const int Usize = 1 << qrM;
	int Hqr[qrM][qrN] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
		0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
		0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
		0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
		0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
		0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
		0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
		0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };

	int hm[qrN] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };

	/*int Hqr[qrM][qrN] = {	1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1 };
	int hm[qrN] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };*/


	double **U = new double*[Usize];
	for (int i = 0; i != Usize; ++i)
	{
		U[i] = new double[qrN + 1];
	}
	for (int i = 0; i != Usize; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}



	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	int qr_new[qrN] = { 0 };
	int qr_old[qrN] = { 0 };
	int qr_most[qrN] = { 0 };
	int qr_v1[qrN] = { 0 };
	int qr_v2[qrN] = { 0 };
	int *codeFromPartialChk = new int[cols];
	int *distri = new int[cols];
	for (int j = 0; j != cols; ++j)
	{
		distri[j] = 0;
	}


	int SuccConceutive = 0;
	int SuccConceutiveMost = -1;

	for (int i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
	for (int i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps->Zmn = Fn[ps->col_num];
			ps = ps->right;
		}
	}

	int judgeFlag = 0;
	int itr_num;

	ofstream os;
	os.open("GLDPC_new_design//information.txt");

	ofstream os1;
	os1.open("GLDPC_new_design//errorBits_information-2.txt");

	int MAXITERNO = _MAXITERNO;

	for (itr_num = 1; itr_num <= MAXITERNO; itr_num++)
	{

		os1 << itr_num << " ********";

		int rightNum = 0;
		for (int i = 0; i < ROWS; i++)
		{


			if (Contained(checkNodeSelected, i))
			{

				ps = M.rhead[i];

				for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < Usize; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}
					ps = ps->right;

				}

				double sum = 0.0;

				qs = M.rhead[i];
				for (int m = 0; m < qrN; m++)
				{
					qs->Lt_extr = 0;
					qs->Lt_extr_str = 0;
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;


						qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

						qs->Lt_extr_str = qs->Lt_extr;

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}

					sum += abs(qs->Lt_extr);

					if (qs->Lt_extr > 0)
						qr_v1[m] = 0;
					else
						qr_v1[m] = 1;

					if (log(qs->outp0 / qs->outp1) > 0)
						qr_v2[m] = 0;
					else
						qr_v2[m] = 1;

					if (qs->Lt_extr > 0)
						codeFromPartialChk[qs->col_num] = 0;
					else
						codeFromPartialChk[qs->col_num] = 1;

					qs = qs->right;
				}

				int diff_1 = 0;
				for (int m = 0; m != qrN; ++m)
				{
					if (qr_v1[m] != qr_v2[m])
						diff_1++;
				}

				int diff_2 = 0;
				ps = M.rhead[i];
				for (int m = 0; m != qrN; ++m)
				{
					if (qr_v1[m] != CodeWord[ps->col_num])
						diff_2++;
					ps = ps->right;
				}

				os1 << " v1<->v2 " << diff_1 << " ";
				os1 << " v1<->CW " << diff_2 << " ";

				int isOri = 0;
				ps = M.rhead[i];
				for (int m = 0; m != qrN; ++m)
				{
					if (qr_v1[m] != CodeWord[ps->col_num])
					{
						isOri++;
					}
					ps = ps->right;
				}

				int isCheck_v1 = 0;
				for (int m = 0; m != qrM; ++m)
				{
					int value = 0;
					for (int n = 0; n != qrN; ++n)
					{
						value ^= qr_v1[n] * Hqr[m][n];
					}
					if (value != 0)
					{
						isCheck_v1++;
					}
				}

				/*int isCheck_v2 = 0;
				for (int m = 0; m != qrM; ++m)
				{
				int value = 0;
				for (int n = 0; n != qrN; ++n)
				{
				value ^= qr_v2[n] * Hqr[m][n];
				}
				if (value != 0)
				{
				isCheck_v2++;
				}
				}*/

				double coefficient[8 + 1] = { 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2 };

				if (isOri == 0)
					rightNum++;

				if (isOri == 0)
					os1 << i << " -> right ";
				else
					os1 << i << " -> wrong ";

				if (isCheck_v1 == 0)
					os1 << i << " ->  partial chk right ";
				else
					os1 << i << " ->  partial chk wrong ";


				if (isCheck_v1 != 0)
				{
					qs = M.rhead[i];
					for (int m = 0; m < qrN; m++)
					{
						qs->Lt_extr *= coefficient[isCheck_v1];
						qs = qs->right;
					}
				}

				if (itr_num % 10 == 0)
				{
					/*if (isCheck_v1 != 0)
					{
					cout << "check v1" << endl;
					qs = M.rhead[i];
					while (qs != NULL)
					{
					cout << qs->col_num << " ";
					qs = qs->right;
					}
					cout << endl;
					}

					if (isCheck_v2 != 0)
					{
					cout << "check v2" << endl;
					qs = M.rhead[i];
					while (qs != NULL)
					{
					cout << qs->col_num << " ";
					qs = qs->right;
					}
					cout << endl;
					}*/


					// cout << "chk " << i << " " << isCheck_v1 << " " << sum << " *******************" << endl;
					if (isCheck_v1 != 0)
					{
						qs = M.rhead[i];
						int m = 0;
						while (qs != NULL)
						{
							// cout << qs->col_num << " ";
							//distri[qs->col_num]++;
							qs = qs->right;
							m++;
						}
						// cout << endl;
					}



					/*if (isCheck_v2 != 0)
					{
					qs = M.rhead[i];
					int m = 0;
					while (qs != NULL)
					{

					cout << qs->col_num << " ";
					qs = qs->right;
					}
					cout << endl;
					}*/


					/*
					if (isCheck_v1 != 0 && diff_1 != 0 && isCheck_v2 == 0)
					{
					qs = M.rhead[i];
					for (int m = 0; m < qrN; m++)
					{
					if (qr_v1[m] != qr_v2[m])
					{
					if (abs(Fn[qs->col_num]) < 5.0)
					{
					if (qr_v2[m] == 0)
					qs->Lt_extr = 0.5*MAXLLR2;
					else
					qs->Lt_extr = -0.5*MAXLLR2;
					}
					}
					qs = qs->right;
					}
					}
					*/
				}

			}
			else
			{

				//2.(SPC) check nodes transmit to variable nodes 

				int chk = 0;

				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					double min = 100000.0;
					ts = M.rhead[i];
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							/*Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;*/
							double temp = abs(ts->Zmn);
							if (temp < min)
							{
								min = temp;
							}
							if (ts->Zmn >= 0)
							{
								Tmn *= 1.0;
							}
							else
							{
								Tmn *= -1.0;
							}
						}
						ts = ts->right;
					}
					Tmn *= min;

					chk ^= (Tmn > 0) ? 0 : 1;

					if (Tmn > MAXLLR2) Tmn = MAXLLR2;
					if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					qs->Lt_extr = Tmn;
					////����������,variable nodes get value from check nodes
					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
					//if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					//if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;
					qs = qs->right;
				}

				/*if (itr_num % 10 == 0)
				{
				if (chk == 1)
				{
				qs = M.rhead[i];
				while (qs != NULL)
				{
				distri[qs->col_num]++;
				qs = qs->right;
				}
				}
				}*/
			}
		}
		os1 << " " << rightNum << " ";
		os1 << endl;

		if (itr_num % 10 == 0)
		{

			// cout << "************************" << endl;
			// for (int j = 0; j != cols; ++j)
			// {
			// 	// if (distri[j] == 1)
			// 		// cout << j << " ";
			// }
			// cout << endl;
		}

		//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (int j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)
			{
				double prod_rmn = 0.0;
				double temp = ps->Zmn;
				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}

				ps->Zmn = Fn[j] + prod_rmn;



				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				// transment ps->Zmn as inputs posibility in the next iteration 
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn <= 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;

				ps = ps->down; //ָ����һ�����  
			}
		}


		os << itr_num << "   ***************** " << endl;
		for (int j = 0; j != cols; ++j)
		{
			if (j == 361)
			{
				ps = M.chead[j];
				os << "Var " << j << " ";
				os << "REAL " << CodeWord[j] << " ";
				os << "Deci " << vhat[j] << " ";
				os << "Zn: " << ps->Zn << " ";
				os << "Fn: " << Fn[ps->col_num] << " ";
				while (ps != NULL)
				{
					os << ps->Lt_extr << " ";
					ps = ps->down;
				}
				os << endl;
			}
		}


		judgeFlag = otherLDPC1.judgeZero(vhat);

		if (itr_num % 10 == 0)
		{



			for (int i = 0; i < ROWS; i++)
			{


				if (Contained(checkNodeSelected, i))
				{

					ps = M.rhead[i];

					for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
					{
						for (int m = 0; m < Usize; m++)
						{
							U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
						}
						ps = ps->right;

					}

					int qr[17] = { 0 };

					qs = M.rhead[i];
					for (int m = 0; m < qrN; m++)
					{
						qr[m] = vhat[qs->col_num];
						qs = qs->right;
					}







					int isCheck = 0;
					for (int m = 0; m != qrM; ++m)
					{
						int value = 0;
						for (int n = 0; n != qrN; ++n)
						{
							value ^= qr[n] * Hqr[m][n];
						}
						if (value != 0)
						{
							isCheck++;
						}
					}

					if (isCheck != 0)
					{
						qs = M.rhead[i];
						while (qs != NULL)
						{
							distri[qs->col_num]++;
							qs = qs->right;
						}
					}



				}
				else
				{

					//2.(SPC) check nodes transmit to variable nodes 

					int chk = 0;

					qs = M.rhead[i];
					while (qs != NULL)
					{
						chk ^= vhat[qs->col_num];
						qs = qs->right;
					}


					if (chk == 1)
					{
						qs = M.rhead[i];
						while (qs != NULL)
						{
							distri[qs->col_num]++;
							qs = qs->right;
						}
					}
				}
			}

			// for (int j = 0; j != cols; ++j)
			// {
			// 	if (distri[j] == 2)
			// 		cout << j << " ";
			// }
			// cout << endl;

			for (int j = 0; j != cols; ++j)
			{
				if (j == 9 ||
					//j == 30 ||
					j == 46 ||
					j == 50 ||
					j == 56 ||
					//j == 62 ||
					j == 70 ||
					//j == 84 ||
					j == 97 ||
					j == 103 ||
					//j == 104 ||
					//j == 131 ||
					j == 161 ||
					//j == 163 ||
					//j == 189 ||
					//j == 196 ||
					//j == 197 ||
					j == 211 ||
					j == 216 ||
					j == 223 ||
					j == 236 ||
					j == 263 ||
					j == 295 ||
					//	j == 296 ||
					//	j == 301 ||
					//	j == 349 ||
					//	j == 355 ||
					//	j == 361 ||
					j == 388 ||
					j == 395 ||
					//	j == 398 ||
					j == 407 ||
					//	j == 419 ||
					j == 421 ||
					j == 426 ||
					j == 428 ||
					j == 440 ||
					j == 454 ||
					j == 468 ||
					j == 473 ||
					j == 481 ||
					//	j == 482 ||
					j == 494 ||
					j == 502 ||
					j == 514 //||
							 //	j == 518
					)
				{
					//if (distri[j] != 0)
					//{
					//if (vhat[j] != CodeWord[j])
					//{
					qs = M.chead[j];
					while (qs != NULL)
					{
						//qs->Zmn = (CodeWord[j] == 0) ? MAXLLR2 : -MAXLLR2;
						qs->Zmn = (vhat[j] == 1) ? MAXLLR2 : -MAXLLR2;
						qs = qs->down;
					}
					//}
					//}
				}
			}

			// cout << endl << endl;
			// for (int j = 0; j != cols; ++j)
			// {
			// 	if (vhat[j] != CodeWord[j])
			// 		cout << j << " ";
			// }
			// cout << endl;
		}


		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
		   //if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/

	delete[]distri;
	delete[]codeFromPartialChk;
	os.close();
	os1.close();

	decodingflag = 1;

	/*

	vector<int> idxMin;
	vector<double> valueMin;
	const int numMin = 9;
	for (int i = 0; i != cols; ++i)
	{
	vector<double>::iterator iter1 = valueMin.begin();
	vector<int>::iterator iter2 = idxMin.begin();
	while (iter1 != valueMin.end() && abs(M.chead[i]->Zn) >= *iter1)
	{
	iter1++;
	iter2++;
	}
	valueMin.insert(iter1, abs(M.chead[i]->Zn));
	idxMin.insert(iter2, i);
	while (idxMin.size() > numMin)
	{
	idxMin.pop_back();
	valueMin.pop_back();
	}
	}
	for (int i = 0; i != idxMin.size(); ++i)
	{
	valueMin[i] = M.chead[idxMin[i]]->Zn;
	}

	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn_pre = ps->Zmn;
	ps->qmn0_pre = ps->qmn0;
	ps->qmn1_pre = ps->qmn1;
	ps = ps->right;
	}
	}



	if (judgeFlag == 0 && SuccConceutive >= itr_num / 5 * 4)
	{
	decodingflag = 2;

	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps = ps->right;
	}
	}




	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (!Contained(checkNodeSelected, i))
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i];
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	//////����������,variable nodes get value from check nodes
	//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	qs = qs->right;
	}

	}
	}



	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	if (!Contained(qrVNode, j))
	{
	ps = M.chead[j];
	while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;


	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;

	ps = ps->down;

	}
	}
	}

	for (int j = 0; j != qrVNode.size(); ++j)
	{
	ps = M.chead[qrVNode[j]];
	vhat[qrVNode[j]] = qr_new[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	if (qr_new[j] == 0)
	ps->Zmn = MAXLLR2;
	else
	ps->Zmn = -MAXLLR2;
	ps = ps->down;
	}
	}


	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	ps->Zmn = -MAXLLR2;
	else
	ps->Zmn = MAXLLR2;

	ps = ps->down;
	}




	judgeFlag = otherLDPC1.judgeZero(vhat);
	//judgeFlag = judgeZero(vhat);



	if (judgeFlag)
	{
	break;
	}

	}
	idxForMin++;
	}
	}

	if (judgeFlag == 0 && (SuccConceutive > 0 && SuccConceutive < itr_num / 5 * 4))
	{
	decodingflag = 3;



	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps = ps->right;
	}
	}



	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (!Contained(checkNodeSelected, i))
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i]; //��������������
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	////����������,variable nodes get value from check nodes
	//qs->Lt_extr = coefficient + (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;


	qs = qs->right;
	}
	}
	}



	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	if (!Contained(qrVNode, j))
	{
	ps = M.chead[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;


	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;

	ps = ps->down; //ָ����һ�����

	}
	}
	}

	for (int j = 0; j != qrVNode.size(); ++j)
	{
	ps = M.chead[qrVNode[j]];
	vhat[qrVNode[j]] = qr_most[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;

	if (qr_most[j] == 0)
	ps->Zmn = MAXLLR2;
	else
	ps->Zmn = -MAXLLR2;

	ps = ps->down; //ָ����һ�����

	}
	}

	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	ps->Zmn = -MAXLLR2;
	else
	ps->Zmn = MAXLLR2;

	ps = ps->down; //ָ����һ�����
	}

	judgeFlag = otherLDPC1.judgeZero(vhat);
	//judgeFlag = judgeZero(vhat);



	if (judgeFlag)
	{
	break;
	}
	}
	idxForMin++;
	}
	}

	if (judgeFlag == 0 && SuccConceutive == 0)
	{
	decodingflag = 4;

	int idxForMin = 0;
	while (judgeFlag == 0 && idxForMin != valueMin.size())
	{
	for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	ps->Zmn = ps->Zmn_pre;
	ps->qmn0 = ps->qmn0_pre;
	ps->qmn1 = ps->qmn1_pre;
	ps = ps->right;
	}
	}



	for (int itr_num = 1; itr_num <= 50; itr_num++)
	{
	for (int i = 0; i < ROWS; i++)
	{
	if (Contained(checkNodeSelected, i))
	{

	ps = M.rhead[i];
	for (int j = 1; j < qrN + 1; j++)  //����״̬�����ȸ�����(s)����(n)
	{
	for (int m = 0; m < Usize; m++)
	{
	U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
	}
	ps = ps->right;
	}

	qs = M.rhead[i];
	for (int m = 0; m < qrN; m++)
	{
	if (qs->qmn0 != qs->qmn1)
	{
	qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][qrN] / U[0][qrN]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
	if (qs->outp0 < MIN) qs->outp0 = MIN;
	qs->outp1 = 1 - qs->outp0;
	if (qs->outp1 < MIN) qs->outp1 = MIN;

	//qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
	qs->Lt_extr = (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

	if (qs->Lt_extr > MAXLLR2)
	qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2)
	qs->Lt_extr = -MAXLLR2;
	}

	qs = qs->right;
	}
	qs = M.rhead[i];
	while (qs != NULL)
	{
	if (qs->col_num < cols - qrM)
	qs->Lt_extr = 0;

	qs = qs->right;
	}
	}
	else
	{
	qs = M.rhead[i];
	while (qs != NULL)
	{
	double Tmn = 1.0;
	double min = 100000.0;
	ts = M.rhead[i];
	while (ts != NULL)
	{
	if ((ts->col_num) != (qs->col_num))
	{
	double temp = abs(ts->Zmn);
	if (temp < min)
	{
	min = temp;
	}
	if (ts->Zmn >= 0)
	{
	Tmn *= 1.0;
	}
	else
	{
	Tmn *= -1.0;
	}
	}
	ts = ts->right;
	}
	Tmn *= min;
	if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	qs->Lt_extr = Tmn;
	//////����������,variable nodes get value from check nodes
	//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	qs = qs->right;
	}
	}

	}
	//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	for (int j = 0; j < cols; j++)
	{
	ps = M.chead[j];
	while (ps != NULL)
	{
	double prod_rmn = 0.0;
	double temp = ps->Zmn;
	qs = M.chead[j];
	while (qs != NULL)
	{
	if ((qs->row_num) != (ps->row_num))
	{
	prod_rmn = prod_rmn + qs->Lt_extr;
	}
	qs = qs->down;
	}


	ps->Zmn = Fn[j] + prod_rmn;

	if (ps->Zmn > MAXLLR2)
	ps->Zmn = MAXLLR2;
	if (ps->Zmn < -MAXLLR2)
	ps->Zmn = -MAXLLR2;

	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;

	ps->Zn = ps->Zmn + ps->Lt_extr;
	if (ps->Zn > MAXLLR2)
	ps->Zn = MAXLLR2;
	if (ps->Zn < -MAXLLR2)
	ps->Zn = -MAXLLR2;

	if (ps->Zn <= 0)
	vhat[j] = 1;
	else
	vhat[j] = 0;




	ps = ps->down; //ָ����һ�����

	}
	}

	ps = M.chead[idxMin[idxForMin]];
	while (ps != NULL)
	{
	if (valueMin[idxForMin] >= 0)
	{
	ps->Zmn = -MAXLLR2;
	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;
	}
	else
	{
	ps->Zmn = MAXLLR2;
	ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	if (ps->qmn1 < MIN)
	ps->qmn1 = MIN;

	ps->qmn0 = 1 - ps->qmn1;
	if (ps->qmn0 < MIN)
	ps->qmn0 = MIN;
	}

	ps = ps->down;
	}

	judgeFlag = otherLDPC1.judgeZero(vhat);

	if (judgeFlag)
	{
	break;
	}
	}
	idxForMin++;
	}
	}

	*/

	//delete[]variableDis;

	for (int i = 0; i != Usize; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}


int qcldpc::log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected, int *QRvalue, int *SPCvalue)
{
	int Hqr[8][17] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
		0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
		0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
		0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
		0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
		0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
		0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
		0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };

	int i, j;
	int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;

	int hm[17] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };
	//double U[256][18] = { { 1.0 } };
	double **U = new double*[256];
	for (int i = 0; i != 256; ++i)
	{
		U[i] = new double[18];
	}
	for (int i = 0; i != 256; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}

	/*int hm[23] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };
	double **U = new double*[2048];
	for (int i = 0; i != 2048; ++i)
	{
	U[i] = new double[24];
	}
	for (int i = 0; i != 2048; ++i)
	{
	if (i == 0)
	U[i][0] = 1.0;
	else
	U[i][0] = 0.0;
	}*/

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	for (i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		//cout<<pn1[i]<<" ";
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}//cout<<endl;


	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  //1. У��ڵ�ĳ�ʼ��
	for (i = 0; i < ROWS; i++)//42  rows/3
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}
	//cout << "iteration begins " << endl;
	int judgeFlag = 0;
	//cout << "ROWS = " << ROWS << endl;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//cout << "itr_num = " << itr_num << endl;
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{

			if (Contained(checkNodeSelected, i))
			{

				/*cout << "row Selected Check Row" << i << endl;

				cout << "**************" << endl;
				ps = M.rhead[i];
				while (ps != NULL)
				{
				if (TotalLLR[ps->col_num] > 0)
				cout << 0 << " ";
				else
				cout << 1 << " ";
				ps = ps->right;
				}
				cout << endl;*/

				ps = M.rhead[i];
				for (int j = 1; j < 18; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < 256; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}

					ps = ps->right;    //ָ����һ�����

				}

				int qr[17] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

				qs = M.rhead[i];
				for (int m = 0; m < 17; m++)
				{
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][17] / U[0][17]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;


						qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));

						//qs->Lt_extr = log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1);


						if (qs->outp1 > qs->outp0)
							qr[m] = 1;

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
						//		int tmpValue = qs->Lt_extr;
						//		QRvalue[tmpValue+30]++;
					}
					//cout << qs->Lt_extr << " ";
					//cout << log(qs->outp0 / qs->outp1) << " ";
					qs = qs->right; //ָ����һ�����
				}//cout<<endl;
				 //cout << endl;

				 //bool isLegal = true;
				 //for (int i = 0; i != 8; ++i)
				 //{
				 //	int sum = 0;
				 //	for (int j = 0; j != 17; ++j)
				 //	{
				 //		sum ^= qr[j] * Hqr[i][j];
				 //	}
				 //	if (sum == 1)
				 //	{
				 //		isLegal = false;
				 //		break;
				 //	}
				 //}
				 //if (isLegal)
				 //{
				 //	//cout << "Legal" << endl;
				 //	qs = M.rhead[i];
				 //	for (int m = 0; m < 17; m++)
				 //	{
				 //	    qs->Lt_extr *= 1;
				 //		qs = qs->right; //ָ����һ�����
				 //	}
				 //}
				 //else
				 //{
				 //	//cout << "nonLegal" << endl;
				 //	qs = M.rhead[i];
				 //	for (int m = 0; m < 17; m++)
				 //	{
				 //		qs->Lt_extr *= 1;
				 //		qs = qs->right; //ָ����һ�����
				 //	}
				 //}
			}
			else
			{
				//cout << "Normal Check Row" << i << endl;
				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					ts = M.rhead[i]; //��������������
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
						}
						ts = ts->right;
					}
					//qs->Lt_extr = 1.15 * (log((1 + Tmn) / (1 - Tmn))) + 5;//����������,variable nodes get value from check nodes
					qs->Lt_extr = 1.15 * (log((1 + Tmn) / (1 - Tmn)));
					if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

					//	int tmp = qs->Lt_extr;
					//	SPCvalue[tmp+30]++;

					//cout << qs->Lt_extr << " ";
					qs = qs->right;
				}//end-while(qs!=NULL)	
				 //cout << endl;
			}

		}//end-for(0;Ldpc_Row)  


		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j < cols; j++)//147
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;




				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)

		judgeFlag = otherLDPC1.judgeZero(vhat);
		//judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
			decodingflag++;
			//if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/

	 /*if (itr_num == 51)
	 {
	 int bitPos[23];
	 int bitVal[23];
	 for (int i = 0; i != 23; ++i)
	 {
	 bitPos[i] = 0 + i * 6;
	 bitVal[i] = vhat[bitPos[i]];
	 }
	 otherLDPC2.spa_decoder_bitPilot(vhat, channelData, 50, sigma, bitPos, bitVal, 23);
	 }*/

	 /*for (int i = 0; i != 23; ++i)
	 {
	 cout << M.chead[0 + i * 6]->Zn << " ";
	 }
	 cout << endl;

	 cout << "iter_num = " << itr_num << endl;*/

	for (int i = 0; i != 256; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}
int qcldpc::log_Decoder_v2(int *vhat, int *code, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected)
{
	int Hqr[8][17] = { 1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,
		0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,
		0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,
		0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,
		0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,
		0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,
		0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,
		0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1 };

	//int i, j;
	//int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;

	int hm[17] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };
	//double U[256][18] = { { 1.0 } };
	double **U = new double*[256];
	for (int i = 0; i != 256; ++i)
	{
		U[i] = new double[18];
	}
	for (int i = 0; i != 256; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}

	/*int hm[23] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };
	double **U = new double*[2048];
	for (int i = 0; i != 2048; ++i)
	{
	U[i] = new double[24];
	}
	for (int i = 0; i != 2048; ++i)
	{
	if (i == 0)
	U[i][0] = 1.0;
	else
	U[i][0] = 0.0;
	}*/

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	int *vhat_pre = new int[cols];

	for (int i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		//cout<<pn1[i]<<" ";
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
		if (pn0[i] > pn1[i])
		{
			vhat_pre[i] = 0;
		}
		else
			vhat_pre[i] = 1;
	}//cout<<endl;

	int *freqChanged = new int[cols];
	for (int i = 0; i != cols; ++i)
	{
		freqChanged[i] = 0;
	}

	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  //1. У��ڵ�ĳ�ʼ��
	for (int i = 0; i < ROWS; i++)//42  rows/3
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}
	//cout << "iteration begins " << endl;
	int judgeFlag = 0;
	//cout << "ROWS = " << ROWS << endl;
	for (int itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//cout << "itr_num = " << itr_num << endl;
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{
			if (Contained(checkNodeSelected, i))
			{

				/*cout << "row Selected Check Row" << i << endl;

				cout << "**************" << endl;
				ps = M.rhead[i];
				while (ps != NULL)
				{
				if (TotalLLR[ps->col_num] > 0)
				cout << 0 << " ";
				else
				cout << 1 << " ";
				ps = ps->right;
				}
				cout << endl;*/

				ps = M.rhead[i];
				for (int j = 1; j < 18; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < 256; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;
					}
					ps = ps->right;    //ָ����һ�����
				}

				int qr[17] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

				qs = M.rhead[i];
				for (int m = 0; m < 17; m++)
				{
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][17] / U[0][17]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;

						qs->Lt_extr = 0.5 * (log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1));
						//qs->Lt_extr = log(qs->outp0 / qs->outp1);


						if (qs->outp1 > qs->outp0)
							qr[m] = 1;

						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}
					//cout << qs->Lt_extr << " ";
					//cout << log(qs->outp0 / qs->outp1) << " ";
					qs = qs->right; //ָ����һ�����
				}//cout<<endl;
				 //cout << endl;

				 //bool isLegal = true;
				 //for (int i = 0; i != 8; ++i)
				 //{
				 //	int sum = 0;
				 //	for (int j = 0; j != 17; ++j)
				 //	{
				 //		sum ^= qr[j] * Hqr[i][j];
				 //	}
				 //	if (sum == 1)
				 //	{
				 //		isLegal = false;
				 //		break;
				 //	}
				 //}
				 //if (isLegal)
				 //{
				 //	//cout << "Legal" << endl;
				 //	qs = M.rhead[i];
				 //	for (int m = 0; m < 17; m++)
				 //	{
				 //	    qs->Lt_extr *= 1;
				 //		qs = qs->right; //ָ����һ�����
				 //	}
				 //}
				 //else
				 //{
				 //	//cout << "nonLegal" << endl;
				 //	qs = M.rhead[i];
				 //	for (int m = 0; m < 17; m++)
				 //	{
				 //		qs->Lt_extr *= 1;
				 //		qs = qs->right; //ָ����һ�����
				 //	}
				 //}
			}
			else
			{
				//cout << "Normal Check Row" << i << endl;
				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					ts = M.rhead[i]; //��������������
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
						}
						ts = ts->right;
					}
					qs->Lt_extr = 1.15 * log((1 + Tmn) / (1 - Tmn));//����������,variable nodes get value from check nodes
					if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

					//cout << qs->Lt_extr << " ";
					qs = qs->right;
				}//end-while(qs!=NULL)	
				 //cout << endl;
			}

		}//end-for(0;Ldpc_Row)  


		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (int j = 0; j < cols; j++)//147
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;




				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)
		// cout << itr_num << ") ";
		// for (int i = 0; i != cols; ++i)
		// {
		// 	if (vhat[i] != code[i])
		// 		cout << i << " ";
		// }
		// cout << endl;

		for (int i = 0; i != cols; ++i)
		{
			if (vhat[i] != vhat_pre[i])
			{
				freqChanged[i] ++;
			}
			vhat_pre[i] = vhat[i];
		}

		/*for (int i = 0; i != cols; ++i)
		{
		if (freqChanged[i] != 0)
		{
		cout << i << "] " << freqChanged[i] << " ";
		}
		}
		cout << endl;*/

		judgeFlag = otherLDPC1.judgeZero(vhat);
		//judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
			decodingflag = 1;
			break;
		}
	}/*********���һ֡������********/

	delete[]freqChanged;
	delete[]vhat_pre;
	for (int i = 0; i != 256; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}
int qcldpc::log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, int *selectRows, int seleRowsN, qcldpc &otherLDPC1, qcldpc &otherLDPC2)
{

	int i, j;
	int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;

	int hm[17] = { 128,64,32,144,200,228,242,121,60,158,79,39,19,9,4,2,1 };
	//double U[256][18] = { { 1.0 } };
	double **U = new double*[256];
	for (int i = 0; i != 256; ++i)
	{
		U[i] = new double[18];
	}
	for (int i = 0; i != 256; ++i)
	{
		if (i == 0)
			U[i][0] = 1.0;
		else
			U[i][0] = 0.0;
	}

	/*int hm[23] = { 1024, 512, 1280, 640, 320, 1184, 592, 296, 1172, 1610, 1829, 1938, 1993, 996, 498, 249, 124, 62, 31, 15, 7, 3, 1 };
	double **U = new double*[2048];
	for (int i = 0; i != 2048; ++i)
	{
	U[i] = new double[24];
	}
	for (int i = 0; i != 2048; ++i)
	{
	if (i == 0)
	U[i][0] = 1.0;
	else
	U[i][0] = 0.0;
	}*/

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	for (i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		//cout<<pn1[i]<<" ";
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}//cout<<endl;


	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  //1. У��ڵ�ĳ�ʼ��
	for (i = 0; i < ROWS; i++)//42  rows/3
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}
	//cout << "iteration begins " << endl;
	int judgeFlag = 0;
	//cout << "ROWS = " << ROWS << endl;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{

			//if (i<ROWS)
			if (Contained(selectRows, i, seleRowsN))
			{
				qs = ps = M.rhead[i];
				for (int j = 1; j < 18; j++)  //����״̬�����ȸ�����(s)����(n)
				{
					for (int m = 0; m < 256; m++)
					{
						U[m][j] = U[m][j - 1] * ps->qmn0 + U[m ^ hm[j - 1]][j - 1] * ps->qmn1;

					}

					ps = ps->right;    //ָ����һ�����

				}



				for (int m = 0; m < 17; m++)
				{
					if (qs->qmn0 != qs->qmn1)
					{
						qs->outp0 = (qs->qmn0 / qs->qmn1 - U[hm[m]][17] / U[0][17]) / ((qs->qmn0 / qs->qmn1) - (qs->qmn1 / qs->qmn0));
						if (qs->outp0 < MIN) qs->outp0 = MIN;
						qs->outp1 = 1 - qs->outp0;
						if (qs->outp1 < MIN) qs->outp1 = MIN;

						qs->Lt_extr = log(qs->outp0 / qs->outp1) - log(qs->qmn0 / qs->qmn1);



						if (qs->Lt_extr > MAXLLR2)
							qs->Lt_extr = MAXLLR2;
						if (qs->Lt_extr < -MAXLLR2)
							qs->Lt_extr = -MAXLLR2;
					}
					//cout << qs->Lt_extr << " ";
					qs = qs->right; //ָ����һ�����
				}//cout<<endl;
				 //cout << endl;
			}
			else
			{
				//2.(SPC) check nodes transmit to variable nodes 
				qs = M.rhead[i];
				while (qs != NULL)
				{
					double Tmn = 1.0;
					ts = M.rhead[i]; //��������������
					while (ts != NULL)
					{
						if ((ts->col_num) != (qs->col_num))
						{
							Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
							if (Tmn > MAXLLR2) Tmn = MAXLLR2;
							if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
						}
						ts = ts->right;
					}
					qs->Lt_extr = log((1 + Tmn) / (1 - Tmn));//����������,variable nodes get value from check nodes
					if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
					if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

					qs = qs->right;
				}//end-while(qs!=NULL)	
			}

		}//end-for(0;Ldpc_Row)  

		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j < cols; j++)//147
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;
				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;
				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)

		judgeFlag = otherLDPC1.judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
			decodingflag = 1;//������ȷ
			break;
		}
	}/*********���һ֡������********/

	 /*if (itr_num == 51)
	 {
	 int bitPos[23];
	 int bitVal[23];
	 for (int i = 0; i != 23; ++i)
	 {
	 bitPos[i] = 0 + i * 6;
	 bitVal[i] = vhat[bitPos[i]];
	 }
	 otherLDPC2.spa_decoder_bitPilot(vhat, channelData, 50, sigma, bitPos, bitVal, 23);
	 }*/

	 /*for (int i = 0; i != 23; ++i)
	 {
	 cout << M.chead[0 + i * 6]->Zn << " ";
	 }
	 cout << endl;

	 cout << "iter_num = " << itr_num << endl;*/

	for (int i = 0; i != 256; ++i)
	{
		delete[]U[i];
	}
	delete[]U;

	return decodingflag;
}
int qcldpc::log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO)
{

	int i, j;
	int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	for (i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		//cout<<pn1[i]<<" ";
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}//cout<<endl;


	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  //1. У��ڵ�ĳ�ʼ��
	for (i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}
	//cout << "iteration begins " << endl;
	int judgeFlag = 0;

	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{
			//2.(SPC) check nodes transmit to variable nodes 
			qs = M.rhead[i];
			while (qs != NULL)
			{
				double Tmn = 1.0;
				ts = M.rhead[i]; //��������������
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
						if (Tmn > MAXLLR2) Tmn = MAXLLR2;
						if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					}
					ts = ts->right;
				}
				qs->Lt_extr = log((1 + Tmn) / (1 - Tmn));//����������,variable nodes get value from check nodes
				if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
				if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

				qs = qs->right;
			}//end-while(qs!=NULL)	


		}//end-for(0;Ldpc_Row)  

		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;
				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;
				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)

		judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
			decodingflag = 1;//������ȷ
			break;
		}
	}/*********���һ֡������********/

	return decodingflag;
}
int qcldpc::log_Decoder(int *vhat, double *inLLR, double *outLLR, int _MAXITERNO)
{

	int i, j;
	int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;

	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	ofstream os;
	os.open("analysis.txt");

	for (i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-inLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = inLLR[i];
	}

	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
										 //1. У��ڵ�ĳ�ʼ��
	for (i = 0; i < ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}
	//cout << "iteration begins " << endl;
	int judgeFlag = 0;

	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{
			//2.(SPC) check nodes transmit to variable nodes 
			qs = M.rhead[i];
			while (qs != NULL)
			{
				
				double Tmn = 1.0;
				ts = M.rhead[i]; //��������������
				//os << "***************" << endl;
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
						//os << ts->col_num << " " << ts->Zmn << " " << exp(ts->Zmn) << endl;
						if (Tmn > MAXLLR2) Tmn = MAXLLR2;
						if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					}
					ts = ts->right;
				}
				//os << "**************" << endl;
				//os << Tmn << endl;
				qs->Lt_extr = log((1 + Tmn) / (1 - Tmn));//����������,variable nodes get value from check nodes
				if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
				if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;
				

				
				

				qs = qs->right;
			}//end-while(qs!=NULL)	


		}//end-for(0;Ldpc_Row)  

		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j < cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;

				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;
				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;

				outLLR[j] = ps->Zn;
				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)

		judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
			decodingflag = 1;//������ȷ
			break;
		}
	}/*********���һ֡������********/

	os.close();

	return decodingflag;
}
int qcldpc::spa_decoder(int *vhat, double *channelData, int _MAXITERNO, double sigma)
{
	int len = cols, N = cols - ROWS;
	int i, j;
	int itr_num;

	double MIN = 1E-10;
	int decodingflag = 0;

	vector<double>  pn0(cols);
	vector<double>  pn1(cols);

	double h = 1.0;
	for (i = 0; i < len; i++)
	{
		pn0[i] = 1 / (1 + exp(-2 * h*channelData[i] / (sigma*sigma)));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�

		if (pn0[i]<MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i]<MIN)
			pn1[i] = MIN;
	}



	/*for(int i = 0; i != 100; ++i)
	{
	cout << pn1[i] << " ";
	}
	cout << endl;*/

	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  // ��ʼ��
	for (i = 0; i<ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps = ps->right;    //ָ����һ�����
		}

	}

	int judgeFlag = 0;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		// ˮƽ����
		// rmn�ĳ�ʼ��
		for (i = 0; i<ROWS; i++)
		{
			qs = ps = M.rhead[i];
			while (ps != NULL)
			{
				//ps->Qmn0 = (ps->qmn0+ ps->qmn1);
				ps->Qmn1 = (1 - 2 * ps->qmn1);
				ps = ps->right;
			}//while(ps!=NULL)	
			 //whileѭ����ͱ���Qmn0~Qmn3������ɣ���ʱ����Fmn0~Fmn3 �ĸ�ֵ
			while (qs != NULL)
			{
				//double Fmn0=1.0;   
				double Fmn = 1.0;
				ts = M.rhead[i];    //������������������Ե�Fmn0~Fmn3.
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						Fmn = Fmn*(ts->Qmn1);
					}
					ts = ts->right;
				}
				qs->rmn0 = (1.0 + Fmn) / 2;
				//cout << qs->rmn0 << " ";
				if (qs->rmn0<MIN) qs->rmn0 = MIN;

				qs->rmn1 = (1.0 - Fmn) / 2;
				if (qs->rmn1<MIN) qs->rmn1 = MIN;

				qs = qs->right;

			}//end-while(qs!=NULL)	
			 //cout << endl;
		}//ˮƽ����for i= 1:rows             


		 //����Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j<cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double   prod_rmn0 = 1.0;//(*(ps+symbol_in_col[i]*rows+j))
				double   prod_rmn1 = 1.0;
				qs = M.chead[j];
				//cout<<"ps->i="<<ps->i<<" ";

				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn0 = prod_rmn0*(qs->rmn0);
						prod_rmn1 = prod_rmn1*(qs->rmn1);
					}
					qs = qs->down;
				}

				double const1 = pn0[j] * prod_rmn0;
				double const2 = pn1[j] * prod_rmn1;

				ps->alphamn = 1 / (const1 + const2); // alphamnΪ��һ������ 

				ps->qmn0 = ps->alphamn*const1;
				ps->qmn1 = ps->alphamn*const2;
				//update pseudo posterior probability     
				//����ǰû�˵�newh(ones_in_col(i),j).rmn0�˽�ȥ
				//------------------------------------���qmn0 qmn1 qmn2 qmn3

				//-------------------------------------

				double const5 = const1*(ps->rmn0);
				double const6 = const2*(ps->rmn1);

				double alpha_n = 1 / (const5 + const6);
				ps->qn0 = alpha_n*const5; // ���qn0��qn1,qn2,qn3
				ps->qn1 = alpha_n*const6;
				//pn0[j] = ps->qn0;
				//pn1[j] = ps->qn1;
				//ps->qmn0 = ps->qn0;
				//ps->qmn1 = ps->qn1;
				//---------------------------tentative decoding
				double temp = 0;
				temp = (ps->qn0>ps->qn1) ? ps->qn0 : ps->qn1;
				if (ps->qn1 == temp)
					vhat[j] = 1;
				else
					vhat[j] = 0;
				ps = ps->down; //ָ����һ����㡣      
			}//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)
		  /*cout << "&&&&&&&&&&&&&" << endl;
		  ps = M.chead[0];
		  while(ps!=NULL)
		  {
		  cout << ps->qmn0 << " ";
		  ps = ps->down;
		  }
		  cout << endl;*/

		judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
			decodingflag = 1;
			break;
		}
	}/*********���һ֡������********/


	return decodingflag;
}
int qcldpc::spa_decoder_bitPilot(int *vhat, double *channelData, int _MAXITERNO, double sigma, int *bitPos, int *bitVal, int bitN)
{
	int len = cols, N = cols - ROWS;
	int i, j;
	int itr_num;

	double MIN = 1E-10;
	int decodingflag = 0;

	vector<double>  pn0(cols);
	vector<double>  pn1(cols);

	double h = 1.0;
	for (i = 0; i < len; i++)
	{
		pn0[i] = 1 / (1 + exp(-2 * h*channelData[i] / (sigma*sigma)));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�

		if (pn0[i]<MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i]<MIN)
			pn1[i] = MIN;
	}

	for (int i = 0; i != bitN; ++i)
	{
		if (bitVal[i] == 0)
		{
			pn0[bitPos[i]] = 1.0;
			pn1[bitPos[i]] = 0.0;
		}
		else
		{
			pn0[bitPos[i]] = 0.0;
			pn1[bitPos[i]] = 1.0;
		}
		vhat[bitPos[i]] = bitVal[i];
	}

	/*for(int i = 0; i != 100; ++i)
	{
	cout << pn1[i] << " ";
	}
	cout << endl;*/

	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  // ��ʼ��
	for (i = 0; i<ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps = ps->right;    //ָ����һ�����
		}

	}

	int judgeFlag = 0;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		// ˮƽ����
		// rmn�ĳ�ʼ��
		for (i = 0; i<ROWS; i++)
		{
			qs = ps = M.rhead[i];
			while (ps != NULL)
			{
				//ps->Qmn0 = (ps->qmn0+ ps->qmn1);
				ps->Qmn1 = (1 - 2 * ps->qmn1);
				ps = ps->right;
			}//while(ps!=NULL)	
			 //whileѭ����ͱ���Qmn0~Qmn3������ɣ���ʱ����Fmn0~Fmn3 �ĸ�ֵ
			while (qs != NULL)
			{
				//double Fmn0=1.0;   
				double Fmn = 1.0;
				ts = M.rhead[i];    //������������������Ե�Fmn0~Fmn3.
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						Fmn = Fmn*(ts->Qmn1);
					}
					ts = ts->right;
				}
				qs->rmn0 = (1.0 + Fmn) / 2;
				//cout << qs->rmn0 << " ";
				if (qs->rmn0<MIN) qs->rmn0 = MIN;

				qs->rmn1 = (1.0 - Fmn) / 2;
				if (qs->rmn1<MIN) qs->rmn1 = MIN;

				qs = qs->right;

			}//end-while(qs!=NULL)	
			 //cout << endl;
		}//ˮƽ����for i= 1:rows             


		 //����Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j<cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double   prod_rmn0 = 1.0;//(*(ps+symbol_in_col[i]*rows+j))
				double   prod_rmn1 = 1.0;
				qs = M.chead[j];
				//cout<<"ps->i="<<ps->i<<" ";

				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn0 = prod_rmn0*(qs->rmn0);
						prod_rmn1 = prod_rmn1*(qs->rmn1);
					}
					qs = qs->down;
				}

				double const1 = pn0[j] * prod_rmn0;
				double const2 = pn1[j] * prod_rmn1;

				ps->alphamn = 1 / (const1 + const2); // alphamnΪ��һ������ 

				ps->qmn0 = ps->alphamn*const1;
				ps->qmn1 = ps->alphamn*const2;

				if (Contained(bitPos, j, bitN))
				{
					if (vhat[j] == 0)
					{
						ps->qmn0 = 1.0;
						ps->qmn1 = 0.0;
					}
					else
					{
						ps->qmn0 = 0.0;
						ps->qmn1 = 1.0;
					}
				}
				else
				{
					double const5 = const1*(ps->rmn0);
					double const6 = const2*(ps->rmn1);

					double alpha_n = 1 / (const5 + const6);
					ps->qn0 = alpha_n*const5; // ���qn0��qn1,qn2,qn3
					ps->qn1 = alpha_n*const6;
					//pn0[j] = ps->qn0;
					//pn1[j] = ps->qn1;
					//ps->qmn0 = ps->qn0;
					//ps->qmn1 = ps->qn1;
					//---------------------------tentative decoding
					double temp = 0;
					temp = (ps->qn0 > ps->qn1) ? ps->qn0 : ps->qn1;
					if (ps->qn1 == temp)
						vhat[j] = 1;
					else
						vhat[j] = 0;
				}
				ps = ps->down; //ָ����һ����㡣      
			}//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)
		  /*cout << "&&&&&&&&&&&&&" << endl;
		  ps = M.chead[0];
		  while(ps!=NULL)
		  {
		  cout << ps->qmn0 << " ";
		  ps = ps->down;
		  }
		  cout << endl;*/

		judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
			decodingflag = 1;
			break;
		}
	}/*********���һ֡������********/
	return decodingflag;
}
int qcldpc::min_sum(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO)
{


	int i, j;
	int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;






	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	for (i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  //1. У��ڵ�ĳ�ʼ��
	for (i = 0; i < ROWS; i++)//42  rows/3
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}



	//cout << "iteration begins " << endl;
	int judgeFlag = 0;
	//cout << "ROWS = " << ROWS << endl;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//cout << "itr_num = " << itr_num << endl;
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{
			//cout << "Normal Check Row" << i << endl;
			//2.(SPC) check nodes transmit to variable nodes 
			qs = M.rhead[i];
			while (qs != NULL)
			{
				double Tmn = 1.0;
				double min = 100000.0;
				ts = M.rhead[i]; //��������������
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						/*Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
						if (Tmn > MAXLLR2) Tmn = MAXLLR2;
						if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;*/
						double temp = abs(ts->Zmn);
						if (temp < min)
						{
							min = temp;
						}
						if (ts->Zmn >= 0)
						{
							Tmn *= 1.0;
						}
						else
						{
							Tmn *= -1.0;
						}
					}
					ts = ts->right;
				}
				Tmn *= min;
				if (Tmn > MAXLLR2) Tmn = MAXLLR2;
				if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
				qs->Lt_extr = Tmn;
				//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));//����������,variable nodes get value from check nodes
				//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
				if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
				if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

				//cout << qs->Lt_extr << " ";
				qs = qs->right;
			}//end-while(qs!=NULL)	
			 //cout << endl;

		}//end-for(0;Ldpc_Row)  

		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j < cols; j++)//147
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;
				double temp = ps->Zmn;
				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;




				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)

		judgeFlag = judgeZero(vhat);
		//judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
			decodingflag++;
			//if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/

	 //vector<int> idxMin;
	 //vector<double> valueMin;
	 //const int numMin = 9;
	 //for (int i = 0; i != cols; ++i)
	 //{
	 //	vector<double>::iterator iter1 = valueMin.begin();
	 //	vector<int>::iterator iter2 = idxMin.begin();
	 //	while (iter1 != valueMin.end() && abs(M.chead[i]->Zn) >= *iter1)
	 //	{
	 //		iter1++;
	 //		iter2++;
	 //	}
	 //	valueMin.insert(iter1, abs(M.chead[i]->Zn));
	 //	idxMin.insert(iter2, i);
	 //	while (idxMin.size() > numMin)
	 //	{
	 //		idxMin.pop_back();
	 //		valueMin.pop_back();
	 //	}
	 //}
	 //for (int i = 0; i != idxMin.size(); ++i)
	 //{
	 //	valueMin[i] = M.chead[idxMin[i]]->Zn;
	 //}

	 //for (int i = 0; i != ROWS; ++i)
	 //{
	 //	ps = M.rhead[i];
	 //	while (ps != NULL)
	 //	{
	 //		ps->Zmn_pre = ps->Zmn;
	 //		ps->qmn0_pre = ps->qmn0;
	 //		ps->qmn1_pre = ps->qmn1;
	 //		ps = ps->right;
	 //	}
	 //}

	 //if (judgeFlag == 0)
	 //{
	 //	decodingflag = 2;

	 //	int idxForMin = 0;
	 //	while (judgeFlag == 0 && idxForMin != valueMin.size())
	 //	{
	 //		for (int i = 0; i != ROWS; ++i)
	 //		{
	 //			ps = M.rhead[i];
	 //			while (ps != NULL)
	 //			{
	 //				ps->Zmn = ps->Zmn_pre;
	 //				ps = ps->right;
	 //			}
	 //		}




	 //		for (int itr_num = 1; itr_num <= 50; itr_num++)
	 //		{
	 //			for (int i = 0; i < ROWS; i++)
	 //			{
	 //				qs = M.rhead[i];
	 //				while (qs != NULL)
	 //				{
	 //					double Tmn = 1.0;
	 //					double min = 100000.0;
	 //					ts = M.rhead[i];
	 //					while (ts != NULL)
	 //					{
	 //						if ((ts->col_num) != (qs->col_num))
	 //						{
	 //							double temp = abs(ts->Zmn);
	 //							if (temp < min)
	 //							{
	 //								min = temp;
	 //							}
	 //							if (ts->Zmn >= 0)
	 //							{
	 //								Tmn *= 1.0;
	 //							}
	 //							else
	 //							{
	 //								Tmn *= -1.0;
	 //							}
	 //						}
	 //						ts = ts->right;
	 //					}
	 //					Tmn *= min;
	 //					if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 //					if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 //					qs->Lt_extr = Tmn;
	 //					//////����������,variable nodes get value from check nodes
	 //					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 //					if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	 //					if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //					qs = qs->right;
	 //				}
	 //			}



	 //			//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 //			for (int j = 0; j < cols; j++)
	 //			{
	 //				
	 //				ps = M.chead[j];
	 //				while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 //				{
	 //					double prod_rmn = 0.0;
	 //					double temp = ps->Zmn;
	 //					qs = M.chead[j];
	 //					while (qs != NULL)
	 //					{
	 //						if ((qs->row_num) != (ps->row_num))
	 //						{
	 //							prod_rmn = prod_rmn + qs->Lt_extr;
	 //						}
	 //						qs = qs->down;
	 //					}


	 //					ps->Zmn = Fn[j] + prod_rmn;

	 //					if (ps->Zmn > MAXLLR2)
	 //						ps->Zmn = MAXLLR2;
	 //					if (ps->Zmn < -MAXLLR2)
	 //						ps->Zmn = -MAXLLR2;


	 //					ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 //					if (ps->qmn1 < MIN)
	 //						ps->qmn1 = MIN;

	 //					ps->qmn0 = 1 - ps->qmn1;
	 //					if (ps->qmn0 < MIN)
	 //						ps->qmn0 = MIN;

	 //					ps->Zn = ps->Zmn + ps->Lt_extr;
	 //					if (ps->Zn > MAXLLR2)
	 //						ps->Zn = MAXLLR2;
	 //					if (ps->Zn < -MAXLLR2)
	 //						ps->Zn = -MAXLLR2;

	 //					if (ps->Zn <= 0)
	 //						vhat[j] = 1;
	 //					else
	 //						vhat[j] = 0;

	 //					ps = ps->down;

	 //					
	 //				}
	 //			}

	 //		
	 //			ps = M.chead[idxMin[idxForMin]];
	 //			while (ps != NULL)
	 //			{
	 //				if (valueMin[idxForMin] >= 0)
	 //					ps->Zmn = -MAXLLR2;
	 //				else
	 //					ps->Zmn = MAXLLR2;

	 //				ps = ps->down;
	 //			}




	 //			judgeFlag = judgeZero(vhat);
	 //			//judgeFlag = judgeZero(vhat);



	 //			if (judgeFlag)
	 //			{
	 //				break;
	 //			}

	 //		}
	 //		idxForMin++;
	 //	}
	 //}




	return decodingflag;
}
int qcldpc::spa_log(int *vhat, double *TotalLLR, double sigma, int _MAXITERNO)
{


	int i, j;
	int itr_num;

	double MIN = 1E-15;
	int decodingflag = 0;






	vector<double>  pn0(cols);  //�ŵ����ݸ������ڵ�ĸ���
	vector<double>  pn1(cols);
	vector<double>  Fn(cols);

	for (i = 0; i < cols; i++)
	{
		pn0[i] = 1.0 / (1.0 + exp(-TotalLLR[i]));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�
		if (pn0[i] < MIN)
			pn0[i] = MIN;
		pn1[i] = 1 - pn0[i];
		if (pn1[i] < MIN)
			pn1[i] = MIN;
		Fn[i] = TotalLLR[i];
	}



	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  //1. У��ڵ�ĳ�ʼ��
	for (i = 0; i < ROWS; i++)//42  rows/3
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			//�����ڵ㴫��У��ڵ�ĳ�ʼ��Ϣq[i][j]=Pi(0)
			ps->qmn0 = pn0[ps->col_num];             //input0[42][7]
													 //cout<<ps->qmn0<<" ";
			ps->qmn1 = pn1[ps->col_num];             //input1[42][7]
			ps->Zmn = Fn[ps->col_num];
			//cout<<Fn[ps->col_num]<<" ";
			ps = ps->right;    //ָ����һ�����
		}//cout<<endl;
	}



	//cout << "iteration begins " << endl;
	int judgeFlag = 0;
	//cout << "ROWS = " << ROWS << endl;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		//cout << "itr_num = " << itr_num << endl;
		//haming check node update
		for (int i = 0; i < ROWS; i++)
		{
			//cout << "Normal Check Row" << i << endl;
			//2.(SPC) check nodes transmit to variable nodes 
			qs = M.rhead[i];
			while (qs != NULL)
			{
				double Tmn = 1.0;
				double min = 100000.0;
				ts = M.rhead[i]; //��������������
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						Tmn = Tmn*(exp(ts->Zmn) - 1) / (1 + exp(ts->Zmn));
						if (Tmn > MAXLLR2) Tmn = MAXLLR2;
						if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
					}
					ts = ts->right;
				}
				//����������,variable nodes get value from check nodes
				qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
				if (qs->Lt_extr>MAXLLR2) qs->Lt_extr = MAXLLR2;
				if (qs->Lt_extr<-MAXLLR2) qs->Lt_extr = -MAXLLR2;

				//cout << qs->Lt_extr << " ";
				qs = qs->right;
			}//end-while(qs!=NULL)	
			 //cout << endl;

		}//end-for(0;Ldpc_Row)  

		 //����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j < cols; j++)//147
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double prod_rmn = 0.0;
				double temp = ps->Zmn;
				qs = M.chead[j];
				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn = prod_rmn + qs->Lt_extr;
					}
					qs = qs->down;
				}
				ps->Zmn = Fn[j] + prod_rmn;
				if (ps->Zmn > MAXLLR2)
					ps->Zmn = MAXLLR2;
				if (ps->Zmn < -MAXLLR2)
					ps->Zmn = -MAXLLR2;

				/** transment ps->Zmn as inputs posibility in the next iteration **/
				ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
				if (ps->qmn1 < MIN)
					ps->qmn1 = MIN;

				ps->qmn0 = 1 - ps->qmn1;
				if (ps->qmn0 < MIN)
					ps->qmn0 = MIN;

				ps->Zn = ps->Zmn + ps->Lt_extr;
				if (ps->Zn > MAXLLR2)
					ps->Zn = MAXLLR2;
				if (ps->Zn < -MAXLLR2)
					ps->Zn = -MAXLLR2;

				if (ps->Zn < 0)
					vhat[j] = 1;
				else
					vhat[j] = 0;




				ps = ps->down; //ָ����һ�����  

			}//cout<<endl;//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)

		judgeFlag = judgeZero(vhat);
		//judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
		   //decodingflag = 1;//������ȷ
			decodingflag++;
			//if(decodingflag==2)
			break;
		}
	}/*********���һ֡������********/

	 //vector<int> idxMin;
	 //vector<double> valueMin;
	 //const int numMin = 9;
	 //for (int i = 0; i != cols; ++i)
	 //{
	 //	vector<double>::iterator iter1 = valueMin.begin();
	 //	vector<int>::iterator iter2 = idxMin.begin();
	 //	while (iter1 != valueMin.end() && abs(M.chead[i]->Zn) >= *iter1)
	 //	{
	 //		iter1++;
	 //		iter2++;
	 //	}
	 //	valueMin.insert(iter1, abs(M.chead[i]->Zn));
	 //	idxMin.insert(iter2, i);
	 //	while (idxMin.size() > numMin)
	 //	{
	 //		idxMin.pop_back();
	 //		valueMin.pop_back();
	 //	}
	 //}
	 //for (int i = 0; i != idxMin.size(); ++i)
	 //{
	 //	valueMin[i] = M.chead[idxMin[i]]->Zn;
	 //}

	 //for (int i = 0; i != ROWS; ++i)
	 //{
	 //	ps = M.rhead[i];
	 //	while (ps != NULL)
	 //	{
	 //		ps->Zmn_pre = ps->Zmn;
	 //		ps->qmn0_pre = ps->qmn0;
	 //		ps->qmn1_pre = ps->qmn1;
	 //		ps = ps->right;
	 //	}
	 //}

	 //if (judgeFlag == 0)
	 //{
	 //	decodingflag = 2;

	 //	int idxForMin = 0;
	 //	while (judgeFlag == 0 && idxForMin != valueMin.size())
	 //	{
	 //		for (int i = 0; i != ROWS; ++i)
	 //		{
	 //			ps = M.rhead[i];
	 //			while (ps != NULL)
	 //			{
	 //				ps->Zmn = ps->Zmn_pre;
	 //				ps = ps->right;
	 //			}
	 //		}




	 //		for (int itr_num = 1; itr_num <= 50; itr_num++)
	 //		{
	 //			for (int i = 0; i < ROWS; i++)
	 //			{
	 //				qs = M.rhead[i];
	 //				while (qs != NULL)
	 //				{
	 //					double Tmn = 1.0;
	 //					double min = 100000.0;
	 //					ts = M.rhead[i];
	 //					while (ts != NULL)
	 //					{
	 //						if ((ts->col_num) != (qs->col_num))
	 //						{
	 //							double temp = abs(ts->Zmn);
	 //							if (temp < min)
	 //							{
	 //								min = temp;
	 //							}
	 //							if (ts->Zmn >= 0)
	 //							{
	 //								Tmn *= 1.0;
	 //							}
	 //							else
	 //							{
	 //								Tmn *= -1.0;
	 //							}
	 //						}
	 //						ts = ts->right;
	 //					}
	 //					Tmn *= min;
	 //					if (Tmn > MAXLLR2) Tmn = MAXLLR2;
	 //					if (Tmn < -MAXLLR2) Tmn = -MAXLLR2;
	 //					qs->Lt_extr = Tmn;
	 //					//////����������,variable nodes get value from check nodes
	 //					//qs->Lt_extr = (log((1 + Tmn) / (1 - Tmn)));
	 //					if (qs->Lt_extr > MAXLLR2) qs->Lt_extr = MAXLLR2;
	 //					if (qs->Lt_extr < -MAXLLR2) qs->Lt_extr = -MAXLLR2;

	 //					qs = qs->right;
	 //				}
	 //			}



	 //			//����2������Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
	 //			for (int j = 0; j < cols; j++)
	 //			{
	 //				
	 //				ps = M.chead[j];
	 //				while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
	 //				{
	 //					double prod_rmn = 0.0;
	 //					double temp = ps->Zmn;
	 //					qs = M.chead[j];
	 //					while (qs != NULL)
	 //					{
	 //						if ((qs->row_num) != (ps->row_num))
	 //						{
	 //							prod_rmn = prod_rmn + qs->Lt_extr;
	 //						}
	 //						qs = qs->down;
	 //					}


	 //					ps->Zmn = Fn[j] + prod_rmn;

	 //					if (ps->Zmn > MAXLLR2)
	 //						ps->Zmn = MAXLLR2;
	 //					if (ps->Zmn < -MAXLLR2)
	 //						ps->Zmn = -MAXLLR2;


	 //					ps->qmn1 = 1.0 / (1.0 + exp(ps->Zmn));
	 //					if (ps->qmn1 < MIN)
	 //						ps->qmn1 = MIN;

	 //					ps->qmn0 = 1 - ps->qmn1;
	 //					if (ps->qmn0 < MIN)
	 //						ps->qmn0 = MIN;

	 //					ps->Zn = ps->Zmn + ps->Lt_extr;
	 //					if (ps->Zn > MAXLLR2)
	 //						ps->Zn = MAXLLR2;
	 //					if (ps->Zn < -MAXLLR2)
	 //						ps->Zn = -MAXLLR2;

	 //					if (ps->Zn <= 0)
	 //						vhat[j] = 1;
	 //					else
	 //						vhat[j] = 0;

	 //					ps = ps->down;

	 //					
	 //				}
	 //			}

	 //		
	 //			ps = M.chead[idxMin[idxForMin]];
	 //			while (ps != NULL)
	 //			{
	 //				if (valueMin[idxForMin] >= 0)
	 //					ps->Zmn = -MAXLLR2;
	 //				else
	 //					ps->Zmn = MAXLLR2;

	 //				ps = ps->down;
	 //			}




	 //			judgeFlag = judgeZero(vhat);
	 //			//judgeFlag = judgeZero(vhat);



	 //			if (judgeFlag)
	 //			{
	 //				break;
	 //			}

	 //		}
	 //		idxForMin++;
	 //	}
	 //}




	return decodingflag;
}


int qcldpc::spa_decoder(int *vhat, double *channelData, int _MAXITERNO, double sigma, int *vs_candidate, int *vs_candidate_hard)
{
	int len = cols, N = cols - ROWS;
	int i, j;
	int itr_num;

	double MIN = 1E-10;
	int decodingflag = 0;

	vector<double>  pn0(cols);
	vector<double>  pn1(cols);

	double h = 1.0;
	for (i = 0; i < len; i++)
	{
		if (!Contained(vs_candidate, i, 23))
		{
			pn0[i] = 1 / (1 + exp(-2 * h*channelData[i] / (sigma*sigma)));// pn:���ʣ��磺pn1(i):��iλΪ1�ĸ��ʡ�

			if (pn0[i]<MIN)
				pn0[i] = MIN;
			pn1[i] = 1 - pn0[i];
			if (pn1[i]<MIN)
				pn1[i] = MIN;
		}
		else
		{
			//pn0[i] = pn1[i] = 0.5;
			//pn0[i] = 1;//bit-pinning
			//pn1[i] = 0;//bit-pinning
			if (vs_candidate_hard[i] == 0)
			{
				pn0[i] = 1; pn1[i] = 0;
			}
			else
			{
				pn0[i] = 0; pn1[i] = 1;
			}
		}
	}

	struct OLNode4BINARY *ps, *qs, *ts;  //����һ��ָ��ʮ��������ָ��
								  // ��ʼ��
	for (i = 0; i<ROWS; i++)
	{
		ps = M.rhead[i];
		while (ps != NULL)
		{
			ps->qmn0 = pn0[ps->col_num];
			ps->qmn1 = pn1[ps->col_num];
			ps = ps->right;    //ָ����һ�����
		}

	}

	int judgeFlag = 0;
	for (itr_num = 1; itr_num <= _MAXITERNO; itr_num++)
	{
		// ˮƽ����
		// rmn�ĳ�ʼ��
		for (i = 0; i<ROWS; i++)
		{
			qs = ps = M.rhead[i];
			while (ps != NULL)
			{
				//ps->Qmn0 = (ps->qmn0+ ps->qmn1);
				ps->Qmn1 = (1 - 2 * ps->qmn1);
				ps = ps->right;
			}//while(ps!=NULL)	
			 //whileѭ����ͱ���Qmn0~Qmn3������ɣ���ʱ����Fmn0~Fmn3 �ĸ�ֵ
			while (qs != NULL)
			{
				//double Fmn0=1.0;   
				double Fmn = 1.0;
				ts = M.rhead[i];    //������������������Ե�Fmn0~Fmn3.
				while (ts != NULL)
				{
					if ((ts->col_num) != (qs->col_num))
					{
						Fmn = Fmn*(ts->Qmn1);
					}
					ts = ts->right;
				}
				qs->rmn0 = (1.0 + Fmn) / 2;
				if (qs->rmn0<MIN) qs->rmn0 = MIN;

				qs->rmn1 = (1.0 - Fmn) / 2;
				if (qs->rmn1<MIN) qs->rmn1 = MIN;

				qs = qs->right;

			}//end-while(qs!=NULL)	

		}//ˮƽ����for i= 1:rows             


		 //����Ϊ��ֱ���� ��ֱ�����ʾ��һ���������ܵ�����У����У�顣
		for (j = 0; j<cols; j++)
		{
			ps = M.chead[j];
			while (ps != NULL)	//�������е�rmn0��rmn3�ó˻���
			{
				double   prod_rmn0 = 1.0;//(*(ps+symbol_in_col[i]*rows+j))
				double   prod_rmn1 = 1.0;
				qs = M.chead[j];
				//cout<<"ps->i="<<ps->i<<" ";

				while (qs != NULL)
				{
					if ((qs->row_num) != (ps->row_num))
					{
						prod_rmn0 = prod_rmn0*(qs->rmn0);
						prod_rmn1 = prod_rmn1*(qs->rmn1);
					}
					qs = qs->down;
				}

				double const1 = pn0[j] * prod_rmn0;
				double const2 = pn1[j] * prod_rmn1;

				ps->alphamn = 1 / (const1 + const2); // alphamnΪ��һ������ 

				ps->qmn0 = ps->alphamn*const1;
				ps->qmn1 = ps->alphamn*const2;
				if (Contained(vs_candidate, j, 23))
				{
					//ps->qmn0=0.5;   
					//ps->qmn1=0.5;
					/*ps->qmn0 = 1;
					ps->qmn1 = 0;*/
					if (vs_candidate_hard[j] == 0)
					{
						ps->qmn0 = 1; ps->qmn1 = 0;
					}
					else
					{
						ps->qmn0 = 0; ps->qmn1 = 1;
					}
				}
				//update pseudo posterior probability     
				//����ǰû�˵�newh(ones_in_col(i),j).rmn0�˽�ȥ
				//------------------------------------���qmn0 qmn1 qmn2 qmn3

				//-------------------------------------

				double const5 = const1*(ps->rmn0);
				double const6 = const2*(ps->rmn1);

				double alpha_n = 1 / (const5 + const6);
				ps->qn0 = alpha_n*const5; // ���qn0��qn1,qn2,qn3
				ps->qn1 = alpha_n*const6;
				//pn0[j] = ps->qn0;
				//pn1[j] = ps->qn1;
				//ps->qmn0 = ps->qn0;
				//ps->qmn1 = ps->qn1;
				//---------------------------tentative decoding
				double temp = 0;
				temp = (ps->qn0>ps->qn1) ? ps->qn0 : ps->qn1;

				if (ps->qn1 == temp)
					vhat[j] = 1;
				else
					vhat[j] = 0;
				ps = ps->down; //ָ����һ����㡣      
			}//end- while(qs!=NULL)
		} // for(j=0;j<cols;j++)
		judgeFlag = judgeZero(vhat);

		if (judgeFlag)
		{  //---����У����˳�,���򷵻ؼ���������
			decodingflag = 1;
			break;
		}
	}/*********��ɵ�������������********/

	 //for(int i=0; i!=cols; ++i)
	 //{
	 //	if(Contained(vs_candidate, i, 23))
	 //	{
	 //		vhat[i]=-1;
	 //	}
	 //}

	return decodingflag;
}
int qcldpc::cpm_min_sum(int *vhat, double *channelData, int _MAXITERNO, double sigma, int Z)
{
	//Initialization
	double *R = new double[cols];
	for (int i = 0; i != cols; ++i)
	{
		R[i] = channelData[i] * 2.0 / sigma / sigma;
	}

	OLNode4BINARY *p, *q;
	for (int i = 0; i != ROWS; ++i)
	{
		p = M.rhead[i];
		while (p != NULL)
		{
			p->L = 0;
			p = p->right;
		}
	}

	//Iteration
	//vector<int> rowIndex(ROWS/Z);
	int *rowIndex = new int[ROWS / Z];
	bool stop = 0;
	for (int i = 0; i != _MAXITERNO; ++i)
	{

		stop = 0;
		for (int k = 0; k != Z; ++k)
		{
			for (int j = 0; j != ROWS / Z; ++j)
			{
				rowIndex[j] = j*Z;
			}
			// first step
			for (int j = 0; j != ROWS / Z; ++j)
			{
				rowIndex[j] += k;
			}
			for (int j = 0; j != cols; ++j)
			{
				double temp = 0.0;
				p = M.chead[j];
				while (p != NULL)
				{
					if (Contained(rowIndex, p->row_num, ROWS / Z))
					{
						temp += p->L;
					}
					p = p->down;
				}
				R[j] -= temp;
			}

			// second step
			for (int t = 0; t != ROWS / Z; ++t)
			{
				p = M.rhead[rowIndex[t]];

				while (p != NULL)
				{
					q = M.rhead[rowIndex[t]];
					double sign = 1.0;
					double min = 1000.0;
					while (q != NULL)
					{
						if (q->col_num != p->col_num)
						{
							if (abs(R[q->col_num])<min)
								min = abs(R[q->col_num]);
							if (R[q->col_num]<0)
								sign = 0 - sign;
						}
						q = q->right;
					}
					p->L = 0.75 * sign * min;
					p = p->right;
				}
			}

			// third step
			for (int j = 0; j != cols; ++j)
			{
				double temp = 0.0;
				p = M.chead[j];
				while (p != NULL)
				{
					if (Contained(rowIndex, p->row_num, ROWS / Z))
					{
						temp += p->L;
					}
					p = p->down;
				}
				R[j] += temp;
			}


		}
		// fourth step
		/*	for(int j = 0; j != cols; ++j)
		{
		if(R[j]>0)
		vhat[j] = 0;
		else
		vhat[j] = 1;
		}
		if(judgeZero(vhat))
		{
		stop = 1;
		break;
		}*/
	}
	for (int j = 0; j != cols; ++j)
	{
		if (R[j]>0)
			vhat[j] = 0;
		else
			vhat[j] = 1;
	}
	delete[]rowIndex;
	delete[]R;
	return 1;
}
int qcldpc::CreatMatrix_OL(char *s)
{
	int numCOL, numROW, RowWeight, numELEMENT;

	fstream in;
	in.open(s);
	in >> numCOL;
	in >> numROW;
	in >> RowWeight;
	in >> numELEMENT;
	numofnonzeroElem = numELEMENT;
	ROWS = numROW;
	cols = numCOL;



	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;

	for (int i = 0; i != ROWS; ++i)
	{
		in >> temp;
		while (temp != numCOL)
		{

			H[i*cols + temp] = 1;
			p1 = new OLNode4BINARY;
			if (p1 == NULL)
			{
				cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
				exit(1);
			}
			//initialize nodes
			//initialize nodes
			p1->row_num = i;
			p1->col_num = temp;
			//p1->e = temp;
			p1->qmn0 = 0.0;
			p1->qmn1 = 0.0;
			//p1->qmn2 = 0.0;
			// p1->qmn3 = 0.0;
			p1->Qmn0 = 0.0;
			p1->Qmn1 = 0.0;
			//p1->Qmn2 = 0.0;
			// p1->Qmn3 = 0.0;
			p1->rmn0 = 0.0;
			p1->rmn1 = 0.0;
			//p1->rmn2 = 0.0;
			//p1->rmn3 = 0.0;
			p1->qn0 = 0.0;
			p1->qn1 = 0.0;
			// p1->qn2 = 0.0;
			// p1->qn3 = 0.0;
			p1->alphamn = 1.0;

			//insert at the direction of row
			if (M.rhead[i] == NULL)
			{
				M.rhead[i] = p1;
				p1->right = NULL;
			}
			else
			{
				p2->right = p1;
				p1->right = NULL;
			}
			p2 = p1;

			//insert at the direction of column
			if (M.chead[temp] == NULL)
			{
				M.chead[temp] = p1;
				p1->down = NULL;
			}
			else
			{
				for (q = M.chead[temp]; q->down != NULL; q = q->down);
				q->down = p1;
				p1->down = NULL;
			}
			in >> temp;
		}

	}

	// cout << ROWS << " " << cols << endl;
	return 1;

}
int qcldpc::CreatMatrix_OL_NB(char *s)
{
	int numCOL, numROW, RowWeight, numELEMENT;

	fstream in;
	in.open(s);
	in >> numCOL;
	in >> numROW;
	in >> RowWeight;
	in >> numELEMENT;
	numofnonzeroElem = numELEMENT;
	ROWS = numROW;
	cols = numCOL;



	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;
	int value = 0;

	for (int i = 0; i != ROWS; ++i)
	{
		in >> temp;
		in >> value;
		// cout << temp << " " << value << " ";
		while (temp != numCOL)
		{

			H[i*cols + temp] = 1;
			p1 = new OLNode4BINARY;
			if (p1 == NULL)
			{
				cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
				exit(1);
			}
			//initialize nodes
			//initialize nodes
			p1->row_num = i;
			p1->col_num = temp;
			//p1->e = temp;
			p1->qmn0 = 0.0;
			p1->qmn1 = 0.0;
			//p1->qmn2 = 0.0;
			// p1->qmn3 = 0.0;
			p1->Qmn0 = 0.0;
			p1->Qmn1 = 0.0;
			//p1->Qmn2 = 0.0;
			// p1->Qmn3 = 0.0;
			p1->rmn0 = 0.0;
			p1->rmn1 = 0.0;
			//p1->rmn2 = 0.0;
			//p1->rmn3 = 0.0;
			p1->qn0 = 0.0;
			p1->qn1 = 0.0;
			// p1->qn2 = 0.0;
			// p1->qn3 = 0.0;
			p1->alphamn = 1.0;

			//insert at the direction of row
			if (M.rhead[i] == NULL)
			{
				M.rhead[i] = p1;
				p1->right = NULL;
			}
			else
			{
				p2->right = p1;
				p1->right = NULL;
			}
			p2 = p1;

			//insert at the direction of column
			if (M.chead[temp] == NULL)
			{
				M.chead[temp] = p1;
				p1->down = NULL;
			}
			else
			{
				for (q = M.chead[temp]; q->down != NULL; q = q->down);
				q->down = p1;
				p1->down = NULL;
			}
			in >> temp;
			in >> value;
			// cout << temp << " " << value << " ";
		}
		// cout << endl;

	}

	// cout << ROWS << " " << cols << endl;
	return 1;
}
int qcldpc::CreatMatrix_OL(char *s, int *CNodeSel, int CNodeN)
{
	int numCOL, numROW, RowWeight, numELEMENT;

	fstream in;
	in.open(s);
	in >> numCOL;
	in >> numROW;
	in >> RowWeight;
	in >> numELEMENT;
	numofnonzeroElem = numELEMENT;
	cols = numCOL;



	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;
	int i = 0;
	for (int r = 0; r != numROW; ++r)
	{
		if (Contained(CNodeSel, r, CNodeN))
		{
			in >> temp;
			while (temp != numCOL)
				in >> temp;
			continue;
		}

		in >> temp;

		while (temp != numCOL)
		{

			H[i*cols + temp] = 1;
			p1 = new OLNode4BINARY;
			if (p1 == NULL)
			{
				cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
				exit(1);
			}
			//initialize nodes
			//initialize nodes
			p1->row_num = r;
			p1->col_num = temp;
			//p1->e = temp;
			p1->qmn0 = 0.0;
			p1->qmn1 = 0.0;
			//p1->qmn2 = 0.0;
			// p1->qmn3 = 0.0;
			p1->Qmn0 = 0.0;
			p1->Qmn1 = 0.0;
			//p1->Qmn2 = 0.0;
			// p1->Qmn3 = 0.0;
			p1->rmn0 = 0.0;
			p1->rmn1 = 0.0;
			//p1->rmn2 = 0.0;
			//p1->rmn3 = 0.0;
			p1->qn0 = 0.0;
			p1->qn1 = 0.0;
			// p1->qn2 = 0.0;
			// p1->qn3 = 0.0;
			p1->alphamn = 1.0;

			//insert at the direction of row
			if (M.rhead[i] == NULL)
			{
				M.rhead[i] = p1;
				p1->right = NULL;
			}
			else
			{
				p2->right = p1;
				p1->right = NULL;
			}
			p2 = p1;

			//insert at the direction of column
			if (M.chead[temp] == NULL)
			{
				M.chead[temp] = p1;
				p1->down = NULL;
			}
			else
			{
				for (q = M.chead[temp]; q->down != NULL; q = q->down);
				q->down = p1;
				p1->down = NULL;
			}
			in >> temp;

		}

		i++;
	}

	OLNode4BINARY *ps;
	/*for (int i = 0; i != ROWS; ++i)
	{
	ps = M.rhead[i];
	while (ps != NULL)
	{
	cout << ps->row_num << " ";
	ps = ps->right;
	}
	cout << endl;

	}*/

	// cout << ROWS << " " << cols << endl;
	return 1;
}
int qcldpc::CreatMatrix_OL_Mackay(char *s)
{
	int numCOL, numROW, maxRowWeight, maxColWeight;
	ifstream in;
	// cout << s << endl;
	in.open(s);

	in >> numCOL;
	in >> numROW;
	ROWS = numROW; cols = numCOL;
	in >> maxColWeight;
	in >> maxRowWeight;
	// cout << numCOL << " " << numROW << endl;
	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[numROW];
	M.chead = new OLNode4BINARY *[numCOL];
	//Initialize the head pointer vector
	for (int i = 0; i != numROW; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != numCOL; ++i)
		M.chead[i] = NULL;

	int *variableNode_degree, *checkNode_degree;
	variableNode_degree = new int[numCOL];
	checkNode_degree = new int[numROW];

	for (int i = 0; i != numCOL; ++i)
	{
		in >> variableNode_degree[i];
	}
	for (int i = 0; i != numROW; ++i)
	{
		in >> checkNode_degree[i];
	}

	OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;

	for (int i = 0; i != numCOL; ++i)
	{
		for (int j = 0; j != maxColWeight; ++j)
		{
			in >> temp;
			if (temp == 0)
				continue;
			p1 = new OLNode4BINARY;
			if (p1 == NULL)
			{
				cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
				exit(1);
			}
			//initialize nodes
			H[(temp - 1)*numCOL + i] = 1;
			p1->row_num = temp - 1;
			p1->col_num = i;
			//insert at the direction of column
			if (M.chead[i] == NULL)
			{
				M.chead[i] = p1;
				p1->down = NULL;
			}
			else
			{
				p2->down = p1;
				p1->down = NULL;
			}
			p2 = p1;

			//insert at the direction of column
			if (M.rhead[temp - 1] == NULL)
			{
				M.rhead[temp - 1] = p1;
				p1->right = NULL;
			}
			else
			{
				for (q = M.rhead[temp - 1]; q->right != NULL; q = q->right);
				q->right = p1;
				p1->right = NULL;
			}
		}
	}
	in.close();



	delete[]variableNode_degree;
	delete[]checkNode_degree;
	return 1;
}
int qcldpc::CreatMatrix_OL_Duan(char *s)
{
	int numCOL, numROW, RowWeight, numELEMENT;

	ifstream in;
	in.open(s);
	in >> numROW;
	in >> numCOL;
	in >> RowWeight;
	in >> numELEMENT;
	ROWS = numROW;
	cols = numCOL;



	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;
	int rowIndex;
	int value;
	for (int i = 0; i != ROWS; ++i)
	{
		for (int j = 0; j != RowWeight; ++j)
		{
			in >> rowIndex;
			in >> temp;
			in >> value;
			H[i*cols + temp] = 1;
			p1 = new OLNode4BINARY;
			if (p1 == NULL)
			{
				cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
				exit(1);
			}
			//initialize nodes
			//initialize nodes
			p1->row_num = i;
			p1->col_num = temp;
			//p1->e = temp;
			p1->qmn0 = 0.0;
			p1->qmn1 = 0.0;
			//p1->qmn2 = 0.0;
			// p1->qmn3 = 0.0;
			p1->Qmn0 = 0.0;
			p1->Qmn1 = 0.0;
			//p1->Qmn2 = 0.0;
			// p1->Qmn3 = 0.0;
			p1->rmn0 = 0.0;
			p1->rmn1 = 0.0;
			//p1->rmn2 = 0.0;
			//p1->rmn3 = 0.0;
			p1->qn0 = 0.0;
			p1->qn1 = 0.0;
			// p1->qn2 = 0.0;
			// p1->qn3 = 0.0;
			p1->alphamn = 1.0;

			//insert at the direction of row
			if (M.rhead[i] == NULL)
			{
				M.rhead[i] = p1;
				p1->right = NULL;
			}
			else
			{
				p2->right = p1;
				p1->right = NULL;
			}
			p2 = p1;

			//insert at the direction of column
			if (M.chead[temp] == NULL)
			{
				M.chead[temp] = p1;
				p1->down = NULL;
			}
			else
			{
				for (q = M.chead[temp]; q->down != NULL; q = q->down);
				q->down = p1;
				p1->down = NULL;
			}
		}
		// cout << endl;
	}

	// cout << ROWS << " " << cols << endl;
	return 1;
}
int qcldpc::CreatMatrix_OL_Li(char *s)
{


	int numCOL, numROW, RowWeight, numELEMENT;

	FILE *chk;
	if ((chk = fopen(s, "r")) == NULL)
	{
		// cout << "cannot open the check matrix !" << endl;
		return 0;
	}
	fscanf(chk, "%d", &numCOL); //the number of columns
	fscanf(chk, "%d", &numROW); //the number of rows
	fscanf(chk, "%d", &RowWeight); //the max row weight 
	fscanf(chk, "%d", &numELEMENT); //the total number of elements
	numofnonzeroElem = numELEMENT;

	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;

	for (int i = 0; i != ROWS; ++i)
	{
		for (int j = 0; j != RowWeight; ++j)
		{
			fscanf(chk, "%d", &temp);
			temp = temp - 1;
			if (temp != -1)
			{
				H[i*cols + temp] = 1;
				p1 = new OLNode4BINARY;
				if (p1 == NULL)
				{
					cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
					exit(1);
				}
				//initialize nodes
				//initialize nodes
				p1->row_num = i;
				p1->col_num = temp;
				//p1->e = temp;
				p1->qmn0 = 0.0;
				p1->qmn1 = 0.0;
				//p1->qmn2 = 0.0;
				// p1->qmn3 = 0.0;
				p1->Qmn0 = 0.0;
				p1->Qmn1 = 0.0;
				//p1->Qmn2 = 0.0;
				// p1->Qmn3 = 0.0;
				p1->rmn0 = 0.0;
				p1->rmn1 = 0.0;
				//p1->rmn2 = 0.0;
				//p1->rmn3 = 0.0;
				p1->qn0 = 0.0;
				p1->qn1 = 0.0;
				// p1->qn2 = 0.0;
				// p1->qn3 = 0.0;
				p1->alphamn = 1.0;

				//insert at the direction of row
				if (M.rhead[i] == NULL)
				{
					M.rhead[i] = p1;
					p1->right = NULL;
				}
				else
				{
					p2->right = p1;
					p1->right = NULL;
				}
				p2 = p1;

				//insert at the direction of column
				if (M.chead[temp] == NULL)
				{
					M.chead[temp] = p1;
					p1->down = NULL;
				}
				else
				{
					for (q = M.chead[temp]; q->down != NULL; q = q->down);
					q->down = p1;
					p1->down = NULL;
				}
			}
		}

	}
	fclose(chk);

	return 1;

}
int qcldpc::CreatMatrix_OL_PEG(char *s)
{
	int numCOL, numROW, RowWeight;

	ifstream in;
	in.open(s);

	in >> numCOL;   //the number of columns
	in >> numROW;   //the number of rows
	in >> RowWeight;//the max row weight 

	int *rowDistri = new int[numROW];
	for (int i = 0; i != numROW; ++i)
	{
		in >> rowDistri[i];
	}



	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;

	for (int i = 0; i != ROWS; ++i)
	{
		for (int j = 0; j != RowWeight; ++j)
		{
			in >> temp;
			if (temp == 0)
				continue;
			temp = temp - 1;
			if (temp != -1)
			{
				H[i*cols + temp] = 1;
				p1 = new OLNode4BINARY;
				if (p1 == NULL)
				{
					cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
					exit(1);
				}
				//initialize nodes
				//initialize nodes
				p1->row_num = i;
				p1->col_num = temp;
				//p1->e = temp;
				p1->qmn0 = 0.0;
				p1->qmn1 = 0.0;
				//p1->qmn2 = 0.0;
				// p1->qmn3 = 0.0;
				p1->Qmn0 = 0.0;
				p1->Qmn1 = 0.0;
				//p1->Qmn2 = 0.0;
				// p1->Qmn3 = 0.0;
				p1->rmn0 = 0.0;
				p1->rmn1 = 0.0;
				//p1->rmn2 = 0.0;
				//p1->rmn3 = 0.0;
				p1->qn0 = 0.0;
				p1->qn1 = 0.0;
				// p1->qn2 = 0.0;
				// p1->qn3 = 0.0;
				p1->alphamn = 1.0;

				//insert at the direction of row
				if (M.rhead[i] == NULL)
				{
					M.rhead[i] = p1;
					p1->right = NULL;
				}
				else
				{
					p2->right = p1;
					p1->right = NULL;
				}
				p2 = p1;

				//insert at the direction of column
				if (M.chead[temp] == NULL)
				{
					M.chead[temp] = p1;
					p1->down = NULL;
				}
				else
				{
					for (q = M.chead[temp]; q->down != NULL; q = q->down);
					q->down = p1;
					p1->down = NULL;
				}
			}
		}

	}
	in.close();
	delete[]rowDistri;
	return 1;

}
int qcldpc::CreatMatrix_OL_Li(char *s, int *CNodeSel, int CNodeN)
{
	// In this function, we will deselect the CNodeN check nodes.
	int numCOL, numROW, RowWeight, numELEMENT;

	FILE *chk;
	if ((chk = fopen(s, "r")) == NULL)
	{
		cout << "cannot open the check matrix !" << endl;
		return 0;
	}
	fscanf(chk, "%d", &numCOL); //the number of columns
	fscanf(chk, "%d", &numROW); //the number of rows
	fscanf(chk, "%d", &RowWeight); //the max row weight 
	fscanf(chk, "%d", &numELEMENT); //the total number of elements
	numofnonzeroElem = numELEMENT;

	/*** construct M chainList for decoding  *****/
	M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
		M.rhead[i] = NULL;
	for (int i = 0; i != cols; ++i)
		M.chead[i] = NULL;

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;
	int i = 0;
	for (int r = 0; r != numROW; ++r)
	{
		if (Contained(CNodeSel, r, CNodeN))
		{
			for (int j = 0; j != RowWeight; ++j)
				fscanf(chk, "%d", &temp);
			continue;
		}

		for (int j = 0; j != RowWeight; ++j)
		{
			fscanf(chk, "%d", &temp);
			temp = temp - 1;
			if (temp != -1)
			{
				H[i*cols + temp] = 1;
				p1 = new OLNode4BINARY;
				if (p1 == NULL)
				{
					cout << "Function CreateMatrix_OL(): Memory exhausted !" << endl;
					exit(1);
				}
				//initialize nodes
				//initialize nodes
				p1->row_num = r;
				p1->col_num = temp;
				//p1->e = temp;
				p1->qmn0 = 0.0;
				p1->qmn1 = 0.0;
				//p1->qmn2 = 0.0;
				// p1->qmn3 = 0.0;
				p1->Qmn0 = 0.0;
				p1->Qmn1 = 0.0;
				//p1->Qmn2 = 0.0;
				// p1->Qmn3 = 0.0;
				p1->rmn0 = 0.0;
				p1->rmn1 = 0.0;
				//p1->rmn2 = 0.0;
				//p1->rmn3 = 0.0;
				p1->qn0 = 0.0;
				p1->qn1 = 0.0;
				// p1->qn2 = 0.0;
				// p1->qn3 = 0.0;
				p1->alphamn = 1.0;

				//insert at the direction of row
				if (M.rhead[i] == NULL)
				{
					M.rhead[i] = p1;
					p1->right = NULL;
				}
				else
				{
					p2->right = p1;
					p1->right = NULL;
				}
				p2 = p1;

				//insert at the direction of column
				if (M.chead[temp] == NULL)
				{
					M.chead[temp] = p1;
					p1->down = NULL;
				}
				else
				{
					for (q = M.chead[temp]; q->down != NULL; q = q->down);
					q->down = p1;
					p1->down = NULL;
				}
			}
		}
		i++;
	}
	fclose(chk);

	return 1;
}
void qcldpc::update()
{
	//if there are some redundant rows, we need update
	//the row and cols and rate after eliminating these rows
	infoLen = cols - now_rows;
	delete[]G;
	G = new int[infoLen*now_rows];
	for (int i = 0; i<infoLen; i++)
		for (int j = 0; j<now_rows; j++)
			G[i*now_rows + j] = 0;
}
bool qcldpc::HasNonZeroInTheRow(int i)
{
	for (int j = 0; j != cols; ++j)
	{
		if (H[i*cols + j] != 0)
		{
			return true;
		}
	}
	return false;
}
void qcldpc::rearrangeXS()
{
	// cout << "8888888888" << endl;
	// cout << ROWS << endl;
	os.open("H_matrix_original.txt");
	for (int i = 0; i != ROWS; ++i)
	{
		for (int j = 0; j != ROWS; ++j)
		{
			os << H[i*cols + j];
		}
		os << endl;
	}
	os.close();
	// cout << "7777777777" << endl;
	/*
	// This is Li's rearrange programming
	int  i,j;
	for(i=0; i<ROWS; i++)
	{
	int k=0;
	if(H[i*cols+i]==0)
	{
	k=i+1;
	if(k<cols-1)
	{
	while(H[i*cols+k]==0)   //��Ϊ��ʱ�˳�ѭ��
	{
	k=k+1;
	if(k==cols-1)
	break;
	}
	for(j = 0; j < ROWS; j ++)	//������Ԫ��(�ѵ�i�к͵�k�н���)
	{
	swap(H[j*cols+i],H[j*cols+k]);
	}
	rearrange[i]=k;  //record the rearrange column number
	}
	}
	else
	rearrange[i]=0;

	if(H[i*cols+i] != 0) // H may be not full-rank.
	{
	//make the A(i,i) become 1
	if(H[i*cols+i]!=1)
	{
	int d=H[i*cols+i];     // record the operation
	for(j = 0; j < cols; j ++)
	{
	H[i*cols+j]=divarray[H[i*cols+j]][d];
	}
	}
	// make the other elements in the i cols become 0 except A(i,i)
	for(j = 0; j < ROWS; j ++)
	{
	if(j != i)
	{  // to make sure the column will be zero except the major element.
	if(H[j*cols+i]!=0)
	{
	int m=H[j*cols+i];//m(i,:)*m(j,i)+m(j,:)---->m(j,:)
	for(int k = 0; k <cols; k ++)
	{
	H[j*cols+k]=add[H[j*cols+k]][mul[H[i*cols+k]][m]];
	//H[j*cols+k] = H[i*cols+k]^H[j*cols+k];
	}
	}
	}
	}
	}
	}
	*/

	// This is my programming 
	int i, j;
	int changed_row = ROWS - 1;
	for (i = 0; i != ROWS;)
	{
		int k;
		if (H[i*cols + i] == 0)
		{
			k = i + 1;
			if (k<cols - 1)
			{
				while (H[i*cols + k] == 0)
				{
					k = k + 1;
					if (k == cols - 1)
						break;
				}
				for (j = 0; j < ROWS; j++)	//������Ԫ��(�ѵ�i�к͵�k�н���)
				{
					swap(H[j*cols + i], H[j*cols + k]);
				}
				rearrange[i] = k;  //record the rearrange column number
				col_permutation4Binary cp;
				cp.col_1 = i;
				cp.col_2 = k;
				permutation_nodes.push_back(cp);
			}
		}
		else
			rearrange[i] = 0;

		//if i-th row has all zero, we need row permutations
		//changed_row = rows-1;
		if (!HasNonZeroInTheRow(i))
		{
			if (changed_row == i) // below i-th row, all has zero
			{
				// cout << "below " << i << "-th row, all has zero" << endl;
				break;
			}
			//cout << i << "-th row has all zero, we need row permutations" << endl;
			for (int k = 0; k != cols; ++k)
			{
				swap(H[changed_row*cols + k], H[i*cols + k]);
			}
			changed_row--;
			//cout << changed_row << " ";
			while (!HasNonZeroInTheRow(i) && changed_row != i)
			{
				//cout << changed_row << " ";
				for (int k = 0; k != cols; ++k)
				{
					swap(H[changed_row*cols + k], H[i*cols + k]);
				}
				changed_row--;
			}
		}


		k = i;
		if (k<cols - 1)
		{
			while (H[i*cols + k] == 0)
			{
				k = k + 1;
				if (k == cols - 1)
					break;
			}
			if (k != i)
			{
				for (j = 0; j < ROWS; j++)	//������Ԫ��(�ѵ�i�к͵�k�н���)
				{
					swap(H[j*cols + i], H[j*cols + k]);
				}
				rearrange[i] = k;  //record the rearrange column number
				col_permutation4Binary cp;
				cp.col_1 = i;
				cp.col_2 = k;
				permutation_nodes.push_back(cp);
			}
		}

		if (H[i*cols + i] != 0) // H may be not full-rank.
		{
			//make the A(i,i) become 1
			if (H[i*cols + i] != 1)
			{
				int d = H[i*cols + i];     // record the operation
				for (j = 0; j < cols; j++)
				{
					H[i*cols + j] = divarray[H[i*cols + j]][d];
				}
			}
			// make the other elements in the i cols become 0 except A(i,i)
			for (j = 0; j < ROWS; j++)
			{
				if (j != i)
				{  // to make sure the column will be zero except the major element.
					if (H[j*cols + i] != 0)
					{
						int m = H[j*cols + i];//m(i,:)*m(j,i)+m(j,:)---->m(j,:)
						for (int k = 0; k <cols; k++)
						{
							//H[j*cols+k]=add[H[j*cols+k]][mul[H[i*cols+k]][m]];
							H[j*cols + k] = H[i*cols + k] ^ H[j*cols + k];
						}
					}
				}
			}
		}
		i++;
	}

	//if there are some redundant rows, we need eliminate them
	now_rows = ROWS;
	while (!HasNonZeroInTheRow(now_rows - 1))
		now_rows--;
	update();

	os.open("H_matrix_rearrange.txt");
	os << now_rows << " " << cols << endl;
	for (int i = 0; i != now_rows; ++i)
	{
		for (int j = 0; j != now_rows; ++j)
		{
			os << H[i*cols + j];
		}
		os << endl;
	}
	os.close();

	int num = 0;
	for (int i = 0; i != now_rows; ++i)
	{
		if (H[i*cols + i] == 1)
			num++;
	}
	//cout << "the ONE number " << endl;
	//cout << num << endl;

	//for (int i = 0; i != permutation_nodes.size(); ++i)
	//{
	//	cout << "(" << permutation_nodes[i].col_1 << " " << permutation_nodes[i].col_2 << ") ";
	//}
	//cout << endl;


	/* ����н�������Ϣ
	cout<<"����:"<<endl;
	for(int ii=0;ii<cols;ii++)
	cout<<rearrange[ii]<<" ";
	cout<<endl<<endl;
	//  rhs=m;// m��[I,D]����ʽ��D��A����B�ĳ˻�
	*/
}
void qcldpc::gen_G()
{

	int i, j;

	int **D = new int *[now_rows];
	for (i = 0; i<now_rows; i++)
		D[i] = new int[cols - now_rows];
	for (i = 0; i<now_rows; i++)
		for (int j = 0; j<cols - now_rows; j++)
			D[i][j] = H[i*cols + j + now_rows];
	//G[j*rows+i]=H[i*cols+j+rows];


	//G[(j-rows)*rows+i]=H[i*cols+j+rows]; //���D��---A������B�ĳ˻�,��G��D��ת��
	for (i = 0; i<cols - now_rows; i++)
	{
		for (j = 0; j<now_rows; j++)
		{
			G[i*now_rows + j] = D[j][i];
			//fprintf(G_FILE,"%2d",G[i*rows+j]);
		}

		//fprintf(G_FILE,"\n");
	}

	ofstream out;
	out.open("G.txt");

	out << cols - now_rows << " " << now_rows << endl;
	for (i = 0; i < cols - now_rows; ++i)
	{
		for (j = 0; j < now_rows; ++j)
		{
			out << G[i*now_rows + j] << " ";
		}
		out << endl;
	}
	out.close();



	for (i = 0; i<now_rows; i++)
	{
		delete[]D[i];
		D[i] = NULL;
	}
	delete[]D;
	D = NULL;
}
void qcldpc::Encoder(int * Mesg, int *CodeWord)
{
	int *Check = new int[codeLen - infoLen];
	//int N=cols-rows;
	int i;
	int CheckLen = now_rows;
	for (i = 0; i < CheckLen; i++)
		Check[i] = 0;
	ArrayMultiply4Binary(Check, Mesg, G, cols - now_rows, now_rows);
	for (i = 0; i < CheckLen; i++)
		CodeWord[i] = Check[i];
	for (i = CheckLen; i < codeLen; i++)
		CodeWord[i] = Mesg[i - CheckLen];
	//ReorderSymbol(CodeWord);
	//Reorder_bits(CodeWord, codeLen);
	Reorder_bits(CodeWord);
	//CheckEncoder(CodeWord);
	delete[]Check;
	Check = NULL;

}

void qcldpc::colSwap(int COL1, int COL2)
{

}

/*
void ldpc::extract_mesg(int *res, int *ref)
{
int i;
for(i = 0; i < rows; i++)
{
if (rearrange[i] != 0)
swap(ref[i], ref[rearrange[i]]);
}
for (i = 0; i < cols-rows; i++)
res[i] = ref[rows + i];
}
*/
void qcldpc::extract_mesg(double *res, double *ref)
{
	for (int i = 0; i != permutation_nodes.size(); ++i)
	{
		int temp = ref[permutation_nodes[i].col_1];
		ref[permutation_nodes[i].col_1] = ref[permutation_nodes[i].col_2];
		ref[permutation_nodes[i].col_2] = temp;
	}
	for (int i = 0; i != cols - now_rows; ++i)
	{
		res[i] = ref[now_rows + i];
	}
}


int qcldpc::Error_Count(int *Mesg, int *Dec_Mesg)
{
	int error_num = 0;
	for (int i = 0; i < cols - now_rows; i++)
	{
		if (Mesg[i] != Dec_Mesg[i])
			error_num++;
	}

	return error_num;
}
int qcldpc::judgeZero(int *vs)
{



	struct OLNode4BINARY *ps;
	int judgeFlag = 1;
	for (int i = 0; i<ROWS; i++)
	{
		ps = M.rhead[i];
		int sum = 0;
		while (ps != NULL)
		{

			if (vs[ps->col_num] != 0)
			{
				sum ^= 1;//add[sum][mul[ps->e][vs[ps->col_num]]];
			}
			ps = ps->right;
		}

		// ��ʱk4��h(i,:)*vhat;
		if (sum != 0)
		{
			judgeFlag = 0;
			break;
		}
	}
	return(judgeFlag);
}

void qcldpc::Reorder_bits(int *u, int len)
{
	for (int i = len; i>0; i--)
	{
		if (rearrange[i - 1] != 0)
			swap(u[i - 1], u[rearrange[i - 1]]);
	}
}
void qcldpc::Reorder_bits(int *u)
{
	for (int i = permutation_nodes.size() - 1; i != -1; --i)
	{
		int temp = u[permutation_nodes[i].col_1];
		u[permutation_nodes[i].col_1] = u[permutation_nodes[i].col_2];
		u[permutation_nodes[i].col_2] = temp;
	}
}
void qcldpc::Reorder_bits(double *u)
{
	for (int i = permutation_nodes.size() - 1; i != -1; --i)
	{
		double temp = u[permutation_nodes[i].col_1];
		u[permutation_nodes[i].col_1] = u[permutation_nodes[i].col_2];
		u[permutation_nodes[i].col_2] = temp;
	}
}


int qcldpc::destroy_crossList()
{
	struct OLNode4BINARY *p1, *p2;
	for (int i = 0; i != ROWS; ++i)
	{
		p1 = M.rhead[i]->right;
		M.rhead[i] = NULL;
		while (p1 != NULL)
		{
			p2 = p1;
			p1 = p1->right;
			delete p2;
		}

	}
	for (int j = 0; j<cols; j++)
	{
		M.chead[j] = NULL;
	}
	delete[] M.rhead;
	delete[] M.chead;

	return 1;
}

void qcldpc::add_i2t(int i, int t)
{
	for (int k = 0; k != cols; ++k)
		H[t*cols + k] = H[t*cols + k] ^ H[i*cols + k];
}



ldpc::ldpc(int rowN, int colN, int q)
{
	infoLen = 1;
	codeLen = cols;
	ROWS = rowN;
	cols = colN;
	now_rows = ROWS;
	qnumber = q;//q-array
}
ldpc::~ldpc()
{
	delete[]numInCol;
	if (destroy_crossList())
	{
		cout << "succssfully delete the crossList !" << endl;
	}
}
int ldpc::CreatMatrix_OL_Li(char *s)
{
	int numCOL, numROW, RowWeight, numELEMENT;

	FILE *chk;
	if ((chk = fopen(s, "r")) == NULL)
	{
		cout << "cannot open the check matrix !" << endl;
		return 0;
	}

	fscanf(chk, "%d", &numCOL); //the number of columns
	fscanf(chk, "%d", &numROW); //the number of rows
	fscanf(chk, "%d", &RowWeight); //the max row weight 
	fscanf(chk, "%d", &numELEMENT); //the total number of elements
	numofnonzeroElem = numELEMENT;



	/*** construct M chainList for decoding  *****/
	/*M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	M.rtail = new OLNode4BINARY *[ROWS];
	M.ctail = new OLNode4BINARY *[cols];*/

	M.rhead = new OLNode4BINARY[ROWS];
	M.chead = new OLNode4BINARY[cols];
	M.rtail = new OLNode4BINARY[ROWS];
	M.ctail = new OLNode4BINARY[cols];

	// cout << "%%%%%%%%%%%%%%%%%%%%%" << endl;

	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
	{
		M.rhead[i].left = NULL;
		M.rhead[i].col_num = -1;
		M.rhead[i].row_num = i;
	}



	for (int i = 0; i != cols; ++i)
	{
		M.chead[i].up = NULL;
		M.chead[i].row_num = -1;
		M.chead[i].col_num = i;
	}
	for (int i = 0; i != ROWS; ++i)
	{
		M.rtail[i].right = NULL;
		M.rtail[i].col_num = cols;
		M.rtail[i].row_num = i;
	}
	for (int i = 0; i != cols; ++i)
	{
		M.ctail[i].down = NULL;
		M.ctail[i].row_num = ROWS;
		M.ctail[i].col_num = i;
	}

	// cout << "&&&&&&&&&&&&&&&&&" << endl;

	numInCol = new int[cols];
	for (int i = 0; i != cols; ++i)
	{
		numInCol[i] = 0;
	}

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;

	for (int i = 0; i != ROWS; ++i)
	{
		for (int j = 0; j != RowWeight; ++j)
		{
			fscanf(chk, "%d", &temp);
			temp = temp - 1;
			if (temp != -1)
			{
				p1 = new OLNode4BINARY;

				p1->row_num = i;
				p1->col_num = temp;


				//insert at the direction of row
				if (j == 0)
				{
					M.rhead[i].right = p1;
					p1->right = NULL;
				}
				else
				{
					p2->right = p1;
					p1->right = NULL;
				}
				p2 = p1;

				if (j == RowWeight - 1)
				{
					p2->right = &M.rtail[i];
				}

				//insert at the direction of column
				if (numInCol[temp] == 0)
				{
					M.chead[temp].down = p1;
					p1->down = NULL;
				}
				else
				{
					for (q = &M.chead[temp]; q->down != NULL; q = q->down);
					q->down = p1;
					p1->down = NULL;
				}
				numInCol[temp]++;
			}
		}

	}
	for (int i = 0; i != cols; ++i)
	{
		int j;
		p1 = &M.chead[i];
		while (p1->down != NULL)
		{
			p1 = p1->down;
		}
		p1->down = &M.ctail[i];
		M.ctail[i].up = p1;
	}

	for (int i = 0; i != ROWS; ++i)
	{
		p1 = &M.rhead[i];
		p2 = p1->right;
		while (p2 != &M.rtail[i])
		{
			p2->left = p1;
			p1 = p2;
			p2 = p2->right;
		}
	}

	for (int i = 0; i != cols; ++i)
	{
		p1 = &M.chead[i];
		p2 = p1->down;
		while (p2 != &M.ctail[i])
		{
			p2->up = p1;
			p1 = p2;
			p2 = p2->down;
		}
	}

	fclose(chk);



	return 1;

}
int ldpc::CreatMatrix_OL_Huang(char *s)
{
	int numCOL, numROW, RowWeight, numELEMENT;

	FILE *chk;
	if ((chk = fopen(s, "r")) == NULL)
	{
		cout << "cannot open the check matrix !" << endl;
		return 0;
	}

	fscanf(chk, "%d", &numCOL); //the number of columns
	fscanf(chk, "%d", &numROW); //the number of rows
	fscanf(chk, "%d", &RowWeight); //the max row weight 
	fscanf(chk, "%d", &numELEMENT); //the total number of elements
	numofnonzeroElem = numELEMENT;



	/*** construct M chainList for decoding  *****/
	/*M.rhead = new OLNode4BINARY *[ROWS];
	M.chead = new OLNode4BINARY *[cols];
	M.rtail = new OLNode4BINARY *[ROWS];
	M.ctail = new OLNode4BINARY *[cols];*/

	M.rhead = new OLNode4BINARY[ROWS];
	M.chead = new OLNode4BINARY[cols];
	M.rtail = new OLNode4BINARY[ROWS];
	M.ctail = new OLNode4BINARY[cols];



	//Initialize the head pointer vector
	for (int i = 0; i != ROWS; ++i)
	{
		M.rhead[i].left = NULL;
		M.rhead[i].right = NULL;
		M.rhead[i].col_num = -1;
		M.rhead[i].row_num = i;

	}
	for (int i = 0; i != cols; ++i)
	{
		M.chead[i].up = NULL;
		M.chead[i].row_num = -1;
		M.chead[i].col_num = i;
	}
	for (int i = 0; i != ROWS; ++i)
	{
		M.rtail[i].right = NULL;
		M.rtail[i].col_num = cols;
		M.rtail[i].row_num = i;
	}
	for (int i = 0; i != cols; ++i)
	{
		M.ctail[i].down = NULL;
		M.ctail[i].row_num = ROWS;
		M.ctail[i].col_num = i;
	}



	numInCol = new int[cols];
	for (int i = 0; i != cols; ++i)
	{
		numInCol[i] = 0;
	}

	struct OLNode4BINARY *p1, *p2, *q;
	p1 = p2 = q = NULL;
	int temp = 0;

	for (int i = 0; i != ROWS; ++i)
	{

		fscanf(chk, "%d", &temp);
		while (temp != numCOL)
		{
			p1 = new OLNode4BINARY;

			//			cout << temp << " ";

			p1->row_num = i;
			p1->col_num = temp;

			//insert at the direction of row
			if (M.rhead[i].right == NULL)
			{
				M.rhead[i].right = p1;
				p1->right = NULL;
			}
			else
			{
				p2->right = p1;
				p1->right = NULL;
			}
			p2 = p1;


			//insert at the direction of column
			if (numInCol[temp] == 0)
			{
				M.chead[temp].down = p1;
				p1->down = NULL;
			}
			else
			{
				for (q = &M.chead[temp]; q->down != NULL; q = q->down);
				q->down = p1;
				p1->down = NULL;
			}
			numInCol[temp]++;

			fscanf(chk, "%d", &temp);
		}
		//		cout << endl;
	}

	for (int i = 0; i != cols; ++i)
	{
		if (numInCol[i] == 0)
		{
			M.chead[i].down = &M.ctail[i];
			M.ctail[i].up = &M.chead[i];
		}
	}

	//I think there is at least one 1 in every row

	for (int i = 0; i != ROWS; ++i)
	{
		p1 = M.rhead[i].right;
		p2 = M.rhead[i].right;
		while (p1 != NULL)
		{
			p2 = p1;
			p1 = p1->right;
		}
		p2->right = &M.rtail[i];
	}

	// cout << "cols = " << cols << endl;

	for (int i = 0; i != cols; ++i)
	{
		int j;

		p1 = &M.chead[i];
		while (p1->down != NULL)
		{
			// if (p1->col_num == cols - 1)
				// cout << p1->row_num << " ";
			p1 = p1->down;
			// if (p1->col_num == cols - 1)
				// cout << " && ";
		}

		// cout << i << " ";
		p1->down = &M.ctail[i];
		M.ctail[i].up = p1;


	}

	// cout << "&&&&&&&&&&&&&" << endl;

	for (int i = 0; i != ROWS; ++i)
	{
		p1 = &M.rhead[i];
		p2 = p1->right;
		while (p2 != &M.rtail[i])
		{
			p2->left = p1;
			p1 = p2;
			p2 = p2->right;
		}
	}

	for (int i = 0; i != cols; ++i)
	{
		p1 = &M.chead[i];
		p2 = p1->down;
		while (p2 != &M.ctail[i])
		{
			p2->up = p1;
			p1 = p2;
			p2 = p2->down;
		}
	}

	fclose(chk);



	return 1;

}
int ldpc::destroy_crossList()
{

	struct OLNode4BINARY *p1, *p2;
	for (int i = 0; i != ROWS; ++i)
	{
		p1 = M.rhead[i].right;
		while (p1 != &M.rtail[i])
		{
			p2 = p1;
			p1 = p1->right;
			delete p2;
		}

	}

	delete[] M.rhead;
	delete[] M.chead;
	delete[] M.rtail;
	delete[] M.ctail;

	return 1;
}
void ldpc::colSwap(int COL1, int COL2)
{
	if (COL1 == COL2)
		return;
	vector<OLNode4BINARY*> nodeArrCol1;
	vector<OLNode4BINARY*> nodeArrCol2;
	OLNode4BINARY *p1, *p2;
	int i, j;



	for (p1 = M.chead[COL1].down; p1 != &M.ctail[COL1]; p1 = p1->down)
	{
		p1->left->right = p1->right;
		p1->right->left = p1->left;
		nodeArrCol1.push_back(p1);
	}
	for (p1 = M.chead[COL2].down; p1 != &M.ctail[COL2]; p1 = p1->down)
	{
		p1->left->right = p1->right;
		p1->right->left = p1->left;
		nodeArrCol2.push_back(p1);
	}
	M.chead[COL1].down = &M.ctail[COL1];
	M.ctail[COL1].up = &M.chead[COL1];
	M.chead[COL2].down = &M.ctail[COL2];
	M.ctail[COL2].up = &M.chead[COL2];



	if (nodeArrCol1.size() != 0)
	{
		nodeArrCol1[0]->up = &M.chead[COL2];
		M.chead[COL2].down = nodeArrCol1[0];
		nodeArrCol1[nodeArrCol1.size() - 1]->down = &M.ctail[COL2];
		M.ctail[COL2].up = nodeArrCol1[nodeArrCol1.size() - 1];

		for (int i = 0; i != nodeArrCol1.size(); ++i)
		{
			int rowIndex = nodeArrCol1[i]->row_num;
			for (p1 = &M.rhead[rowIndex], p2 = M.rhead[rowIndex].right; p2 != NULL; p1 = p1->right, p2 = p2->right)
			{
				if (p1->col_num < COL2 && p2->col_num > COL2)
				{
					nodeArrCol1[i]->col_num = COL2;
					nodeArrCol1[i]->left = p1;
					nodeArrCol1[i]->right = p2;
					p1->right = nodeArrCol1[i];
					p2->left = nodeArrCol1[i];
					break;
				}
			}
		}
	}

	if (nodeArrCol2.size() != 0)
	{
		nodeArrCol2[0]->up = &M.chead[COL1];
		M.chead[COL1].down = nodeArrCol2[0];
		nodeArrCol2[nodeArrCol2.size() - 1]->down = &M.ctail[COL1];
		M.ctail[COL1].up = nodeArrCol2[nodeArrCol2.size() - 1];

		for (int i = 0; i != nodeArrCol2.size(); ++i)
		{
			int rowIndex = nodeArrCol2[i]->row_num;
			for (p1 = &M.rhead[rowIndex], p2 = M.rhead[rowIndex].right; p2 != NULL; p1 = p1->right, p2 = p2->right)
			{
				if (p1->col_num < COL1 && p2->col_num > COL1)
				{
					nodeArrCol2[i]->col_num = COL1;
					nodeArrCol2[i]->left = p1;
					nodeArrCol2[i]->right = p2;
					p1->right = nodeArrCol2[i];
					p2->left = nodeArrCol2[i];
					break;
				}
			}
		}
	}

}
void ldpc::rowSwap(int ROW1, int ROW2)
{
	if (ROW1 == ROW2)
		return;
	vector<OLNode4BINARY*> nodeArrRow1;
	vector<OLNode4BINARY*> nodeArrRow2;
	OLNode4BINARY *p1, *p2;
	for (p1 = M.rhead[ROW1].right; p1 != &M.rtail[ROW1]; p1 = p1->right)
	{
		p1->up->down = p1->down;
		p1->down->up = p1->up;
		nodeArrRow1.push_back(p1);
	}
	for (p1 = M.rhead[ROW2].right; p1 != &M.rtail[ROW2]; p1 = p1->right)
	{
		p1->up->down = p1->down;
		p1->down->up = p1->up;
		nodeArrRow2.push_back(p1);
	}

	M.rhead[ROW1].right = &M.rtail[ROW1];
	M.rtail[ROW1].left = &M.rhead[ROW1];
	M.rhead[ROW2].right = &M.rtail[ROW2];
	M.rtail[ROW2].left = &M.rhead[ROW2];



	if (nodeArrRow1.size() != 0)
	{
		nodeArrRow1[0]->left = &M.rhead[ROW2];
		M.rhead[ROW2].right = nodeArrRow1[0];
		nodeArrRow1[nodeArrRow1.size() - 1]->right = &M.rtail[ROW2];
		M.rtail[ROW2].left = nodeArrRow1[nodeArrRow1.size() - 1];

		for (int i = 0; i != nodeArrRow1.size(); ++i)
		{
			int colIndex = nodeArrRow1[i]->col_num;
			for (p1 = &M.chead[colIndex], p2 = M.chead[colIndex].down; p2 != NULL; p1 = p1->down, p2 = p2->down)
			{
				if (p1->row_num<ROW2 && p2->row_num>ROW2)
				{
					nodeArrRow1[i]->row_num = ROW2;
					nodeArrRow1[i]->up = p1;
					nodeArrRow1[i]->down = p2;
					p1->down = nodeArrRow1[i];
					p2->up = nodeArrRow1[i];
					break;
				}
			}
		}
	}



	if (nodeArrRow2.size() != 0)
	{
		nodeArrRow2[0]->left = &M.rhead[ROW1];
		M.rhead[ROW1].right = nodeArrRow2[0];
		nodeArrRow2[nodeArrRow2.size() - 1]->right = &M.rtail[ROW1];
		M.rtail[ROW1].left = nodeArrRow2[nodeArrRow2.size() - 1];



		for (int i = 0; i != nodeArrRow2.size(); ++i)
		{
			int colIndex = nodeArrRow2[i]->col_num;
			for (p1 = &M.chead[colIndex], p2 = M.chead[colIndex].down; p2 != NULL; p1 = p1->down, p2 = p2->down)
			{
				if (p1->row_num<ROW1 && p2->row_num>ROW1)
				{
					nodeArrRow2[i]->row_num = ROW1;
					nodeArrRow2[i]->up = p1;
					nodeArrRow2[i]->down = p2;
					p1->down = nodeArrRow2[i];
					p2->up = nodeArrRow2[i];
					break;
				}
			}
		}
	}



}
bool ldpc::RU_PREPROCESS(int colErased, int MaxRow)
{

	OLNode4BINARY *p;

	//	cout << "the number of columns earsed: " << colErased << endl;
	//	int ColErased = colErased;
	//	int ColErased = colErased - 1;
	int rowNum = 0;
	int startRow = 0;
	int startCol = colErased + 1;
	vector<int> Omega;
	for (int i = startRow; i != MaxRow; ++i)
	{
		int numOne = 0;
		p = M.rhead[i].right;
		while (p != &M.rtail[i])
		{
			if (p->col_num >= startCol)
				numOne++;
			p = p->right;
		}
		if (numOne == 1)
			Omega.push_back(i);
	}
	int N = 0;
	while (Omega.size() != 0)
	{
		N++;
		//cout << " Omega size = " << Omega.size() << endl;
		vector<int> colFirst(Omega.size(), -1);
		for (int i = 0; i != Omega.size(); ++i)
		{
			p = M.rhead[Omega[i]].right;
			while (p != &M.rtail[Omega[i]])
			{
				if (p->col_num >= startCol)
				{
					colFirst[i] = p->col_num;
					break;
				}
				p = p->right;
			}
		}
		vector<int> OmegaTemp;
		vector<int> colFirstTemp;
		for (int i = 0; i != Omega.size(); ++i)
		{
			bool isThere = false;
			for (int j = 0; j != colFirstTemp.size(); ++j)
			{
				if (colFirstTemp[j] == colFirst[i])
				{
					isThere = true;
					break;
				}
			}
			if (!isThere)
			{
				colFirstTemp.push_back(colFirst[i]);
				OmegaTemp.push_back(Omega[i]);
			}
		}
		Omega = OmegaTemp;

		rowNum += Omega.size();

		/*cout << "Omega ..." << endl;
		for (int i = 0; i != Omega.size(); ++i)
		{
		cout << Omega[i] << " ";
		}
		cout << endl;
		cout << "colFirstTemp ..." << endl;
		for (int i = 0; i != colFirstTemp.size(); ++i)
		{
		cout << colFirstTemp[i] << " ";
		}
		cout << endl;*/
		for (int i = 0; i != Omega.size(); ++i)
		{
			rowSwap(startRow + i, Omega[i]);
		}

		/*cout << "col swap" << endl;
		for (int i = 0; i != Omega.size(); ++i)
		{
		cout << "(" << startCol + i << "," << colFirstTemp[i] << ") ";
		colSwap(startCol + i, colFirstTemp[i]);
		}*/

		for (int i = 0; i != Omega.size(); ++i)
		{
			int colIndex;
			p = M.rhead[startRow + i].right;
			while (p != &M.rtail[startRow + i])
			{
				if (p->col_num >= startCol)
				{
					colIndex = p->col_num;
					break;
				}
				p = p->right;
			}
			colSwap(startCol + i, colIndex);
		}

		startCol += Omega.size();
		startRow += Omega.size();

		Omega.clear();
		for (int i = startRow; i < MaxRow; ++i)
		{
			int numOne = 0;
			p = M.rhead[i].right;
			while (p != &M.rtail[i])
			{
				if (p->col_num >= startCol)
					numOne++;
				p = p->right;
			}
			if (numOne == 1)
				Omega.push_back(i);
		}

	}





	int m_g = rowNum;
	for (int col = cols - 1, j = rowNum - 1; col != cols - rowNum - 1; --col, --j)
	{
		vector<int> rowNeedChanged;
		p = M.chead[col].down;
		while (p != &M.ctail[col])
		{
			if (p->row_num >= rowNum)
			{
				rowNeedChanged.push_back(p->row_num);
			}
			p = p->down;
		}
		for (int i = 0; i != rowNeedChanged.size(); ++i)
		{
			addRow1ToRow2(j, rowNeedChanged[i]);
		}
	}

	colSwap(42, cols - rowNum - 1);

	ofstream os;
	os.open("H_RU.txt");

	int *rowValue = new int[cols];

	for (int i = 0; i != ROWS; ++i)
	{
		if (i == rowNum)
			os << endl;
		for (int j = 0; j != cols; ++j)
		{
			rowValue[j] = 0;
		}
		p = M.rhead[i].right;
		while (p != &M.rtail[i])
		{
			rowValue[p->col_num] = 1;
			p = p->right;
		}
		for (int j = 0; j != cols; ++j)
		{
			if (j == cols - ROWS)
				os << "   " << rowValue[j];
			else if (j == cols - rowNum)
				os << "   " << rowValue[j];
			else
				os << rowValue[j];

		}
		os << endl;
	}

	delete[]rowValue;

	os.close();


	int g = ROWS - rowNum;
	int *g_rowValue = new int[g];

	// cout << "g" << g << endl;

	os.open("H_g.txt");
	os << g << endl;
	os << g << endl;
	os << g << endl;
	os << g << endl;

	for (int i = rowNum; i != ROWS; ++i)
	{
		for (int j = 0; j != g; ++j)
		{
			g_rowValue[j] = 0;
		}
		p = M.rhead[i].right;
		while (p != &M.rtail[i])
		{
			if (p->col_num >= cols - ROWS && p->col_num <= cols - ROWS + g - 1)
			{
				g_rowValue[p->col_num - (cols - ROWS)] = 1;
			}
			p = p->right;
		}

		for (int j = 0; j != g; ++j)
		{
			if (g_rowValue[j] == 1)
				os << j << " ";
		}
		os << g << endl;
	}

	delete[]g_rowValue;
	os.close();



	ldpc gH(g, g, 2);
	gH.CreatMatrix_OL_Huang("H_g.txt");

	// cout << "g = " << g << endl;

	bool isSingular = false;
	for (int i = 0; i != g; ++i)
	{
		int rowIndex = -1;
		p = gH.M.chead[i].down;
		while (p != &gH.M.ctail[i])
		{
			if (p->row_num >= i)
			{
				rowIndex = p->row_num;
				break;
			}
			p = p->down;
		}
		if (rowIndex == -1)
		{
			isSingular = true;
			break;
		}



		gH.rowSwap(i, rowIndex);

		p = gH.M.chead[i].down;
		vector<int> needChange;
		while (p != &gH.M.ctail[i])
		{
			if (p->row_num > i)
			{
				needChange.push_back(p->row_num);
			}
			p = p->down;
		}

		for (int j = 0; j != needChange.size(); ++j)
		{
			gH.addRow1ToRow2(i, needChange[j]);
		}



	}

	ofstream os1;
	os1.open("H_gT.txt");
	int *g_rowValueT = new int[g];
	for (int j = 0; j != g; ++j)
	{
		for (int k = 0; k != g; ++k)
		{
			g_rowValueT[k] = 0;
		}
		p = gH.M.rhead[j].right;
		while (p != &gH.M.rtail[j])
		{
			g_rowValueT[p->col_num] = 1;
			p = p->right;
		}

		for (int k = 0; k != g; ++k)
		{
			os1 << g_rowValueT[k] << " ";
		}
		os1 << endl;
	}
	delete[]g_rowValueT;
	os1.close();

	return isSingular;
}
int ldpc::getElement(int row, int col)
{
	int value = 0;

	OLNode4BINARY *p;
	p = M.rhead[row].right;
	while (p != &M.rtail[row])
	{
		if (p->col_num == col)
		{
			value = 1;
			break;
		}
		p = p->right;
	}

	return value;
}
void ldpc::addRow1ToRow2(int row1, int row2)
{
	if (row1 == row2)
		return;
	OLNode4BINARY *p1, *p2, *p3;
	p1 = M.rhead[row1].right;

	while (p1 != &M.rtail[row1])
	{
		p2 = &M.rhead[row2];
		p3 = M.rhead[row2].right;
		while (p3 != NULL)
		{

			if (p2->col_num < p1->col_num && p3->col_num > p1->col_num)
			{
				OLNode4BINARY *p = new OLNode4BINARY;
				p->col_num = p1->col_num;
				p->row_num = row2;
				p2->right = p;
				p->left = p2;
				p3->left = p;
				p->right = p3;
				OLNode4BINARY *q1, *q2;
				q1 = &M.chead[p->col_num];
				q2 = M.chead[p->col_num].down;
				while (q2 != NULL)
				{

					if (q1->row_num < p->row_num && q2->row_num > p->row_num)
					{
						q1->down = p;
						p->up = q1;
						q2->up = p;
						p->down = q2;
						break;
					}
					q1 = q1->down;
					q2 = q2->down;
				}
				break;
			}
			if (p2->col_num < p1->col_num && p3->col_num == p1->col_num)
			{
				OLNode4BINARY *p = p3;
				p->left->right = p->right;
				p->right->left = p->left;
				p->up->down = p->down;
				p->down->up = p->up;
				delete p;
				break;
			}
			p2 = p2->right;
			p3 = p3->right;
		}
		p1 = p1->right;
	}
}