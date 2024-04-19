#ifndef H_bLDPC
#define H_bLDPC
#include <vector>
#include <fstream>
using namespace std;
#include "CrossList4Binary.h"


const double MAXLLR2 = 40;

struct col_permutation4Binary
{
	int col_1;
	int col_2;
};

class qcldpc
{
public:
	//input: base matrix row number,
	//       base matrix column number, 
	//       CPM size and q-ary
	qcldpc(int, int, int, int);
	//delete the memory
	~qcldpc();
	//according the current H matrix, 
	//to compute the infoLen and codeLen and rate
	//and to resize the G matrix
	void update();
	// to creat the cross List
	int CreatMatrix_OL(char *s);
	int CreatMatrix_OL_NB(char *s);
	int CreatMatrix_OL(char *s, int *CNodeSel, int CNodeN);
	int CreatMatrix_OL_Li(char *s);
	int CreatMatrix_OL_Li(char *s, int *CNodeSel, int CNodeN);
	int CreatMatrix_OL_Mackay(char *s);
	int CreatMatrix_OL_Duan(char *s);
	int CreatMatrix_OL_PEG(char *s);
	// to create the G matrix
	void gen_G();
	void add_i2t(int, int);
	// using Gaussian Elimination to transform the H matrix
	void rearrangeXS();
	// check whether the i-th row has non-zero elements
	bool HasNonZeroInTheRow(int i);

	void colSwap(int COL1, int COL2);
	void rowSwap(int ROW1, int ROW2);


	int get_infoLen() { return infoLen; }
	int get_codeLen() { return codeLen; }
	int get_rows() { return base_rows*CPM_size; } // the row number of H matrix
	int get_cols() { return base_cols*CPM_size; } // the column number of H matrix

												  // transforming the meessage
	void Encoder(int*, int*);
	// to check whether the codeword is a legal codeword
	int judgeZero(int *);
	// getting the message from the codeword after decoded.
	void extract_mesg(double *, double *);
	// computing the difference between orignial message and the message after decoded
	int Error_Count(int*, int*);
	// changing the order of codeword
	void Reorder_bits(int*, int);
	// changing the order of codeword using structure I defined
	void Reorder_bits(int*);
	void Reorder_bits(double*);
	// decoding algorithm: sum-product algorithm 
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO);
	int log_Decoder(int *vhat, double *inLLR,double *outLLR, int _MAXITERNO);
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, int *selectRows, int seleRowsN, qcldpc &otherLDPC1, qcldpc &otherLDPC2);
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected);
	int min_sum(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected);
	int min_sum_v2(int *CodeWord, int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected);
	int min_sum_v3(int *CodeWord, int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected);
	int min_sum_v4(int *CodeWord, int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected);


	int min_sum(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO);
	int spa_log(int *vhat, double *TotalLLR, double sigma, int _MAXITERNO);
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected, int *QRvalue, int *SPCvalue);
	int log_Decoder_v2(int *vhat, int *code, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, vector<int> &checkNodeSelected);
	int spa_decoder(int *vhat, double *channelData, int _MAXITERNO, double sigma);
	int spa_decoder_bitPilot(int *vhat, double *channelData, int _MAXITERNO, double sigma, int *bitPos, int *bitVal, int bitN);
	int spa_decoder(int *vhat, double *channelData, int _MAXITERNO, double sigma, int *vs_condidate, int *vs_candi_hard);

	int cpm_min_sum(int *vhat, double *chanelData, int _MAXITERNO, double sigma, int Z);

	bool Contained(int *vector, int ele, int length);
	bool Contained(vector<int> &vector, int ele);
	// the crossList which is to store the sparse H matrix
	CrossList4BINARY  M;
	int destroy_crossList(); // delete the cross List
							 //private:
public:
	// file stream operator os
	ofstream os;
	// the message length and the code length
	// they are the final version after eliminating 
	// the redundant rows
	int  infoLen, codeLen;
	// the base matrix row number and column number and cpm size
	int  base_rows, base_cols, CPM_size;
	// the original H matrix' rows and column
	int ROWS, cols;
	// the rows of H matrix after eliminating redundant rows
	int now_rows;
	//int  rows, original_rows, cols, num;
	// number of non-zero elements in H matrix
	int  numofnonzeroElem;
	// q-ary
	int  qnumber;
	int  *G;
	int  *H;
	int  *rearrange;

	vector<col_permutation4Binary> permutation_nodes;
};

struct unit
{
	unit()
	{
		colIndex = -1;
		a = vector<int>();
	}
	int colIndex;
	vector<int> a;
};

class ldpc
{
public:
	//input: base matrix row number,
	//       base matrix column number, 
	//       CPM size and q-ary
	ldpc(int, int, int);
	//delete the memory
	~ldpc();
	//according the current H matrix, 
	//to compute the infoLen and codeLen and rate
	//and to resize the G matrix
	void update();
	// to creat the cross List
	int CreatMatrix_OL(char *s);
	int CreatMatrix_OL(char *s, int *CNodeSel, int CNodeN);
	int CreatMatrix_OL_Li(char *s);
	int CreatMatrix_OL_Huang(char *s);
	int CreatMatrix_OL_Li(char *s, int *CNodeSel, int CNodeN);
	int CreatMatrix_OL_Mackay(char *s);
	int CreatMatrix_OL_Duan(char *s);
	int CreatMatrix_OL_PEG(char *s);
	// to create the G matrix
	void gen_G();
	void add_i2t(int, int);
	// using Gaussian Elimination to transform the H matrix
	void rearrangeXS();
	// check whether the i-th row has non-zero elements
	bool HasNonZeroInTheRow(int i);

	void colSwap(int COL1, int COL2);
	void rowSwap(int ROW1, int ROW2);
	void addRow1ToRow2(int row1, int row2);
	bool RU_PREPROCESS(int, int);
	int getElement(int row, int col);
	bool isSingular();

	int get_infoLen() { return infoLen; }
	int get_codeLen() { return codeLen; }
	int get_rows() { return base_rows*CPM_size; } // the row number of H matrix
	int get_cols() { return base_cols*CPM_size; } // the column number of H matrix

												  // transforming the meessage
	void Encoder(int*, int*);
	// to check whether the codeword is a legal codeword
	int judgeZero(int *);
	// getting the message from the codeword after decoded.
	void extract_mesg(double *, double *);
	// computing the difference between orignial message and the message after decoded
	int Error_Count(int*, int*);
	// changing the order of codeword
	void Reorder_bits(int*, int);
	// changing the order of codeword using structure I defined
	void Reorder_bits(int*);
	void Reorder_bits(double*);
	// decoding algorithm: sum-product algorithm 
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO);
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, int *selectRows, int seleRowsN, qcldpc &otherLDPC1, qcldpc &otherLDPC2);
	int log_Decoder(int *vhat, double *TotalLLR, double *channelData, double sigma, int _MAXITERNO, qcldpc &otherLDPC1, qcldpc &otherLDPC2);
	int spa_decoder(int *vhat, double *channelData, int _MAXITERNO, double sigma);
	int spa_decoder_bitPilot(int *vhat, double *channelData, int _MAXITERNO, double sigma, int *bitPos, int *bitVal, int bitN);
	int spa_decoder(int *vhat, double *channelData, int _MAXITERNO, double sigma, int *vs_condidate, int *vs_candi_hard);

	int cpm_min_sum(int *vhat, double *chanelData, int _MAXITERNO, double sigma, int Z);

	bool Contained(int *vector, int ele, int length);
	// the crossList which is to store the sparse H matrix
	CrossListV24BINARY  M;
	int destroy_crossList(); // delete the cross List
							 //private:
public:
	// file stream operator os
	ofstream os;
	// the message length and the code length
	// they are the final version after eliminating 
	// the redundant rows
	int  infoLen, codeLen;
	// the base matrix row number and column number and cpm size
	int  base_rows, base_cols, CPM_size;
	// the original H matrix' rows and column
	int ROWS, cols;
	// the rows of H matrix after eliminating redundant rows
	int now_rows;
	//int  rows, original_rows, cols, num;
	// number of non-zero elements in H matrix
	int  numofnonzeroElem;
	// q-ary
	int  qnumber;
	int  *G;
	int  *H;
	int  *rearrange;
	int *numInCol;
	vector<col_permutation4Binary> permutation_nodes;
};

static int divarray[2][2] = { { 0,0 },{ 0,1 } };
static int add[2][2] = { { 0,1 },{ 1,0 } };
static int mul[2][2] = { { 0,0 },{ 0,1 } };

#endif

