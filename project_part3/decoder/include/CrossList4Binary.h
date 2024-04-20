#ifndef   H_CROSSLIST4BINARY
#define   H_CROSSLIST4BINARY
//# include <time.h>
struct OLNode4BINARY {
	int row_num, col_num, e;       //the subcript of row and col of non_zero elements

	double pr;
	double lr;


	double qmn0;
	double qmn1;
	double Zmn;//log(qmn0/qmn1)
	double Zmn_pre;
	double qmn0_pre;
	double qmn1_pre;

	double outp0;
	double outp1;
	double Lt_extr;
	double temp;
	double Lt_extr_str;

	double Qmn0;
	double Qmn1;
	//double Qmn2;
	//double Qmn3;
	double rmn0;
	double rmn1;
	// double Lmn;  //log(rmn0/rmn1)
	//double rmn2;
	//double rmn3;

	double qn0;
	double qn1;
	double Zn;
	//double qn2;
	//double qn3;

	//float R; // for CPM-RMSA
	double L; // for CPM-RMSA

	double alphamn; //.
	double wmn; // for WBF decoding
	int flag; //check whether the row_num'th check node be satisfied or not
	struct OLNode4BINARY   *right, *down;  //the pointers of row and col
	struct OLNode4BINARY   *left, *up;     //the 

};//OLNode;*Olink;

typedef struct CrossList4BINARY {
	struct OLNode4BINARY     **rhead, **chead;    //the head pointer vector of rows and cols
	struct OLNode4BINARY     **rtail, **ctail;    //the head pointer vector of rows and cols
										   //int          mu, nu, tu;       //
}CrossList4BINARY;                   //CrossList

typedef struct CrossListV24BINARY {
	struct OLNode4BINARY     *rhead, *chead;    //the head pointer vector of rows and cols
	struct OLNode4BINARY     *rtail, *ctail;    //the head pointer vector of rows and cols
										 //int          mu, nu, tu;       //
}CrossListV24BINARY;                   //CrossList





#endif


