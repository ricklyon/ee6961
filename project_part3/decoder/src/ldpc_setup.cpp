#include "ldpc_setup.h"
#include <math.h>
qcldpc* ldpc;

uint64_t ldpc_initial(uint32_t checkLen, uint32_t codeLen, char* input_buf) {

	ldpc = new qcldpc(checkLen, codeLen, 1, 2);
	(*ldpc).CreatMatrix_OL_NB(input_buf);
	(*ldpc).rearrangeXS();
	(*ldpc).gen_G();
	uint64_t ptr = reinterpret_cast<uint64_t>(ldpc);

	return ptr;

}

uint64_t get_ldpc_infoLen(uint64_t& ptr) {
	qcldpc* ldpc_ptr = reinterpret_cast<qcldpc*>(ptr);
	uint64_t infoLen = (*ldpc).infoLen;
	return infoLen;
}


void ldpc_clear(uint64_t &ptr) {
	qcldpc* ldpc_ptr = reinterpret_cast<qcldpc*>(ptr);
	delete ldpc_ptr;
}



vector<uint32_t> ldpc_encoder(uint64_t &ptr, vector<uint32_t> doubleMesg, uint32_t mesgLen, uint32_t codeLen) {
	qcldpc* ldpc_ptr = reinterpret_cast<qcldpc*>(ptr);

	int* mesg = new int[mesgLen];
	int* code = new int[codeLen];
	for (int i = 0; i < mesgLen; ++i)
	{
		mesg[i] = (int)doubleMesg[i];
	}
	(*ldpc_ptr).Encoder(mesg, code);

	vector<uint32_t> output_p(codeLen);
	for (int i = 0; i < codeLen; ++i)
		output_p[i] = code[i];

	delete[]mesg;
	delete[]code;
	return output_p;
}



vector<double> ldpc_decoder(uint64_t& ptr, vector<double> LeBits, uint32_t codeLen_bit) {
	qcldpc* ldpc_ptr = reinterpret_cast<qcldpc*>(ptr);

	int bitsGF = 1;
	int infoLen_bit = (*ldpc_ptr).infoLen * bitsGF;

	double* inLLR = new double[codeLen_bit];
	double* outLLR = new double[codeLen_bit];
	int* vhat = new int[codeLen_bit];
	for (int i = 0; i < codeLen_bit; ++i) inLLR[i] = LeBits[i];
	(*ldpc_ptr).log_Decoder(vhat, inLLR, outLLR, 100);


	vector<double> output_APP_c(codeLen_bit);
	for (int i = 0; i < codeLen_bit; ++i)
		output_APP_c[i] = outLLR[i];


	delete[]inLLR;
	delete[]outLLR;
	delete[]vhat;
	return output_APP_c;
}

vector<double> ldpc_Ext(uint64_t& ptr, vector<double> LeBits) {
	qcldpc* ldpc_ptr = reinterpret_cast<qcldpc*>(ptr);

	int bitsGF = 1;
	int infoLen_bit = (*ldpc_ptr).infoLen * bitsGF;
	int codeLen_bit = (*ldpc_ptr).codeLen * bitsGF;

	double* inLLR = new double[codeLen_bit];
	double* outLLR = new double[infoLen_bit];
	
	for (int i = 0; i < codeLen_bit; ++i) inLLR[i] = LeBits[i];
	(*ldpc_ptr).extract_mesg(outLLR, inLLR);


	vector<double> output_APP_m(infoLen_bit);
	for (int i = 0; i < infoLen_bit; ++i)
		output_APP_m[i] = outLLR[i];


	delete[]inLLR;
	delete[]outLLR;
	return output_APP_m;
}
