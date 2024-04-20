#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>,
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "bLDPC.h"

//PYBIND11_MAKE_OPAQUE(vector<uint8_t>);
//PYBIND11_MAKE_OPAQUE(vector<uint32_t>);
//PYBIND11_MAKE_OPAQUE(vector<int>);
//PYBIND11_MAKE_OPAQUE(vector<double>);

/* library of functions */


uint64_t ldpc_initial(uint32_t checkLen, uint32_t codeLen, char* input_buf);
uint64_t get_ldpc_infoLen(uint64_t& ptr);
void ldpc_clear(uint64_t& ptr);
vector<uint32_t> ldpc_encoder(uint64_t& ptr, vector<uint32_t> doubleMesg, uint32_t mesgLen, uint32_t codeLen);
vector<double> ldpc_decoder(uint64_t& ptr, vector<double> LeBits, uint32_t codeLen_bit);
vector<double> ldpc_Ext(uint64_t& ptr, vector<double> LeBits);

namespace py = pybind11;

PYBIND11_MODULE(pybind_11_ldpc_setup, mod) {
	mod.def("ldpc_initial", &ldpc_initial, "ldpc initialization.");
	mod.def("get_ldpc_infoLen", &get_ldpc_infoLen, "get_ldpc_infoLen");
	mod.def("ldpc_clear", &ldpc_clear, "ldpc clear");
	mod.def("ldpc_encoder", &ldpc_encoder, "ldpc_encoder");
	mod.def("ldpc_decoder", &ldpc_decoder, "ldpc_decoder");
	mod.def("ldpc_Ext", &ldpc_Ext, "ldpc_Ext");
}