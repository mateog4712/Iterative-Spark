/**
* @mainpage
*
* Space-efficient sparse variant of an RNA (loop-based) free energy
* minimization algorithm (RNA folding equivalent to the Zuker
* algorithm).
*
* The results are equivalent to RNAfold.
*/
#define NDEBUG
#include "base_types.hh"
#include "cmdline.hh"
#include "Spark.cc"
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <cstring>
#include <string>
#include <cassert>

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}


// std::string remove_structure_intersection(std::string restricted, std::string structure){
// 	cand_pos_t length = structure.length();
// 	for(cand_pos_t i=0; i< length; ++i){
// 		if(restricted[i] == '(' || restricted[i] == ')') structure[i] = '.';
		
// 		if (output[i] == '['){
// 			output[i] = '(';
// 		}
// 		if (output[i] == ']'){
// 			output[i] = ')';
// 		}
// 	}
// 	return output;
// }
// /**
//  * @brief returns a vector of pairs which represent the start and end indices for each disjoint substructure in the structure
//  * 
//  * @param CL_ Candidate list
//  * @return total number of candidates
//  */
// std::string find_disjoint_substructure(std::string structure){
// 	cand_pos_t length = structure.length();
// 	string restricted= std::string (n,'.')
// 	for(cand_pos_t i = 0; i< length;++i){
// 		// if()
// 	}
// }

// std::string obtainRelaxedStems(std::string restricted, std::string pkfree_structure, sparse_tree tree){
// 	cand_pos_t length = restricted.length();

// 	//Gresult <- G1
// 	std::string relaxed = pkfree_structure;

// 	cand_pos_t i = 0;
// 	cand_pos_t j = 0;

// 	for(cand_pos_t k=0;k<length;k++){
// 		if(tree.tree[k] > -1){
// 			i = k;
// 			j = G2_pair[k];
// 			if(i < j){ //for each ij in G2
// 				if( (G1[i] != G2[i]) && (G1[j] != G2[j]) ){//if ij not in G1
// 					//include bulges of size 1
// 					if(paired_structure(i-1,j+1,G1_pair,length) || paired_structure(i+1,j-1,G1_pair,length) ){
// 						Gresult[i] = G2[i];
// 						Gresult[j] = G2[j];
// 					//include loops of size 1x1
// 					}else if( paired_structure(i-2,j+1,G1_pair,length) || paired_structure(i-1,j+2,G1_pair,length) || \
// 							paired_structure(i+1,j-2,G1_pair,length) || paired_structure(i+2,j-1,G1_pair,length) ){
// 						Gresult[i] = G2[i];
// 						Gresult[j] = G2[j];
// 					//include loops of size 1x2 or 2x1
// 					}else if( paired_structure(i-2,j+2,G1_pair,length) || paired_structure(i+2,j-2,G1_pair,length) ){
// 						Gresult[i] = G2[i];
// 						Gresult[j] = G2[j];
// 					}else if( paired_structure(i-3,j+2,G1_pair,length) || paired_structure(i-2,j+3,G1_pair,length) || \
// 							paired_structure(i+2,j-3,G1_pair,length) || paired_structure(i+3,j-2,G1_pair,length) ){

// 						Gresult[i] = G2[i];
// 						Gresult[j] = G2[j];
// 					}
// 				}
// 			}
// 		}
// 	}
// }


/**
* @brief Simple driver for @see Spark.
*
* Reads sequence from command line or stdin and calls folding and
* trace-back methods of SparseMFEFold.
*/
int main(int argc,char **argv) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
	std::getline(std::cin,seq);
	}
	cand_pos_t n = seq.length();

	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (n,'.');

	if(restricted != "" && restricted.length() != n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		std::cout << seq << std::endl;
		std::cout << restricted << std::endl;
		exit(0);
	}

	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}

	

	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;
	// sparse_tree tree(restricted,n);

	// noGU = args_info.noGU_given;
	// seqtoRNA(seq);

	// SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted);

	// if(args_info.dangles_given) sparsemfefold.params_->model_details.dangles = dangle_model;
	// pseudoknot = ~args_info.pseudoknot_free_given;
	// pk_only = args_info.pk_only_given;  
	
	// cmdline_parser_free(&args_info);

	// int count = 0;
	// for(int i = 1;i<=n;++i){
	// 	if(tree.tree[i].pair > i || (tree.tree[i].pair < i && tree.tree[i].pair > 0)) count = 4;
	// 	if(tree.tree[i].pair < 0 && count > 0){
	// 		sparsemfefold.WI_Bbp[i] = (5 - count)*PUP_penalty;
	// 		// sparsemfefold.WIP_Bbp[i] = (5 - count)*PUP_penalty;
	// 		count--;
	// 	}
	// }
	// count = 0;
	// for(int i = n;i>0;--i){
	// 	if(tree.tree[i].pair > i || (tree.tree[i].pair < i && tree.tree[i].pair > 0)) count = 4;
	// 	if(tree.tree[i].pair < 0 && count > 0){
	// 		sparsemfefold.WI_Bp[i] = (5 - count)*PUP_penalty;
	// 		count--;
	// 	}
	// }
	
	// energy_t mfe = fold(sparsemfefold.seq_, tree,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.CLVP_,sparsemfefold.CLBE_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.taVP_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.VP_, sparsemfefold.WVe_, sparsemfefold.WV_, sparsemfefold.dwvp_, sparsemfefold.WMB_,sparsemfefold.dwmbi_,sparsemfefold.WMBP_,sparsemfefold.WMBA_,sparsemfefold.WI_,sparsemfefold.dwi1_,sparsemfefold.WIP_,sparsemfefold.dwip1_,sparsemfefold.WI_Bbp,sparsemfefold.WIP_Bbp,sparsemfefold.WIP_Bp,sparsemfefold.WI_Bp,sparsemfefold.n_,sparsemfefold.garbage_collect_);		
	// std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.CLBE_,sparsemfefold.CLVP_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.taVP_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.WI_Bbp,sparsemfefold.WIP_Bbp,sparsemfefold.WIP_Bp,sparsemfefold.WI_Bp,sparsemfefold.n_,tree, mark_candidates);
	
	// std::ostringstream smfe;
	// smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;
	// std::cout << seq << std::endl;

	// std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

	return 0;
}