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
#include "Spark.hh"
#include "sparse_tree.cc"
#include <sstream>
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


std::string remove_structure_intersection(std::string restricted, std::string structure){
	cand_pos_t length = structure.length();
	for(cand_pos_t i=0; i< length; ++i){
		if(restricted[i] == '(' || restricted[i] == ')') structure[i] = '.';
		
		if (structure[i] == '['){
			structure[i] = '(';
		}
		if (structure[i] == ']'){
			structure[i] = ')';
		}
	}
	return structure;
}
// /**
//  * @brief returns a vector of pairs which represent the start and end indices for each disjoint substructure in the structure
//  * 
//  * @param CL_ Candidate list
//  * @return total number of candidates
//  */
std::string find_disjoint_substructure(std::string structure, std::vector<cand_pos_t>){
// 	cand_pos_t length = structure.length();
// 	string restricted= std::string (n,'.')
// 	for(cand_pos_t i = 0; i< length;++i){
// 		// if()
// 	}
}

std::string obtainRelaxedStems(std::string restricted, std::string pkfree_structure, sparse_tree &tree){
	cand_pos_t length = restricted.length();

	//Gresult <- G1
	std::string relaxed = pkfree_structure;

	cand_pos_t i = 0;
	cand_pos_t j = 0;

	// for(cand_pos_t k=0;k<length;k++){
	// 	if(tree.tree[k] > -1){
	// 		i = k;
	// 		j = tree.tree[k].pair;
	// 		if(i < j){ //for each ij in G2
	// 			if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
	// 				//include bulges of size 1
	// 				if(paired_structure(i-1,j+1,G1_pair,length) || paired_structure(i+1,j-1,G1_pair,length) ){
	// 					Gresult[i] = G2[i];
	// 					Gresult[j] = G2[j];
	// 				//include loops of size 1x1
	// 				}else if( paired_structure(i-2,j+1,G1_pair,length) || paired_structure(i-1,j+2,G1_pair,length) || \
	// 						paired_structure(i+1,j-2,G1_pair,length) || paired_structure(i+2,j-1,G1_pair,length) ){
	// 					Gresult[i] = G2[i];
	// 					Gresult[j] = G2[j];
	// 				//include loops of size 1x2 or 2x1
	// 				}else if( paired_structure(i-2,j+2,G1_pair,length) || paired_structure(i+2,j-2,G1_pair,length) ){
	// 					Gresult[i] = G2[i];
	// 					Gresult[j] = G2[j];
	// 				}else if( paired_structure(i-3,j+2,G1_pair,length) || paired_structure(i-2,j+3,G1_pair,length) || \
	// 						paired_structure(i+2,j-3,G1_pair,length) || paired_structure(i+3,j-2,G1_pair,length) ){

	// 					Gresult[i] = G2[i];
	// 					Gresult[j] = G2[j];
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	return relaxed;
}
void seqtoRNA(std::string &sequence){
	bool DNA = false;
    for (char &c : sequence) {
      	if (c == 'T' || c == 't') {
			c = 'U';
			DNA = true;
		}
    }
	noGU = DNA;
}

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

	seqtoRNA(seq);

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
	Dangle dangle = 2;
	if(args_info.dangles_given) dangle = dangles;

	sparse_tree tree(restricted,n);
	

	bool verbose;
	verbose = args_info.verbose_given;

	energy_t method1_energy = INF;
	energy_t method2_energy = INF;
	energy_t method3_energy = INF;
	energy_t method4_energy = INF;

	//Method1
	std::string method1_structure = Spark(seq,restricted,method1_energy,dangle,false,true,file);


	//Method2
	std::string temp = Spark(seq,restricted,method2_energy,dangle,true,true,file);
	temp = remove_structure_intersection(restricted,temp);
	std::string method2_structure = Spark(seq,temp,method2_energy,dangle,false,true,file);

	//Method3
	std::string temp2 = Spark(seq,restricted,method3_energy,dangle,false,false,file);
	std::string relaxed = obtainRelaxedStems(restricted,temp2,tree);
	std::string pk_structure = Spark(seq,relaxed,method3_energy,dangle,true,true,file);
	temp2 = remove_structure_intersection(relaxed,pk_structure);
	std::string method3_structure = Spark(seq,temp2,method3_energy,dangle,false,true,file);



	//Method4
	std::vector< std::pair<int,int> > disjoint_substructure_index;



	// energy_t energy = std::min(std::min(std::min(method1_energy,method2_energy),method3_energy),method4_energy);
	std::cout << seq << std::endl;
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method1_energy/100.0 ;

	std::ostringstream smfe2;
	smfe2 << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method2_energy/100.0 ;

	std::ostringstream smfe3;
	smfe3 << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method3_energy/100.0 ;

	std::cout << "Method1: " << method1_structure << " ("<<smfe.str()<<")"<<std::endl;
	std::cout << "Method2: " << method2_structure << " ("<<smfe2.str()<<")"<<std::endl;
	std::cout << "Method3: " << method3_structure << " ("<<smfe3.str()<<")"<<std::endl;

	return 0;
}