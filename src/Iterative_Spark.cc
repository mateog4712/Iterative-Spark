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
		
		if (structure[i] == '[') structure[i] = '(';
		
		if (structure[i] == ']') structure[i] = ')';
	}
	return structure;
}
// /**
//  * @brief returns a vector of pairs which represent the start and end indices for each disjoint substructure in the structure
//  * 
//  * @param CL_ Candidate list
//  * @return total number of candidates
//  */
void find_disjoint_substructure(std::string structure, std::vector< std::pair<int,int> > &pair_vector){
	cand_pos_t n = structure.length();
	cand_pos_t count = 0;
	cand_pos_t i = 0;
	for(cand_pos_t k=0; k<n;++k){
		if(structure[k] == '('){
			if(count == 0) i = k;
			count++;

		}else if(structure[k] == ')'){
			count--;
			if(count == 0){
				std::pair <int,int> ij_pair (i,k);
				pair_vector.push_back(ij_pair);
			}
		}
	}
}
/**
 * @brief Fills the pair array
 * p_table will contain the index of each base pair
 * X or x tells the program the base cannot pair and . sets it as unpaired but can pair
 * @param structure Input structure
 * @param p_table Restricted array
 */
void detect_pairs(const std::string &structure, std::vector<cand_pos_t> &p_table){
	cand_pos_t i, j, count = 0, length = structure.length(),last_j=length;
	std::vector<cand_pos_t>  pairs;
	pairs.push_back(length);

	for (i=length-1; i >=0; --i){
		if ((structure[i] == 'x') || (structure[i] == 'X'))
			p_table[i] = -1;
		else if (structure[i] == '.')
			p_table[i] = -2;
		if (structure[i] == ')'){
			pairs.push_back(i);
			count++;
		}
		if (structure[i] == '('){
			j = pairs[pairs.size()-1];
			pairs.erase(pairs.end()-1);
			p_table[i] = j;
			p_table[j] = i;
			count--;
		}
	}
	pairs.pop_back();
	if (pairs.size() != 0) std::cout << pairs[0] << std::endl << structure << std::endl;
	if (pairs.size() != 0)
	{
		fprintf (stderr, "The given structure is not valid: more left parentheses than right parentheses: \n");
		exit (1);
	}
}

cand_pos_t paired_structure(cand_pos_t i, cand_pos_t j, std::vector<cand_pos_t> &pair_index){
	cand_pos_t n = pair_index.size();
	return (i >= 0 && j < n && (pair_index[i] == j));
}

std::string obtainRelaxedStems(std::string restricted, std::string pkfree_structure){
	cand_pos_t n = restricted.length();

	//Gresult <- G1
	std::string relaxed = restricted;

	cand_pos_t i = 0;
	cand_pos_t j = 0;
	
	std::vector<cand_pos_t> G1_pair;
	std::vector<cand_pos_t> G2_pair;
	G1_pair.resize(n,-2);
	G2_pair.resize(n,-2);
	detect_pairs(restricted,G1_pair);
	detect_pairs(pkfree_structure,G2_pair);

	
	for(int k=0;k<n;++k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i-1,j+1,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include bulges of size 1
					}else if(paired_structure(i-2,j+1,G1_pair) || paired_structure(i-1,j+2,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i-2,j+2,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-3,j+2,G1_pair) || paired_structure(i-2,j+3,G1_pair)){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}

	for(int k=n-1;k>=0;--k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (restricted[i] != pkfree_structure[i])){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i+1,j-1,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
						
					//include bulges of size 1
					}else if(paired_structure(i+1,j-2,G1_pair) || paired_structure(i+2,j-1,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i+2,j-2,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if(paired_structure(i+2,j-3,G1_pair) || paired_structure(i+3,j-2,G1_pair) ){
						relaxed[i] = pkfree_structure[i];
						relaxed[j] = pkfree_structure[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}
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
std::string remove_x(std::string structure){
	for(char &c: structure) {if(c == 'x') c = '.';}
	return structure; 
}

std::string method2(std::string &seq, std::string &restricted, energy_t &method2_energy, Dangle dangle, std::string file){

	std::string pk_only_output = Spark(seq,restricted,method2_energy,dangle,true,true,file);
	std::string pk_free_removed = remove_structure_intersection(restricted,pk_only_output);
	std::string no_x_restricted = remove_x(restricted);

	if(pk_only_output != no_x_restricted) return Spark(seq,pk_free_removed,method2_energy,dangle,false,true,file);
	else return pk_only_output;
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
	if(args_info.dangles_given) dangle = dangle_model;

	sparse_tree tree(restricted,n);
	

	bool verbose;
	verbose = args_info.verbose_given;

	energy_t method1_energy = INF;
	energy_t method2_energy = INF;
	energy_t method3_energy = INF;
	energy_t method4_energy = INF;
	energy_t final_energy = INF;

	std::string final_structure;

	//Method1
	std::string method1_structure = Spark(seq,restricted,method1_energy,dangle,false,true,file);
	if(method1_energy < final_energy){
		final_energy = method1_energy;
		final_structure=method1_structure;
	}

	//Method2
	std::string method2_structure = method2(seq,restricted,method2_energy,dangle,file);
	if(method2_energy < final_energy){
		final_energy = method2_energy;
		final_structure=method2_structure;
	}
	//Method3
	std::string pk_free = Spark(seq,restricted,method3_energy,dangle,false,false,file);
	std::string relaxed = obtainRelaxedStems(restricted,pk_free);
	std::string method3_structure = method2(seq,relaxed,method3_energy,dangle,file);
	if(method3_energy < final_energy){
		final_energy = method3_energy;
		final_structure=method3_structure;
	}



	//Method4
	std::vector< std::pair<int,int> > disjoint_substructure_index;
	find_disjoint_substructure(restricted,disjoint_substructure_index);
	std::string disjoint_structure = restricted;
	for(auto current_substructure_index : disjoint_substructure_index){
		cand_pos_t i = current_substructure_index.first;
		cand_pos_t j = current_substructure_index.second;
		energy_t energy = INF;

		std::string subsequence = seq.substr(i,j-i+1);
		std::string substructure = restricted.substr(i,j-i+1);

		std::string pk_free = Spark(subsequence,substructure,energy,dangle,false,false,file);
		std::string relaxed = obtainRelaxedStems(substructure,pk_free);
		disjoint_structure.replace(i,j-i+1,relaxed);
	}
	std::string method4_structure = method2(seq,disjoint_structure,method4_energy,dangle,file);
	if(method4_energy < final_energy){
		final_energy = method4_energy;
		final_structure=method4_structure;
	}



	std::cout << seq << std::endl;
	
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << final_energy/100.0 ;
	std::cout << final_structure << " ("<<smfe.str()<<")"<<std::endl;
	
	
	if(args_info.verbose_given){
	std::ostringstream smfe1;
	smfe1 << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method1_energy/100.0 ;

	std::ostringstream smfe2;
	smfe2 << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method2_energy/100.0 ;

	std::ostringstream smfe3;
	smfe3 << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method3_energy/100.0 ;

	std::ostringstream smfe4;
	smfe4 << std::setiosflags(std::ios::fixed) << std::setprecision(2) << method4_energy/100.0 ;

	std::cout << "Method1: " << method1_structure << " ("<<smfe1.str()<<")"<<std::endl;
	std::cout << "Method2: " << method2_structure << " ("<<smfe2.str()<<")"<<std::endl;
	std::cout << "Method3: " << method3_structure << " ("<<smfe3.str()<<")"<<std::endl;
	std::cout << "Method4: " << method4_structure << " ("<<smfe4.str()<<")"<<std::endl;
	}
	
	cmdline_parser_free(&args_info);

	return 0;
}