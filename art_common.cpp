#include <iostream>
#include <vector>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include "art_common.h"
#include "io_funct.h"

using namespace std;

/**
 * define not_found as maximum value of the unsigned long
 * */
#define not_found  ULONG_MAX 

/**
 * core ART algorithm and functions
 */


void ART::data_normalize(ART::a_lines& lines, ART::minmax& mm){
	//declaring iterator
	list<ART::iline>::iterator itline = lines.begin();
	
	//all lines has to have the same number of columns
	const unsigned long lsize = lines.front().val.size();
	
	//initialize minmax
	//it is a variable for storing max and min values for every column
	ART::mmax m;
	m.min=FLT_MAX;
	m.max=FLT_MIN;
	mm.assign(lsize,m);
	
        //line size 
	const unsigned long linesize = lines.size();
	//searching max and min through lines 
	for(itline=lines.begin();itline!=lines.end();itline++){
	   ART::iline line = *itline;
	   
	   for(unsigned int i=0;i<lsize;i++){
	   	if(line.val[i]>mm[i].max) mm[i].max = line.val[i];
	   	if(line.val[i]<mm[i].min) mm[i].min = line.val[i];
	   	}
	   }
	   

	   //data normalizing
	   for(itline=lines.begin();itline!=lines.end();itline++){
	   	ART::iline line = *itline;
	   	for(unsigned int i=0;i<lsize;i++){
	   		line.val[i] = (double)(line.val[i]-mm[i].min)/(mm[i].max-mm[i].min);
	   	}
   		
  	   
  	   	//insert line with normalized values
	   	lines.push_front(line);
	   }

    //deleting old lines
    for(unsigned long i=0;i<linesize;i++) lines.pop_back();
}


unsigned long ART::instance_in_sequence(ART::prototype_sequence& prot_seq, unsigned long iinst){
	unsigned long prot_number=not_found;
	for(unsigned long i=0;i<prot_seq.size();i++){
		//every prototype is compound of different number of instances 
		for(unsigned long j=0;j<prot_seq[i].size();j++){
			if(prot_seq[i][j]==iinst){
				prot_number=i;
				break;
				}
		}
		if(prot_number!=not_found) break;
	}
	return prot_number;
}





/**
 * find a particular instance(example Ek) or prototype in a sequence
 * and give as a result the index of instance(prototype) in the sequence
 * item means instance or prototype.
 * If must_find is set up as true and the example (item) has not been found -- write error 
 * message and stop the program.
 * It is used by ART 1, ART 2A and ART 2A-C algorithms
 **/
unsigned long ART::find_item(ART::a_instance& sample, ART::instance& inst, bool must_find){
	unsigned long index=not_found;
	
	for(unsigned long i=0;i<sample.size();i++){
	if(sample[i]==inst){
		index=i;
		break;
	}
	}
		
	if(!must_find)return index;
	else{
		if(index != not_found) return index;
		else{
			cerr << "\n\tfind_item(): Error: sample '" << inst << "' was not find in sequence.\n" << endl;
			     exit(22);
			     }
	}
}



// It converts the part of input lines which is going to be processed by the ART algorithm  into vectors.
// Vectors are simpler to work with
ART::a_instance ART::toVectors(ART::a_lines& lines){
	ART::a_instance res;
	
	list<ART::iline>::iterator itline;
	for(itline=lines.begin();itline!=lines.end();itline++){
		ART::iline line = *itline;
		res.push_back(line.val);
	}
	return res;
}


