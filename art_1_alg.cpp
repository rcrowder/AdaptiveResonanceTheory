
#include <iostream>
#include <vector>
#include <climits>
#include <cstdlib>
#include "art_common.h"


//core algorithm and specific functions for art_1


using namespace std;

//note: sample[i]=instance=example=Ek
//note: cluster is represented by a prototype (prototype is a centroid of the cluster)

/**
 * define not_found as maximum value of the unsigned long
 * */
#define not_found  ULONG_MAX 

/**
 * input parameters - will be initialize  in art_distance function;
 * */
in_param param;


/* *
 * Add a instance (Ek) to a particular cluster (it will make the prototype more similar)
 * The prototype will be more similar with an instance.
 * */
void add_instance_bin(instance& inst, prototype& prot){
	//'1' which are not in a inst[iinst] (but they are in prototype)
	//are nullifying in the prototype
	for(unsigned int i=0;i<inst.size();i++){
		if((inst[i]==0)&&(prot[i]==1))prot[i]=0;
	}
}



/* *
 * Removing an instance(Ek) from the particular prototype.
 * Remove an instance with index 'iinst' in 'sample' from prototype
 * with index 'iprot' in 'prot'. But also remove particular index 
 * from prototype sequence.
 * */
void remove_instance(a_instance& sample, unsigned long iinst, a_prototype& prot,unsigned long iprot, prototype_sequence& seq){
	//find and erase in the prototype sequence the instance which should be deleted
	for(unsigned int i=0;i<seq[iprot].size();i++) if(seq[iprot][i]==iinst){
		seq[iprot].erase(seq[iprot].begin()+i);
		break;
	}

	//if the particular prototype is empty now - delete whole prototype
	//delete also line (prototype) in prototype sequence
	if(seq[iprot].empty()){
		prot.erase(prot.begin()+iprot);
		seq.erase(seq.begin()+iprot);
	}

	//if it is not empty - re-create it from the rest instances
	else{
		//build prototype but without instance which should be deleted
		//at first -- prototype is the first item in the prototype sequence
		prot[iprot]=sample[seq[iprot][0]];
		//it started from 2nd member because the first is already in (1st sample = prototype
		//continually add others samples
		for(unsigned int i=1;i<seq[iprot].size();i++) add_instance_bin(sample[seq[iprot][i]],prot[iprot]);
	}
}



/* *
 * It create a cluster by creating of new prototype. And it also create a new sequence in prot_seq
 * One line in prot_seq = one cluster represented by a prototype
 * */
void create_cluster(a_instance& inst, unsigned long iinst, a_prototype& prot, prototype_sequence& prot_seq){
	//check if the instance is not presented in some other prototype
	//if the instance is already included in some prototype - recreate the prototype
	//and exclude this instance
	unsigned long prot_number=instance_in_sequence(prot_seq,iinst);
	if(prot_number!=not_found) remove_instance(inst,iinst,prot,prot_number,prot_seq);
	
	// create a new prototype
	prot.push_back(inst[iinst]);
	// create a new prototype sequence and insert the first index of instance
	vector <unsigned long> new_seq;
	new_seq.push_back(iinst);
	prot_seq.push_back(new_seq);
}




/** 
 * Return number of ones in an instance(Ek) or prototype.
 **/
inline long ones(instance& inst){
	unsigned long number=0;
	for(unsigned int i=0;i<inst.size();i++) if(inst[i]==1) number++;
	return number;
}

/** 
 * It returns the number of common 1's in both instances (prototypes).
 * Common ones means 1's at the same position in inst1 and inst2
 **/
inline long common_ones(instance& inst1,instance& inst2){
	unsigned long number=0;
	//size of both instances have to be the same
	for(unsigned int i=0;i<inst1.size();i++) if((inst1[i]==1)&&(inst2[i]==1)) number++;
	return number;
}


/** 
 * It returns the prototype with highest score (which was not used yet)
 * The score is counted for a current instance (Ek) 
 * If it is returned empty prototype -- was not possible (for any reason) to find it.
 * @param inst Ek
 * @param prot set of prototypes
 * @param used set of already tested prototypes
 * @param beta parameter set by an user
 **/
prototype best_prototype_bin(instance& inst,a_prototype& prot,a_prototype& used,float beta){
	//prototypes with the same score
	vector < prototype > same_score;
	const unsigned long usize=used.size();
	const unsigned long psize=prot.size();
	
	prototype empty;
	
	//if the number of already used prototypes and the number of
	// prototypes are the same return empty prototype. (no best prototype.)
	if(used.size()==prot.size()) return empty;

	double *score = new double[psize];
	//setting initial value for scoring prototypes
	for(unsigned int i=0;i<psize;i++)score[i]=-1;

	//set score for every prototype
        for(unsigned int i=0;i<psize;i++){
	    //search if prototype is not among already used prototypes
	    bool usedb=false;
            for(unsigned int j=0;j<usize;j++){
	       	if(prot[i]==used[j]){
		   usedb=true;
		   break;
	        }
	    }
	    //is proto[i] among the used ??
	    if(usedb) continue;
	    //if not count it's score
	    else score[i]= double(common_ones(prot[i],inst))/(beta+ones(prot[i]));

      }
	
     //find prototype with highest score
     double higher=-1;
     for(unsigned int i=0;i<psize;i++){
   	if(score[i]==higher){
   	same_score.push_back(prot[i]);
     }
     else{
   	 if(score[i]>higher){
	    //erase the old list
	    if(same_score.size()>0) same_score.erase(same_score.begin(),same_score.end());
	    same_score.push_back(prot[i]);
	    higher=score[i];
   	 }
      }
   }

   delete(score);

   //the result is an empty prototype
   if(same_score.size()==0) return empty;
   
   //the result is the only one possible best prototype 
   if(same_score.size()==1) return same_score[0];
   
   //if there is more best prototypes with the same score -- random choosing
   unsigned long index = 0;
   index = (rand() % same_score.size());
   return same_score[index];
}





/**
 * It provides the binary ART itself.
 * How exactly it is working can be found at www.fi.muni.cz/~xhudik/art
 * @param sample -- set of input examples (Ek's)
 * @param par -- all input parameters needed by ART 1
 **/
Clust art_1(a_instance& sample, in_param& par) {

    Clust results;

    //initializing of global variable 
    param = par;

    //prototype with highest score
    prototype P;

    //list of all prototypes
    vector < prototype > prot;

    //the best set of prototypes in all passes
    vector < prototype > prot_best;

    //sequences of samples Ek from which prototype has been created
    //it is possible to reconstruct a prototype from the sequence
    //defined in art_common.h
    prototype_sequence prot_seq;

    //the best sequence of prototypes of the whole history
    prototype_sequence prot_seq_best;

    //list of prototypes which were used already 
    vector < prototype > used;
    
    //vector of all samples which were changed in a previous run
    //(run is crossing through all examples)
    //From this variable it is counted a fluctuation
    vector < bool > changed(sample.size(),true);
	
   //this is used for counting how many instances were transformed
   //into another prototype (in %)
   float fluctuation=100;

   //the lowest error of the whole history
   //it is initialized as some impossible number (higher than 100% can't be), to avoid problems with first iteration
   float fluctuation_best=120;
	
   //how many times it run throughout the samples
   unsigned long pass = 0;
	
   //initialization of random numbers
   srand((unsigned int)time(0)); 

   //do cycle while error is higher than parameter -e or number of passes parameter -E
    while((pass<param.pass)&&(fluctuation>param.error)) {
        //nullifying changed values
        changed.assign(changed.size(),false);
		  
        //cycle for instances
	//sample[i] = Ek; sample - all examples
        for(unsigned int i=0;i<sample.size();i++) {

            //zeroing 'used' prototypes 
            used.erase(used.begin(),used.end());

            //declaring similarity 
            double similarity;

            do {
	        //find the best prototype for current Ek
                P=best_prototype_bin(sample[i],prot,used,param.beta);

                //if there is no best prototype (list of prototypes is empty,
                //or number of used and number of prototypes are the same 
                if(P.empty()){
		   
		   //check if the instance is not included already in some other prototype
                   const unsigned long prot_index=instance_in_sequence(prot_seq,i);
                   if(prot_index!=not_found){ 
 
		     //if so, remove it (recreate prototype--without the instance)
                     remove_instance(sample,i,prot,prot_index,prot_seq);
		  }
		     
		   create_cluster(sample,i,prot, prot_seq);
		     changed[i]=true;
		     break;
                }

		used.push_back(P);
		
		//measure1 = (how many 1 have P and sample[i](Ek) in common) / (beta+how many 1 has P)
		//to have 1's in common = have 1's on the same positions
		const double measure1 = double(common_ones(P,sample[i])) / (par.beta+ones(P));

		//measure2 = (how many 1 has Ek) / (beta + number of dimensions )
		const double measure2 = double(ones(sample[i])) / (par.beta+sample[i].size());

		//if P and sample[i] are similar enough
		if(measure1>=measure2) {
			
		   //similarity = (how many 1 have Ek and P in common) / (how many 1's has Ek)		
		   //to have 1's in common = have 1's on the same positions
		   similarity=double(common_ones(P,sample[i])) / ones(sample[i]);

		   //if the similarity is sufficient  -- than sample[i] going to be a member of P
      		   if(similarity >= par.vigilance) {

			   //if the instance is already included in some  prototype -- find it
			   const unsigned long prot_index=instance_in_sequence(prot_seq,i);
			   if(prot_index!=not_found) {
                  		//test if founded prototype is not actual one in that case go for another Ek
                                if(prot[prot_index]==P) {
			          	break;
                                } else{
				    //re-build prototype - without the sample (Ek)
				    remove_instance(sample,i,prot,prot_index,prot_seq);
				    }
                            }
                            
                            //index of P(current prototype) in prototypes
                            const unsigned long Pindex = find_item(prot,P,true);

                            //add instance to the current prototype
                            add_instance_bin(sample[i],prot[Pindex]);
                            prot_seq[Pindex].push_back(i);
                            changed[i]=true;
                            break;
                        }
                        // try other best P
                        else {
                           continue;
                        }
                    
                    } //measure1>=measure2
                    
                    //if prototype is not enough similar to instance(sample[i]) then create a new prototype
                    else {

		        //check if the instance is not already in some other prototype
			const unsigned long prot_index=instance_in_sequence(prot_seq,i);
                    	if(prot_index!=not_found){ 
               		   //if so, remove it (recreate prototype--without the instance)
                           remove_instance(sample,i,prot,prot_index,prot_seq);
			}
			create_cluster(sample,i,prot, prot_seq);
			changed[i]=true;
			break;
                  }

            } while(prot.size()!=sample.size());

        } // for sample
      
       //count statistics for this pass
       unsigned long number_changed=0;
       for(unsigned long c=0;c<changed.size();c++) if(changed[c]) number_changed++;
       fluctuation = ((float)number_changed/sample.size())*100;
       pass++;
       cout.precision(3);
       cout << "Pass: " << pass <<", fluctuation: " << fluctuation << "%" << ", clusters: " << prot.size() << endl;
		
       //test if this iteration has not lower error
       if(fluctuation < fluctuation_best){
	    //if it is so - assign the new best results
	    prot_best = prot;
	    prot_seq_best = prot_seq;
	    fluctuation_best = fluctuation;
	}

    } // while
    //create results
    results.proto = prot_best;
    results.proto_seq = prot_seq_best;
    results.fluctuation  = fluctuation_best;

    return results;
}

