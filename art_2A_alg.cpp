#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <climits>
#include "art_common.h"

//core algorithm and specific functions for art_2A


using namespace std;

/**
 * Define not_found as maximum value of the unsigned long
 * */
#define not_found  ULONG_MAX 

/**
 * Input parameters - will be initialize in art_2A function;
 * */
in_param param;

/**
 * Vector normalizing.
 * The results will be stored in the parameter lines (the parameter will be changed)
 * @param norm_index tells if it is first normalization (0), or the second (1).
 * Which line.norm will be filled up (0 or 1)
 * @param lines input examples
 **/
void vec_normalize(a_lines& lines, const unsigned int norm_index){
	//declaring iterator
	list<iline>::iterator itline = lines.begin();
	
	//all lines has to have the same number of columns (dimensions)
	//So, we take the first line
	const unsigned long lsize = lines.front().val.size();

        //line size 
	const unsigned long linesize = lines.size();

	//passing through lines 
	for(itline=lines.begin();itline!=lines.end();itline++){
	   iline line = *itline;
	   
	   double sum = 0;
	   for(unsigned long i=0;i<lsize;i++) sum += line.val[i] * line.val[i];
	   double norm = sqrt(sum);
	   for(unsigned long i=0;i<lsize;i++) line.val[i] = float (line.val[i]) / norm;
	   
	   //storing a normalization parameter
	   line.norm[norm_index]=norm;
	   //insert line with normalized values
	   lines.push_front(line);
	}

       //deleting old lines
       for(unsigned long i=0;i<linesize;i++) lines.pop_back();

}


/**
 * Nullifying values (in examples) which are < theta
 * The results will be stored in the parameter: line
 * @param theta  it is given by an user or default value
 * @param lines input examples
 **/
unsigned long denoising(a_lines& lines, const float theta){
   
   //take dimension
   const unsigned long lsize = lines.front().val.size();
	
   //how many examples were changed because of theta
   unsigned long changed = 0;
   bool changedb;
	
   //line size 
   const unsigned long linesize = lines.size();
   //declaring iterator
   list<iline>::iterator itline = lines.begin();

   for(itline=lines.begin();itline!=lines.end();itline++){
	iline line = *itline;
	changedb = false;
	for(unsigned int i=0;i<lsize;i++) 
	   if((line.val[i]>0)&&(abs(line.val[i])<=theta)){
	       line.val[i]=0;
	       changedb = true;
	    }
	    if(changedb) changed++;
	    lines.push_front(line);
	}
	
	//deleting old lines
	for(unsigned long i=0;i<linesize;i++) lines.pop_back();
	return changed;
}


/**
 * Count score (similarity) how similar inst and prot are.
 * The output (similarity) will be their dot product 
 * @param inst it is some instance (example, Ek)
 * @param prot it is some prototype
 **/
inline float countScore(const prototype& prot, const instance& inst){
	float score=0;
	for(unsigned int i=0;i<prot.size();i++) score += inst[i]*prot[i];
	return score;
}


/**
 * Add an example (Ek) to a particular cluster.
 * It means that it moves the prototype toward the example.
 * The prototype will be more similar with the example.
 * P'=sqrt( sum_i((1-beta)*Pi + beta*Eki)^2 )
 * @param inst Ek
 * @param prot some prototype
 * @param beta it is given by an user
 * */
void add_instance(instance& inst, prototype& prot, const float beta){
	prototype tmp;
	double norm = 0;
	//make vector  tmp=(1-beta)*P + beta*Ek
	for(unsigned int i=0;i<inst.size();i++)tmp.push_back((double)(1-beta)*prot[i]+beta*inst[i]);
	//count vector norm semi = sqrt(tmp^2)
	for(unsigned int i=0;i<inst.size();i++)norm += tmp[i]*tmp[i];
	norm = sqrt(norm);
	//count prototype
	prototype pr_tmp;
	for(unsigned int i=0;i<inst.size();i++)pr_tmp.push_back((double)(1/norm)*tmp[i]);
	prot=pr_tmp;
}


/**
 * Removing an instance(Ek) from the particular prototype.
 * Remove the instance with index 'iinst' in 'sample' from prototype
 * with index 'iprot' in 'prot'. But also remove particular index 
 * from prototype sequence
 **/
void remove_instance(a_instance& sample, unsigned long iinst, a_prototype& prot,unsigned long iprot, prototype_sequence& seq, const float beta , float& vigilance){
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
	//if it is not empty - re-create it from the rest examples
	else{
		//build prototype but without instance which should be deleted
		//at first -- prototype is the first item in the prototype sequence
		prot[iprot]=sample[seq[iprot][0]];
			
		//if PE < vigilance -- it won't stop (infinite looping)
		double score=countScore(sample[seq[iprot][0]],sample[seq[iprot][0]]);
		if(score < vigilance){
			float tmpv= vigilance;
			vigilance = score;
			cerr << "\nWARNING: vigilance is too high (" << tmpv << "). What means infinite looping!!!\n";
			cerr << "Vigilance was decreased: vigilance=" << vigilance << endl;
		}

		//continually add others examples
		//it started from 2nd member because the first is already in
		for(unsigned int i=1;i<seq[iprot].size();i++){
		add_instance(sample[seq[iprot][i]],prot[iprot],beta);
		}
		
	}
}




/**
 * Create a new prototype and also create a new sequence in prot_seq
 * One line in prot_seq = one cluster represented by a prototype
 * @param inst set of all examples
 * @param iinst ID of particular Ek
 * @param prot set of prototypes
 * @param prot_seq set of all prototypes with indexes of member's Ek
 * @param vigilance it is set by an user
 * */
void create_prototype(a_instance& inst, unsigned long iinst, a_prototype& prot, prototype_sequence& prot_seq, float& vigilance){
	
	//if PE < vigilance -- it won't stop
	double score=countScore(inst[iinst],inst[iinst]);
	if(score < vigilance){
	float tmpv;
	if((score-vigilance)<0.0001) {
		tmpv= vigilance;
		vigilance= vigilance - (0.0001 + 0.0001*vigilance);
	}else{ 
		tmpv= vigilance;
		vigilance = score;
	}

		cerr << "\nWARNING: vigilance is too high (" << tmpv << "). What means infinite looping!!!\n";
		cerr << "Vigilance was decreased: vigilance=" << vigilance << endl;
	}
	
	// create a new prototype
	prot.push_back(inst[iinst]);
	// create a new prototype sequence and insert the first index of instance
	vector <unsigned long> new_seq;
	new_seq.push_back(iinst);
	prot_seq.push_back(new_seq);
}



/** 
 * Returns a prototype with highest similarity (score) -- which was not used yet.
 * The score is counted for a particular instance Ek and all the prototypes.
 * If it is returned empty prototype -- was not possible (for some reason) to find the best
 * @param inst example Ek
 * @param prot set of prototypes
 * @param used set of already tested prototypes
 **/
prototype best_prototype_2A(instance& inst,a_prototype& prot,a_prototype& used){

   //prototypes with the same score
   vector < prototype > same_score;
   const unsigned long usize=used.size();
   const unsigned long psize=prot.size();
	
   prototype empty;

   //if the number of already used prototypes and the number of
   // prototypes are the same return empty protot. (no best protot.)
   if(used.size()==prot.size()) return empty;

   double *score = new double[psize];
   //setting initial value(the minimum for type double for this particular architecture) for scoring prototypes
   for(unsigned int i=0;i<psize;i++)score[i]=DBL_MIN;

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
        else score[i]=countScore(prot[i],inst);
   }

   //find prototype with highest score
   double higher=DBL_MIN;
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
 * ART 2A algorithm, inputs: examples and input parameters given by an user
 * How exactly it is working can be found at www.fi.muni.cz/~xhudik/art/drafts
 * @param sample  set if input examples (Eks)
 * @param par all input parameters set by an user or default
 **/
Clust art_2A(a_instance& sample, in_param& par) {

    Clust results;

    //initializing of global variable 
    param = par;

    //prototype with highest score
    prototype P;

    //list of all prototypes
    vector < prototype > prot;

    //the best representation of the prototypes of the whole history
    vector < prototype > prot_best;

    //sequences of samples Ek from which prototype has been created
    //it is possible to reconstruct a prototype from the sequence
    //defined in art_common.h
    prototype_sequence prot_seq;

    //the best representation of the prototypes of the whole history
    prototype_sequence prot_seq_best;

    //list of prototypes which were used already
    vector < prototype > used;
    
   //initialization of random numbers
   srand((unsigned int)time(0)); 
	
   float fluctuation = 100;

   //the lowest error of the whole history
   //it is initialized as some impossible number(higher than 100% can't be), to avoid problems with first iteration
   float fluctuation_best=120;  

   //how many times it run throughout the samples
   unsigned long pass = 0;

   //how many Ek's has been reassign to other cluster (prototype) in a previous pass (run)
   vector < bool > changed(sample.size(),true);


   //do cycle while error is higher than the parameter -e  or a number of passes is lower than the parameter -E
    while((pass<param.pass)&&(fluctuation>param.error)) {
        //nullifying changed values
        changed.assign(sample.size(),false);

	//cycle for instances
        for(unsigned long i=0;i<sample.size();i++) {

            //zeroing 'used' prototypes 
            used.erase(used.begin(),used.end());
            
            do {
	        //find the best prototype for current Ek
                P=best_prototype_2A(sample[i],prot,used);

                //if there is no best prototype 
                if(P.empty()){
		     //check if the instance is not included already in some other prototype
		     const unsigned long prot_index=instance_in_sequence(prot_seq,i);
		     if(prot_index!=not_found){ 
			   //if so, remove it (recreate prototype--without the instance)
			   remove_instance(sample,i,prot,prot_index,prot_seq,param.beta,param.vigilance);
		     }
                	 create_prototype(sample,i,prot, prot_seq,param.vigilance);
                	 changed[i]=true;
                	 break;
                }

	        //add P among 'used' 
		used.push_back(P);
		
		//count similarity between P and Ek (it is called "score") and alpha*sum_i Eki
		const double score = countScore(P,sample[i]);
		double alphaSum=0;
		for(unsigned long j=0;j<sample[i].size();j++) alphaSum += param.alpha*sample[i][j];

	        //if similarity is sufficient -- sample[i] is member of the P
		if(score>=alphaSum){
		     if(score >= param.vigilance) {
		     
			//if the example Ek is already included in some prototype -- find it
                  	const unsigned long prot_index=instance_in_sequence(prot_seq,i);
                  	if(prot_index!=not_found){

			        //test if the found prototype is not actual one (P) in that case try - go for another Ek
                                if(prot[prot_index]==P) {
				    break;
                                } else{
				       //re-build prototype - without the sample
                                       remove_instance(sample,i,prot,prot_index,prot_seq,param.beta,param.vigilance);
				   }
                        }
                            
                        // find an index of P in prototypes
                        const unsigned long Pindex = find_item(prot,P,true);

                        //add instance to the current prototype
                        add_instance(sample[i],prot[Pindex],param.beta);
                        prot_seq[Pindex].push_back(i);
                        changed[i]=true;
                        break;
		     }

		     // try other best P
                     else {
                        continue;
                     }
		} //score=>alphaSize
                    
                //if prototype is not enough similar to the example(sample[i]) then create a new prototype
                else {
		   //check if the instance is not already in some other prototype
                   const unsigned long prot_index=instance_in_sequence(prot_seq,i);
                   if(prot_index!=not_found){ 
		      //if so, remove it (recreate prototype--without the instance)
                      remove_instance(sample,i,prot,prot_index,prot_seq,param.beta,param.vigilance);
		   }
		   create_prototype(sample,i,prot, prot_seq,param.vigilance);
		   changed[i] = true;
		   break;
		}
	    } while(prot.size()!=sample.size());

        } // for sample

	//count statistics for this pass
	unsigned long number_changed=0;
	for(unsigned long j=0;j<changed.size();j++) if(changed[j]) number_changed++;
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

