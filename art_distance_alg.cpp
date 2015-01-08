#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "art_common.h"
#include "art_distance_alg.h"
#include <time.h>

//core algorithm and specific functions for art_distance

using namespace std;
using namespace ART;


/* *
 * define not_found as maximum value of the unsigned long
 * */
#define not_found  ULONG_MAX 


/**
 * The dimensionality of the task. It is Initialized in art_distance()
 * Since all the examples (Ek) have the same dimensionality so we can take it as sample[0].size()
 * It will be initialized in art_distance
 * */
unsigned long dim = -1;

/**
 * LU decomposition of a covariance matrix is stored in mat. Every cluster has it's own structure lu.
 * All of them are stored as a item in LU vector. 
 * lu structure contains covariance matrix mean for every variable (dimension), and gsl_permutation p
 * which is important for solving Ax=b equation.
 * Every prototype has it's own lu structure because the computation is much faster
 * This is for Mahalanobis distance only. It is not functional yet!
 * */
struct lu{
   /**
    * matrix after LU decomposition. Important for gsl_linalg_LU_solve
    * */
   gsl_matrix *mat;

   /**
    *arithmetic mean for every column of the matrix
    * */
   gsl_vector *mean;

   /**
    *important for gsl_linalg_LU_solve
    * */
   gsl_permutation *p;

   /**
    * number of examples in a cluster
    * */
   unsigned long size;
};
vector < lu > LU;


/**
 * Structure which stored the best prototype and it's index for some actual example (Ek)
 * */
struct bestP{
      /**
       * best prototype
       * */
      prototype prot;
	
      /**
       * index of the best prototype (prototype's ID). It is the number of some line in Clust::proto_seq
       * */
      unsigned long iprot;
	};

//input parameters - will be set up in art_distance function;
in_param param;

/**
 * Input: a prototype and instance (example). Output: their correlation distance 
 * Equation D(x,y) = 0.5 + ((x-x_mean)*(y-y_mean)')/2*(sqrt((x-x_mean)(x-x_mean)')*sqrt((y-y_mean)*(y-y_mean)'))
 * It is a bit different than original equation because we need it to be in range [0,1]. And result 1 is if the similarity is
 * the biggest. The parameter dimen is number of dimensions.
 * It is working for ART distance (with parameter -d 4)
 * @param dimen it is number of dimensions
 * @param inst example Ek
 * @param prot prototype 
 **/
double ART::correl(const prototype& prot, const instance& inst, unsigned long dimen){
	double *inarr = new double[dimen];
	double *prarr = new double[dimen];
	
	//inst
	copy(inst.begin(),inst.end(),inarr);
	double inst_mean = gsl_stats_mean(inarr,1,dimen);
	gsl_vector_view inst_view = gsl_vector_view_array(inarr,dimen);
	
	//prot
	copy(prot.begin(),prot.end(),prarr);
	double prot_mean = gsl_stats_mean(prarr,1,dimen);
	gsl_vector_view prot_view = gsl_vector_view_array(prarr,dimen);
	
	//count inst - inst_mean
	gsl_vector_add_constant(&inst_view.vector,-inst_mean);

	//count prot - prot_mean
	gsl_vector_add_constant(&prot_view.vector,-prot_mean);

	//count sqrt(inst_v * inst_v') * sqrt(prot_v * prot_v')
	double bottom = 2* (gsl_blas_dnrm2(&inst_view.vector) * gsl_blas_dnrm2(&prot_view.vector));
	
	//count inst_v * prot_v
	double upper;
	gsl_blas_ddot(&inst_view.vector,&prot_view.vector,&upper);

	//if bottom == 0 (it would be dividing by 0, change it into a small positive number)
	if(bottom==0)bottom = DBL_MIN;
	
	//correlation equation = 1- (upper/bottom). It is a bit changed to range become [0,1].
	//If the result is 1 it means the similarity is highest (x=y)
	double result = 0.5 +  upper/bottom;
	
	delete(inarr);
	delete(prarr);

	return result;
}


/**
 * Input: a prototype and instance (example). Output: their dot product (euclidean distance) 
 * Equation D(x,y) = 1 - sqrt( (sum(x[i]-y[i])^2)/dimen). 
 * It is working for ART distance (with parameter -d 1) 
 * @param dimen it is number of dimensions
 * @param prot prototype
 * @param inst example Ek
 **/
double ART::euclidean(const prototype& prot, const instance& inst, unsigned long dimen){
	double tmp=0;
	for(unsigned int i=0;i<dimen;i++) tmp += (inst[i] - prot[i]) * (inst[i] - prot[i]);
	return 1.0 - (sqrt(tmp/(double)dimen));
}



/**
 * Input: a prototype and instance (example). Output: dot product (modified euclidean distance) 
 * Equation D(x,y) = log(dimen) - sqrt(sum(x[i]-y[i])^2). 
 * It is working for ART distance (with parameter -d 2)
 * IT IS NOT WORKING PROPERLY YET!
 * @param dimen it is number of dimensions
 * @param prot prototype
 * @param inst example Ek
 **/
double ART::euclidean_m(const prototype& prot, const instance& inst, unsigned long dimen){
	double tmp=0;
	for(unsigned int i=0;i<dimen;i++) tmp += (inst[i] - prot[i]) * (inst[i] - prot[i]);
	return log10((double)dimen*(double)dimen) - sqrt(tmp);
}




/**
 * Input: a prototype and instance (example). Output: their mahalanobis distance .
 * It is working for ART distance (with parameter -d 6). The parameter dimen is number of dimensions.
 * STILL UNDER DEVELOPMENT - NOT WORKING YET
 * @param dimen it is number of dimensions
 * @param iprot is ID of the prot in prototype_sequence
 * @param prot prototype
 * @param inst example Ek 
 **/
inline double mahalanobis(const prototype& prot, const instance& inst, const unsigned long iprot, unsigned long dimen){

	if(LU[iprot].size<2) 
		return ART::euclidean(prot, inst, dimen);

	double result;
	double *exarr = new double[dimen];
	gsl_vector *x = gsl_vector_alloc(dimen);
	copy(inst.begin(),inst.end(),exarr);
	gsl_vector_view example = gsl_vector_view_array(exarr,dimen);
	
	//count C^{-1}(x-mean); where C^{-1} is inverse covariance matrix, result is in x
	gsl_linalg_LU_solve(LU[iprot].mat, LU[iprot].p,&example.vector,x);

	gsl_blas_ddot(&example.vector,x,&result);

	//return sqrt( (Ekj - mean)^t * C^{-1} * (Ekj - mean);
	// 1- and division by dimen is just addition (like euclidean distance) - the more similar the closer to 1
	//double res = 1 - sqrt(result/dimen);
	double res = result/dimen;

	gsl_vector_free(x);
	delete(exarr);

	return res;
}

/**
 * Input: a prototype and instance (example). Output: their manhattan (block) distance.
 * Equation: D(x,y) = 1 - (sum( fabs(x[i]-y[i]))/dimen. The parameter dimen is number of dimensions.
 * It is working for ART distance (with parameter -d 3) 
 * @param dimen it is number of dimensions
 * @param prot prototype
 * @param inst example Ek 
 **/
double ART::manhattan(const prototype& prot, const instance& inst, unsigned long dimen){
	double tmp = 0;
	for(unsigned int i=0;i<prot.size();i++) tmp += fabs(inst[i] - prot[i]);
	return 1 - (tmp / dimen);
}


/**
 * Input: a prototype and instance (example). Output: their minkowski distance.
 * Equation: D(x,y) = 1 - (sum( fabs(x[i]-y[i])^p)/dimen )^1/p. The parameter dimen is number of dimensions.
 * It is working for ART distance (with parameter -d 5). The parameter p is given by an user 
 * @param dimen it is number of dimensions
 * @param prot prototype
 * @param inst example Ek 
 * @param power power is given by an user or it is taken the default value (3). Minkowski with power 1 is 
 * Manhattan distance, with power 2 it is Euclidean distance
 **/
double ART::minkowski(const prototype& prot, const instance& inst, int power, unsigned long dimen){
	double res = 0;
	double tmp;
	//count res = sum(fabs(x-y)^p);
	for(unsigned int i=0;i<dimen;i++){ 
		tmp = fabs(inst[i] - prot[i]);
		res += gsl_pow_int(tmp,power);
		}
	//divide by number of dimensions
	res /= dimen;
	
	//count res_comp = res^{1/p}; unfortunatelly it is necessary to do it in complex numbers
	gsl_complex res_comp;
	GSL_SET_COMPLEX(&res_comp,res,0);
	double pow = 1.0 / power;
	res_comp = gsl_complex_pow_real(res_comp,pow);	
	
	//res= 1- (res)^{1/p}
	res = 1 - GSL_REAL(res_comp);
	if((GSL_IMAG(res_comp)!=0)&&(GSL_IMAG(res_comp)!=GSL_NAN)) {
		cerr << "ERROR: minkowski function: a complex number has arisen:"<<GSL_REAL(res_comp)<<"+"<<GSL_IMAG(res_comp)<< "i\n\t...Exiting" <<endl;
		exit(50);
		}	

	return res;
}


/**
 * Count score (similarity) how similar inst and prot are.
 * The output depends on distance measure which set an user up
 * @param inst it is some instance (example, Ek)
 * @param P it is the best prototype
 **/
inline double countScore(const bestP& P, const instance& inst){
	double result=0;
	//parameter dim is already set up in art_distance()
	switch(param.distance){
		case 1:
			result=euclidean(P.prot, inst, dim);
			break;
		case 2:
			result=euclidean_m(P.prot, inst, dim);
			break;
		case 3:
			result=manhattan(P.prot,inst, dim);
			break;
		case 4:
			result=correl(P.prot,inst,dim);
			break;
		case 5:
			//param.power is a parameter p;
			result=minkowski(P.prot,inst, param.power, dim);
			break;
		case 6:
			result=mahalanobis(P.prot,inst, P.iprot, dim);
			break;

	}
	
	return result;
}


/**
 * It counts LU decomposition of covariance matrix (CM). Only for Mahalanobis distance
 **/
void countCM(const a_instance& sample, const prototype_sequence& seq, unsigned long iprot){
	//number of examples in the cluster
	unsigned long size = seq[iprot].size();

	lu res;// = LU[iprot];
	res.mean = gsl_vector_alloc(dim);
	res.mat = gsl_matrix_alloc(dim,dim);
	res.p = gsl_permutation_alloc(dim);
	res.size = size;

	//make transpose matrix to samples from cluster
	gsl_matrix *trans = gsl_matrix_alloc(dim,size);
	double *arr= new double[dim];
	for(unsigned int i=0;i<size;i++){
		//arr ;
		copy(sample[seq[iprot][i]].begin(),sample[seq[iprot][i]].end(),arr);
		gsl_vector_view vec = gsl_vector_view_array(arr,dim);
		gsl_matrix_set_col(trans,i,&vec.vector);
	}

	//count mean for every dimension
	double mean;
	for(unsigned long i=0;i<dim;i++){
		mean=0;
		for(unsigned long j=0;j<size;j++){
			mean+=gsl_matrix_get(trans,i,j);
		}
		mean /= size;
		gsl_vector_set(res.mean,i,mean);
	}


	//count covariance matrix from trans
	double cov;
	gsl_vector *veci= gsl_vector_alloc(size);
	gsl_vector *vecj= gsl_vector_alloc(size);
	for(unsigned long i=0;i<dim;i++){
		for(unsigned long j=0;j<=i;j++){
			gsl_matrix_get_row(veci,trans,i);
			gsl_matrix_get_row(vecj,trans,j);

			cov = gsl_stats_covariance_m(veci->data,veci->stride,vecj->data,vecj->stride,size,gsl_vector_get(res.mean,i),gsl_vector_get(res.mean,j));

			//covar. matrix is symmetric
			gsl_matrix_set(res.mat,i,j,cov);
			//if it is not an diagonal element then set the symmetric element up
			if(i!=j)gsl_matrix_set(res.mat,j,i,cov);
		}
	}

	int signum;
	//make cholesky decomposition of the covariance matrix
	gsl_linalg_LU_decomp(res.mat,res.p,&signum);

	
	
	//give results to appropriate LU structure
	LU[iprot] = res;

	//freeing temporary variables
	gsl_matrix_free(trans);
	gsl_vector_free(veci);
	gsl_vector_free(vecj);
}


/**
 * Add an example (Ek) to a particular cluster.
 * It means that it moves the prototype toward the example.
 * The prototype will be more similar with the example.
 * P'= sum_i(1-beta*Pi + beta*Eki)
 * @param inst Ek
 * @param prot some prototype
 * @param beta it is given by an user
 * */
inline void add_instance(const instance& inst, prototype& prot, const float beta){
	//make vector  prot=(1-beta)*P + beta*E
	for(unsigned int i=0;i<inst.size();i++)prot[i]=((double)(1-beta)*prot[i]+beta*inst[i]);
	
}





/**
 * Removing an instance(Ek) from the particular prototype.
 * Remove an instance with index 'iinst' in 'sample' from prototype
 * with index 'iprot' in 'prot'. But also remove particular index 
 * from prototype sequence.
 * @param sample all input examples
 * @param iinst ID of the instance
 * @param prot set of prototypes
 * @param iprot ID of the prototype
 * @param seq the struct with all prototypes and example IDs
 * @param P best (the most similar) prototype for particular 'sample[iinst]'
 * @param beta parameter set by an user
 **/
void remove_instance(const a_instance& sample, const unsigned long iinst, a_prototype& prot,const unsigned long iprot, prototype_sequence& seq, bestP& P, const float beta){

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
		
		//if P is not empty, and iprot (removed prototype index) is lower than P.iprot (best prototype index) 
		//then decrease P.iprot number (because there is 1 prototype less)
		if(!P.prot.empty()) if(iprot< P.iprot){ 
			if(P.iprot <=0) cerr << "\nError occurred: remove_instance():P.iprot = 0;\nIt can not be decreased!" << endl;
			P.iprot--;
		}
		
		//if distance is mahalanobis delete particular CM
		if(param.distance==6){ 
			gsl_matrix_free(LU[iprot].mat);
			gsl_permutation_free(LU[iprot].p);
			delete(LU[iprot].mean);
			LU.erase(LU.begin()+iprot);
			}
	}
	//if it is not empty - re-create it from the rest instances
	else{
		//build prototype but without instance which should be deleted
		//at first -- prototype is the first item in the prototype sequence
		prot[iprot]=sample[seq[iprot][0]];
			
		//continually add others samples
		//it started from 2nd member because the first is already in
		for(unsigned int i=1;i<seq[iprot].size();i++){
		add_instance(sample[seq[iprot][i]],prot[iprot],beta);
		}
		
		//if the distance measure is mahalanobis, count LU covariance matrix decomposition
		if((param.distance==6)&&(seq[iprot].size()>1)) countCM(sample,seq,iprot);
	}
}



/**
 * Create a new prototype and also create a new sequence in prot_seq
 * One line in prot_seq = one cluster represented by a prototype
 * @param inst set of all examples
 * @param iinst ID of particular Ek
 * @param prot set of prototypes
 * @param prot_seq set of all prototypes with indexes of member's Ek
 * */
void create_cluster(const a_instance& inst, const unsigned long iinst, a_prototype& prot, prototype_sequence& prot_seq){
		
	// create a new prototype from the instance
	prot.push_back(inst[iinst]);
	
	// create a new prototype sequence and insert the first index of instance
	vector <unsigned long> new_seq;
	new_seq.push_back(iinst);

	prot_seq.push_back(new_seq);

	//if the measure is mahalanobis - create new CM
	if(param.distance==6){
		lu newlu;
		newlu.mat = gsl_matrix_alloc(dim,dim);
		newlu.mean = gsl_vector_alloc(dim);
		newlu.p = gsl_permutation_alloc(dim);
		newlu.size = 1;
		LU.push_back(newlu);
	}

}




/** 
 * Returns a prototype with highest similarity (score) -- which was not used yet.
 * The score is counted for a particular instance Ek and all the prototypes.
 * If it is returned empty prototype -- was not possible (for some reason) to find the best
 * @param inst example Ek
 * @param prot set of prototypes
 * @param used set of already tested prototypes
 **/
bestP best_prototype(const instance& inst,const a_prototype& prot,a_prototype& used){
	//prototypes with the same score and their indexes
	vector < bestP > same_score;

	const unsigned long usize=used.size();
	const unsigned long psize=prot.size();
	
	bestP empty;
	//temporary variable
	bestP P;

	//if the number of already used prototypes and the number of
	// prototypes are the same return empty protot. (no best protot.)
	if(used.size()==prot.size()) return empty;

	double *score = new double[psize];
	//setting initial value(the minimum for type double for this particular architecture) for scoring prototypes
	for(unsigned int i=0;i<psize;i++)score[i]=DBL_MIN;

	//set score for every prototype
        for(unsigned int i=0;i<psize;i++){
        
	   //search if a prototype is not among already used prototypes
	   bool usedb=false;
	       for(unsigned int j=0;j<usize;j++){
		  if(prot[i]==used[j]){
		     usedb=true;
		     break;
		  }
	       }
	    //is proto[i] among the used ??
	    if(usedb) continue;
	    //if it is not count it's score
	    else {
	       P.prot=prot[i];
	       P.iprot=i;
	       score[i] = countScore(P,inst);
	    }
      }

      //find prototype with highest score
      double higher=DBL_MIN;
      bestP same;
      for(unsigned int i=0;i<psize;i++){
	 if(score[i]==higher){
		same.prot = prot[i];
		same.iprot = i;
   		same_score.push_back(same);
   	}
   	else{
   		if(score[i]>higher){
   		
		     //erase the old list
		     if(same_score.size()>0) same_score.erase(same_score.begin(),same_score.end());
		     same.prot = prot[i];
		     same.iprot = i;
		     same_score.push_back(same);
		     higher=score[i];
   		}
   	 }
   }
   
   delete(score);
   
   // random choosing from the best prototypes
   unsigned long index = 0;
	
   //the result is an empty prototype
   if(same_score.size()==0){ 
	return empty;
   }
   
   //the result is the only one possible best prototype 
   if(same_score.size()==1){
	return same_score[0];
   }
   
   //if there is more best prototypes with the same score -- random choosing
   index = (rand() % same_score.size());

   return same_score[index];
}




/**
 * ART distance algorithm, inputs: examples and input parameters given by an user
 * How exactly it is working can be found at www.fi.muni.cz/~xhudik/art/drafts
 * @param sample  set if input examples (Eks)
 * @param par all input parameters set by an user or default
 **/
ART::Clust ART::art_distance(const ART::a_instance& sample, ART::in_param& par) {
    Clust results;

    //P will store the prototype with highest score and it's index
    bestP P;

    //list of all prototypes
    vector < prototype > prot;

    //the best representation of the prototypes of the whole history
    vector < prototype > prot_best;
    
    //sequences of samples from which prototype has been created
    //it is possible to reconstruct a prototype from the sequence
    prototype_sequence prot_seq;

    //the best representation of the prototype sequence of the whole history
    prototype_sequence prot_seq_best;

    //list of prototypes which were used already
    vector < prototype > used;
    
    //initializing of global variable 
    param = par;

    //initializing of dimensionality
    dim = sample[0].size();


   //initialization of random numbers
   srand((unsigned int)time(NULL)); 
	
   float fluctuation = 100;

   //the best error of the whole history
   //it is initialized as some impossible number, to avoid problems with first iteration
   float fluctuation_best=120;
   unsigned long pass = 0;

   //how many Ek's has been reassign to other cluster (prototype) in a previous pass (run)
   vector < bool > changed(sample.size(),true);

   //do cycle while error is higher the parameter -e or number of passes is lower than the 
   //the parameter -E
   while((pass<param.pass)&&(fluctuation>param.error)) {
        
        //nullifying 'changed'
        changed.assign(sample.size(),false);

	//cycle for examples Ek
        for(unsigned long i=0;i<sample.size();i++) {

            //nullifying 'used' prototypes 
            used.erase(used.begin(),used.end());
            
            do {
	       //find the best prototype for this Ek
	       P=best_prototype(sample[i],prot,used);

               //if there is no best prototype 
               if(P.prot.empty()){
                  
		  //check if the instance is not included already in some other prototype
                  const unsigned long prot_index=ART::instance_in_sequence(prot_seq,i);
                  if(prot_index!=not_found){ 
                    	
		     //if so, remove it (recreate prototype--without the instance)
                     remove_instance(sample,i,prot,prot_index,prot_seq,P,param.beta);
		  }
		  create_cluster(sample,i,prot, prot_seq);
		  changed[i]=true;
		  break;
               }
	
	       //add P among 'used' 
	       used.push_back(P.prot);
	       const double score = countScore(P,sample[i]);
	       double alphaSum=0;
	       for(unsigned long j=0;j<sample[i].size();j++) alphaSum += param.alpha*sample[i][j];

	       if(score>=alphaSum) {
	
		  //if similarity is sufficient -- sample[i] is member of the P
		  if(score >= param.vigilance) {
                  
		     //if the instance is already included in some other prototype -- find it
                     const unsigned long prot_index=ART::instance_in_sequence(prot_seq,i);
                     if(prot_index!=not_found) {
                  		
			//test if founded prototype is not actual one(P) in that case try another instance
                        if(prot_index==P.iprot) {
                           //example is already included into the actual prototype - do not change
                           break;
                           } else{
                                //re-build prototype - without the sample
                                remove_instance(sample,i,prot,prot_index,prot_seq,P,param.beta);
			   }
                     }
                            
                     //add instance to the current prototype
                     add_instance(sample[i],P.prot,param.beta);


                     prot_seq[P.iprot].push_back(i);
                     if(param.distance==6) countCM(sample,prot_seq,P.iprot);
                     changed[i]=true;
                     break;
                  }
                  
		  // try other best P.prot
                  else {
                     continue;
                  }
                    
	       } //score=>alphaSize 
                  
	       //if prototype is not enough similar to instance(sample[i]) then create a new prototype
               else {
		  
		  //check if the instance is not already in some other prototype
                  const unsigned long prot_index=ART::instance_in_sequence(prot_seq,i);
                  if(prot_index!=not_found){ 
                  
		     //if so, remove it (recreate prototype--without the instance)
                     remove_instance(sample,i,prot,prot_index,prot_seq,P,param.beta);
		  }
		  create_cluster(sample,i,prot,prot_seq);
		  changed[i] = true;
		  break;
	       }
            } while(prot.size()!=sample.size());

        } // for sample
		
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

    //giving back (maybe) changed input parameters during computing
    par = param;
    
    //create results
    results.proto = prot_best;
    results.proto_seq = prot_seq_best;
    results.fluctuation  = fluctuation_best;
    return results;
}

