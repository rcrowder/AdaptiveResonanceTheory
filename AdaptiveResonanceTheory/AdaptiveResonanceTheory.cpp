/*
 *  ART distance (Adaptive Resonance Theory) implementation
 *
 *  xhudik@gmail.com
 *
 */

#include <ctime>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <list>
#include <sstream>
#include <vector>
#include "..\get_options.h"
#include "..\io_funct.h"
#include "..\art_common.h"
#include "..\art_distance_alg.h"

using namespace std;
using namespace ART;

/* *
 * All input parameters for ART algorithm are stored in the structure
 * */
in_param par;

/**
 * Output file name -- the name can be changed by user
 * */
const char* o_file="art";


/**
 * Output file stream 
 * */
ofstream out;

/**
 * input file name and file -- the name can be changed by user
 * */
const char* i_file;

/**
 * input file stream 
 * */
ifstream in;


/**
 * how many input examples contain only 0's
 * They are not processed by the ART algorithm
 * */
unsigned long zero_all = 0 ,zero_real = 0;

/**
 * The structure of the input: array of array (matrix)
 * It is the set of input examples (Ek's)
 * */
a_instance sample;

/**
 * Input lines are stored in the a_lines
 * They consist of Ek + columns which are not going to processed
 * */
a_lines lines;

/**
 * Variable for storing min and max value for every column (feature)
 * */
ART::minmax mm;

/**
 * It reads user's input arguments and change them
 **/
void read_args(int argc, char* argv[]){
   int opt;
   bool	iflg=false,
	errflg=false;
	float tmpf;
	//bool tmpi;
	int tmpint;


   while ((opt = GetOpt::getopt(argc, argv, "i:o:b:v:a:e:s:E:d:p:")) != EOF) {
    switch (opt) {
    case 'b':
      tmpf=(float)atof(GetOpt::optarg);
      if((tmpf<=1)&&(tmpf>=0)) par.beta=tmpf;
      	else errflg=true;
      break;
      
    case 's':
    	tmpint=atoi(GetOpt::optarg);
      if(tmpint>=0){ par.skip=tmpint;
	cout<<"par.skip=========="<<par.skip<<" a optarg="<<GetOpt::optarg<<" tmpi="<<tmpint<<endl;}
      else errflg=true;
      break;
	
	case 'p':
    	tmpint=atoi(GetOpt::optarg);
      if(tmpint>=0) par.power=tmpint;
      else errflg=true;
      break;

    case 'v':
      tmpf=(float)atof(GetOpt::optarg);
      if(((tmpf==0)||(tmpf>0))&&(tmpf<=1)) par.vigilance=tmpf;
      else errflg=true;
      break;

    
    case 'a':
    	tmpf=(float)atof(GetOpt::optarg);
    	if((tmpf>0)||(tmpf==0)) par.alpha=tmpf;
    	else errflg=true;
      break;
    
    case 'o':
      o_file=GetOpt::optarg;
      break;

    case 'e':
    	par.error=(float)atof(GetOpt::optarg);
    	if((par.error<0)||(par.error>100)) errflg = true;
      break;
      
    
	 case 'd':
    	tmpint=atoi(GetOpt::optarg);
    	switch(tmpint){
    		case 1:
    			par.distance=1;
    			break;
    		case 2:
    			par.distance=2;
    			break;
    		case 3:
    			par.distance=3;
    			break;
    		case 4:
    			par.distance=4;
    			break;
    		case 5:
    			par.distance=5;
    			break;
    		case 6:
    			par.distance=6;
    			break;

    		default:
    			errflg=true;
    			break;
        	}
      break;

    case 'E':
    	par.pass=atol(GetOpt::optarg);
      break;


   //a must argument
    case 'i':
        if (GetOpt::optarg == NULL) errflg = true;
        else{ 
	   i_file = GetOpt::optarg;
	   iflg=true;
	}
      break;

    case '?': errflg = true; // Notice that '?' indicates a bad option
//    default:  errflg = true;
    } //switch
  } //while

   //necessary arguments have to be filled otherwise error
   if(!iflg)errflg=true;


  if (errflg) {
    cerr << "usage: art_distance [-o -s -b -v -a -e -E] -i input_file \n\n";
    cerr << "The options mean:\n";
    cerr << "\t-i (a must!) input_file\n";
    cerr << "\t-o prefix for output files\n";
    cerr << "\t-b beta (learning rate, [0,1]). Default 0.5\n";
    cerr << "\t-v vigilance ( 0<= vigilance <=1 )\n";
    cerr << "\t-a alpha ( alpha <= #columns^-0.5)\n";
    cerr << "\t-d distance measure: \n";
    cerr << "\t\t1 euclidean: 1 - sqrt( Sum((x-y)^2)/#dimen )\n";
    cerr << "\t\t2 modified euclidean: log(#dimen^2) - sqrt( Sum((x-y)^2) ) NOT WORKING YET\n";
    cerr << "\t\t3 manhattan: 1 - Sum(|x - y|)/#dimen\n";
    cerr << "\t\t4 correlation distance: 0.5 + (x-Xmean)*(y-Ymean)'/(2*sqrt((x-Xmean)(x-Xmean)' * sqrt((y-Ymean)(y-Ymean)')))\n";
    cerr << "\t\t5 minkowski distance: 1 - (Sum((|x - y|)^p)/#dimen)^(1/p); set up -p\n";
    cerr << "\t\t6 mahalanobis: 1 - sqrt( (x-mean) * C^-1 * (x-mean)') NOT WORKING YET\n";
    cerr << "\t-p x (x is positive integer); it is a power for the minkowski distance\n";
    cerr << "\t-s x (x is positive integer) - skip last x columns\n";
    cerr << "\t-e fluctuation: % of examples which are re-assigned\n";
    cerr << "\t-E integer: maximum number of passes through the input examples\n\n\n";
    cerr << "art distance is based on:\nhttp://www.fi.muni.cz/~xhudik/art\n\nFrank, Kraiss, Kuhlen,\"Comparative analysis of Fuzzy ART and ART-2A network clustering performance\"";
    cerr << " and  Recognition\", Neural Networks, vol. 9, pp. 544--559, 1998\n\n\n" << flush;
    exit(1);
  }
  cout << flush;

}



/**
 * It counts the original values of all examples (before normalizations).
 * And the examples are converted into output lines.
 * The type of results will be a vector (not a list)
 **/
vec_lines createVecLines(a_lines& lines,const ART::minmax mm){
	vec_lines res;
	list<iline>::iterator itline;

	itline=lines.begin();
	iline line = *itline;
	
	const unsigned long valsize = line.val.size();
	
	for(itline=lines.begin();itline!=lines.end();itline++){
		line = *itline;

		//count original(from input example) value of every feature in particular line
		for(unsigned long i=0;i<valsize;i++){ 
		line.val[i]= (double) line.val[i]*(mm[i].max - mm[i].min) + mm[i].min;
		}
		res.push_back(line);
	}
	return res;
}



/**
 * closing input and output files
 * */
void clean(){
   in.close();
   out.close();
}



/**
 *  The main function: main()
 *
 *  Initialize, perform the ART distance algorithm
 *  @param argc the number of input arguments
 *  @param argv the list of input arguments
 **/
int main(int argc, char* argv[])
{
   //initializing default values of input parameters
   par.beta = 0.5f;			// learning parameter [0,1] 
   par.vigilance = 0.1f;	// 0 <= vigilance <= 1 
   par.alpha = FLT_MAX;		//range: alpha <= #columns^-0.5
   par.distance=1;			// 1=Euclidean, 2=mahalanobis, 3=manhattan (block), 4=correlation, 5=minkowski
   par.power = 3;			//only for minkowski measure 
   par.error = 5;			//error rate - parameter "-e"
   par.pass = 100;			//maximum pass through - parameter "-E"
   par.skip = 0;			//number of skipped columns (-s)


   cout << "ART distance (Adaptive Resonance Theory) -- clustering algorithm.\n";
   cout << "Copyright Tomas Hudik, Jan Zizka\nContact: xhudik@fi.muni.cz\n\n" << endl;
   read_args(argc,argv);
  
  //print all option and input parameter
  cout << "Input file: " << i_file << endl;
  stringstream outname;
  outname << o_file << "_results";
  cout << "Output file: " << outname.str() << endl;
  if(par.skip>0) cout << "Number of skipped columns: " << par.skip << endl;
  	
  cout << "Max # of passes: " << par.pass << endl;
  cout << "Maximal allowed fluctuation: " << par.error << "%" << endl;

  cout << "distance: ";
  if(par.distance==1) cout << " euclidean" << endl;
  if(par.distance==2) cout << " modified euclidean" << endl;
  if(par.distance==3) cout << " manhattan (block)" << endl;
  if(par.distance==4) cout << " correlation" << endl;
  if(par.distance==5) cout << " minkowski ; power = " << par.power << endl;
  if(par.distance==6){ cout << " mahalanobis NOT WORKING. MODIFIED EUCLIDEAN DISTANCE USED INSTEAD" << endl;
    par.distance =2;
    }

  cout << "beta: " << par.beta << endl;
  cout << "vigilance: " << par.vigilance << "\n";
  if(par.skip>0) cout << "Skip the last " << par.skip << " columns from input examples\n\n";

  openFile(in,i_file,true);
  openFile(out,outname.str().c_str(),true);
  cout << "Reading input file: " << i_file << flush;
  readFile(in,lines,par.skip);
  cout << " ...OK (read " << lines.size() << " examples)" << endl;

   //check if a number of input instances is at least 1, otherwise it doesn't make sense
   if(lines.size()<2){
      cerr << "ERROR -- size of input data is " << lines.size() 
	   << "\nThere has to be at least 2 input instances!\n "
	   << "-- exiting" << endl;
      exit(1);
   }
   
   //removing samples which contains only 0's
   //return the number of 0's samples
   cout << "\nRemoving examples containing only 0's ..." << flush;
	
   zero_all = zero_free(lines);
   if(zero_all>0){ 
      cerr << "\nWARNING - in the input examples are " << zero_all ;
      cerr << " samples which contain only 0's. They were removed!\n";
   }
   else cout << "OK" << endl;
	
   list<iline>::iterator itline;
   itline = lines.begin();
   iline line = *itline;
   	
   float possiblemax = float(1.0/sqrt((double)line.val.size()));

   //if the alpha value is not set - use default and write the values
   cout << "alpha";
   if(par.alpha==FLT_MAX) {
      cout << " - default value ((#columns^-0.5)/2)";
      par.alpha = float(1.0/sqrt((double)line.val.size())) / 2.f;
   }
   
   cout << ": " << par.alpha << endl;
   if(par.alpha > possiblemax) {
      cerr << "ERROR alpha =" << par.alpha << ", but alpha for this task should be: alpha <= " << possiblemax;
      cerr <<")\n\n\tExiting..." << endl;
      exit(50);
   }

   cout << "\nInput samples normalizing ..." << flush;
   	
   //normalize data and create minmax values for all columns
   data_normalize(lines,mm);
   cout << "OK" << endl << flush;
  	
   //check if a number of samples is at least 1, otherwise ART doesn't make sense
   if(lines.size()<2){
      cerr << "ERROR -- size of prepared data is " << lines.size() 
	   << "\nThere has to be at least 2 instances!\n "
	   << "-- exiting" << endl;
      exit(1);
   }
	
   cout << "\nThe number of entering examples into ART: " << lines.size() << endl;
	
   //converting lines into the vectors
   sample = toVectors(lines);

   //controlling if the vigilance is set up correctly	
   //More precisely: if the score of an instance with the same instance is < vigilance
   //If it is so, we have to decrease vigilance otherwise it is infinite loop
   cout << "\nControlling correctness of vigilance against examples ..." << flush;
   double score=0;
   float tmpv;
   
   //in the following procedures will be needed variable dim (number of dimensions)
   unsigned long dim=sample[0].size();
   for(unsigned long i = 0;i<sample.size();i++){
      if(par.distance==1) score = euclidean(sample[i],sample[i],dim);
      if(par.distance==2) score = euclidean_m(sample[i],sample[i],dim);
      if(par.distance==3) score = manhattan(sample[i],sample[i],dim);
      if(par.distance==4) score = correl(sample[i],sample[i],dim);
      if(par.distance==5) score = minkowski(sample[i],sample[i],par.power,dim);
      //because menahttan point versus point is euclidean - we use euclidean measure instead of manhattan
      if(par.distance==6) score = euclidean(sample[i],sample[i],dim);

      if(score < par.vigilance){
	    if((score-par.vigilance)<0.0001) {
	    	tmpv= par.vigilance;
	        par.vigilance = par.vigilance - (0.0001f + 0.0001f*par.vigilance);
	    }else{ 
		tmpv = par.vigilance;
		par.vigilance = (float)score;
	    }
	    cerr << "\nWARNING: vigilance is too high (" << tmpv << "). What would mean infinite looping!!!\n";
	    cerr << "Vigilance was decreased: vigilance=" << par.vigilance << endl;
      }
   }
   cout << "OK" << endl << flush;


    
   cout << "\nComputing started\n" << endl;
   
   //counting
   //vigtmp stores the value given by user  -- it can be changed in art_distance()
   float vigtmp = par.vigilance;
   Clust clusters;
   clusters = art_distance(sample,par);

   cout << "\nComputing finished" << endl;
   cout << "Writing results ..." << flush;
	
   //write results
   out << "Adaptive Resonance Theory -- (ART distance)\n\n\n";
   out << "Input parameters:\n";
   if(par.vigilance!=vigtmp) out << "WARNING: the original vigilance value had to be decreased (to avoid infinite looping)\n";
   out << "\tvigilance: " << par.vigilance;
   out << "\n\tbeta: " << par.beta << "\n";
   out << "\talpha: " << par.alpha << "\n";
   out << "\tdistance: ";
   switch (par.distance) {
	case 1:
		out << "euclidean";
		break;
	case 2:
		out << "modified euclidean";
		break;
	case 3:
		out << "manhattan";
		break;
	case 4:
		out << "correlation";
		break;
	case 5:
		out << "minkowski ; power = " << par.power;
		break;
	case 6:
		out << "mahalanobis" << par.power;
		break;

   }
    
   if(zero_real!=0){ 
      out << "\nWARNING: After normalizing remained " << zero_real;
      out << " examples which are not included in the next results.";
   }
   out << "\n";
	
   out << "\tMax # of passes: " << par.pass << endl;
   out << "\tMaximal allowed fluctuation: " << par.error << "%" << endl;

   out << "\nThe number of examples entering into ART distance: " << sample.size() << endl;
   out << "\n\n\nResults:" << endl;
	
   //convert lines into vectors
   vec_lines res_lines;
   res_lines = createVecLines(lines,mm); 

   //write results
   //note: examples in sample are not changes in any way, therefore we can
   //substitute them for lines
   write_results(out,o_file,res_lines,clusters,zero_all+zero_real);
   cout << "OK\n" << endl;

   clean();
   return 0;
}

