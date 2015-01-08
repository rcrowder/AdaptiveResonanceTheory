/*
 *
 *  ART 1 (Adaptive Resonance Theory) implementation
 *  (binary input only)
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
#include "get_options.h"
#include "io_funct.h"
#include "art_common.h"
#include "art_1_alg.h"


using namespace std;
		
/**
 * Structure for storing all the input parameters needed by the ART algorithm
 * */
in_param par;

/**
 * Output file name -- the name can be changed by an user
 * */
const char* o_file="art";

/**
 * Output file stream
 * */
ofstream out;

/**
 * input file name and file -- the name can be changed by an user
 * */
const char* i_file;

/**
 * Input file stream
 * */
ifstream in;


/**
 * how many input examples contain only 0's
 * They are not processed by the ART algorithm
 * */
unsigned long zero_all = 0 ,zero_real = 0;

/**
 * the structure of input : array of array (matrix)
 * It is the set of input examples (Ek's)
 * */
a_instance sample;

/**
 * Input lines are stored in the a_lines
 * They consist of Ek + columns which are not going to processed
 * */
a_lines lines;

/**
 * read user's input arguments and change some default values
 **/
void read_args(int argc, char** argv){
   int opt;
   bool	iflg=false,
	errflg=false;
	int tmpi;

   while ((opt = GetOpt::getopt(argc, argv, "i:o:b:v:e:E:s:")) != EOF) {
    switch (opt) {
    case 'b':
      tmpi=atoi(GetOpt::optarg);
      par.beta=tmpi;
      break;
      
    case 's':
    	tmpi=atoi(GetOpt::optarg);
      if(tmpi>=0) par.skip=tmpi;
      else errflg=true;
      break;

    case 'v':
      float tmp;
      tmp=atof(GetOpt::optarg);
      if(((tmp==0)||(tmp>0))&&(tmp<=1)) par.vigilance=tmp;
      else errflg=true;
      break;

    case 'o':
      o_file=GetOpt::optarg;
      break;

    case 'e':
    	par.error=atof(GetOpt::optarg);
    	if((par.error<0)||(par.error>100)) errflg = true;
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
    cerr << "usage: art [-o -s -b -v -e -E] -i input_file \n\n";
    cerr << "The options mean:\n";
    cerr << "\t-i (a must!) input_file\n";
    cerr << "\t-o prefix for output files\n";
    cerr << "\t-b beta (small integer). Default 1\n";
    cerr << "\t-v vigilance ( 0<= vigilance <=1 )\n";
    cerr << "\t-s x (x is positive integer) - skip last x columns\n";
    cerr << "\t-e fluctuation: % of examples which are re-assigned\n";
    cerr << "\t-E integer: maximum number of passes through the input examples \n\n\n";
    cerr << "art 1 is based on:\nhttp://www.fi.muni.cz/~xhudik/art\n\nCarpenter and Grossberg, The ART of Adaptive Pattern Recognition by a Self-Organizing Neural Network, Computer, 21, 3, 1988\n\n" << flush;
    exit(1);
  }
  cout << flush;

}



/* 
 * Converting input lines into output lines (vector type)
 * The type of results will be a vector (not a list)
 */
vec_lines createVecLines(a_lines& lines){
	vec_lines res;
	list<iline>::iterator itline;
	for(itline=lines.begin();itline!=lines.end();itline++){
		res.push_back(*itline);
	}
	return res;
}




/**
 * Closing input and output files
 * */
void clean(){
   in.close();
   out.close();
}



/**
 *  The main function: main()
 *
 *  Initialize, perform the ART 1 algorithm
 *  @param argc the number of input arguments
 *  @param argv the list of input arguments
 **/
int main(int argc, char** argv)
{	
   //initializing default values of input parameters
   par.beta = 1;		// a small integer  !!! be aware this is different than other ART where it is a learning constant
   par.vigilance = 0.1;	// 0 <= vigilance <= 1 
   par.error = 5; //error rate - parameter "-e"
   par.pass = 100; //maximum pass through - parameter "-E"


   cout << "ART 1 (Adaptive Resonance Theory - for binary input) -- clustering algorithms.\n";
   cout << "Copyright Tomas Hudik, Jan Zizka\nContact: xhudik@fi.muni.cz\n\n" << endl;
   read_args(argc,argv);
  
  //print all options and input parameters
  cout << "Input file: " << i_file << endl;
  stringstream outname;
  outname << o_file << "_results";
  cout << "Output file: " << outname.str() << endl;
  if(par.skip>0) cout << "Number of skipped columns: " << par.skip << endl;
  	
  //error rate (or number of pass)
  if(par.error!=5){ cout << "Continuously fluctuating percentage of clustered instances: " << par.error << "%" << endl;
  }
  if(par.pass!=LONG_MAX){ cout << "Max # of passing through instances: " << par.pass << endl;
  	}

  cout << "beta: " << par.beta << endl;
  cout << "vigilance: " << par.vigilance << "\n";
  if(par.skip>0) cout << "Skip the last " << par.skip << " columns from input examples\n\n";

  openFile(in,i_file,true);
  openFile(out,outname.str().c_str(),true);
  cout << "Reading input file: " << i_file << flush;
  readFile(in,lines,par.skip);
  check_binary(lines);
  cout << " ...OK (read " << lines.size() << " examples)" << endl;


   
   //check if a number of input samples is at least 1, otherwise it doesn't make sense
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
   if(zero_all>0){ cerr << "\nWARNING - in the input examples are " << zero_all ;
	cerr << " samples which contain only 0's. They were removed!\n";}
   else cout << "OK" << endl;
	
   
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
    
   cout << "\nComputing started\n" << endl;
   
   
   //vigtmp stores the vigilance value given by user  -- it can be changed in art_1()
   float vigtmp = par.vigilance;
   Clust clusters;   
   
   //the main function of the algorithm
   clusters = art_1(sample,par);
	
   cout << "\nComputing finished" << endl;
   cout << "Writing results ..." << flush;
	
   //write results
   out << "Adaptive Resonance Theory (binary form)\n";
   out << "Input parameters:\n";
   if(par.vigilance!=vigtmp) out << "WARNING: the original vigilance value had to be decreased (to avoid infinite looping)\n";
   out << "\tvigilance: " << par.vigilance;
   out << "\n\tbeta: " << par.beta << "\n";
	
   out << "Continuously fluctuating percentage of clustered instances: " << par.error << "%" << endl;
   if(par.pass!=(LONG_MAX-1)) out << "Max # of passing through instances: " << par.pass << endl;
   out << "\nThe number of examples entering into ART: " << sample.size() << endl;
   out << "\n\n\nResults:" << endl;
	
   //convert lines into vectors
   vec_lines res_lines;
   res_lines = createVecLines(lines); 

   //write results
   //note: examples in sample are not changes in any way, therefore we can
   //substitute them for lines
   write_results(out,o_file,res_lines,clusters,zero_all+zero_real);
   cout << "OK\n" << endl;

   clean();
   return 0;
}



