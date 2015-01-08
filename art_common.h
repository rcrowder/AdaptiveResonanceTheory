#ifndef _ART_COMMON_
#define _ART_COMMON_

#include <vector>
#include <list>
#include <string>


/**
 * @authors &copy;Tomas Hudik, Jan Zizka<br>Contact: <A HREF="mailto:xhudik at gmail.com">xhudik\htmlonly &#64; \endhtmlonly gmail.com</A><br>
 * The software package is released under GNU GLP license
 * */


/**
 * @mainpage Adaptive Resonance Theory
 *
 * <h3> Basic description</h3>
 *
 * <em>Adaptive Resonance Theory</em> was developed by psychologists few decades ago. 
 * Later on, it was adapted to machine learning. Nowadays it is a whole framework
 * consisted of a bunch of algorithms ranging from supervised to unsupervised learning.
 *
 *
 *
 * This manual documents the software package that includes the ART algorithms for
 * unsupervised learning only. It is a family of four programs based on 
 * different ART algorithms. All of them are clustering algorithms and they are command-line programs. 
 * All of them are able to process only <b>numerical continuous values</b> (with the
 * exception of ART 1) and cannot handle <b>missing values</b>. Another important feature is 
 * that ART clustering algorithms do not need to know a number of how many clusters should be created 
 * the same is true for, e.g.: EM and XMeans). A list of implemented algorithms is in the following lines:
 * <br/>
 *
 * <em>art_1</em> -- it is the original binary ART algorithm. Input examples (in literature it is called <em>Ek</em>) 
 * can only be vectors containing 0's or 1's. The program consists of:
 * <em>art_1.cpp, art_1_alg.cpp, art_1_alg.h, get_options.cpp, get_options.h, art_common.cpp, art_common.h, 
 * io_funct.cpp, io_funct.h</em>  
 * <br/>
 *
 *
 * <em>art_2A</em> -- it is the basic ART algorithm for the real numbers input. It consists
 * of: <em>art_2A.cpp, art_2A_alg.cpp, art_2A_alg.h, get_options.cpp, get_options.h, art_common.cpp, art_common.h,
 * io_funct.cpp, io_funct.h</em>
 * <br/>
 *
 * <em>art_2A-C</em> -- it is improved art_2A. It is also a clustering algorithm for real numbers
 * using a complement coding which can obtain better results but it takes a
 * longer time. It consists of: <em>art_2A-C.cpp, art_2A-C_alg.cpp, art_2A-C_alg.h, get_options.cpp, get_options.h, 
 * art_common.cpp, art_common.h, io_funct.cpp, io_funct.h</em>
 * <br/>
 *
 * <em>art_distance</em> -- it is based on ART 2A-E (Euclidean version of ART). It is 
 * a clustering algorithm for real numbers. It is using various distance measures (not just Euclidean). 
 * The different distance measures can be set by the parameter -d. It consists of <em>art_distance.cpp, 
 * art_distance_alg.cpp, art_distance_alg.h, get_options.cpp, get_options.h, art_common.cpp, art_common.h,
 * io_funct.cpp, io_funct.h</em><br>
 * 
 * In the following figure, there are differences among various ART versions: <img src="art_versions.jpg">
 *  
 * 
 * 
 * <br/>
 * <br/>
 *
 * <h3>Usage</h3>
 * After choosing the right ART program (it is recommended to use art_distance - because it usually 
 * gives the best results), we run the program, e.g.:<br>
 * <code>art_distance -i iris.csv -o output -s 1 -v 0.1</code>
 *
 * A monitor output should look like: 
 * <pre>
 *
 * 
 *
 * Input file: iris.csv
 * Output file: output_results
 * Number of skipped columns: 1
 * Max # of passes: 100
 * Maximal allowed fluctuation: 5%
 * distance:  euclidean
 * beta: 0.5
 * vigilance: 0.1
 * Skip the last 1 columns from input examples
 *
 * Opening file: "iris.csv" ...OK
 * Opening file: "output_results" ...OK
 * Reading input file: iris.csv ...OK (read 150 examples)
 *
 * Removing examples containing only 0's ...OK
 * alpha - default value ((#columns^-0.5)/2): 0.25
 *
 * Input samples normalizing ...OK
 * The number of entering examples into ART: 150
 * Controlling correctness of vigilance against examples ...OK
 * 
 * Computing started
 * Pass: 1, fluctuation: 100%, clusters: 3
 * Pass: 2, fluctuation: 38.7%, clusters: 3
 * Pass: 3, fluctuation: 3.33%, clusters: 3
 *
 * Computing finished
 * Writing results ...OK
 *
 * </pre>
 * 
 * This runs art_distance with input <em>iris.csv</em>. In the first paragraph of monitor output's, 
 * there is information about input parameters. Because this dataset contains also
 * information about a classification class (it was in the last column) we skipped this last one column
 * by setting the <em>-s</em> parameter to 1. Later on, the program removed all
 * examples that contained only 0's. Otherwise, the program would make a never-ending loop. After that, 
 * <em>sample normalization</em> takes a place. This converts all values of every
 * column into the range [0,1]. Then a computation starts. In the first pass,
 * every example was assigned to a new prototype. Therefore fluctuation was 100%.
 * In the second run 38% of examples were re-assigned to another cluster (different
 * than the one to which they were assigned in the previous run). In the final
 * (third) run just 3% of them were re-assigned. Because the maximal allowed
 * fluctuation is 5% (it can be changed by the <em>-e</em> parameter) and the actual fluctuation
 * is 3%, the program stops. 5% is the default value. If the fluctuation after many runs (passes) is
 * still very high the program stops after reaching a maximum number of passes 
 * (<em>-E</em> parameter). The default value is 100 passes.<br>
 * The output is written into <em>output_results</em> file. It is a text file which contains all information
 * about the clustering. There is written almost the same information
 * as which a user can see on the monitor. In the files <em>output_clust_0</em>,...,<em>output_clust_3</em>, there are created clusters
 * (like EM, XMeans).
 * The files contains the examples which belong to the particular cluster
 * (prototype).
 * <br>
 * 
 * An input for the all programs has to be in the csv (comma separated values) format. If the first character
 * on a line is '#' it is considered as a comment. There can be also
 * columns (features) which are not used for computation. For example, if a classification class needs 
 * to be hidden. Exactly it was the case in the previous example. In this case it is necessary to
 * use the parameter <em>-s some_integer</em> which tells how many last columns should not be used
 * for the computation (they will be skipped).<br>
 * 
 * ART based algorithms are sensitive to the example order. It means that for the same dataset with
 * changed example order it can be achieved a different result. Therefore it is a good practice to 
 * change the example order several times and then average the results. This sensitivity to the example order can
 * be compared with using different random seeds in another algorithms. They also need to be run more times with 
 * different random seeds.<br>
 * More examples and comparison with some other clustering methods can be found 
 * <a href="examples.html"> here</a>.
 *
 *
 * <br>
 * <br>
 *
 * <h3>The description of parameters</h3> 
 * <ul>
 * <dt><b>-i</b> specifies the input csv file (dataset) with a set of input examples (sometimes called instances, samples or Ek).
 *  It is the only one mandatory parameter in every program.
 *
 * <dt><b>-e</b> it is a stop criterion. It specifies how many % of input examples can be assigned to another cluster than 
 * in the previous run (pass): how stable the created clusters are. This parameter is called fluctuation. Program will stop if
 * fluctuation is lower than it was defined in the parameter -e. The default value is 5%.
 *
 * <dt><b>-E</b> it is another stop criterion. It defines a maximum number of passes through the input examples.
 * If the algorithm did not find a correct solution (fluctuation is lower than -e) in <em>n</em> passes (runs) through 
 * the examples, it stops. This number <em>n</em> is specified by -E parameter. Then the result is not the last run but
 * the best one. It means that if -E is 50 but the lowest fluctuation was in the 34th run, the result is taken from 34th run.
 * The default value is 100 passes.<br>
 * We usually combine both parameters. For example: try to find a solution with fluctuation less than 10% in maximum
 * 25 passes through examples. In that case the stop criterions look like: <code>-e 10 -E 25</code>
 *
 * <dt><b>-s</b> An integer which specifies how many last columns should not be skipped (not processed by the algorithm).
 * For example:
 * <pre>
 * 5.3,3.7,1.5,0.2,Iris-setosa
 * 5.0,3.3,1.4,0.2,Iris-setosa
 * 7.0,3.2,4.7,1.4,Iris-versicolor
 * 6.4,3.2,4.5,1.5,Iris-versicolor
 * </pre>
 * The last column is some information which can be interesting for a user (therefore it is good to leave it in the dataset)
 * but it is not possible to feed the algorithm with it. Therefore we use <code>-s 1</code> and the column will be skipped during
 * the computation. The default value is 0.
 *
 *<dt><b>-o</b> specifies a prefix for output files. If <code>-o out</code> then <em>out_results</em> will contain every information
 * about the computation and the clusters will be written into <em>out_clust_X</em> where X is the number of a particular cluster. These cluster
 * files contain prototype (the cluster's centroid) and a list of examples which belong to the cluster
 * 
 * <dt><b>-b</b> <em>beta</em> in ART 1 is a small positive integer. It influences a number of created clusters. The higher value the 
 * higher number of created clusters. The default is 1. <br>
 * For all other ART implementations (based on the real-value input), <em>beta</em> is a learning constant. It has a range [0, 1]. The default
 * value is 0.5.
 * 
 * <dt><b>-v</b> <em>vigilance</em> together with alpha (ART for real numbers input) or beta (ART 1) sets up the similarity threshold. 
 * This threshold influences a minimum similarity under which an example Ek will be accepted by the prototype. The higher value 
 * the higher number of clusters. It has a range [0,1]. The default is 0.1. For a deeper understanding take a look at 
 * <a href="art1_draft.jpg">ART 1</a> (<a href="art2_draft.jpg">ART 2</a> respectively) and source codes.
 *
 * <dt><b>-a</b> <em>alpha</em> it is used by the real-value ART algorithms. Together with vigilance it sets up a similarity threshold. 
 * This threshold influences a minimum similarity which is necessary for the example Ek to be accepted by the prototype. 
 * The range: [0,1/sqrt(dim)] where dim is a number of dimensions. The default is 1/sqrt(dim) * 1/2. For a deeper understanding 
 * take a look at  <a href="art2_draft.jpg">ART 2</a> and source codes
 * 
 * <dt><b>-t</b> <em>theta</em>  it is used only by ART 2A. It is a denoising parameter. If a value Ek_i (Ek is an example 
 * and its i-th column) is lower than theta then the number will be changed to 0. Its range: [0,1/dim^-0.5]. The default is 0.00001.
 * 
 * <dt><b>-d</b> <em>distance</em> it is used only by the art_distance. It sets up a distance measure:
    *   <ol><li> Euclidean distance
    *    <li> Manhattan distance
    *    <li> Correlation distance
    *    <li> Minkowski distance
    *   </ol>
 *
 * <dt><b>-p</b> <em>power</em> it is used only for the Minkowski distance in art_distance. It sets up the power for the Minkowski 
 * distance measure. The default value is 3. Minkowski with the power 1 is the Manhattan distance. Minkowski with power 2 is 
 * the Euclidean distance.
 *
 *
 * 
 * </ul>
 *
 *
 *
 * 
 *
 * <br><br>
 * <h3>Compilation</h3>
 * <p>Required software:<br>
 * <ol>
 * <li> A C++ compiler with standard STL library.
 * (the best possibility is GNU gcc 4.1.2 or newer)
 * <li> The mathematical library GNU GSL. It is needed only by
 * art_distance. If you do not need this program you don't need the GNU GSL
 * library. 
 * </ol><br></p>
 * 
 * <p><b>note 1: </b>both recommended packages are freely downloadable and distributed under
 * GNU GPL license. <br>
 * If you would like to make yourself sure that the GSL library is installed correctly:
 * <ol><li> download <a href="http://users.visualserver.org/xhudik/art/gsl_test.cpp"> 
 * http://users.visualserver.org/xhudik/art/gsl_test.cpp </a>
 * <li> run command: <code>g++ gsl_test.cpp -o gsl_test -lgsl -lgslcblas -lm</code>
 * <li> run <code>./gsl_test</code> -- if the result is: J0(5) = -0.177597 -- the GNU GSL software
 * is installed correctly</ol>
 * 
 * Hint: maybe you will have to add some paths pointing to where your GSL
 * library is installed and to it's source codes. In this case the compilation command would looks like:<br> 
 * <code>g++ gsl_test.cpp -o gsl_test -L/usr/local/lib -I/usr/local/include -lgsla -lgslcblas -lm</code><br>
 * (just change the paths according to your own environment). If it is working add the same paths to 
 * the <em>Makefile</em> (line:<br>
 * <code> \$(CC) \$(CFLAGS) \$(COMMON) \$(OBJS4) -o \$(PROG4) -L/usr/local/lib -I/usr/local/include 
 * -lgsl -lgslcblas -lm -lstdc++  -lc  </code><br>
 *
 *
 * 
 * <b>note 2:</b> MS Windows users do not have to compile sources - it is possible to
 * download binaries from <a href="http://users.visualserver.org/xhudik/art">users.visualserver.org/xhudik/art</a>
 * <br/>
 *
 * When the required software is installed, we unpack ART package into a
 * directory. In this directory, run the command: <em>make</em>. After this,
 * everything should be successfully compiled and in the directory there should appear the
 * programs: <em> art_1, art_2A, art_2A-C, art_distance</em><br>
 * A bit deeper description of the compilation process can be found in the package's file: <em>README</em>
 * </p>
 * <br>
 * <br>
 *
 *
 * <h3>A deeper description of source files</h3>
 *
 * <p><em>art_common.cpp, art_common.h</em>: They define and implement functions and
 * types which are used by more than one implemented ART algorithm. (e.g.:
 * Clust -- result's representation, iline -- the represenatation of some line of the input,... )</p>
 * <br/>
 *
 * <p><em>io_funct.cpp, io_funct.h</em>: These files define and implement
 * everything about input and output</p>
 * 
 * <p><em>get_options.cpp, get_options.h</em>: It is Martin Roesch's program (released
 * under GNU GPL license) which parse command lines' input arguments</p>
 *
 * <p>All the other files are about some particular ART clustering algorithm.<br>
 * Note: The documentation can be sometimes confusing. Therefore we recommend
 * to look also into the source codes.<br/>
 *
 * How all the files are interconnected, it is shown in the following picture:
 * <img src="files.jpg">
 * Files in gray squares create the standalone applications.
 *
 * A deeper and more understandable description of ART 1 algorithm itself
 * can be found in <a href="art1_draft.jpg">ART 1</a>. And all other algorithms (ART 2A, ART 2A-C, ART Distance)
 * is shown in <a href="art2_draft.jpg">ART 2</a>. In fact, both of them are very similar.
 *</p>
 *
 * 
 *
 */

namespace ART
{

/**
 * types and functions common for more ART programs
 **/

/**
 * Definition of a simple instance (instance=example (Ek))
 **/
typedef std::vector<double> instance;

/**
 * Definition of a simple prototype. In fact, it is the same as instance.
 * It is easier to read code if it is divided. 
 **/
typedef std::vector<double> prototype;

/**
 * Definition of whole set of  instances
 * */
typedef std::vector< instance > a_instance;

/**
 * Definition of whole set of prototypes
 * */
typedef std::vector< prototype > a_prototype;


/**
 * This struct represents an input line (from csv input file)
 **/
struct iline{

   /**
    * Key are columns of input lines which should not be processed by ART algorithm.
    * The number of unprocessed columns is set by the parameter -s (skip last n columns)
    * */
	std::vector< std::string > key;

   /**
    * It represents all columns which should be processed. In fact, it is 1 .. col - key<br>
    * Where: col - number of all columns<br>
    * key - iline::key
    * */
   std::vector< double > val;

   /**
    * Norms are 2 numbers by which was an example (iline::val) normalized.
    * There are 2 values - because 2x normalization(ART 2A), and max min.
    * It is used only by ART 2A. It is useless in other ART algorithms (norms are set to 1)
    **/
	double norm[2]; 
};


/**
 * All ilines are stored in this type
 * */
typedef std::list< iline > a_lines;

/**
 * The results will be represented in the following form
 * It's a vector not a list!!! So, it is not the same as a_lines
 * */
typedef std::vector< iline > vec_lines;


/**
 * prototype_sequence is a matrix where lines are prototypes (clusters) and
 * rows are ID of instances (examples Ek)
 * */
typedef std::vector < std::vector < unsigned long> > prototype_sequence;


/**
 * In the structure Clust are stored all the results.<br>
 * More specifically: prototypes, fluctuation (error) and sequence
 * of all examples for each prototype
 **/
struct Clust{

/** 
  * proto is a set of created prototypes
  * */
   a_prototype proto;

   /** 
    * proto_seq it is a sequence of sequences (a matrix). Where each line
    * represents one prototype and each column in the line represents some
    * example's ID.<br>
    * Example:<br>
    * 1 2 4<br>
    * 7 3 5
    *<br><br>
    * The first prototype consists of the ID's examples: 1, 2 and 3<br>
    * The second cluster consist of the following examples: 7, 3 and 5<br>
    * An example with ID 5 is a vector. In fact, it is an input line
    * */
      prototype_sequence proto_seq;
 
   /** 
    * How many examples
    * were re-assign (they are in a different cluster then they were before) 
    * */
     float fluctuation;
   };
   

/**
 * In this structure are stored all input parameters important for run every art algorithm. 
 * They can be given by user or they are set as default
 **/
struct in_param{
   /** 
    * Input parameter beta (-b) in ART 1 is a small positive integer. It influences a number of created clusters. The higher value the 
    * higher number of created clusters. The default is 1.<br>
    * For another ART implementations (based on real value input) <em>beta</em> is a learning constant. It has a range [0, 1]. The default
    * is 0.5
    * */
   float beta;

  /**
   * positive integer or 0 - skip the last n columns (in input examples) -- default value:0
   * */
   unsigned long skip;


   /**
    * Input parameter vigilance (-v) together with alpha (ART for real numbers input) or beta (ART 1) set up a similarity threshold. 
    * This threshold influence a minimum similarity under which an example Ek will be accepted by prototype. The higher value 
    * the higher number of clusters. It has a range [0,1]. The default is 0.1.
    * */
   float vigilance;

   /** 
    * Input parameter theta (-t) denoising parameter. If a value Ek_i (Ek is a example and
    * it's ith column) is lower than theta then the number will be changed to 0. It is used only by ART 2A. It's range:
    * [0,1/dim^-0.5]. Default 0.00001
    * */
   float theta;

   /** 
    * Input parameter alpha (-a)  is used by real value ART algorithms. Together with vigilance set up a similarity threshold. 
    * This threshold influences a minimum similarity which is necessary for the example Ek to be accepted by the prototype. 
    * The range: [0,1/sqrt(dim)] where dim is a number of dimensions. The default is 1/sqrt(dim) * 1/2
    * */
   float alpha;

   /** 
    * Input parameter distance (-d) set up a distance measure:
    *   <ol><li> Euclidean distance
    *    <li>Modified Euclidean distance -- it is in a testing mode. Euclidean distance use equation 1 - E/dim where E is Euclidean distance
    *    and dim is a number of dimensions. Modified Euclidean use equation log(dim^2) - E. This distance in some cases can achieve a better
    *    performance than standard Euclidean distance. However, it is recommended to use standard Euclidean distance. DO NOT USE IT
    *    <li> Manhattan distance
    *    <li> Correlation distance
    *    <li> Minkowski distance
    *   </ol>
    *    It works only for art_distance. Default Euclidean distance measure
    **/ 
   int distance;

   /** 
    * Input parameter power (-p) it is used only for Minkowski distance in art_distance. It set up the power for Minkowski 
    * distance measure. The default is 3. Minkowski with the power 1 is Manhattan distance. Minkowski with power 2 is 
    * Euclidean distance
    * */
   int power;

   /**
    * An input parameter --  a number of passes (-E), it is a maximum number of how many times an example Ek 
    * can be re-assigned. If it reach this number the program will stop. The default is 100
    * */
   unsigned long pass;

   /** 
    * An input parameter -- fluctuation (-e), it is a highest possible error rate (%). It means a maximum
    * number (in %) of how many instances can be re-assign. If the real fluctuatio is lower than -e
    * then program will stop. Default is 5% examples.
    * */
   float error;
};
   
   
typedef Clust Clusters;
   

/**
 * This structure stores a minimum and maximum for some feature (column). 
 * Not just from one line but from all of them. Therefore it can't be in 
 * i_line structure.  This is useful only for ART 2A-C and ART Distance
 **/
struct mmax{

   /**
    * A minimum value of all examples (Ek) for particular column
    * */
	double min;

   /**
    * A maximum value of all examples (Ek) for particular column
    * */
   double max;
	};
typedef std::vector < mmax > minmax;


/**
 * Data normalizing (it is NOT the same as vector normalizing)
 * It is used in ART 2A-C and ART Distance
 * */
void data_normalize(a_lines& lines, minmax& mm);


/**
 * Is the example (Ek) included in some prototype already?
 * 
 * instance_in_sequence searching through the prototype sequence (prot_seq).
 * The goal is to find out the index of specified (iinst) sample(example Ek). This index is given 
 * as a result. If the example is not contained in some existing prototype the result will be "not_found"
 * It is used by all ART algorithms
 * */
unsigned long instance_in_sequence(prototype_sequence& prot_seq, unsigned long iinst);


/**
 * Find a particular instance or prototype (inst) in a sequence (sample). The result is the index 
 * of the instance or prototype in the sequence. The program will stop if must_find=true and the 
 * algorithm don't find the particular instance 
 * */
unsigned long find_item(a_instance& sample, instance& inst, bool must_find);

/**
 * It converts the part of input lines which is going to be processed by the ART algorithm  into vectors.
 * Vectors are simpler to work with
 **/
a_instance toVectors(a_lines& lines);

}; //namespace ART

#endif

