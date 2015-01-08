
#ifndef _IOFUNCT_
#define _IOFUNCT_

#include <fstream>
#include <list>
#include <iostream>
#include "art_common.h"

/**
 * Reading and writing files -- all io functions for ART
 **/
 

/**
 * A generic function for opening an input fstream (file)
 * @param name name of the file
 * @param writeInfo  after opening should be written " ... OK" ?
 * @param file input  file stream
 * Output: error state
 * */
int openFile(std::ifstream & file, const char* name, bool writeInfo=false);


/**
 * A generic function for opening an output fstream (file)
 * @param name name of the file
 * @param writeInfo  after opening should be written " ... OK" ? 
 * @param file input file stream
 * Output: error state
 * */
int openFile(std::ofstream & file, const char* name, bool writeInfo=false);

/**
 * This function read an input csv file.
 * @param file input file stream
 * @param lines read lines of the input
 * @param skip how many columns should not be processed
 * if skip=4 -- last 4 columns won't be processed
 * */
void readFile(std::ifstream& file,ART::a_lines& lines, unsigned long skip);


/**
 * Check if all the input examples are binary (only 0 or 1)
 * Don't consider skip columns.
 * Used by art_1
 * */
void check_binary(ART::a_lines& s);

/**
 * Check if every column of every input example is bigger than 0 (positive values are allowed only)
 * Used by ART 2A
 * */
void check_positiveness(ART::a_lines& s);

/**
 * Write output results to the output file
 * @param out output file stream
 * @param clust_prefix which prefix output files should have. An user can set this by 
 * the parameter -o
 * @param lines examples Ek
 * @param clust created clusters (prototypes and ordered examples for every cluster)
 * @param zero how many input examples were composed just out of 0's (they can't be processed)
 * */
void write_results(std::ofstream& out, const char* clust_prefix,ART::vec_lines& lines, ART::Clusters&  clust,unsigned long zero);


/**
 * Remove examples composed just from 0's.
 * The output is number of these examples
 * */
unsigned long zero_free(ART::a_lines& lines);

/**
 * Overloaded operator for printing out instances and prototypes
 * */
std::ostream& operator<<(std::ostream& out, const ART::instance& inst);

/**
 * Overloaded operator for printing out whole lines. Also with columns which were not processed
 * (parameter -s)
 * */
std::ostream& operator<<(std::ostream& out, const ART::iline line);

#endif

