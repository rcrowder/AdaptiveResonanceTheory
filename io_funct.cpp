#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <climits>
#include <list>
#include "io_funct.h"
//#include "art_1_alg.h"

using namespace std;
using namespace ART;

/**
 * The number of columns in the first line.
 * The same number have to have all others input examples
 * */
unsigned long columns=0;

/**
 * It tells whether or not the first line has been read already. 
 * If it so - variable columns (dimensionality) was set also
 * */
bool first_read=false;

int openFile(ifstream& file, const char* name, bool writeInfo){
	if(writeInfo) cout << "Opening file: \"" << name << "\"" << flush;
	file.open(name);
	if(file.fail()) {
	   cerr << "\n\n\n!!!ERROR: I can not open file: \"" << name << "\"\n\t\aExiting ...\n" << endl;
	   exit(10);
	}
	if(writeInfo) cout << " ...OK\n" << flush;
	return file.fail();
}

int openFile(ofstream& file, const char* name, bool writeInfo){
	file.open(name);
	if(writeInfo) cout << "Opening file: \"" << name << "\"" << flush;
	if(file.fail()) {
	   cerr << "\n\n\n!!!ERROR: I can not open file: \"" << name << "\"\n\t\aExiting ...\n" << endl;
	   exit(11);
	}
	if(writeInfo) cout << " ...OK\n" << flush;
	return file.fail();
}



ostream& operator<<(ostream& out, const instance& inst){
   for(unsigned long i=0;i<(inst.size()-1);i++) out << inst[i] << ",";
   out << inst[inst.size()-1];
   return out;
}


ostream& operator<<(ostream& out, const iline line){
	const unsigned long keysize = line.key.size();
	const unsigned long valsize = line.val.size();
	
	for(unsigned long i=0;i<(valsize-1);i++)out << line.val[i] << ",";
	out <<line.val[valsize-1];
	if(keysize>0){
		out << ",";
		for(unsigned long i=0;i<(keysize-1);i++)out << line.key[i] << ",";
		out << line.key[keysize-1];
		}
   return out;
}



/**
 * The predicate function for the function zero_free.
 * If the particular line is composed just from 0's then true
 * otherwise false
 * */
inline bool zero(iline& line){ 
	unsigned long size = line.val.size();
	instance tmpvect;
	for(unsigned long i=0;i<size;i++) tmpvect.push_back(0);
	if(line.val==tmpvect)return true;
		else return false;
	}


unsigned long zero_free(a_lines& lines){
	//the size of all lines before nullifying
	unsigned long beforesize = lines.size();
	
	//remove samples where every member of val part is 0 (keys can be whatever)
	lines.remove_if(zero);
	
	//the size of all lines after nullifying
	unsigned long aftersize = lines.size();	
//return the number of instances (lines), where all the values were 0 (not keys)
return beforesize - aftersize;
}


void check_binary(a_lines& s){
   //check every number
   list <iline>::iterator itline;
   itline=s.begin();
   iline line = *itline;
   const unsigned long size = line.val.size();
   unsigned long index = 0;
   
   for(itline=s.begin();itline!=s.end();itline++){
   	line=*itline;
      for(unsigned long j=0;j<size;j++)
		if((line.val[j]!=0)&&(line.val[j]!=1)){
	       cerr << "\n\nERROR value at the line " << index+1 << ",column " << j+1;
	       cerr << " is "<<line.val[j]<<". It can be 1 or 0 (ART 1 is just for binary input)\nWhole data are in data_tmp.";
	       cerr << "\n\n Exiting...\n" << endl;
	       exit(22);
	       }
	   index++;
   }

} 

void check_positiveness(a_lines& s){
   //check every number
   list <iline>::iterator itline;
   itline=s.begin();
   iline line = *itline;
   const unsigned long size = line.val.size();
   unsigned long index = 0;
   
   for(itline=s.begin();itline!=s.end();itline++){
   	line=*itline;
      for(unsigned long j=0;j<size;j++)
		if(line.val[j]<0){
	       cerr << "\n\nERROR: value at the line " << index+1 << ",column " << j+1;
	       cerr << " is "<<line.val[j]<<". Only positive values are allowed for ART 2A\nWhole data are in data_tmp.";
	       cerr << "\n\n Exiting...\n" << endl;
	       exit(22);
	       }
	   index++;
   }

}


/**
 *Process some input line. Read it from input stream s. The output is read line
 * */
inline iline processLine(istream& s, const unsigned long skip ){
   float number;
   iline res;
   unsigned long index = 0;
   const unsigned long skip_from = columns - skip;
   while(!s.eof()){
      string snumber;
      getline(s,snumber,',');
      //insert string if the column should be skipped
      if(index>skip_from) res.key.push_back(snumber);
      	else {
      	//this row won't be skipped - transform on float
      	number=(float)atof(snumber.c_str());
      	res.val.push_back(number);
      	}
      index++;
   }
   //initial set up of normalization multilayer
   res.norm[0]=1;
   res.norm[1]=1;
   return res;
}


void readFile(ifstream& file,a_lines& lines, unsigned long skip){
   string readline;
   //initialization of line size
   unsigned long line_size=0;
   //number of lines (vectors, samples)
   unsigned long nlines=1;
   bool first_run=true;
   
   while(!file.eof()){
      getline(file,readline,'\n');
      //treating the last line and empty and comment lines
      if((readline.size()<1)||(readline[0]=='#')) continue;
      
      stringstream sline(stringstream::in|stringstream::out);
      sline.str(readline);
      //if the first line has not been read yet - find out a number of columns
   	if(first_read==false){
   		string str=string(sline.str());
   		while(str.find(",")<str.size()){
   			str=str.substr(str.find(",")+1);
   			columns++;
   		}
   		first_read=true;
   	}
      iline inst;
      inst=processLine(sline,skip);
      //size of the current line
      unsigned long actsize = inst.key.size()+inst.val.size();
      if(!first_run) if(line_size==actsize) lines.push_back(inst);
		     else{
			if(actsize>0){
			cerr << "\n\n!!!ERROR the line " << nlines << " has ";
			cerr << actsize << " columns. Previous lines had ";
			cerr << line_size <<" columns ...\n\nExiting...\n\n";
			exit(11);
			}else cerr << "\nline " << nlines << "is empty - skipping " << endl;
		     }
		  else{
		  	  if(first_run){
		     	first_run=false;
		     	line_size=inst.key.size()+inst.val.size();
		     	}
		     lines.push_back(inst);
		  }
     nlines++;
     //if number of samples is too big
     if(nlines==ULONG_MAX) {
     	cerr << "The number of samples (" << nlines << ")reached the maximum for your architecture.\n";
     	cerr << "You should decrease the number of samples, or change computer"; 
     	cerr << "(with higher number of bits)\n\nExiting ... \n" << endl;
     	exit(5);
     }
   }

}



void write_results(std::ofstream& out, const char* clust_prefix, vec_lines& lines, Clusters&  clust, unsigned long zero){
	out << "Number of clusters: " << clust.proto.size() << endl;
	
	out << "Fluctuation: " << clust.fluctuation << endl;
	out << "\n\nPrototypes and cluster files\n";
	out <<"(First column is a file name. Second column is a prototype)"<< endl;
	for(unsigned long i=0;i<clust.proto.size();i++){
		out << clust_prefix << "_clust_" << i << ".csv" << "\t" << clust.proto[i] << "\n";
	}
	
	for(unsigned long i=0;i<clust.proto_seq.size();i++){
		ofstream of;
		stringstream name;
		name << clust_prefix << "_clust_" << i << ".csv";
		openFile(of,name.str().c_str(),false);
		of << "# The list of instances for the prototype: " << clust.proto[i] << "\n\n";
		for(unsigned long j=0;j<(clust.proto_seq[i].size()-1);j++){
			of << lines[clust.proto_seq[i][j]] << "\n";
		}
		of << lines[clust.proto_seq[i][clust.proto_seq[i].size()-1]] << "\n";
		of.close();
	} 
	//if in the input were examples which contained only 0's - it is a special cluster
	if(zero>0){
		out << "A special cluster -- were are the examples contained only 0's\n";
		//out << "(these examples are not listed in the other parts of the results)\n";
		
		//write down all 0's samples
		for(unsigned long i=0;i<zero;i++){
			for(unsigned long j=0;j<(lines[0].val.size()-1);j++) out << "0,";
			if(i<(zero-1))out << "0;\n";
			else out << "0\n";
		}
	}
	out << flush;
}


