#ifndef OPTIONS_H
#define OPTIONS_H

#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/typesizes.h>
#include <bits/types.h>
#include <getopt.h>
using namespace std;

class Options
{
public:
	Options(){};
	Options(string s1, string s2, string s3, int s4, double s5, double s6, double s7, string s8, string s9, bool s10)
	{
		bamfile = s1;
		outdir = s2;
		readtype = s3;
		N_cores = s4;
		conf = s5;
		mu = s6;
		sigma = s7;
		cempath = s8;
		samtoolspath = s9;
		output_all = s10;
	};

	string bamfile;
	string outdir;
	string readtype;
	int N_cores;
	double conf;
	double mu;
	double sigma;
	string cempath;
	string samtoolspath;
	bool output_all;
	void assign(string s1, string s2, string s3, int s4, double s5,  double s6, double s7, string s8, string s9, bool s10);
	int parse_options(int argc, char* argv[]);
	string usage();
};


extern string OPT_bamfile;
extern string OPT_outdir;
extern string OPT_readtype;
extern int OPT_Ncores;
extern double OPT_conf;
extern bool OPT_help;
extern double OPT_mu;
extern double OPT_sigma;
extern string OPT_cempath;
extern string OPT_samtoolspath;
extern bool OPT_output_all

static const char *short_options = "b:r:o:p:c:h:m:v";

#define OPT_CEMPATH 301
#define OPT_SAMPATH 302
#define OPT_OUTPUT_ALL 303


static struct option long_options[] = {
  {"bam",                 required_argument,      0,      'b'},
  {"help",                no_argument,            0,      'h'},
  {"out-dir",			  required_argument,      0,      'o'},
  {"readtype",            required_argument,      0,      'r'},
  {"threads",             required_argument,      0,      'p'},
  {"conf",                required_argument,      0,      'c'},
  {"mean",                required_argument,      0,      'm'},
  {"stdvar",                required_argument,      0,      'v'},
  {"cempath",                required_argument,      0,      OPT_CEMPATH},
  {"sampath",                required_argument,      0,      OPT_SAMPATH},
  {"outputall",                no_argument,            0,      OPT_OUTPUT_ALL},
//...
  {0,0,0,0} // terminator

};

#endif

