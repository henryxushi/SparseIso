/*
 *  common.h
 */
#include "options.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <libgen.h>
#include <iostream>
#include <sstream>

string OPT_bamfile = "NaN";
string OPT_outdir = "NaN";
string OPT_readtype = "NaN";
int OPT_Ncores = 1;
double OPT_conf = 0.5;
bool OPT_help = false;
double OPT_mu = 200;
double OPT_sigma = 20;
string OPT_cempath = "";
string OPT_samtoolspath = "";
bool OPT_output_all = false;

void Options::assign(string s1, string s2, string s3, int s4, double s5, double s6, double s7, string s8, string s9, bool s10)
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
}

int Options::parse_options(int argc, char* argv[]) {

  int option_index = 0;
  int next_option;

  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
    case -1:     /* Done with options. */
      break;
    case 'b':
      OPT_bamfile = optarg;
      break;
    case 'o':
      OPT_outdir = optarg;
      break;
    case 'r':
      OPT_readtype = optarg;
      break;
    case 'h':
      OPT_help = true;
      break;
    case 'p':
      OPT_Ncores = atoi(optarg);
      break;
    case 'c':
      OPT_conf = atof(optarg);
      break;
    case 'm':
      OPT_mu = atof(optarg);
      break;
    case 'v':
      OPT_sigma = atof(optarg);
      break;
    case OPT_CEMPATH:
      OPT_cempath = optarg;
      break;
    case OPT_SAMPATH:
      OPT_samtoolspath = optarg;
      break;
    case OPT_OUTPUT_ALL:
      OPT_output_all = true;
      break;
    default:
      std::cout << usage();
      exit(1);
    }
  } while(next_option != -1);

  if (OPT_help) {
    std::cout << usage() ;
    exit(1);
  }
  
  //check input argument
  
  if (OPT_bamfile == "" || OPT_bamfile == "NaN") {
  	std::cerr << "Please check bam file option\n";
    std::cout << usage() ;
    exit(1);
  }
  
  if (OPT_readtype != "p" && OPT_readtype != "s")
  {
  	std::cerr << "Please check read type option\n";
  	std::cout << usage() ;
    exit(1);
  }
  
  if (OPT_outdir == "" || OPT_outdir == "NaN") {
    std::cerr << "Please check output file directory option\n";
    std::cout << usage() ;
    exit(1);
  }

  if (OPT_Ncores < 0)
  {
  	std::cerr << "Please check threads option\n";
  	std::cout << usage() ;
    exit(1);
  }
  
  if (OPT_conf < 0 || OPT_conf > 1)
  {
  	std::cerr << "Please check confidence option\n";
  	std::cout << usage() ;
    exit(1);
  }

  if (OPT_mu < 0)
  {
    std::cerr << "Please check the mean of fragment length distribution\n";
    std::cout << usage() ;
    exit(1);
  }

  if (OPT_sigma < 0)
  {
    std::cerr << "Please check the standard deviation of fragment length distribution\n";
    std::cout << usage() ;
    exit(1);
  }
  
  assign(OPT_bamfile, OPT_outdir, OPT_readtype, OPT_Ncores, OPT_conf, OPT_mu, OPT_sigma, OPT_cempath, OPT_samtoolspath, OPT_output_all);

  return 0;
}


string Options::usage () {

  std::stringstream usage_info;
  usage_info
    << std::endl
    << "===============================================================================" << std::endl
    << " Usage: SparseIso [--bam/b] <filename>  [opts] " << std::endl
    << "===============================================================================" << std::endl
    << " **Required :" << std::endl
    << " --bam/-b <string>         " << ": the name of bam file" << std::endl
    << " --out-dir/-o <string>     " << ": the directory stored all output files" << std::endl
    << " --readtype/-r <p/s>       " << ": the type of reads paired-end(p) or single-end(s)" << std::endl;

  usage_info
    << std::endl
    << " ** Prerequisite program (not needed to set if in system path) :" << std::endl
    << " --cempath <string>      " << ": the path to processsam." << std::endl
    << " --sampath <string>      " << ": the path to samtools." << std::endl;

  usage_info
    << std::endl
    << " ** Optional :" << std::endl
    << " --threads/-p <int>        " << ": the number of threads to be used, default: 1." << std::endl
    << " --outputall               " << ": output all isoforms enumberated from the graph." << std::endl
    << " --conf/-c <double>        " << ": the confidence level for isoform identification, default: 0.5." << std::endl
    << " --mean/-m <double>        " << ": the mean of fragment length distribution, default: 200." << std::endl
    << " --stdvar/-v <double>      " << ": the variance of fragment length distribution, default: 20." << std::endl
    << " --help/-h                 " << ": display the help information."<< std::endl
    << std::endl;
    
  usage_info
    << "================================================================================" << std::endl
    << std::endl;

  return (usage_info.str());
}

