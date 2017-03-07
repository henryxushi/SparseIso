#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <stdio.h>
#include <list>
#include <math.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <set>
#include <map>
#include <limits.h>
#include "utility.h"
#include "readinstance.h"
#include "options.h"
#include "Info.h"
#include "filter_bam.h"
#include "getinsertsize.h"
#include <boost/threadpool.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
using namespace std;

class Job
{
public:
	Job(Info* info_in, int totalReads_in, int idx_in) 
	{
		info = info_in;
		totalReads = totalReads_in;
		idx = idx_in;
	}
	~Job() {}
	void run()
	{
		//cout << "Sampling\t" << idx << endl;
		info->infoidx = idx;
		info->totalReads = totalReads;
		if (info->single)
		{
			info->run_single();
		}
		else
		{
			info->run();
		}
	}
private:
	Info* info;
	int totalReads;
	int idx;
};

int main(int argc, char* argv[])
{
	Options opt;
	if (opt.parse_options(argc,argv))
		exit(-1);

	if (!boost::filesystem::exists(opt.outdir))
		boost::filesystem::create_directory(opt.outdir);

	string outbamfile = opt.bamfile;
	// filter reads for paired end data
	stringstream samtools_stream, instance_stream, outinstance_stream;
	stringstream samfile;
	stringstream instlog_stream;
	if (opt.readtype == "p")
	{
		cout << "Input paired-end data" << endl;
		cout << "Filtering bam file..." << endl;
		filter_bam(opt.bamfile,opt.outdir,outbamfile);
	}
	else
		cout << "Input single-end data" << endl;

	// run processsam
	outinstance_stream << opt.outdir << "/sparseiso_inst";
	samfile << opt.outdir << "/sparseiso_sam.sam";
	instlog_stream << opt.outdir << "/inst_log.txt";
	
	if (opt.samtoolspath.size()>0)
	{
		if (opt.samtoolspath.at(opt.samtoolspath.size()-1) != '/')
			opt.samtoolspath.append("/");

	}

	if (opt.cempath.size()>0)
	{
		if (opt.cempath.at(opt.cempath.size()-1) != '/')
			opt.cempath.append("/");

	}

	samtools_stream << opt.samtoolspath << "samtools view " << outbamfile << " > " << samfile.str();
	if (system(NULL) == 0)
	{
		cerr << "Processor not available" << endl;
		exit(-1);
	}
	else
		assert(system(samtools_stream.str().c_str())==0);
	double mu = 200, sigma = 20;
	if (opt.readtype == "p")
	{
		//estimate fragment length distribution
		getinsertsize(samfile.str(),mu,sigma);
		cout << "Estimated mean and stdvar: " << mu << " " << sigma << endl; 
	}
	instance_stream << "nice " << opt.cempath << "processsam --no-coverage -c 1 -g 1 -u 0.05 -d . -o " << outinstance_stream.str() << " " << samfile.str() << " > " << instlog_stream.str() << " 2>&1 "; 
	if (system(NULL) == 0)
	{
		cerr << "Processor not available" << endl;
		exit(-1);
	}
	else
		assert(system(instance_stream.str().c_str())==0);
	outinstance_stream << ".instance";
	string outinstance = outinstance_stream.str();
	//assert(system(instance_stream.str().c_str()) == 0);

	//estimate fragment length distribution and estimate bias
	//we need first read the instance file then we can get the fragment length and bias
	//because we need to get the genes with single isoform
	/*vector<double> bias(100,1);
	//read bias from ground truth
	string biasfile = "/media/73a039e0-7805-4ca6-90d0-b63ed265dcfe/Henry/RNA-seq/SparseIso_sim/RNASeqReadSimulator-master/demo/input/sampleposbias1.txt";
	ifstream biasreader;
	biasreader.open(biasfile.c_str());
	string line;
	int biasidx = 0;
	if (biasreader.is_open())
	{
		while (biasreader.good())
		{
			getline(biasreader,line);
			if (!line.empty())
			{
				bias[biasidx++] = atof(line.c_str());
			}
		}
	}*/


	//get info list in read instance class

	/*string filename1 = argv[1];
	string outputfile = argv[2];
	double frag_mu = atof(argv[3]);
	double frag_sigma = atof(argv[4]);
	string readtype = argv[5];
	string Numcores = argv[7];

	Options opt;
	opt.readtype = readtype;
	opt.outdir = argv[6];
	opt.N_cores = atoi(Numcores.c_str());*/
	if (opt.N_cores > 1)
		cout << "Multi_core enabled " << endl;

	cout << "Processing instance..." << endl;
	readinstance loadData(opt);
	loadData.TOTALNUMREADS = 0;
	int totalReads = 0;
	if (opt.readtype == "p")
	{
		loadData.readinstance_p(outinstance, mu, sigma);
		totalReads = round(loadData.TOTALNUMREADS/2);
	}
	else
	{
		loadData.readinstance_se(outinstance);
		totalReads = round(loadData.TOTALNUMREADS);
	}
	
	//vector<Info> infolist;
	//loadData.getInfoList(infolist);

	cout << "Start sampling..." << endl;
	if (opt.N_cores > 1)
	{
		boost::threadpool::pool multi_tp(opt.N_cores);
		for (int i = 0; i < loadData.infolist.size(); i++)
		{
			if (!loadData.infolist[i].valid)
				continue;
			if (loadData.infolist[i].valid)
			{
				//cout << "Sampling\t" << i << endl;
				boost::shared_ptr<Job> job(new Job(&loadData.infolist[i],totalReads,i));
				multi_tp.schedule(boost::bind(&Job::run,job));
			}
			
		}
		multi_tp.wait();
	}
	else
	{
		for (int i = 0; i < loadData.infolist.size(); i++)
		{
			if (!loadData.infolist[i].valid)
				continue;
			//cout << "Sampling\t" << i << endl;
			loadData.infolist[i].infoidx = i;
			loadData.infolist[i].totalReads = totalReads;
			if (loadData.infolist[i].single)
			{
				loadData.infolist[i].run_single();
			}
			else
			{
				loadData.infolist[i].run();
			}
		}
		

		
	}
	cout << "Start writing..." << endl;
	for (int i = 0; i < loadData.infolist.size(); i++)
	{
		if (!loadData.infolist[i].valid)
			continue;
		if (loadData.infolist[i].single)
		{
			if (loadData.infolist[i].output_bool)
			{
				//cout << "Writing\t" << i << endl;
				loadData.infolist[i].write(opt.outdir,i+1,opt.output_all);
			}
		}
		else
		{
			//cout << "Writing\t" << i << endl;
			loadData.infolist[i].write(opt.outdir,i+1,opt.output_all);
		}
		
	}
	/*cout << "Finished reading" << endl;
	ofstream outfile1;
	outfile1.open(outputfile.c_str(),std::ios_base::app);
	outfile1 << "totalNumReads:\t" << round(loadData.TOTALNUMREADS/2) << endl;*/
	

	
	return 0;
}
