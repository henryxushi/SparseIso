#ifndef READINSTANCE_H
#define READINSTANCE_H

#include "utility.h"
#include "options.h"
using namespace std;



class readinstance
{
public:
	readinstance(){};
	readinstance(Options &in_opt)
	{
		opt = in_opt;
	};
	Options opt;
	vector<Info> infolist;
	void readinstance_se(string inputfile);
	void readinstance_p(string inputfile, double frag_mu, double frag_sigma);
	void getInfoList(vector<Info>& outinfolist);
	int totalNumofPaths;
	double TOTALNUMREADS;
	void runSparseIso_Threaded();
	
};


#endif
