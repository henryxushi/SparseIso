#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "Info.h"
using namespace std;

const int MAX_NUM_GENE = 100000;
const int MAX_NUM_PATH = 100;
const double MAX_INTRON_COV = 0.95;

struct Seg
{
	int node;
	vector<int> visited;
};

class Instance
{
public:
	Instance(){}
	string location;
	string strand;
	vector<vector<int> > exonbound;
	vector<double> exonAbundance;
	vector<vector<int> > SegConf;
	vector<vector<int> > SegConfPET;
	int ReadLen;
	int NumofExon;
	int NumofSegs;
	double maxExonCov;
	double frag_mu;
	double frag_sigma;

	vector<int> SegAbundance;
	vector<int> connSeg; // find SGTypes with connection, count number of exons
	vector<vector<int> > PETconnect;
	vector<int> PET_switch;
	vector<vector<vector<int> > > PETreads;
	vector<int> SegAbundancePET;
	vector<set<int> > inc_exon;
	
	int NumofPaths;
	
	void process_inst(Info& info);
	void process_inst_se(Info& info);
	void correct_y(vector<vector<int> > &paths,Info &info);
	void remove_low();
	void getmap(map<int,set<int> > &exonmap, double &th);
	int enumerate_path(map<int,set<int> > &exonmap, int not_source_list[], int not_sink_list[], vector<vector<int> > &paths);
	void mod_X(vector<int> &path, vector<double> &proc_abun,vector<int> &segpath);
	void mod_y(map<int, vector<int> > &iscontained, vector<double> &SegAbundance1, int &output, double &th1, double &min_abun);
	void mod_y_se(map<int, vector<int> > &iscontained, vector<double> &SegAbundance1, vector<double> &SegReads, vector<double> &SegLength, int &output, double &th1, double &min_abun);
	void findcontain(map<int, vector<int> > &iscontained); // contain: index start from 1
	int calintronlength(int &s1, int &s2, vector<int> &path, int &l1, int &l2, int &l3, int label, int &bias);
	void intronlength_helper(int &pos1, int &pos2, vector<int> &path, int &intron_l, int &l3, int &bias);
	void exonlength_plus(int &pos1, int &pos2, vector<int> &path, int &intron_l,int &l3);//calculate exon length between two points
	void exonlength_reverse(int &pos1, int &pos2, vector<int> &path, int &intron_l, int &l3);
	int segStart(int &s1, int &s2);
	void segStart(int &s1, int &s2, set<int> &s3);
	int segStart(int &s1);
};

#endif
