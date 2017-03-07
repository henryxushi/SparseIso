#ifndef INFO_H
#define INFO_H

#include <vector>
#include <string>
#include <set>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

class Info
{
public:
	Info(){
		single = true;
	}
	MatrixXd X;
	MatrixXd X_m; // X with length; preference of select: long segment to short segment
	VectorXd sumX_m; // sum X horizontally
	MatrixXd X_exon;
	VectorXd y;
	VectorXd y_m;
	VectorXd R;
	VectorXd L;
	VectorXd bias;
	VectorXd transL;
	VectorXd exonL;
	VectorXd lw; //length weight of segment
	vector<vector<int> > exonbound;
	string label;
	vector<int> final_isoidx;
	vector<double> final_FPKM;
	MatrixXd beta_frac_low_high;
	int totalReads;
	int infoidx;
	void run();
	void run_single();
	void clear();
	void clear_info();
	bool single;
	bool bias_flag;
	void preprocess();
	void modX_exon(vector<vector<int> > &paths);
	void modX(vector<vector<int> > &paths);
	void mody(vector<double> &s_y);
	void addy(double s_y);
	void modR(vector<double> &s_R);
	void addR(double s_R);
	void modL(vector<double> &s_L);
	void addL(double s_L);
	void modexonbound(vector<vector<int> > &s_bound);
	void modbias(vector<double> &s_bias);
	void correct_bias();
	void getIsoIdxForX(vector<int> &current_idx, vector<int> &new_idx); //choose the isoform to explain all reads
	void calFPKM(VectorXd s_beta, VectorXd& FPKM, VectorXd& s_reads);
	void write(string outdir, int idx, bool output_all_bool);
	void getExonRegion(int isoidx, set<set<int> >& exon_region); //get the exon boundaries for isoform isoidx
	double get_muth(double p, int idx, VectorXd s_beta);
	bool output_bool;
	bool valid;
	//bool output_all_bool;
};

#endif
