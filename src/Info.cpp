#include "Info.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <sstream>
#include "sampling.h"
#include "basic_math.h"
//#include <boost/date_time/posix_time/ptime.hpp>
//#include <boost/date_time/microsec_time_clock.hpp>
using namespace std;

void print_v(vector<double> v);
void print_vi(vector<int> v);
void print_vv(vector<vector<double> > v);
bool vector_matchX(MatrixXd X, int a, int b);
void vectorcp(vector<double> &a, vector<double> &b); // cp a to b
void matcp(vector<vector<double> > &a, vector<vector<double> > &b);
void vectorcptoMat(vector<vector<double> > &a, MatrixXd &b);
void connect_region_recur(set<set<int> > region_in, set<set<int> >& region_out);
void connect_region(set<set<int> > region_in, set<set<int> >& region_out);


void Info::run()
{
	//cout << "Start sampling" << endl;
	int N = X.rows();
	//assert(N > 0);
	int M = X.cols();

	//cout << X << endl;
	//cout << y << endl;
	//cout << "---------------------------" << endl;
	if (M > 1 and N > 1)
	{
		bool verbose = true;
		preprocess();
		N = X.rows();
		M = X.cols();
		if (true)
		{
			/*cout << "----------y----------" << endl;
			cout << y << endl;
			cout << "----------X----------" << endl;
			cout << X << endl;
			cout << "---------------------" << endl;*/
		}
		//cout << "Finish preprocessing " << endl;
		double v0 = pow(10,-6);
		vector<int> current_idx;
		vector<int> new_idx;
		getIsoIdxForX(current_idx,new_idx);
		int niso = current_idx.size() + new_idx.size();

		Sampler sample;
		//initialize parameters
		//cout << "Initializing parameter" << endl;
		double a1 = 3, a2 = 4;
		double a_sigma0 = pow(10,2), b_sigma0 = pow(10,-2);
		int iter = 500;
		//int nburn = 300;

		VectorXd s_rho;
		s_rho = y;
		VectorXd s_tau = VectorXd::Ones(M);
		VectorXd s_beta = VectorXd::Zero(M);
		VectorXd s_gamma = VectorXd::Ones(M);
		VectorXd s_z = v0 * VectorXd::Ones(M);
		/*for (int i = 0; i < M; i++)
			s_z[i] = v0;
		*/
		for (int i = 0; i < new_idx.size(); i++)
		{
			//cout << new_idx[i] << endl;
			s_z(new_idx[i]) = 1;
		}
		
		//s_gamma = s_z.array() * s_tau.array(); 
		double s_sigma2 = 1e-5;
		double s_w = niso / M;
		vector<VectorXd > samples_beta(iter,VectorXd(M));
		vector<VectorXd > samples_z(iter,VectorXd(M));

		MatrixXd XTX;
		MatrixXd XT = X.transpose();
		XTX = XT * X;
		//MatrixXf XT;
		
		//transpose(X,XT);
		//matMul(XT,X,XTX);

		//cout << "Start sampling iteration" << endl;
		for (int i = 0; i < iter; i++)
		{
			//cout << "Sampling iteration " << i << endl;
			MatrixXd A;
			MatrixXd diag_gamma;

			vector<double> diag_v(M,0);
			for (int k = 0; k < M; k++)
				diag_v[k] = s_sigma2/s_gamma(k);
			diag(diag_v,diag_gamma);
			/*if (infoidx == 5 and verbose)
			{	
				cout << "----------diag_gamma----------" << endl;
				cout << diag_gamma << endl;
			}*/
			//cout << "Finish_1 " << i << endl;
			//matSum(XTX,diag_gamma,A);
			A = XTX + diag_gamma;
			MatrixXd invA;
			invA = A.inverse();
			
			MatrixXd mu_betaM = invA * XT * s_rho;
			assert(mu_betaM.cols()==1);
			VectorXd mu_beta(Map<VectorXd>(mu_betaM.data(),mu_betaM.rows()*mu_betaM.cols()));

			MatrixXd invsigma_beta = 1/s_sigma2 * A;
			
			/*if (infoidx == 5 and verbose)
			{	
				cout << "----------XTX----------" << endl;
				cout << XTX << endl;
				cout << "----------diag_gamma----------" << endl;
				cout << diag_gamma << endl;
				cout << "----------A----------" << endl;
				cout << A << endl;
				cout << "----------invA----------" << endl;
				cout << invA << endl;
				cout << "----------mu_beta----------" << endl;
				cout << mu_beta << endl;
				cout << "----------invsigma_beta----------" << endl;
				cout << invsigma_beta << endl;
			}*/
			vector<int> beta_idx(M,0);
			for (int k = 0; k < M; k++)
				beta_idx[k] = k;
			//random_shuffle(beta_idx.begin(),beta_idx.end());
			double p = 0;
			for (int k = 0; k < 10; k++)
			{
				//cout << "trnormrnd\n";
				for (int ii = 0; ii < M; ii++)
				{
					/*if (k == 0 and infoidx == 5 and verbose)
					{
						cout << "s_beta " << s_beta << endl;
						cout << "---------------------" << endl;
					}*/
					int idx = beta_idx[ii];
					//idx = ii;
					VectorXd s_betatemp = s_beta;
					s_betatemp(idx) = 0;
					double vi2 = 1/invsigma_beta(idx,idx);
					VectorXd btemp,ctemp;
					btemp = mu_beta - s_beta;
					btemp(idx) = 0;
					mattoVecR(invsigma_beta,ctemp,idx);
					ctemp(idx) = 0;
					double ui = mu_beta(idx) + vi2 * btemp.dot(ctemp);
					double mu_th = get_muth(p, idx, s_beta);
					//sample truncated normal distribution
					/*if (k == 0 and infoidx == 5 and verbose)
					{
						cout << "mu_beta(idx) " << mu_beta(idx) << endl;
						cout << "---------------------" << endl;
						cout << btemp << endl;
						cout << "---------------------" << endl;
						cout << ctemp << endl;
						cout << "---------------------" << endl;
						cout << "s_beta(idx) " << s_beta(idx) << " mu: " << ui << " vi2: " << vi2 << " muth: " << mu_th << endl;
					}*/			
					double betai = sample.trnormrnd(s_beta(idx),ui,vi2,mu_th);
					s_beta(idx) = betai;
				}
			}
			samples_beta[i] = s_beta;
			//cout << "Finish beta sampling " << i << endl;
			
			//cout << "s_beta\n";
			//print_v(s_beta);
			
			//cout << "prob_z\n";
			MatrixXf p_all(M,2);
			for (int m = 0; m < M; m++)
			{
				if (s_beta(m) > pow(10,-4))
				{
					//double s_p0 = log(1-s_w) + log(normcdf(s_beta[m],0,sqrt(v0*s_tau[m])) - normcdf(s_beta[m]-pow(10,-3),0,sqrt(v0*s_tau[m])));
					//double s_p1 = log(s_w) + log(normcdf(s_beta[m],0,sqrt(s_tau[m])) - normcdf(s_beta[m]-pow(10,-3),0,sqrt(s_tau[m])));
					//cout << "beta: " << s_beta(m) << " tau: " << s_tau(m) << endl;

					double s_p0 = log(1-s_w) + lognormpdf(s_beta(m),0,sqrt(v0*s_tau(m)));
					if (isinf(s_p0) or s_p0 < -pow(10,6))
						s_p0 = -pow(10,6);
					double s_p1 = log(s_w) + lognormpdf(s_beta(m),0,sqrt(s_tau(m)));
					if (isinf(s_p1) or s_p1 > pow(10,6))
						s_p1 = pow(10,6);
					p_all(m,0) = lognormpdf(s_beta(m),0,sqrt(v0*s_tau(m)));
					p_all(m,1) = lognormpdf(s_beta(m),0,sqrt(s_tau(m)));
					//cout << "p0: " << lognormpdf(s_beta(m),0,sqrt(v0*s_tau(m))) << " p1: " << lognormpdf(s_beta(m),0,sqrt(s_tau(m))) << endl;
					//int aa;
					//cin >> aa;
					double pw = exp(s_p0-s_p1);
					if (pw > 10000)
						pw = 10000;
					double p0 = pw / (pw+1);
					double p1 = 1 - p0;
					if (p0 < 0.05)
						p0 = 0.05;
					double randnum = (double) rand() / (RAND_MAX);
					//cout << m << " " << p0 << " " << randnum << " " << s_p0 << " " << s_p1 << endl;
					if (randnum < p0)
						s_z(m) = v0;
					else
						s_z(m) = 1;
				}
				else
					s_z(m) = v0;
			}
			//samples_z[i] = s_z;
			
			//cout << "s_z\n";
			//print_v(s_z);
			
			//cout << "Finish z sampling calculation " << i << endl;

			s_sigma2 = 1e-5; // implement sigma sampling later
	
			for (int m = 0; m < M; m++)
			{
				//cout << "beta: " << s_beta(m) << " " << s_z(m) << endl;
				double sampled_tau = sample.invgamrnd(a1+0.5,a2+pow(s_beta(m),2)/2/s_z(m));
				if (sampled_tau > 10)
					sampled_tau = 10;
				s_tau(m) = sampled_tau;
				s_gamma(m) = s_tau(m) * s_z(m);
			}

			vector<int> current_idx1;
			for (int m = 0; m < M; m++)
			{
				if (s_z(m) == 1)
					current_idx1.push_back(m);
			}
			//cout << "Adding new isoforms" << endl;
			vector<int> new_idx1;
			/*if (infoidx == 5 and 1)
			{
				cout << "=====================" << endl;
				
				cout << "----------beta----------" << endl;
				cout << s_beta << endl;
				cout << "----------tau----------" << endl;
				cout << s_tau << endl;
				cout << "----------p----------" << endl;
				cout << p_all << endl;
				cout << "----------z----------" << endl;
				cout << s_z << endl;
				cout << "--------------------" << endl;
			}*/
			
			getIsoIdxForX(current_idx1,new_idx1);
			for (int m = 0; m < new_idx1.size(); m++)
			{
				s_z(new_idx1[m]) = 1;
				s_tau(new_idx1[m]) = 0.5;
				s_gamma(new_idx1[m]) = 0.5;
			}
			/*if (infoidx == 5 and 1)
			{
				cout << "----------zmod----------" << endl;
				cout << s_z << endl;
				cout << "=====================" << endl;
				//int aa;
				//cin >> aa;
			}*/
			//cout << "Finish adding new isoforms" << endl;
			samples_z[i] = s_z;
			//s_w = 1e-3; // implement sigma sampling later

		}
		
		//cout << "Finished sampling" << endl;
		
		int burnin = iter - 200;
		vector<double> final_z(M,0);
		vector<double> final_beta(M,0);
		beta_frac_low_high = MatrixXd::Zero(M,2);
		set<double> beta_helper;
		for (int ii = 0; ii < M; ii++)
		{
			beta_helper.clear();
			double sumz = 0;
			for (int i = burnin; i < iter; i++)
			{
				beta_helper.insert(samples_beta[i](ii));
				sumz+= samples_z[i](ii);
			}
			if (beta_helper.size()%2 == 1)
			{
				set<double>::iterator it = beta_helper.begin();
				std::advance(it,(beta_helper.size()-1)/2);
				final_beta[ii] = *it;
			}
			else
			{
				set<double>::iterator it = beta_helper.begin();
				std::advance(it,(beta_helper.size())/2 - 1);
				double beta1 = *it;
				std::advance(it,1);
				double beta2 = *it;
				final_beta[ii] = (beta1+beta2)/2;
			}
			set<double>::iterator it = beta_helper.begin();
			std::advance(it,round(0.025*beta_helper.size()));
			beta_frac_low_high(ii,0) = (*it)/final_beta[ii];
			it = beta_helper.begin();
			std::advance(it,beta_helper.size()-1-round(0.025*beta_helper.size()));
			beta_frac_low_high(ii,1) = (*it)/final_beta[ii];
			final_z[ii] = (double)sumz / (double) (iter-burnin);
		}
		VectorXd final_beta_v = VectorXd::Map(&final_beta[0],final_beta.size());
		VectorXd FPKM, s_reads;
		calFPKM(final_beta_v,FPKM,s_reads);
		for (int ii = 0; ii < FPKM.size(); ii++)
		{
			bool select = final_z[ii] >= 0.5 and s_reads(ii) > 12;
			//bool select = final_z[ii] >= 0.5 and s_beta(ii) > 0.02 and s_reads(ii) > 12;
			//bool select = final_z[ii] >= 0.5;
			if (select)
				final_isoidx.push_back(1);
			else
				final_isoidx.push_back(0);
			final_FPKM.push_back(FPKM(ii));
		}
		/*
		cout << "----------------" << endl;
		cout << "final_z:" << endl;
		//cout << X << endl;
		print_v(final_z);
		cout << "----------------" << endl;
		cout << "final_beta:" << endl;
		print_v(final_beta);
		*/
		//int aa;
		//cin >> aa;
	}

}

void Info::preprocess()
{
	vector<int> X_label;
	vector<double> y1;
	vector<vector<double> > X1;
	vector<double> L1;
	vector<double> R1;
	vector<double> lwv;
	X_label.resize(X.rows(),0);

	for (int i = 0; i < X.rows(); i++)
	{
		double count = 0;
		for (int j = 0; j < X.cols(); j++)
		{
			count+=X(i,j);
		}
		if (count == 0)
			X_label[i] = 1;
	}
	for (int i = 0; i < X.rows(); i++)
	{
		if (X_label[i] == 0)			
		{
			X_label[i] = 1;
			double sumR = R(i);
			double sumL = L(i);
			for (int j = i+1; j < X.rows(); j++)
			{
				if (X_label[j] == 0)
				{
					if (vector_matchX(X,i,j))
					{
						X_label[j] = 1;
						sumR+=R(j);
						sumL+=L(j);
					}
					
				}
			}
			vector<double> s_X1(X.cols(),0);
			for (int j = 0; j < s_X1.size();j++)
				s_X1[j] = X(i,j);
			//print_v(s_X1);
			//int aa;
			//cin >> aa;
			X1.push_back(s_X1);
			L1.push_back(sumL);
			R1.push_back(sumR);
			y1.push_back(sumR/sumL);
			lwv.push_back(sumL/100);
		}
		
	}
	
	//X = MatrixXd::Map(&X1[0][0],X1.size(),X1[0].size());
	vectorcptoMat(X1,X);
	y = VectorXd::Map(&y1[0],y1.size());
	R = VectorXd::Map(&R1[0],R1.size());
	L = VectorXd::Map(&L1[0],L1.size());
	lw = VectorXd::Map(&lwv[0],lwv.size());
	X_m = MatrixXd::Zero(X.rows(),X.cols());
	y_m = VectorXd::Zero(y.size());
	sumX_m = VectorXd::Zero(X.rows());
	//cout << "mapping 1" << endl;
	for (int i = 0; i < X.rows(); i++)
	{
		y_m(i) = lw(i) * y(i);
		y(i) = lw(i) * y(i);//////////////////
		for (int j = 0; j < X.cols(); j++)
		{
			X_m(i,j) = lw(i) * X(i,j);
			X(i,j) = lw(i) * X(i,j); //////////////////////////////
			sumX_m(i) += X_m(i,j);
		}
	}
	//cout << "mapping 2" << endl;
	//get transL
	transL = VectorXd::Zero(X.cols());
	for (int i = 0; i < X.rows(); i++)
	{
		for (int j = 0; j < X.cols(); j++)
		{
			if (X(i,j) > 0)
				transL(j) += L(i);
		}
	}

	//cout << "mapping 3" << endl;
}

void Info::getIsoIdxForX(vector<int> &current_idx, vector<int> &new_idx)
{
	vector<int> current_exon_idx(y.size(),0);
	for (int i = 0; i < current_idx.size(); i++)
	{
		for (int j = 0; j < X.rows(); j++)
		{
			if (X(j,current_idx[i]) > 0)
				current_exon_idx[j] = 1;
		}
	}
	
	set<int> slt_exon_idx;
	for (int i = 0; i < y.size(); i++)
	{
		if (current_exon_idx[i] == 0 and sumX_m(i) > 0)
			slt_exon_idx.insert(i);
	}
	set<int> slt_exon_idx_temp;
	for (int i = 0; i < X.cols(); i++)
		slt_exon_idx_temp.insert(i);
	while (!slt_exon_idx.empty())
	{
		vector<double> sumX_slt;
		matSum(X_m,slt_exon_idx,slt_exon_idx_temp,1,sumX_slt);
		//cout << "----------" << endl;
		//print_v(sumX_slt);
		//cout << "----------" << endl;
		int maxidx = maxidxV(sumX_slt);
		new_idx.push_back(maxidx);
		for (int i = 0; i < y.size(); i++)
		{
			if (X(i,maxidx) > 0)
				slt_exon_idx.erase(i);
		}
		//cout << slt_exon_idx.size() << " " << slt_exon_idx_temp.size() << " " << y.size() << " "<< maxidx << endl;
	}


}

void Info::modX_exon(vector<vector<int> > &paths)
{
	//cout << "Start modifying X_exon" << endl;
	assert(paths.size()>0);
	X_exon = MatrixXd::Zero(paths[0].size(),paths.size());
	for (int i = 0; i < paths.size(); i++)
	{
		for (int j = 0; j < paths[i].size(); j++)
		{
			X_exon(j,i) = paths[i][j];
		}
	}
	//cout << "Finish modifying X_exon" << endl;
}

void Info::modX(vector<vector<int> > &paths)
{
	//cout << "Start modifying X" << endl;
	assert(paths.size()>0);
	if (paths.size()>1)
		single = false;
	X = MatrixXd::Zero(paths[0].size(),paths.size());
	for (int i = 0; i < paths.size(); i++)
	{
		for (int j = 0; j < paths[i].size(); j++)
		{
			//cout << paths[i][j] << " ";
			X(j,i) = paths[i][j];
			if (paths[i][j] == 0)
				X(j,i) = 0;
			//cout << X(j,i) << " ";
		}
	}
	//cout << X << endl;
	//cout << "-----------------------" << endl;
}

void Info::mody(vector<double> &s_y)
{
	y = VectorXd::Map(&s_y[0],s_y.size());
}

/*void Info::addy(double s_y)
{
	y.push_back(s_y);
}*/

void Info::modR(vector<double> &s_R)
{
	R = VectorXd::Map(&s_R[0],s_R.size());
}

/*void Info::addR(double s_R)
{
	R.push_back(s_R);
}*/

void Info::modL(vector<double> &s_L)
{
	L = VectorXd::Map(&s_L[0],s_L.size());
}

/*void Info::addL(double s_L)
{
	L.push_back(s_L);
}*/

void Info::clear()
{
	X.resize(0,0);
	X_m.resize(0,0); 
	sumX_m.resize(0); 
	X_exon.resize(0,0);
	y.resize(0);
	R.resize(0);
	L.resize(0);
	bias.resize(0);
	transL.resize(0);
	exonL.resize(0);
	lw.resize(0);
	label = "";
	single = true;
}

void Info::clear_info()
{
	X.resize(0,0);
	X_m.resize(0,0); 
	sumX_m.resize(0); 
	y.resize(0);
	R.resize(0);
	L.resize(0);
	bias.resize(0);
	transL.resize(0);
	exonL.resize(0);
	lw.resize(0);
	single = true;
}

void Info::modexonbound(vector<vector<int> > &s_bound)
{
	exonbound = s_bound;
}

void Info::modbias(vector<double> &s_bias)
{
	bias = VectorXd::Map(&s_bias[0],s_bias.size());
}

void print_v(vector<double> v)
{
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

void print_vi(vector<int> v)
{
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

void Info::correct_bias()
{

}

void print_vv(vector<vector<double> > v)
{
	for (int i = 0; i < v.size(); i++)
		print_v(v[i]);
}

bool vector_matchX(MatrixXd X, int a, int b)
{
	for (int i = 0; i < X.cols(); i++)
	{
		if (X(a,i) != X(b,i))
			return false;
	}
	return true;
}

void vectorcp(vector<double> &a, vector<double> &b)
{
	b.resize(a.size(),0);
	for (int i = 0; i < a.size(); i++)
		b[i] = a[i];
}
void matcp(vector<vector<double> > &a, vector<vector<double> > &b)
{
	b.resize(a.size(),vector<double> (a[0].size(),0));
	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[i].size(); j++)
		{
			b[i][j] = a[i][j];
		}
	}
}

void vectorcptoMat(vector<vector<double> > &a, MatrixXd &b)
{
	b = MatrixXd::Zero(a.size(),a[0].size());
	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[0].size(); j++)
		{
			b(i,j) = a[i][j];
		}
	}
}

void Info::calFPKM(VectorXd s_beta, VectorXd &FPKM, VectorXd &s_reads)
{
	s_reads = VectorXd::Zero(s_beta.size());
	for (int i = 0; i < R.size(); i++)
	{
		double sum_fpkm = 0;
		for (int j = 0; j < s_reads.size(); j++)
		{
			if (X(i,j) > 0)
			{
				sum_fpkm += s_beta(j);
			}
		}

		for (int j = 0; j < s_reads.size(); j++)
		{
			if (X(i,j) > 0)
			{
				s_reads(j) += s_beta(j) / sum_fpkm * R(i);
			}
		}
	}
	FPKM = s_reads.array() / transL.array() * pow(10,9) / totalReads;
}

void Info::write(string outdir, int gene_idx, bool output_all_bool)
{
	stringstream output_all, output;
	output << outdir << "/SparseIso.gtf";
	output_all << outdir << "/SparseIso_enumerate.gtf";
	string output_str, output_all_str;
	output_str = output.str();
	output_all_str = output_all.str();
	//ofstream outfile, outfile_all;
	//outfile.open(output_str.c_str());
	//outfile_all.open(output_all_str.c_str());

	FILE *outfile;
	FILE *outfile_all;
	if (gene_idx == 1)
	{
		outfile = fopen(output_str.c_str(),"w");
		if (output_all_bool)
			outfile_all = fopen(output_all_str.c_str(),"w");
	}
	else
	{
		outfile = fopen(output_str.c_str(),"a+");
		if (output_all_bool)
			outfile_all = fopen(output_all_str.c_str(),"a+");
	}
	if (false) // or single
	{

	}
	else
	{
		int idx = label.find_first_of("\t");
		string chr = label.substr(0,idx);
		int idx1 = label.find_last_of("\t");
		string strand = label.substr(idx1+1);
		stringstream outgene_stream;
		outgene_stream << "CUFF." << gene_idx;
		string outgene = outgene_stream.str();

		for (int i = 0; i < final_isoidx.size(); i++)
		{
			//cout << "final_isoidx.size(): " << final_isoidx.size() << " X: " << X.cols() << endl;

			set<set<int> > exon_region;
			//cout << "check 1" << endl;
			getExonRegion(i,exon_region);
			//cout << "check 2" << endl;
			set<set<int> >::iterator it = exon_region.begin();
			set<int>::iterator iti = it->begin();
			int e_start = *iti;
			it = exon_region.end();
			std::advance(it,-1);
			iti = it->begin();
			std::advance(iti,1);
			int e_end = *iti;
			stringstream trans_stream;
			trans_stream << outgene << "." << i;
			string outtrans = trans_stream.str();
			/*if (output_all_bool)
				fprintf(outfile_all, "%s\tSparseIso\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; frac \"%.6f\"; conf_lo \"%.6f\"; conf_hi \"%.6f\"; cov \"20\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i], 1.0, final_FPKM[i] * beta_frac_low_high(i,0), final_FPKM[i] * beta_frac_low_high(i,1));
			if (final_isoidx[i] == 1)
				fprintf(outfile, "%s\tSparseIso\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; frac \"%.6f\"; conf_lo \"%.6f\"; conf_hi \"%.6f\"; cov \"20\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i], 1.0, final_FPKM[i] * beta_frac_low_high(i,0), final_FPKM[i] * beta_frac_low_high(i,1));
			*/
			bool writevalid = true;
			if (exon_region.size()==1)
				if (final_FPKM[i] < 5)
					writevalid = false;

			if (output_all_bool)
				fprintf(outfile_all, "%s\tSparseIso\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; \n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i]);
			if (final_isoidx[i] == 1 and writevalid)
				fprintf(outfile, "%s\tSparseIso\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i]);


			//cout << "check 3" << endl;
			int exonNum = 1;
			for (it = exon_region.begin(); it != exon_region.end(); it++)
			{
				iti = it->begin();
				int exon_start = *iti;
				iti = it->begin();
				std::advance(iti,1);
				int exon_end = *iti;
				/*if (output_all_bool)
					fprintf(outfile_all, "%s\tSparseIso\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\"; frac \"%.6f\"; conf_lo \"%.6f\"; conf_hi \"%.6f\"; cov \"20\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i], 1.0, final_FPKM[i] * beta_frac_low_high(i,0), final_FPKM[i] * beta_frac_low_high(i,1));
				if (final_isoidx[i] == 1)
					fprintf(outfile, "%s\tSparseIso\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\"; frac \"%.6f\"; conf_lo \"%.6f\"; conf_hi \"%.6f\"; cov \"20\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i], 1.0, final_FPKM[i] * beta_frac_low_high(i,0), final_FPKM[i] * beta_frac_low_high(i,1));
				*/
				if (output_all_bool)
					fprintf(outfile_all, "%s\tSparseIso\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i]);
				if (final_isoidx[i] == 1 and writevalid)
					fprintf(outfile, "%s\tSparseIso\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i]);

				exonNum++;
			}
			//cout << "check 4" << endl;
			//fprintf(outfile_all, '%s\tSparseIso\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id "%s"; transcript_id "%s"; FPKM "%.6f"; frac "%.6f"; conf_lo "900"; conf_hi "1100"; cov "20";\n',chr.c_str(), e_start, e_end, strand.c_str(); ,outgene, outtrans, FPKM_temp(idx_ident(j)), frac_temp(idx_ident(j)));
		}
	}
	if (output_all_bool)
		fclose(outfile_all);
	fclose(outfile);
}

void Info::getExonRegion(int isoidx, set<set<int> >& exon_region)
{
	set<set<int> > exon_temp;
	//cout << "check 1.1" << endl;
	for (int i = 0; i < X_exon.rows(); i++)
	{
		//cout << "i: " << i << " " << X_exon.rows() << " isoidx: " << isoidx << " " << X_exon.cols() << endl;
		if (X_exon(i,isoidx) > 0)
		{
			set<int> s_exon_temp;
			s_exon_temp.insert(exonbound[i][0]);
			s_exon_temp.insert(exonbound[i][1]);
			exon_temp.insert(s_exon_temp);
		}
	}
	connect_region_recur(exon_temp, exon_region);


}

void connect_region_recur(set<set<int> > region_in, set<set<int> >& region_out)
{
	int quit = 0;
	set<set<int> > temp = region_in;
	while (quit == 0)
	{
		connect_region(temp,region_out);
		if (region_out.size() == temp.size())
			quit = 1;
		else
			temp = region_out;
	}
}


void connect_region(set<set<int> > region_in, set<set<int> >& region_out)
{
	int count = 1, jump = 0;
	vector<set<int> > region_in_v;
	for (set<set<int> >::iterator it = region_in.begin(); it != region_in.end(); it++)
		region_in_v.push_back(*it);
	set<int> zeros;
	zeros.insert(0);
	zeros.insert(0);
	region_in_v.push_back(zeros);
	for (int i = 0; i < region_in_v.size()-1; i++)
	{
		set<int>::iterator iti = region_in_v[i].begin();
		std::advance(iti,1);
		int e_end = *iti;
		int e_start = 0;
		
		set<int>::iterator it1i = region_in_v[i+1].begin();
		e_start = *it1i;
		
		if (jump == 0)
		{
			if (e_start - e_end == 1)
			{
				iti = region_in_v[i].begin();
				int temp1 = *iti;
				std::advance(iti,1);
				int temp2 = *iti;
				iti = region_in_v[i+1].begin();
				if (temp1 > *iti)
					temp1 = *iti;
				std::advance(iti,1);
				if (temp2 < *iti)
					temp2 = *iti;

				set<int> temp;
				temp.insert(temp1);
				temp.insert(temp2);
				region_out.insert(temp);
				region_in_v[i+1] = temp;
				jump = 1;
			}
			else
			{
				region_out.insert(region_in_v[i]);
			}
		}
		else
			jump = 0;
	}
}
double Info::get_muth(double p, int idx, VectorXd s_beta)
{
	MatrixXd mu_tempM = X * s_beta;
	VectorXd mu_temp;
	mattoVecC(mu_tempM,mu_temp,0);
	double min_v = 0;
	for (int i = 0; i < mu_temp.size(); i++)
	{
		if (X(i,idx) > 0)
		{
			if (min_v > mu_temp(i)-p)
				min_v = mu_temp(i)-p;
		}
	}
	return -min_v;
}

void Info::run_single()
{
	double reads = R.sum();
	double length = L.sum();
	if (reads > 12)
	{
		output_bool = true;
		final_FPKM.push_back(reads / length * pow(10,9) / totalReads);
		final_isoidx.push_back(1);
		beta_frac_low_high = MatrixXd::Zero(1,2);
		beta_frac_low_high(0,0) = 0.5;
		beta_frac_low_high(0,1) = 1.5;
	}
	else
	{
		output_bool = false;
		final_isoidx.push_back(0);
	}

}