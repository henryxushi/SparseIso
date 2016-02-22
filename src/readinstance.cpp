#include "readinstance.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <list>
#include <math.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <set>
#include <map>
#include <limits.h>


void readinstance::readinstance_p(string inputfile, double frag_mu, double frag_sigma)
{
	//cout << "Processing\t" << inputfile << "\toutput\t" << outputfile << endl;
	ifstream infile1;
	infile1.open(inputfile.c_str());
	
	string line;
	Instance * inst = new Instance();
	Info info;
	inst->frag_mu = frag_mu; // set the mean and std of the fragment length distribution
	inst->frag_sigma = frag_sigma;
	inst->maxExonCov = -1;
	int instcount = 0;
	int outexon = 0, outSGT = 0, outPET = 0;
	vector<vector<int> > ExonBound;
	vector<double> exonAbundance;
	vector<vector<int> > SegConf;
	vector<int> SegAbundance;
	int totalNumOfReads = 0;
	int sumReadstemp = 0;
	int totalNumOfPaths = 0;
	int totalNumOfGenes = 0;
	int countlinetest = 0;
	if (infile1.is_open())
	{
		while (infile1.good())
		{
			getline(infile1,line);
			//cout << "line: " << ++countlinetest << endl;
			int idx; int a;
			idx = line.find_first_of("\t");
			string substring = line.substr(0,idx);
			vector<int> singleExonBound;
			vector<int> singleSegConf;
			vector<vector<int> > singlePETreads1;
			vector<int> singlePETreads;
			vector<int> singlePETconn;
			
			if (!substring.compare("Instance") || line.empty())
			{
				if (instcount > 0)
				{
					//process....
					if (instcount % 1 == 0)
						cout << "process\t" << instcount << "\t" << inst->location << endl;//"\t" << inst->NumofSegs << endl; 

					if (!inst->strand.compare("-"))
					{
						std::reverse(inst->exonbound.begin(), inst->exonbound.end());
						std::reverse(inst->exonAbundance.begin(), inst->exonAbundance.end());
						for (int i = 0; i < inst->SegConf.size(); i++)
							std::reverse(inst->SegConf[i].begin(), inst->SegConf[i].end());
						std::reverse(inst->SegConf.begin(), inst->SegConf.end());
						std::reverse(inst->SegAbundance.begin(), inst->SegAbundance.end());
						std::reverse(inst->connSeg.begin(), inst->connSeg.end());
						for (int i = 0; i < inst->PETconnect.size(); i++)
						{
							int temp = 0, temp1 = 0;
							temp = inst->NumofSegs - inst->PETconnect[i][0] + 1;
							temp1 = inst->NumofSegs - inst->PETconnect[i][1] + 1;
							inst->PETconnect[i][0] = temp1;
							inst->PETconnect[i][1] = temp;
						}
					}
					//if (inst->NumofSegs < 200 && inst->strand.compare("."))
					if (inst->strand.compare("."))
					{
						cout << "Start processing" << endl;
						if (opt.readtype == "p")
							inst->process_inst(info);
						else if (opt.readtype == "s")
							inst->process_inst_se(info);
						else
							exit(-1);
						if (info.valid)
							infolist.push_back(info);
						//cout << infolist[0].X << endl;
						//cout << "---------------------------" << endl;
						//process_inst(inst,outputfile, totalNumOfPaths);
						cout << "Current Num of Genes: " << ++totalNumOfGenes << endl;
						totalNumOfReads += sumReadstemp;
						sumReadstemp = 0;
					}
					else
					{
						cout << "pass\t" << inst->NumofSegs << endl;
					}
					
					//cout << "fail check1" << endl;
					for (int i = 0; i < inst->exonbound.size(); i++)
					{
						vector<int>().swap(inst->exonbound[i]);
					}
					for (int i = 0; i < inst->SegConf.size(); i++)
					{
						vector<int>().swap(inst->SegConf[i]);
						inst->SegConf[i].clear();
					}
					vector<vector<int> >().swap(inst->exonbound);
					vector<vector<int> >().swap(inst->SegConf);
					vector<int>().swap(inst->SegAbundance);
					vector<int>().swap(inst->connSeg);
					vector<double>().swap(inst->exonAbundance);
					inst->exonbound.clear();
					inst->SegConf.clear();
					inst->SegAbundance.clear();
					inst->connSeg.clear();
					inst->exonAbundance.clear();
					delete inst;
					inst = NULL;
					inst = new Instance();
					Info info;
					inst->frag_mu = frag_mu;
					inst->frag_sigma = frag_sigma;
					outexon = 0; outexon = 1;
					for (int i = 0; i < ExonBound.size(); i++)
					{
						vector<int>().swap(ExonBound[i]);
					}
					for (int i = 0; i < SegConf.size(); i++)
					{
						vector<int>().swap(SegConf[i]);
					}
					vector<vector<int> >().swap(ExonBound);
					ExonBound.clear();
					vector<vector<int> >().swap(SegConf);
					SegConf.clear();
					vector<int>().swap(SegAbundance);
					SegAbundance.clear();
					vector<double>().swap(exonAbundance);
					exonAbundance.clear();
					//if (instcount > 16897)
					//{
						//int aa;
						//cin >> aa;
					//}
					cout << "Finished instance" << endl;
				}
				instcount = instcount + 1;
			}
			else if (!substring.compare("Boundary"))
			{
				inst->location = line.substr(idx+1);
				int idx1 = line.find_last_of("\t");
				inst->strand = line.substr(idx1+1);
			}
			else if (!substring.compare("ReadLen"))
			{
				inst->ReadLen = atoi(line.substr(idx+1).c_str());
			}
			else if (!substring.compare("Segs"))
			{
				inst->NumofExon = atoi(line.substr(idx+1).c_str());
				outexon = 1;
			}
			else if (!substring.compare("Reads"))
			{
				sumReadstemp = atoi(line.substr(idx+1).c_str());
				TOTALNUMREADS += sumReadstemp;
			}
			else if (outexon == 1)
			{
				if (!substring.compare("Refs"))
				{
					inst->exonbound = ExonBound;
					inst->exonAbundance = exonAbundance;
					outexon = 0;
				}
				else
				{
					int idx1 = line.find_last_of("\t");
					int idx2 = line.find_last_of("\t",idx1-1);

					//cout << atoi(line.substr(0,idx).c_str()) << "\t" << atoi(line.substr(idx+1).c_str()) << endl;
					//cout << line << endl << "aaaa\t" << atof(line.substr(idx1+1).c_str()) << endl;
					double exon_abun1 = (double)(atof(line.substr(idx1).c_str()));
					double exon_abun = (double)(atof(line.substr(idx2,idx1-idx2).c_str()));
					if (exon_abun1 > inst->maxExonCov)
					{
						inst->maxExonCov = exon_abun1;
					}
					exonAbundance.push_back(exon_abun);
					singleExonBound.push_back(atoi(line.substr(0,idx).c_str()));
					singleExonBound.push_back(atoi(line.substr(idx+1).c_str()));
					ExonBound.push_back(singleExonBound);
					singleExonBound.clear();
				}

			}
			else if (!substring.compare("SGTypes"))
			{
				inst->NumofSegs = atoi(line.substr(idx+1).c_str());
				outSGT = 1;
			}
			else if (outPET == 1)
			{
				if (!substring.compare("Coverage"))
				{
					//inst->PETconnect = PETconnect;
					//inst->PETreads = PETreads;
					outPET = 0;
				}
				else
				{
					int idx1 = line.find_first_of("\t");
					singlePETconn.push_back(atoi(line.substr(0,idx1).c_str()));
					int idx2 = line.find_first_of("\t",idx1+1);
					singlePETconn.push_back(atoi(line.substr(idx1+1,idx2-idx1).c_str()));
					idx1 = line.find_first_of("\t",idx2+1);
					singlePETconn.push_back(atoi(line.substr(idx2+1,idx1-idx2).c_str()));
					inst->PETconnect.push_back(singlePETconn);
					singlePETconn.clear();
					//cout << singlePETconn[0] << " " << singlePETconn[1] << " " << singlePETconn[2] << endl;

					getline(infile1,line);
					//cout << "line: " << ++countlinetest << endl;
					idx1 = line.find_first_of(" ");
					int preidx = -1;
					int count_exon = 0;
					while (idx1 != -1)
					{
						string seg = line.substr(preidx+1,idx1-preidx);
						idx2 = seg.find_first_of(":");
						singlePETreads.push_back(atoi(seg.substr(0,idx2).c_str()));
						singlePETreads.push_back(atoi(seg.substr(idx2+1).c_str()));
						//cout << singlePETreads[0] << " " << singlePETreads[1] << endl;
						singlePETreads1.push_back(singlePETreads);
						singlePETreads.clear();
						preidx = idx1;
						idx1 = line.find_first_of(" ",preidx+1);

					}
					//int aa; cin >> aa;
					inst->PETreads.push_back(singlePETreads1);
					singlePETreads1.clear();
				}
			}
			else if (outSGT == 1)
			{
				if (!substring.compare("PETypes"))
				{
					inst->SegConf = SegConf;
					inst->SegAbundance = SegAbundance;
					outSGT = 0;
					outPET = 1;
				}
				else
				{
					int idx_last = line.find_last_of(" ",idx+1);
					SegAbundance.push_back(atoi(line.substr(idx_last+1,idx-idx_last).c_str()));
					string confstring = substring;
					//cout << line << endl;
					//cout << line.substr(idx_last+1,idx-idx_last).c_str() << endl;
					int idx1 = confstring.find_first_of(" ");
					int preidx = -1;
					int count_exon = 0;
					while (idx1 != -1)
					{
						//cout << confstring.substr(preidx+1,idx1-preidx) << "\t";
						singleSegConf.push_back(atoi(confstring.substr(preidx+1,idx1-preidx).c_str()));
						count_exon += atoi(confstring.substr(preidx+1,idx1-preidx).c_str());
						preidx = idx1;
						idx1 = confstring.find_first_of(" ",idx1+1);
					}
					SegConf.push_back(singleSegConf);
					inst->connSeg.push_back(count_exon);
				}
			}


		}
		
					
	}
}
