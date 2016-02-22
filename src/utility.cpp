#include "utility.h"
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
using namespace std;

bool ifcontain_vec(vector<int> &a, vector<int> &b);
bool ifcontain_arr(int a[], vector<int> &b);
double normpdf(double x, double mu, double sigma); // not exactly the probability, it is relative one
double phi(double x);
double phi_range(double x, double y);
double effectL(int l1, int l2, int l3, int r, double mu, double sigma, int bias);

void Instance::process_inst(Info& info)
{
	//enumberate path
	info.clear();
	//cout << "Finish clearing info" << endl;
	if (maxExonCov > -1) //modified 06/15/2015, let low FPKM pass, filter in analysis part
	{
		int output = 0;
		remove_low(); //modified 06/15/2015 rely on prunning to remove false segments

		vector<double> proc_abun;
		map<int, vector<int> > iscontained;
		findcontain(iscontained);
		//outfile1 << "Instance\t" << location << endl;
		info.label = location;
		//outfile1 << "ExonBound\t" << exonbound.size() << endl;
		//for (int i = 0; i < exonbound.size(); i++)
		//{
		//	outfile1 << exonbound[i][0] << "\t" << exonbound[i][1] << endl;
		//}
		info.modexonbound(exonbound);
		
		double th1 = 0;
		double min_th1 = 0;

		mod_y(iscontained, proc_abun, output, th1, min_th1);

		map<int,set<int> > exonmap; 
		getmap(exonmap, th1); // exonmap minvalue 1
		for (map<int,set<int> >::iterator mapit = exonmap.begin();mapit!=exonmap.end();++mapit)
		{
			//cout << mapit->first << ": ";
			//for (set<int>::iterator setit = mapit->second.begin(); setit != mapit->second.end();setit++)
				//cout << *setit <<" ";
			//cout << endl;
		}

		//if(!exonmap.empty())
		//{
			set<int>::iterator it;
			int not_source_list[NumofExon];
			int not_sink_list[NumofExon];
			for (int i = 0; i < NumofExon; i++)
			{
				not_source_list[i] = 0;
				not_sink_list[i] = 0;
			}

			for (map<int,set<int> >::iterator ii = exonmap.begin(); ii != exonmap.end(); ++ii)
			{
				int key = (*ii).first;
				set<int> temp = (*ii).second;
				if (!temp.empty())
				{
					not_sink_list[key - 1] = 1;
					for (it = temp.begin(); it != temp.end(); ++it)
					{
						not_source_list[*it - 1] = 1; // not a source node
					}
				}
			}
			
			vector<vector<int> > paths;
			double mod_y_th = 0;
			int quit = enumerate_path(exonmap, not_source_list, not_sink_list, paths);
			while(quit == -1)
			{
				for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
				{
					vector<int>().swap(it->second);
					//exonmap.erase(it);
				}

				map<int,vector<int> >().swap(iscontained);
				iscontained.clear();
				findcontain(iscontained);
				th1 = min_th1;
				if (th1 == 10000)
				{
					break;
				}
				//cout << "new threshold\t" << th1 << endl;
				output = 0;
				proc_abun.clear();
				vector<double>().swap(proc_abun);
				mod_y(iscontained, proc_abun, output, th1, min_th1);
				//cout << "current number of segs in inst\t" << SegConf.size() << "\t" << NumofSegs << "\t" << SegAbundance.size() << endl;
				for (map<int,set<int> >::iterator it = exonmap.begin();it!=exonmap.end();++it)
				{
					set<int>().swap(it->second);
					//exonmap.erase(it);
				}

				map<int,set<int> >().swap(exonmap);
				exonmap.clear();
				map<int,set<int> > exonmap; 
				getmap(exonmap, th1); // exonmap minvalue 1
				set<int>::iterator it;
				int not_source_list1[NumofExon];
				int not_sink_list1[NumofExon];
				for (int i = 0; i < NumofExon; i++)
				{
					not_source_list1[i] = 0;
					not_sink_list1[i] = 0;
				}
				for (map<int,set<int> >::iterator ii = exonmap.begin(); ii != exonmap.end(); ++ii)
				{
					int key = (*ii).first;
					set<int> temp = (*ii).second;
					if (!temp.empty())
					{
						////cout << inst->NumofExon << endl;
						////cout << key << endl;
						not_sink_list1[key - 1] = 1;
						for (it = temp.begin(); it != temp.end(); ++it)
						{
							not_source_list1[*it - 1] = 1; // not a source node
						}
					}
				}
				quit = enumerate_path(exonmap, not_source_list1, not_sink_list1, paths);
			}
			//cout << "final paths\t" << paths.size() << endl;
			for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
			{
				vector<int>().swap(it->second);
			}
			map<int,vector<int> >().swap(iscontained);
			iscontained.clear();
			findcontain(iscontained);
			/*output = 1;
			if (paths.size() == 0)
			{
				outfile1 << -1 << "\t" <<-1 << "\t" << -1 << endl;
			}
			else
			{
				proc_abun.clear();
				vector<double>().swap(proc_abun);
				mod_y(iscontained, inst, proc_abun, outfile1, output, th1, min_th1);
			}*/
			
			vector<vector<int> > paths1;
			for (int i = 0; i < paths.size(); i++)
			{
				if (paths[i].size() > 2)
				{
					vector<int> singlepath1;
					int idxpath = 1;
					for (int j = 0; j < NumofExon; j++)
					{
						if (j+1 == paths[i][idxpath])
						{
							singlepath1.push_back(1);
							idxpath++;
						}
						else
							singlepath1.push_back(0);
					}
					paths1.push_back(singlepath1);
				}

			}
			

			if (paths.size() == 0)
			{
				info.valid = false;
				//outfile1 << "Abundance" << endl << -1 << "\t" << -1 <<"\t" << -1 << endl;
				//outfile1 << "Paths" << endl;
				//outfile1 << -1 << endl << -1 << endl;
			}
			else
			{
				correct_y(paths1,info);
				/*for (int i = 0; i < paths.size(); i++)
				{
					int countpath1 = 0;
					for (int ii = 0; ii < paths[i].size(); ii++)
						countpath1++;
					if (countpath1 > 3)
					{
						mod_X(paths[i],inst->SegConf,proc_abun,inst->NumofSegs,inst->NumofExon,outfile1);
						totalNumOfPaths++;
					}
					
				}*/
			}

			//cout << "Current Num of Paths: " << totalNumOfPaths << endl; 

			for (int i = 0; i < paths.size(); i++)
			{
				paths[i].clear();
				vector<int>().swap(paths[i]);
			}
			paths.clear();
			vector<vector<int> >().swap(paths);
			
		//}
		for (map<int,set<int> >::iterator it = exonmap.begin();it!=exonmap.end();++it)
		{
			set<int>().swap(it->second);
			//exonmap.erase(it);
		}

		map<int,set<int> >().swap(exonmap);
		exonmap.clear();

		for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
		{
			vector<int>().swap(it->second);
			//exonmap.erase(it);
		}

		map<int,vector<int> >().swap(iscontained);
		iscontained.clear();
	}
}

void Instance::process_inst_se(Info& info)
{
	//enumberate path
	info.clear();
	//cout << "Finish clearing info" << endl;
	if (maxExonCov > -1) 
	{
		int output = 0;
		//ofstream outfile1;
		//outfile1.open(outputfile.c_str(),std::ios_base::app);
		//remove_low();

		vector<double> proc_abun;
		vector<double> proc_R;
		vector<double> proc_L;
		map<int, vector<int> > iscontained;
		findcontain(iscontained);
		//outfile1 << "Instance\t" << location << endl;
		info.label = location;
		//outfile1 << "ExonBound\t" << exonbound.size() << endl;
		//for (int i = 0; i < exonbound.size(); i++)
		//{
		//	outfile1 << exonbound[i][0] << "\t" << exonbound[i][1] << endl;
		//}
		info.modexonbound(exonbound);
		
		double th1 = 0;
		double min_th1 = 0;

		mod_y_se(iscontained, proc_abun, proc_R, proc_L, output, th1, min_th1);

		map<int,set<int> > exonmap; 
		getmap(exonmap, th1); // exonmap minvalue 1
		//for (map<int,set<int> >::iterator mapit = exonmap.begin();mapit!=exonmap.end();++mapit)
		//{
			//cout << mapit->first << ": ";
			//for (set<int>::iterator setit = mapit->second.begin(); setit != mapit->second.end();setit++)
				//cout << *setit <<" ";
			//cout << endl;
		//}

		//if(!exonmap.empty())
		//{
			set<int>::iterator it;
			int not_source_list[NumofExon];
			int not_sink_list[NumofExon];
			for (int i = 0; i < NumofExon; i++)
			{
				not_source_list[i] = 0;
				not_sink_list[i] = 0;
			}

			for (map<int,set<int> >::iterator ii = exonmap.begin(); ii != exonmap.end(); ++ii)
			{
				int key = (*ii).first;
				set<int> temp = (*ii).second;
				if (!temp.empty())
				{
					not_sink_list[key - 1] = 1;
					for (it = temp.begin(); it != temp.end(); ++it)
					{
						not_source_list[*it - 1] = 1; // not a source node
					}
				}
			}
			
			vector<vector<int> > paths;
			double mod_y_th = 0;
			int quit = enumerate_path(exonmap, not_source_list, not_sink_list, paths);
			while(quit == -1)
			{
				for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
				{
					vector<int>().swap(it->second);
					//exonmap.erase(it);
				}

				map<int,vector<int> >().swap(iscontained);
				iscontained.clear();
				findcontain(iscontained);
				th1 = min_th1;
				if (th1 == 10000)
				{
					break;
				}
				//cout << "new threshold\t" << th1 << endl;
				output = 0;
				proc_abun.clear();
				vector<double>().swap(proc_abun);
				proc_R.clear();
				vector<double>().swap(proc_R);
				proc_L.clear();
				vector<double>().swap(proc_L);
				mod_y_se(iscontained, proc_abun, proc_R, proc_L, output, th1, min_th1);
				//cout << "current number of segs in inst\t" << SegConf.size() << "\t" << NumofSegs << "\t" << SegAbundance.size() << endl;
				for (map<int,set<int> >::iterator it = exonmap.begin();it!=exonmap.end();++it)
				{
					set<int>().swap(it->second);
					//exonmap.erase(it);
				}

				map<int,set<int> >().swap(exonmap);
				exonmap.clear();
				map<int,set<int> > exonmap; 
				getmap(exonmap, th1); // exonmap minvalue 1
				set<int>::iterator it;
				int not_source_list1[NumofExon];
				int not_sink_list1[NumofExon];
				for (int i = 0; i < NumofExon; i++)
				{
					not_source_list1[i] = 0;
					not_sink_list1[i] = 0;
				}
				for (map<int,set<int> >::iterator ii = exonmap.begin(); ii != exonmap.end(); ++ii)
				{
					int key = (*ii).first;
					set<int> temp = (*ii).second;
					if (!temp.empty())
					{
						//cout << inst->NumofExon << endl;
						//cout << key << endl;
						not_sink_list1[key - 1] = 1;
						for (it = temp.begin(); it != temp.end(); ++it)
						{
							not_source_list1[*it - 1] = 1; // not a source node
						}
					}
				}
				quit = enumerate_path(exonmap, not_source_list1, not_sink_list1, paths);
			}
			//cout << "final paths\t" << paths.size() << endl;
			for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
			{
				vector<int>().swap(it->second);
			}
			map<int,vector<int> >().swap(iscontained);
			iscontained.clear();
			findcontain(iscontained);
			
			vector<vector<int> > paths1;
			for (int i = 0; i < paths.size(); i++)
			{
				if (paths[i].size() > 2)
				{
					vector<int> singlepath1;
					int idxpath = 1;
					for (int j = 0; j < NumofExon; j++)
					{
						if (j+1 == paths[i][idxpath])
						{
							singlepath1.push_back(1);
							idxpath++;
						}
						else
							singlepath1.push_back(0);
					}
					paths1.push_back(singlepath1);
				}

			}
			

			if (paths.size() == 0)
			{
				//outfile1 << "Abundance" << endl << -1 << "\t" << -1 <<"\t" << -1 << endl;
				//outfile1 << "Paths" << endl;
				//outfile1 << -1 << endl << -1 << endl;
				info.valid = false;
			}
			else
			{
				info.valid = true;
				info.mody(proc_abun);
				info.modR(proc_R);
				info.modL(proc_L);
				//outfile1 << "Abundance" << endl;
				//for (int i = 0; i < proc_abun.size(); i++)
				//	outfile1 << proc_abun[i] << "\t" << proc_R[i] << "\t" << proc_L[i] << endl;

				//outfile1 << "Paths" << endl;
				vector<vector<int> > segpath_all;
				for (int i = 0; i < paths.size(); i++)
				{
					vector<int> segpath;
					int countpath1 = 0;
					for (int ii = 0; ii < paths[i].size(); ii++)
						countpath1++;
					if (countpath1 > 0) //may skip single exon transcripts
					{
						mod_X(paths[i],proc_abun, segpath);
					}
					if (segpath.size() > 0)
						segpath_all.push_back(segpath);
					
				}
				if (segpath_all.size() == 0)
					info.valid = false;
				else
				{
					info.modX(segpath_all);
					info.modX_exon(paths);
				}
			}
			

			for (int i = 0; i < paths.size(); i++)
			{
				paths[i].clear();
				vector<int>().swap(paths[i]);
			}
			paths.clear();
			vector<vector<int> >().swap(paths);
		//}
		for (map<int,set<int> >::iterator it = exonmap.begin();it!=exonmap.end();++it)
		{
			set<int>().swap(it->second);
			//exonmap.erase(it);
		}

		map<int,set<int> >().swap(exonmap);
		exonmap.clear();

		for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
		{
			vector<int>().swap(it->second);
		}

		map<int,vector<int> >().swap(iscontained);
		iscontained.clear();
	}
}

void Instance::correct_y(vector<vector<int> > &paths, Info &info)
{
	info.valid = false;
	int NofPET = PETconnect.size();
	int intron_l[exonbound.size()-1]; //calculate intron length between exons
	for (int i = 0; i < exonbound.size()-1; i++)
	{
		if (!strand.compare("-"))
			intron_l[i] = exonbound[i][0] - exonbound[i+1][1] - 1;
		else
			intron_l[i] = exonbound[i+1][0] - exonbound[i][1] - 1;
	}
	////cout << "---------------------------------" <<endl;
	////cout << "check connseg" << endl;
	//for (int i = 0; i < inst->connSeg.size(); i++)
		////cout << inst->connSeg[i] << endl;
	//cout << "CatofPETreads " << PETconnect.size() << endl;
	vector<vector<int> > new_paths;
	vector<int> new_Nreads;
	vector<double> new_length;
	vector<int> start_nodes;
	for (int i = 0; i < paths.size(); i++)
	{
		vector<int> temp;
		new_paths.push_back(temp);
	}
	
	for (int i = 0; i < NofPET; i++)
	{
		int s1 = PETconnect[i][0]; // connect between SGTypes
		int s2 = PETconnect[i][1];

		if (true)//s1 != s2)
		{
			//cout << "processing " << s1 << " " << s2 << endl;
			// collect all paths that contain both s1 and s2
			vector<int> containboth;
			vector<int> containbool;
			containboth.clear();
			for (int ii = 0; ii < paths.size(); ii++)
			{
				int* temppath = &paths[ii][0];
				if (ifcontain_arr(temppath, SegConf[s1-1]) && ifcontain_arr(temppath, SegConf[s2-1]))
				{
					containboth.push_back(ii); // contain idx - 1
					containbool.push_back(1);
				}
				else
					containbool.push_back(0);
				temppath = NULL;
					
			}
			//cout << containboth.size() << endl;
			//if the isoform conatin the two segments, test the length from s1 to s2
			set<int> catofIntron;
			vector<int> catofIntron_vec;
			vector<int> effectiveL;
			map<int,vector<double> > intron2effL;
			for (int ii = 0; ii < containboth.size(); ii++)
			{
				vector<int> lengthv;
				int l1 = 0, l2 = 0, l3 = 0, bias = 0;
				int intron_temp = calintronlength(s1,s2,paths[containboth[ii]],l1,l2,l3,0,bias);
				//cout << "Check Length: " << s1 << " " << s2 << endl;
				//cout << l1 << " " << l2 << " " << l3 << " "  << bias << endl;
				if (intron2effL.find(intron_temp) == intron2effL.end())
				{
					//cout << "find intron_length " << intron_temp << endl;
					intron2effL[intron_temp].push_back(l1);
					intron2effL[intron_temp].push_back(l2);
					intron2effL[intron_temp].push_back(l3);
					intron2effL[intron_temp].push_back(bias);
				}
				int ltemp1 = 0;
				effectiveL.push_back(ltemp1);
				catofIntron_vec.push_back(intron_temp);
				catofIntron.insert(intron_temp);
			}

			// calculate effective length for each 
			for (set<int>::iterator it = catofIntron.begin(); it!=catofIntron.end();++it)
			{
				int l1 = intron2effL[*it][0];
				int l2 = intron2effL[*it][1];
				int l3 = intron2effL[*it][2];
				int bias = intron2effL[*it][3];
				intron2effL[*it].clear();
				intron2effL[*it].push_back(effectL(l1, l2, l3, ReadLen, frag_mu, frag_sigma,bias));
			}
			

			int NumofCat = catofIntron.size(); // for this kind of reads, we need to split into NumofCat parts
			
			map<int,int> readsmap; 
		
			//for each read, assign to different categories
			int readsAssigned[NumofCat];
			for (int jj = 0; jj < NumofCat; jj++)
			{
				readsAssigned[jj] = 0;
			}	
			
			for (int j = 0; j < PETreads[i].size(); j++)
			{
				double readlen = (double) PETreads[i][j][0];
				readlen += 2 * ReadLen;
				//assign reads to most probabale one
				double pdf_all = 0;
				double pdfReads[NumofCat];
				for (int jj = 0; jj < NumofCat; jj++)
				{
					pdfReads[jj] = 0;
				}

				int idxjj = 0;
				for (set<int>::iterator it = catofIntron.begin(); it!=catofIntron.end();++it)
				{
					double frag_length = readlen - (double)*it;
					double pdf = normpdf(frag_length, frag_mu, frag_sigma);
					pdf_all += pdf;
					pdfReads[idxjj] = pdf;
					idxjj++;
				}

				vector<int> assigned_node;
				vector<double> assigned_weight;
				double total_weight = 0;

				int total_reads = 0, total_assigned = 0;
				for (int jj = 0; jj < NumofCat; jj++)
				{
					if (pdfReads[jj] / pdf_all > 0.1) // use 0.2 to make the reads assigned to different categories
					{
						assigned_node.push_back(jj);
						assigned_weight.push_back(pdfReads[jj] / pdf_all);
						total_weight += pdfReads[jj] / pdf_all;
						
					}					
				}

				for (int jj = 0; jj < assigned_node.size(); jj++)
				{
					double assigned_reads = PETreads[i][j][1] * assigned_weight[jj] / total_weight;
					readsAssigned[assigned_node[jj]] += assigned_reads; //similar reads can be removed here
					total_assigned += assigned_reads;
					total_reads += assigned_reads;
				}

			}

			// find whether this kind of reads matched to paths
			int idxjj = 0;
			for (set<int>::iterator it = catofIntron.begin(); it!=catofIntron.end();++it)
			{
				//cout << "Assign all: " << readsAssigned[idxjj] << " effectL " <<  intron2effL[*it][0] << endl;
				if(readsAssigned[idxjj] > 0 && intron2effL[*it][0] > 1)
				{
					//cout << "check Here " << *it << endl;
					int idxii = 0; int idxkk = 0;
					for (int ii = 0; ii < paths.size(); ii++)
					{
						if(containbool[ii] == 1)
						{
							//cout << "111\t" << catofIntron_vec[idxkk] << " " << *it << endl;
							if (catofIntron_vec[idxkk] == *it)
							{
								new_paths[ii].push_back(1);
								//cout << "+1" << endl;
								idxii++;
							}
							else
								new_paths[ii].push_back(0);
							idxkk++;
						}
						else
							new_paths[ii].push_back(0);
					}
					new_Nreads.push_back(readsAssigned[idxjj]);
					//cout << "check assigned reads: " << readsAssigned[idxjj] << endl;
					new_length.push_back(intron2effL[*it][0]);
					//cout << "check new length: " << intron2effL[*it][0] << endl;
					//cout << "Nreads: " << readsAssigned[idxjj] << endl;
					start_nodes.push_back(segStart(s1));
					set<int> single_inc_exon;
					segStart(s1,s2,single_inc_exon);
					inc_exon.push_back(single_inc_exon);
					//cout << new_Nreads.size() << "," << new_length.size() << "," << inc_exon.size() << endl;
				}
				idxjj++;
			}
			
			set<int>().swap(catofIntron);
			catofIntron.clear();
			vector<int>().swap(catofIntron_vec);
			catofIntron_vec.clear();
		}
	}



	//cout << "---------------------------------" <<endl;
	//cout << new_paths.size() << "\t" << new_Nreads.size() << endl;
	vector<double> infoy,infoR,infoL;
	//outfile << "Abundance" << endl;
	for (int i = 0; i < new_Nreads.size(); i++)
	{
		//outfile << new_Nreads[i] / new_length[i] << "\t" << new_Nreads[i] << "\t" << new_length[i] << "\t" << start_nodes[i] << "\t";
		//for(set<int>::iterator it = inc_exon[i].begin(); it != inc_exon[i].end(); ++it)
		//{
		//	outfile << *it << ",";
		//}
		//outfile << endl;
		infoy.push_back(new_Nreads[i] / new_length[i]);
		infoR.push_back(new_Nreads[i]);
		infoL.push_back(new_length[i]);
		//cout << new_Nreads[i] << " ";
	}
	info.mody(infoy);
	info.modR(infoR);
	info.modL(infoL);
	//cout << endl;
	//for (int i = 0; i < new_length.size(); i++)
	//{
		//cout << new_length[i] << " ";
	//}
	//cout << endl;

	//for (int i = 0; i < new_length.size(); i++)
	//{
		//cout << new_Nreads[i] / new_length[i] << " ";
	//}
	//cout << endl;

	//outfile << "Paths" << endl;
	//for (int i = 0; i < new_paths.size(); i++)
	//{
		
		//for (int j = 0; j < paths[i].size(); j++)
			//outfile << paths[i][j] << " ";
		//outfile << endl;
		
		//for (int j = 0; j < new_paths[i].size(); j++)
		//	outfile << new_paths[i][j] << " ";
		//outfile << endl;
		
		//for (int j = 0; j < new_paths[i].size(); j++)
			//cout << new_paths[i][j] << " ";
		//cout << endl;
	//}
	info.modX_exon(paths);
	info.modX(new_paths);
	if (paths.size() > 0 and new_paths.size() > 0)
		info.valid = true;
	//cout << "=================================" << endl;

}

void Instance::remove_low()
{
	//remove exons
	//remove here will increase false positive (observed from real data)
	double th = 0.01;
	double abun_th = 15;
	double abun_th_low = 2;
	vector<int> removed_exon;
	//int count = 0;
	for (int i = 0; i < NumofExon; i++)
	{
		if (exonAbundance[i] >= MAX_INTRON_COV)
			removed_exon.push_back(i+1);
	}

	for (int i = 0; i < removed_exon.size(); i++)
	{
		//cout << "remove exon\t" << i + 1 << endl;
		vector<int>().swap(exonbound[removed_exon[i]-1-i]);
		exonbound.erase(exonbound.begin()+removed_exon[i]-1-i);
		exonAbundance.erase(exonAbundance.begin()+removed_exon[i] - 1 -i);
		for (int j = 0; j < NumofSegs; j++)
		{
			SegConf[j].erase(SegConf[j].begin()+removed_exon[i] - 1 - i);
		}
	}
	NumofExon = NumofExon - removed_exon.size();

	set<int> removed_seg;
	set<int> removed_pet;
	for (int i = 0; i < SegConf.size(); i++)
	{
		for (int j = i + 1; j < SegConf.size(); j++)
		{
			int countsame = 0;
			for (int ii = 0; ii < SegConf[i].size(); ii++)
			{
				countsame += abs(SegConf[i][ii] - SegConf[j][ii]);
			}
			if (countsame == 0)
			{
				SegAbundance[i] += SegAbundance[j];
				removed_seg.insert(removed_seg.end(),j+1);
			}
		}
	}
	
	for (int j = 0; j < SegConf.size(); j++)
	{
		int count = 0;
		for (int ii = 0; ii < SegConf[j].size(); ii++)
		{
			if (SegConf[j][ii] == 1)
				count++;
		}
		if (count == 0)
		{
			removed_seg.insert(removed_seg.end(),j+1);
		}

	}
	
	
	int removetemp = 0;
	for(set<int>::iterator it = removed_seg.begin(); it != removed_seg.end(); ++it)
	{	
		for (int i = 0; i < PETconnect.size(); i++)
		{
			int s1 = PETconnect[i][0];
			int s2 = PETconnect[i][1];
			if (s1 == *it || s2 == *it)
				removed_pet.insert(i+1);
		}
		vector<int>().swap(SegConf[*it-1-removetemp]);
		SegConf.erase(SegConf.begin() + *it - 1 - removetemp);
		SegAbundance.erase(SegAbundance.begin() + *it - 1 - removetemp);
		removetemp++;
	}

	
	connSeg.clear();
	connSeg.resize(SegConf.size(),0);
	for (int j = 0; j < SegConf.size(); j++)
	{
		
		for (int ii = 0; ii < SegConf[j].size(); ii++)
		{
			if (SegConf[j][ii] == 1)
			{
				connSeg[j]++;
			}
		}
	}
	
	removetemp = 0;
	for (set<int>::iterator it = removed_pet.begin(); it != removed_pet.end(); ++it)
	{
		////cout << inst->PETconnect[*it-1-removetemp][0] << " " << inst->PETconnect[*it-1-removetemp][1] << endl;
		vector<int>().swap(PETconnect[*it-1-removetemp]);
		PETconnect.erase(PETconnect.begin() + *it - 1 - removetemp);
		vector<vector<int> >().swap(PETreads[*it-1-removetemp]);
		PETreads.erase(PETreads.begin() + *it - 1 - removetemp);

		removetemp++;
	}

	for (int i = 0; i < PETconnect.size(); i++)
	{
		int s1_decrease = 0, s2_decrease = 0;
		int s1 = PETconnect[i][0];
		int s2 = PETconnect[i][1];
		for (set<int>::iterator it = removed_seg.begin(); it != removed_seg.end(); ++it)
		{
			if (*it < s1)
				s1_decrease++;
			if (*it < s2)
				s2_decrease++;
		}
		////cout << "s1: " << s1 << "->" << s1-s1_decrease << " ";
		////cout << "s2: " << s2 << "->" << s2-s2_decrease << " " << inst->SegConf.size() << endl;
		PETconnect[i][0] = s1-s1_decrease;
		PETconnect[i][1] = s2-s2_decrease;
	}

	NumofSegs = connSeg.size();
	set<int>().swap(removed_seg);
}

void Instance::getmap(map<int,set<int> > &exonmap, double &th)
{
	for(int i = 0; i < connSeg.size(); i++)
	{
		if (connSeg[i] > 1)
		{
			vector<int> singleSegConf = SegConf[i];
			for(int ii = 0; ii < singleSegConf.size(); ii++)
			{
				for(int j = ii + 1; j < singleSegConf.size(); j++)
				{
					if(singleSegConf[ii] + singleSegConf[j] == 2)
					{
						set<int>::iterator it = exonmap[ii+1].end();
						exonmap[ii+1].insert(it,j + 1);
						ii = j;
					}
				}
			}
			vector<int>().swap(singleSegConf);
			singleSegConf.clear();
		}
	}
	//use pair information to connect reads
	/*
	for (int i = 0; i < PETconnect.size(); i++)
	{
		//get pair of Segment
		int s1 = PETconnect[i][0];
		int s2 = PETconnect[i][1];
		//cout << "s1 " << s1 << " s2 " << s2 << " size: " << inst->SegConf.size() << endl;
		//test if connect directly
		int l1 = 0, l2 = 0, l3 = 0, bias = 0;
		vector<int> temppath;

		for (int j = 0; j < SegConf[0].size(); j++)
		{
			temppath.push_back((SegConf[s1-1][j] || SegConf[s2-1][j]));
		}
	
		int intron_temp = calintronlength(s1,s2,temppath,l1,l2,l3,0,bias);
		double eL = effectL(l1, l2, l3, ReadLen, frag_mu, frag_sigma,bias);
		if ((double) PETconnect[i][2] / eL > th+0.01)
		{
			int min_s1 = 1000000, min_s2 = 1000000;
			int max_s1 = -1, max_s2 = -1;
			for (int j = 0; j < SegConf[s1-1].size(); j++)
			{
				if (SegConf[s1-1][j] == 1)
				{
					if (min_s1 > j+1) 
						min_s1 = j+1;
					if (max_s1 < j+1)
						max_s1 = j+1;
				}
				if (SegConf[s2-1][j] == 1)
				{
					if (min_s2 > j+1) 
						min_s2 = j+1;
					if (max_s2 < j+1)
						max_s2 = j+1;
				}

			}

			int s_start = 0, s_end = 0;
			if (!strand.compare("-"))
			{
				s_start = max_s2;
				s_end = min_s1;
				
			}
			else if (!strand.compare("+"))
			{
				s_start = max_s1;
				s_end = min_s2;
			}
			if (s_start < s_end)
			{
				set<int>::iterator it = exonmap[s_start].end();
				set<int>::iterator it1 = exonmap[s_start].find(s_end);
				if (it1 == exonmap[s_start].end())
				{
					exonmap[s_start].insert(it,s_end);
					cout << "add edge: " << s_start << " " << s_end << " with el " << eL << " and reads " << PETconnect[i][2] << endl;
					cout << "eL calculated from " << l1 << " " << l2 << " " << l3 << " " << bias << endl;
				}
			}
		}
	}*/
}

int Instance::enumerate_path(map<int,set<int> > &exonmap, int not_source_list[], int not_sink_list[], vector<vector<int> > &paths)
{
	for (int i = 0; i < NumofExon; i++)
	{
		if (not_source_list[i] == 0)
		{
			//cout << "source " << i + 1 << endl;
			set<int>::iterator it = exonmap[0].end();
			exonmap[0].insert(it,i+1);
		}
	}
	for (int i = 0; i < NumofExon; i++)
	{
		//cout << not_sink_list[i] << endl;
		if (not_sink_list[i] == 0)
		{
			//cout << "sink " << i + 1 << endl;
			set<int>::iterator it = exonmap[i+1].end();
			exonmap[i+1].insert(it,NumofExon+1);
		}

	}

	//vector<vector<int> > paths;
	//vector<int> single_path;
	//apply depth first search
	list<Seg> stack;
	Seg* seg = new Seg();
	seg->node = 0;
	seg->visited.push_back(0);
	int goal = NumofExon + 1;
	stack.push_back(*seg);
	vector<int>().swap(seg->visited);
	delete seg;
	seg = NULL;
	int countpath = 0;
	while(!stack.empty())
	{
		Seg seg1 = stack.back();
		stack.pop_back();
		if (seg1.node == goal)
		{
			if (seg1.visited.size() > 3) // path with seg larger than 1
			{
				paths.push_back(seg1.visited);
				//mod_X(seg1.visited, SegConf, NumofSegs, NumofExon, outfile);
				countpath++;
				//path_all.push_back(seg1.visited);
			}
		}
		else
		{
			set<int> conn = exonmap[seg1.node];
			for (set<int>::iterator it = conn.begin(); it != conn.end(); ++it)
			{
				Seg* push_seg = new Seg();
				push_seg->node = *it;
				push_seg->visited = seg1.visited;
				push_seg->visited.push_back(*it);
				stack.push_back(*push_seg);
				vector<int>().swap(push_seg->visited);
				delete push_seg;
				push_seg = NULL;

			}
		}
		if (countpath > MAX_NUM_PATH)
		{
			for (int i = 0; i < paths.size(); i++)
			{
				paths[i].clear();
				vector<int>().swap(paths[i]);
			}
			paths.clear();
			vector<vector<int> >().swap(paths);
			return -1;
		}
			
	}
	if (countpath <= MAX_NUM_PATH)
	{
		/*for (int i = 0; i < paths.size(); i++)
			mod_X(paths[i],SegConf,NumofSegs,NumofExon,outfile);
		for (int i = 0; i < paths.size(); i++)
		{
			paths[i].clear();
			vector<int>().swap(paths[i]);
		}
		paths.clear();
		vector<vector<int> >().swap(paths);	*/
		return 0;
	}
	else
	{
		for (int i = 0; i < paths.size(); i++)
		{
			paths[i].clear();
			vector<int>().swap(paths[i]);
		}
		paths.clear();
		vector<vector<int> >().swap(paths);
		return -1;
	}

	
} 

void Instance::mod_X(vector<int> &path, vector<double> &proc_abun, vector<int> &segpath)
{
	//transpose of final X !!!!!!!!!!!
	//cout << NumofExon << endl;
	int path_vec[NumofExon];
	memset(path_vec,0,sizeof(path_vec));
	for (int i = 1; i < path.size() - 1; i++)
	{
		path_vec[path[i]-1] = 1;
		//cout << path[i] << " ";
	}
	//cout << endl;

	int seg_map[NumofSegs];
	memset(seg_map,0,sizeof(seg_map));
	int countmap = 0; double abun = 0;
	for (int i = 0; i < NumofSegs; i++)
	{
		if(ifcontain_arr(path_vec,SegConf[i]))
		{
			seg_map[i] = 1;
			countmap++;
			abun += proc_abun[i];
		}
		else
		{
			seg_map[i] = 0;
		}
	}
	double path_mean_cov = abun / (double) countmap;
	if (countmap > 0 && abun > 0.1)
	{
		for (int i = 0; i < NumofSegs; i++)
		{
			segpath.push_back(seg_map[i]);
		}
	//	for (int i = 0; i < NumofExon; i++)
	//	{
	//		outfile << path_vec[i] << " ";
	//	}
	//	outfile << endl;
	//	for (int i = 0; i < NumofSegs; i++)
	//	{
	//		outfile << seg_map[i] << " ";
	//	}
	//	outfile << endl;
	}
	else
	{
		segpath.clear();
	}
	//
}

void Instance::mod_y(map<int, vector<int> > &iscontained, vector<double> &SegAbundance1, int &output, double &th1, double &min_abun)
{
	//cout << "here" << endl;
	//vector<vector<int> > &SegConf, vector<int> &SegAbundance, vector<vector<int> > &ExonBound, vector<int> &connSeg, vector<double> &SegAbundance1, int &NumofSegs, int NumofExon, int READLEN
	vector<double> SegReads;
	int READLEN = ReadLen;
	

	int L[NumofExon];
	//vector<double> SegAbundance1;
	double max_abun = -1.0;
	double th = -1;
	int counttemp = 0;
	vector<double> sumy_v;
	vector<double> adjusted_length_v;
	for (int i = 0; i < NumofExon; i++)
	{
		L[i] = (exonbound[i])[1] - (exonbound[i][0]) + 1;
		//cout << L[i] << endl;
	}
	for (int i = 0; i < NumofSegs; i++)
	{
		double sumy = SegAbundance[i];
		double abun = 0;
		int adjusted_length = 1;
		map<int,vector<int> >::iterator it;
		it = iscontained.find(i+1);
		if (it != iscontained.end())
		{
			for (int ii = 0; ii < iscontained[i+1].size(); ii++)
			{
				sumy += SegAbundance[(iscontained[i+1])[ii]-1];
			}
			if (connSeg[i] == 1)
			{
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						adjusted_length = L[j];
						break;
					}
				}
				
			}
			else
			{
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						adjusted_length = L[j];
						break;
					}
				}
				adjusted_length = min(adjusted_length, READLEN);
				
			}
		}
		else
		{
			if (connSeg[i] == 1)
			{
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						adjusted_length = L[j];
						break;
					}
				}
			}
			else
			{
				int sumlength = 0; int count = 0; int init_length = 0; int final_length = 0;
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						final_length = L[j];
						if (count > 0)
							sumlength += L[j];
						else
							init_length = L[j];
						count++;
					}
				}
				adjusted_length = READLEN - (sumlength - final_length);
				adjusted_length = min(init_length, adjusted_length);
				adjusted_length = min(final_length, adjusted_length);
			}
		}
		if (adjusted_length <= 0)
			adjusted_length = 1;
		abun = sumy / ((double)adjusted_length);
		SegReads.push_back(sumy);
		SegAbundance1.push_back(abun);
		adjusted_length_v.push_back((double)adjusted_length);
		sumy_v.push_back(sumy);
		if (abun > max_abun)
			max_abun = abun;
		//cout << abun << "\t" << sumy << "\t" << adjusted_length << "\t" << connSeg[i] << endl;
		if (output == 1)
		{
			//outfile << "aaa\t" << ++counttemp << endl;
			//cout << abun << "\t" << sumy << "\t" << adjusted_length << endl;
			//outfile << abun << "\t" << sumy << "\t" << adjusted_length << endl;
			//cout << "Seg\t" << i << endl;
			/*for(int j = 0; j < SegConf[i].size(); j++)
				cout << SegConf[i][j] << " ";
			cout << endl;*/
		}
	}

	if (output == 0)
	{
		vector<int> removed_seg;
		set<int> removed_pet;
		min_abun  = 10000;
		for (int i = 0; i < NumofSegs; i++)
		{
			/*
			if ((connSeg[i] > 1 && SegAbundance1[i] < th * max_abun) || SegAbundance1[i] < th1+0.001)
			{
				removed_seg.push_back(i+1);
			}
			else
			{
				if(min_abun>SegAbundance1[i])
					min_abun = SegAbundance1[i];
			}*/

			if (connSeg[i] > 1 and SegReads[i] < th1+1)
			{
				removed_seg.push_back(i+1);
			}
			else
			{
				if(min_abun>SegReads[i] & connSeg[i] > 1)
					min_abun = SegReads[i];
			}
		}

		//cout << "remove ";
		for (int j = 0; j < removed_seg.size(); j++)
		{
			//cout << removed_seg[j] <<" ";
			for (int i = 0; i < PETconnect.size(); i++)
			{
				int s1 = PETconnect[i][0];
				int s2 = PETconnect[i][1];
				if (s1 == removed_seg[j] || s2 == removed_seg[j])
					removed_pet.insert(i+1);
			}
			vector<int>().swap(SegConf[removed_seg[j]-1-j]);
			SegConf.erase(SegConf.begin() + removed_seg[j] - 1 - j);
			connSeg.erase(connSeg.begin() + removed_seg[j] - 1 - j);
			SegAbundance.erase(SegAbundance.begin() + removed_seg[j] - 1 - j);
			SegAbundance1.erase(SegAbundance1.begin() + removed_seg[j] - 1 - j);
			//cout << "remove\t" << removed_seg[j] << endl;
		}
		//cout << endl;
		//cout << "remove\t" << removed_seg.size() << "\tsegs" << endl;
		NumofSegs = NumofSegs - removed_seg.size();
		//cout << "current number of segs\t" << NumofSegs <<endl;
		
		//cout << "removePET " << removed_pet.size() << endl;
		int removetemp = 0;
		for (set<int>::iterator it = removed_pet.begin(); it != removed_pet.end(); ++it)
		{
			//cout << inst->PETconnect[*it-1-removetemp][0] << " " << inst->PETconnect[*it-1-removetemp][1] << endl;
			vector<int>().swap(PETconnect[*it-1-removetemp]);
			PETconnect.erase(PETconnect.begin() + *it - 1 - removetemp);
			vector<vector<int> >().swap(PETreads[*it-1-removetemp]);
			PETreads.erase(PETreads.begin() + *it - 1 - removetemp);

			removetemp++;
		}

		for (int i = 0; i < PETconnect.size(); i++)
		{
			int s1_decrease = 0, s2_decrease = 0;
			int s1 = PETconnect[i][0];
			int s2 = PETconnect[i][1];
			for (int j = 0; j < removed_seg.size(); j++)
			{
				if (removed_seg[j] < s1)
					s1_decrease++;
				if (removed_seg[j] < s2)
					s2_decrease++;
			}
			//cout << "s1: " << s1 << "->" << s1-s1_decrease << " ";
			//cout << "s2: " << s2 << "->" << s2-s2_decrease << " " << inst->SegConf.size() << endl;
			PETconnect[i][0] = s1-s1_decrease;
			PETconnect[i][1] = s2-s2_decrease;
		}

		removed_seg.clear();
		vector<int>().swap(removed_seg);
		/*for (int i = 0; i < inst->SegConf.size(); i++)
		{
			cout << "Seg " << i << ": ";
			for (int j = 0; j < inst->SegConf[i].size(); j++)
			{
				cout << inst->SegConf[i][j] << " ";
			}
			cout << endl;
		}

		vector<int> removed_exon;
		//check if only single exon left
		for (int i = 0; i < inst->NumofExon; i++)
		{
			int countconn = 0;
			for (int j = 0; j < inst->SegConf.size(); j++)
			{
				if (inst->SegConf[j][i] == 1)
				{
					int countconn1 = 0;
					for (int jj = 0; jj < inst->SegConf[j].size();jj++)
					{
						countconn1 += inst->SegConf[j][jj];			
					}
					if (countconn1 > 1)
						countconn++;
				}

			}

			if (countconn == 0)
				removed_exon.push_back(i+1);
		}

		for (int i = 0; i < removed_exon.size(); i++)
		{
			cout << "remove single exon chain: " << removed_exon[i] << endl;
			vector<int>().swap(inst->exonbound[removed_exon[i]-1-i]);
			inst->exonbound.erase(inst->exonbound.begin()+removed_exon[i]-1-i);
			inst->exonAbundance.erase(inst->exonAbundance.begin()+removed_exon[i] - 1 -i);
			for (int j = 0; j < inst->NumofSegs; j++)
			{
				inst->SegConf[j].erase(inst->SegConf[j].begin()+removed_exon[i] - 1 - i);
			}
		}
		inst->NumofExon = inst->NumofExon - removed_exon.size();

		for (int j = 0; j < inst->SegConf.size(); j++)
		{
			int countconn1 = 0;
			for (int jj = 0; jj < inst->SegConf[j].size();jj++)
			{
				countconn1 += inst->SegConf[j][jj];			
			}
			if (countconn1 == 0)
				removed_seg.push_back(j+1);
			}


		for (int j = 0; j < removed_seg.size(); j++)
		{
			vector<int>().swap(inst->SegConf[removed_seg[j]-1-j]);
			inst->SegConf.erase(inst->SegConf.begin() + removed_seg[j] - 1 - j);
			inst->connSeg.erase(inst->connSeg.begin() + removed_seg[j] - 1 - j);
			inst->SegAbundance.erase(inst->SegAbundance.begin() + removed_seg[j] - 1 - j);
			//cout << "remove\t" << removed_seg[j] << endl;
		}
		cout << "remove\t" << removed_seg.size() << "\tsegs" << endl;
		inst->NumofSegs = NumofSegs - removed_seg.size();
		cout << "current number of segs\t" << inst->NumofSegs <<endl;

		for (int i = 0; i < inst->SegConf.size(); i++)
		{
			cout << "Seg " << i << ": ";
			for (int j = 0; j < inst->SegConf[i].size(); j++)
			{
				cout << inst->SegConf[i][j] << " ";
			}
			cout << endl;
		}*/

	}

}

void Instance::mod_y_se(map<int, vector<int> > &iscontained, vector<double> &SegAbundance1, vector<double> &SegReads, vector<double> &SegLength, int &output, double &th1, double &min_abun)
{
	//cout << "here" << endl;
	//vector<vector<int> > &SegConf, vector<int> &SegAbundance, vector<vector<int> > &ExonBound, vector<int> &connSeg, vector<double> &SegAbundance1, int &NumofSegs, int NumofExon, int READLEN
	int READLEN = ReadLen;
	

	int L[NumofExon];
	//vector<double> SegAbundance1;
	double max_abun = -1.0;
	double th = -1;
	int counttemp = 0;
	vector<double> sumy_v;
	vector<double> adjusted_length_v;
	for (int i = 0; i < NumofExon; i++)
	{
		L[i] = (exonbound[i])[1] - (exonbound[i][0]) + 1;
		//cout << L[i] << endl;
	}
	for (int i = 0; i < NumofSegs; i++)
	{
		double sumy = SegAbundance[i];
		double abun = 0;
		int adjusted_length = 1;
		map<int,vector<int> >::iterator it;
		it = iscontained.find(i+1);
		if (it != iscontained.end())
		{
			for (int ii = 0; ii < iscontained[i+1].size(); ii++)
			{
				sumy += SegAbundance[(iscontained[i+1])[ii]-1];
			}
			if (connSeg[i] == 1)
			{
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						adjusted_length = L[j];
						break;
					}
				}
				
			}
			else
			{
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						adjusted_length = L[j];
						break;
					}
				}
				adjusted_length = min(adjusted_length, READLEN);
				
			}
		}
		else
		{
			if (connSeg[i] == 1)
			{
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						adjusted_length = L[j];
						break;
					}
				}
			}
			else
			{
				int sumlength = 0; int count = 0; int init_length = 0; int final_length = 0;
				for(int j = 0; j < NumofExon; j++)
				{
					if((SegConf[i])[j] == 1)
					{
						final_length = L[j];
						if (count > 0)
							sumlength += L[j];
						else
							init_length = L[j];
						count++;
					}
				}
				adjusted_length = READLEN - (sumlength - final_length);
				adjusted_length = min(init_length, adjusted_length);
				adjusted_length = min(final_length, adjusted_length);
			}
		}
		if (adjusted_length <= 0)
			adjusted_length = 1;
		abun = sumy / ((double)adjusted_length);
		SegAbundance1.push_back(abun);
		SegReads.push_back(sumy);
		SegLength.push_back(adjusted_length);
		adjusted_length_v.push_back((double)adjusted_length);
		sumy_v.push_back(sumy);
		if (abun > max_abun)
			max_abun = abun;
		//cout << abun << "\t" << sumy << "\t" << adjusted_length << "\t" << connSeg[i] << endl;
		if (output == 1)
		{
			//outfile << "aaa\t" << ++counttemp << endl;
			//cout << abun << "\t" << sumy << "\t" << adjusted_length << endl;
			//outfile << abun << "\t" << sumy << "\t" << adjusted_length << endl;
			//cout << "Seg\t" << i << endl;
			/*for(int j = 0; j < SegConf[i].size(); j++)
				cout << SegConf[i][j] << " ";
			cout << endl;*/
		}
	}

	if (output == 0)
	{
		vector<int> removed_seg;
		set<int> removed_pet;
		min_abun  = 10000;
		for (int i = 0; i < NumofSegs; i++)
		{
			//cout << SegAbundance1[i] - th * max_abun << endl;
			if ((connSeg[i] > 1 && SegAbundance1[i] < th * max_abun) || SegAbundance1[i] < th1+0.001)
			{
				removed_seg.push_back(i+1);
			}
			else
			{
				if(min_abun>SegAbundance1[i])
					min_abun = SegAbundance1[i];
			}
		}

		//cout << "remove ";
		for (int j = 0; j < removed_seg.size(); j++)
		{
			//cout << removed_seg[j] <<" ";
			for (int i = 0; i < PETconnect.size(); i++)
			{
				int s1 = PETconnect[i][0];
				int s2 = PETconnect[i][1];
				if (s1 == removed_seg[j] || s2 == removed_seg[j])
					removed_pet.insert(i+1);
			}
			vector<int>().swap(SegConf[removed_seg[j]-1-j]);
			SegConf.erase(SegConf.begin() + removed_seg[j] - 1 - j);
			connSeg.erase(connSeg.begin() + removed_seg[j] - 1 - j);
			SegAbundance.erase(SegAbundance.begin() + removed_seg[j] - 1 - j);
			SegAbundance1.erase(SegAbundance1.begin() + removed_seg[j] - 1 - j);
		}
		//cout << "remove\t" << removed_seg.size() << "\tsegs" << endl;
		NumofSegs = NumofSegs - removed_seg.size();
		//cout << "current number of segs\t" << NumofSegs <<endl;
		
		int removetemp = 0;
		for (set<int>::iterator it = removed_pet.begin(); it != removed_pet.end(); ++it)
		{
			//cout << inst->PETconnect[*it-1-removetemp][0] << " " << inst->PETconnect[*it-1-removetemp][1] << endl;
			vector<int>().swap(PETconnect[*it-1-removetemp]);
			PETconnect.erase(PETconnect.begin() + *it - 1 - removetemp);
			vector<vector<int> >().swap(PETreads[*it-1-removetemp]);
			PETreads.erase(PETreads.begin() + *it - 1 - removetemp);

			removetemp++;
		}

		for (int i = 0; i < PETconnect.size(); i++)
		{
			int s1_decrease = 0, s2_decrease = 0;
			int s1 = PETconnect[i][0];
			int s2 = PETconnect[i][1];
			for (int j = 0; j < removed_seg.size(); j++)
			{
				if (removed_seg[j] < s1)
					s1_decrease++;
				if (removed_seg[j] < s2)
					s2_decrease++;
			}
			//cout << "s1: " << s1 << "->" << s1-s1_decrease << " ";
			//cout << "s2: " << s2 << "->" << s2-s2_decrease << " " << inst->SegConf.size() << endl;
			PETconnect[i][0] = s1-s1_decrease;
			PETconnect[i][1] = s2-s2_decrease;
		}

		removed_seg.clear();
		vector<int>().swap(removed_seg);

	}

}

void Instance::findcontain(map<int, vector<int> > &iscontained) // contain: index start from 1
{
	for (int i = 0; i < SegConf.size(); i++)
	{
		for (int j = 0; j < SegConf.size(); j++)
		{
			if (ifcontain_vec(SegConf[i],SegConf[j]) && i != j)
			{
				iscontained[j+1].push_back(i+1);
				//cout << j+1 << "\tcontain\t" << i+1 << endl;
			}
		}
	}
}

bool ifcontain_vec(vector<int> &a, vector<int> &b)
{
	int tempmax = INT_MAX;
	int ifcontain = 0; int idx1 = tempmax; int countj = 0;
	for (int ii = 0; ii < a.size(); ii++)
	{
		if (a[ii] - b[ii] == -1)
		{
			ifcontain = -1;
			return false;
		}
		else if (a[ii] - b[ii] == 1)
		{
			if (idx1 > ii)
			{
				idx1 = ii;
			}
		} 
		if (b[ii] == 1 && countj >= 0)
		{
			//cout << "here\t" << idx1 << "\t" << ii << endl;
			if (idx1 < ii)
			{
				ifcontain = -1;
				return false;
			}
			countj += 1;
		}
		
		ifcontain += a[ii] - b[ii];
	}
	if (ifcontain > 0)
	{
		return true;
	}
}

bool ifcontain_arr(int a[], vector<int> &b) //path to segs
{
	int tempmin = INT_MAX;
	int tempmax = -1;
	//int ifcontain = 0; int idx1 = tempmax; int countj = 0;
	for (int ii = 0; ii < b.size(); ii++)
	{
		if (a[ii] - b[ii] == -1)
			return false;
		if (b[ii] == 1)
		{
			//countj++;
			if (ii < tempmin)
				tempmin = ii;
			if (ii > tempmax)
				tempmax = ii;
		}
	}
	if (tempmin == tempmax)
		return true;
	else
	{
		for (int i = tempmin+1; i <=tempmax-1; i++)
		{
			if (a[i] - b[i] == 1)
				return false;
		}	
	}

	return true;
}

void print_vector(vector<int> aa)
{
	for (int i = 0; i < aa.size(); i++)
	{
		cout << aa[i] << " ";
	}
	cout << endl;
}

int Instance::calintronlength(int &s1, int &s2, vector<int> &path, int &l1, int &l2, int &l3, int label, int &bias)
{
	int exonLength[NumofExon];
	for (int i = 0; i < NumofExon; i++)
	{
		exonLength[i] = exonbound[i][1] - exonbound[i][0] + 1;
	}
	int s1_interl = 0, s2_interl = 0;
	if (label == 0)
	{
		if(connSeg[s1-1] == 1)
		{
			for (int i = 0; i < SegConf[s1-1].size(); i++)
			{
				if (SegConf[s1-1][i] == 1)
					l1 = exonLength[i]-ReadLen+1;
			}
		}
		else
		{
			int sumlength = 0; int count = 0; int init_length = 0; int final_length = 0;
			for(int j = 0; j < NumofExon; j++)
			{
				if((SegConf[s1-1])[j] == 1)
				{
					final_length = exonLength[j];
					if (count > 0)
						sumlength += exonLength[j];
					else
						init_length = exonLength[j];
					count++;
				}
			}
			s1_interl = sumlength - final_length;
			l1 = ReadLen - (sumlength - final_length) - 1;
			l1 = min(init_length, l1);
			l1 = min(final_length, l1);
		}

		
		if(connSeg[s2-1] == 1)
		{
			for (int i = 0; i < SegConf[s2-1].size(); i++)
			{
				if (SegConf[s2-1][i] == 1)
					l2 = exonLength[i]-ReadLen+1;
			}
		}
		else
		{
			int sumlength = 0; int count = 0; int init_length = 0; int final_length = 0;
			for(int j = 0; j < NumofExon; j++)
			{
				if((SegConf[s2-1])[j] == 1)
				{
					final_length = exonLength[j];
					if (count > 0)
						sumlength += exonLength[j];
					else
						init_length = exonLength[j];
					count++;
				}
			}
			s2_interl = sumlength - final_length;
			l2 = ReadLen - (sumlength - final_length) - 1;
			l2 = min(init_length, l2);
			l2 = min(final_length, l2);
		}
	}
	
	l3 = 0;


	int s1_end = 0, s2_start = 0;
	for(s1_end = SegConf[s1-1].size()-1; s1_end>= 0; s1_end--)
	{
		if (SegConf[s1-1][s1_end] == 1)
			break;
	}
	for(s2_start = 0; s2_start < SegConf[s2-1].size(); s2_start++)
	{
		if (SegConf[s2-1][s2_start] == 1)
			break;
	}
	int intron_in = 0; // the method to calculate intron_in now is wrong
	int temps1 = 0, temps2 = 0;
	int tempbias1 = 0, tempbias2 = 0;
	if (connSeg[s1-1] > 1)
	{
		if (!strand.compare("-"))
		{
		    //tempbias1 = inst->exonbound[s1_end][1] - inst->ReadLen + s1_interl;
		    tempbias1 = exonbound[s1_end][1] - min(exonLength[s1_end]-1, ReadLen - s1_interl - 1 - 1);
		    //tempbias1 = inst->exonbound[s1_end][1] - l1;
		}
		else
		{
			tempbias1 = exonbound[s1_end][0] + min(exonLength[s1_end]-1, ReadLen - s1_interl - 1 - 1);
			//tempbias1 = inst->exonbound[s1_end][0] + l1;
		}
		int templ3_1 = 0;
		for (int i = 0; i < SegConf[s1-1].size(); i++)
		{
			if (SegConf[s1-1][i] == 1 && i != s1_end)
				templ3_1 += exonLength[i];
		}
		if (templ3_1 > ReadLen)
			l3 += exonLength[s1_end] - l1;
		else
			l3 += exonLength[s1_end] - (ReadLen - templ3_1) - l1;
		l3 = max(0,l3);
	}
	else 
	{
		if (!strand.compare("-"))
		    tempbias1 = exonbound[s1_end][1] - l1 - ReadLen + 1 + 1;		
		else
			tempbias1 = exonbound[s1_end][0] + l1 + ReadLen - 1 - 1;
	}


	if (!strand.compare("-"))
	    temps1 = exonbound[s1_end][0];
	else
		temps1 = exonbound[s1_end][1];
	


	if (connSeg[s2-1] > 1)
	{
		if (!strand.compare("-"))
		{
		    //l3 += exonLength[s2_start] - l2;
		    //tempbias2 = inst->exonbound[s2_start][0] + inst->ReadLen - s2_interl;
		    tempbias2 = exonbound[s2_start][0] + min(exonLength[s2_start]-1, ReadLen - s2_interl - 1 - 1);
		    //tempbias2 = inst->exonbound[s2_start][0] + l2;
		}
		
		else
		{
			//l3 += exonLength[s2_start] - l2;
			//tempbias2 = inst->exonbound[s2_start][1] - inst->ReadLen + s2_interl;
			tempbias2 = exonbound[s2_start][1] - min(exonLength[s2_start]-1, ReadLen - s2_interl - 1 - 1);
			//tempbias2 = inst->exonbound[s2_start][1] - l2;
		}
		int templ3_1 = 0;
		for (int i = 0; i < SegConf[s2-1].size(); i++)
		{
			if (SegConf[s2-1][i] == 1 && i != s2_start)
				templ3_1 += exonLength[i];
		}
		if (templ3_1 > ReadLen)
			l3 += exonLength[s2_start] - l2;
		else
			l3 += exonLength[s2_start] - (ReadLen - templ3_1) - l2;
		l3 = max(0,l3);

	}
	else 
	{
		if (!strand.compare("-"))
		    tempbias2 = exonbound[s2_start][0] + l2 + ReadLen - 1 - 1;
		else
			tempbias2 = exonbound[s2_start][1] - l2 + 1 - ReadLen + 1;
	}


	if (!strand.compare("-"))
	    temps2 = exonbound[s2_start][1];
	else
		temps2 = exonbound[s2_start][0];
	
	//calculate bias (the length of two ends may overlap)
	if (!strand.compare("-"))
	{
	    bias = tempbias1 - tempbias2 - 1;
	    //intron_in = tempbias1 - tempbias2 - 1;
	    bias = min(bias,0);
	}
	else
	{
		bias = tempbias2 - tempbias1 - 1;
	//	intron_in = tempbias2 - tempbias1 - 1;
		bias = min(bias,0);
	}

	//cout << "tempbias1 " << tempbias1 << " tempbias2 " << tempbias2 << " " << strand << " " << " " << strand.size() << " " << !strand.compare("-") << endl;
	intronlength_helper(tempbias1, tempbias2, path, intron_in, l3, bias);
	//cout << "intron_in " << l3 << " bias " << bias << endl;
	//l3 = max(l3,0);
	//bias = min(bias,0);
	/*
	if (s1_end >= s2_start)
	{
		intron_in = 0;
		//cout << "l3 calculate " << l3 << " " << inst->connSeg[s1-1] << " " << inst->connSeg[s2-1] << endl;
		if (connSeg[s1-1] == 1 || connSeg[s2-1] == 1)
			l3 = 0;
		else
		{
			if (s1_end > s2_start)
				l3 = 0;
			else
			{
				if (bias < 0)
					l3 = 0;
				else
					l3 = exonLength[s2_start] - (ReadLen-s1_interl) - (ReadLen-s2_interl);
				l3 = max(0,l3);
			}
			//l3 = l3 / 2 - inst->ReadLen;
			 
		}
		//cout << "l3 calculate1 " << l3 << endl;
	}
	else
	{

		if (!strand.compare("-"))
			intron_in = temps1 - temps2 - 1;
		else
			intron_in = temps2 - temps1 - 1;

		for (int i = s1_end + 1; i < s2_start; i++)
		{
			if (path[i] == 1)
			{
				intron_in -= exonLength[i];
				l3 += exonLength[i];
			}
		}
	}
	*/

	if (s1 == s2)
	{
		l2 = l1;
		//l3 = 0;
		if (connSeg[s1-1] == 1)
		{
			for (int i = 0; i < SegConf[s1-1].size(); i++)
			{
				if (SegConf[s1-1][i] == 1)
				{
					//bias = -exonLength[i];
					break;
				}
			}
		}
		else
		{
			//bias = -2*ReadLen+2;
		}
		//if (connSeg[s1-1] > 1)
		//	intron_in = temps1 - temps2 - 1;
	}


	return intron_in;

}

int Instance::segStart(int &s1)
{

	int s1_start = 0;
	for(s1_start = 0; s1_start < SegConf[s1-1].size(); s1_start++)
	{
		if (SegConf[s1-1][s1_start] == 1)
			break;
	}
	return s1_start; // should be s1_start+1

}

void Instance::segStart(int &s1, int &s2, set<int> &s3)
{
	for(int i = 0; i < SegConf[s1-1].size(); i++)
	{
		if (SegConf[s1-1][i] == 1)
			s3.insert(i+1);
	}

	for(int i = 0; i < SegConf[s2-1].size(); i++)
	{
		if (SegConf[s2-1][i] == 1)
			s3.insert(i+1);
	}
}

void Instance::intronlength_helper(int &pos1, int &pos2, vector<int> &path, int &intron_l, int &l3, int &bias)
{
	intron_l = 0;
	l3 = 0;
	bias = 0;
	if (!strand.compare("-"))
	{
		//cout << "come neg\n";
		if (pos1 < pos2)
		{
			intron_l = 0;
			l3 = 0;
			int l3_temp = 0;
			exonlength_reverse(pos2,pos1,path,intron_l,bias);
			intron_l = -intron_l;
			bias += 1;
			bias = -bias;
		}
		else
		{
			//check from pos1 to pos2 the length of intron and exon
			exonlength_reverse(pos1,pos2,path,intron_l,l3);
		}
	}
	else
	{
		//cout << "come pos\n";
		if (pos1 > pos2)
		{
			intron_l = 0;
			l3 = 0;
			int l3_temp = 0;
			exonlength_plus(pos2,pos1,path,intron_l,bias);
			intron_l = -intron_l;
			bias += 1;
			bias = -bias;
		}
		else
		{
			//check from pos1 to pos2 the length of intron and exon
			exonlength_plus(pos1,pos2,path,intron_l,l3);
		}
	}
	l3 = max(l3,0);
	bias = min(bias,0);
}

void Instance::exonlength_plus(int &pos1, int &pos2, vector<int> &path, int &intron_l,int &l3)
{
	intron_l = 0;
	l3 = 0;
	int matched_in = 0;
	int preidx = -1;
	for (int i = 0; i < path.size(); i++)
	{
		if (path[i] == 1)
		{
			if (matched_in == 1)
			{
				intron_l += exonbound[i][0] - exonbound[preidx][1] - 1;
				if (pos2 <= exonbound[i][1])
				{
					l3 += pos2 - exonbound[i][1];
					break;
				}

				preidx = i;
			}
			if (pos1 >= exonbound[i][0] and pos1 <= exonbound[i][1])
			{
				preidx = i;
				matched_in = 1;
				if (pos2 <= exonbound[i][1])
				{
					l3 = pos2 - pos1 - 1;
					//cout << "comehere " << pos2 << " " << pos1 << " " << l3 << endl;
					intron_l = 0;
					break;
				}
				l3 += exonbound[i][1] - pos1;
			}
		}
	}
}

void Instance::exonlength_reverse(int &pos1, int &pos2, vector<int> &path, int &intron_l, int &l3)
{
	/*intron_l = 0;
	l3 = 0;
	int matched_in = 0;
	int preidx = -1;
	for (int i = 0; i < path.size(); i++)
	{
		if (path[i] == 1)
		{
			if (matched_in == 1)
			{
				intron_l += exonbound[preidx][0] - exonbound[i][1] - 1;
				if (pos1 <= exonbound[i][1])
				{
					l3 += pos1 - exonbound[i][1];
					break;
				}

				preidx = i;
			}
			if (pos2 >= exonbound[i][0] & pos2 <= exonbound[i][1])
			{
				preidx = i;
				matched_in = 1;
				if (pos1 <= exonbound[i][1])
				{
					l3 = pos1 - pos2 - 1;
					intron_l = 0;
					break;
				}
				l3 += exonbound[i][1] - pos2 - 1;
			}
		}
	}*/
	intron_l = 0;
	l3 = 0;
	int matched_in = 0;
	int preidx = -1;
	for (int i = 0; i < path.size(); i++)
	{
		if (path[i] == 1)
		{
			if (matched_in == 1)
			{
				intron_l += exonbound[preidx][0] - exonbound[i][1] - 1;
				if (pos2 >= exonbound[i][0])
				{
					l3 += exonbound[i][1] - pos2;
					break;
				}

				preidx = i;
			}
			if (pos1 >= exonbound[i][0] & pos1 <= exonbound[i][1])
			{
				preidx = i;
				matched_in = 1;
				if (pos2 >= exonbound[i][0])
				{
					l3 = pos1 - pos2 - 1;
					intron_l = 0;
					break;
				}
				l3 += pos1 - exonbound[i][0];
			}
		}
	}
}


double normpdf(double x, double mu, double sigma)
{
	return exp(-(double)1/2/pow(sigma,2)*pow(x-mu,2));
}


double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}

double phi_range(double x, double y)
{
	return (phi(y)-phi(x));
}

double effectL(int l1, int l2, int l3, int r, double mu, double sigma, int bias)
{
	int minl = min(l1,l2);
	int maxl = max(l1,l2);
	double l = 0;
	for (int i = bias; i < bias+minl; i++)
	{
		double rangel = i + l3 + 2 * r;
		double ranger = rangel + maxl;
		rangel = (rangel - mu - 1) / sigma;
		ranger = (ranger - mu) / sigma;
		l += phi_range(rangel,ranger);
	}
	return l;
}

