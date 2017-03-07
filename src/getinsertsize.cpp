#include "getinsertsize.h"

using namespace std;

double calmedian(vector<int> &fraglen);

void getinsertsize(string filename, double &mu, double &stdvar)
{
	int N = 10000000;
	int count = 1;
	vector<int> fraglen;
	ifstream infile;
	infile.open(filename.c_str());
	string line;
	if (infile.is_open())
	{
		while (infile.good())
		{
			getline(infile,line);
			istringstream ss(line);
			string s1,s2,s3,s4,s5,s6,s7,s8,s10;
			int s9;
			ss >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9 >> s10;
	        if (s7 != "=")
	        	continue;
	        int idx = s6.find_first_of("N");
	        if (idx > 0)
	        	continue;
	        if (s9>0)
	        	fraglen.push_back(s9);
	        if (count > N)
	        	break;
	        count++;
		}
    }
    mu = calmedian(fraglen);
    vector<int> fraglendiff;
    for (int i = 0; i < fraglen.size(); i++)
    {
    	if (abs(fraglen[i]-round(mu)) > 3*round(mu))
    		continue;
    	fraglendiff.push_back(abs(fraglen[i]-round(mu)));
    }
    double MAD = calmedian(fraglendiff);
    stdvar = 1.4826*MAD;
}

double calmedian(vector<int> &fraglen)
{
	sort(fraglen.begin(),fraglen.end());
	int N = fraglen.size();
	if (N%2==1)
		return (double)fraglen[(N+1)/2-1];
	else
		return ((double)(fraglen[N/2-1]+fraglen[N/2]))/2;
}

