 /*
	fqStatN outputs statistic of ambiguous reference characters ‘N’ in fq-file.
	Sample output:

	'N' POSITION STATISTICS
	pos     count   % of total 'N'
	-----------------------------
	 0      29588   71.1%
	 3      4314    10.4%
	20      69      0.166%
	...
	READ TEMPLATE STATISTICS
	position  10        20        30        40          patterns count
	01234567890123456789012345678901234567890123456789
	----------------------------------------------------------------------
	N.................................................    29566     0.172%
	...N..............................................     4313     0.0251%
	...............................N..NN.............N      213     0.00124%
	...
	'N' relative to the total number of nucleotides: 0.00483%
	Reads that include 'N' relative to the total number of reads: 0.211%


	Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)

	This program is free software. It is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the	GNU General Public License for more details.
  */

#include "fqStatN.h"
#include <fstream>
#include <algorithm>    // std::sort

using namespace std;

const string Product::Title = "fqStatN";
const string Product::Version = "1.0";
const string Product::Descr = "Statistics of 'N' in fq-file";

const string OutFile = string(Product::Title) +  "_out.txt";
const string HelpOutFile = "duplicate standart output to " + OutFile + " file.";

enum eOptGroup	{ oOPTION };	// oOTHER should be the last
const char* Options::_OptGroups [] = { NULL };
const BYTE	Options::_GroupCount = oOPTION + 1;

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::_Options [] = {
	{ 'o', "out",	0,	tENUM,	oOPTION, FALSE,	vUNDEF, 2, NULL, HelpOutFile.c_str(), NULL },
	{ 't', "time",	0,	tENUM,	oOPTION, FALSE,	vUNDEF, 2, NULL, "print run time", NULL },
	{ 'h', "help",	0,	tHELP,	oOPTION, vUNDEF, vUNDEF, 0, NULL, "print usage information", NULL }
};

const BYTE Options::_OptCount = oHELP + 1;
const BYTE Options::_UsageCount = 1;
const Options::Usage Options::_Usages[] = {
	{ vUNDEF, "file.fq", true, NULL }
};

ofstream outfile;					// file ostream duplicated cout; inizialised by file in code
dostream dout(cout, outfile);	// stream's duplicator


int main(int argc, char* argv[])
{
	if (argc < 2)	return Options::PrintUsage(false);			// output tip
	int fileInd = Options::Tokenize(argc, argv, "FQ file");
	if( fileInd < 0 )	return 1;								// wrong option

	int ret = 0;					// main() return code
	if( Options::GetBVal(oOUTPUT) )	outfile.open( OutFile.c_str() );
	//setlocale(LC_ALL, "");

	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		const char* fName = FS::CheckedFileName(argv[fileInd]);	// alignment name
		FT::CheckType(fName, FT::FQ);
		dout << fName << SepCl;
		fflush(stdout);
		FqFile fq(fName);
		StatN::Scan(fq);
	}
	catch(const Err &e)			{ ret = 1;	dout << e.what() << EOL; }
	catch(const exception &e)	{ ret = 1;	dout << e.what() << EOL; }
	catch(...)					{ ret = 1;	dout << "Unregistered error\n"; }
	timer.Stop(true);
	if( outfile.is_open() )		outfile.close();	// in case of holding execution by user
	return ret;
}

/************************  class StatN ************************/

void StatN::Scan	(FqFile & fqFile)
{
	UINT i, n;	// n used as N counter in Read in 1st part (counting), and as counter in 2th (output)
	ULONG cntTotalN = 0, cntTotalReads = 0;
	readlen k;
	bool insert;
	const char *read;
	vector<TemplN> templs;

	templs.reserve(20);
	fqFile.GetSequence();
	readlen rLen = fqFile.ReadLength();
	Array<char>	buf(rLen+1);
	Array<ULONG> distr(rLen);	// array of 'N' frequencies

	ULONG cnt = 0;
	// GET OCCURENCES
	do {
		read = fqFile.GetCurrRead();
		buf.Clear();
		for(n=0, i=0; i<rLen; i++)
			if( read[i] == cN ) {
				buf[n++] = i+1;
				distr[i]++;
			}
		// insert new line in statistic
		if(n > 0) {
			insert = false;
			for(i=0; i<templs.size(); i++)
				if( templs[i].Count == n
				&& !strcmp(buf.Data(), templs[i].Pos.Data()) ) {
					templs[i].CountRead++;
					cntTotalReads++;
					insert = true;
					break;
				}
			if( !insert )
			//{	cout << (++cnt) << EOL;
				templs.push_back(TemplN(n, buf));
			//}
			cntTotalN += n;
		}
	}
	while( fqFile.GetSequence() );

	// OUTPUT RESULT
	dout << fqFile.Count() << " reads\n";
	if( templs.size() ) {
		dout<< "'N' POSITION STATISTICS\n"
			<< "pos\tcount\t% of total 'N'\n"
			<< "------------------------------\n";
		for(k=0; k<rLen; k++)
			if( distr[k] )
				dout << setw(2) << int(k) << TAB << distr[k] << TAB
						<< sPercent(Percent(distr[k], cntTotalN), 3, 0, false) << EOL;

		dout << "\nREAD PATTERN STATISTICS\n";
		//sort(templs.begin(), templs.end(), Compare);
		// OUTPUT HEADER
		dout << "position";
		// output tens
		for(k=1, n=6; k<=rLen/10; n=0, k++) {
			for(; n<8; n++) 	dout << ' ';
			dout << k*10;
		}
		dout << "\tcount\n";
		// output units
		for(k=0; k<rLen/10; k++)
			for(n=0; n<10; n++)		dout << n;
		// output rest of units
		k = rLen%10;
		for(n=0; n<k; n++)		dout << n;
		dout << EOL;
		for(k=0; k<rLen+2*8+4; k++)	dout << HPH;
		dout << EOL;

		// OUTPUT ENTRIES
		for(i=0; i<templs.size(); i++) {
			for(n=k=0; k<rLen; k++)
				if(templs[i].Pos[n]-1 == k) {	dout << cN; n++; }
				else							dout << DOT;
			dout << setw(9) << templs[i].CountRead;
			dout << TAB << sPercent(Percent(templs[i].CountRead, fqFile.Count()), 3, 0, false) << EOL;
		}
	
		dout << "\n'N' relative to the total number of nucleotides: " 
			<< sPercent(Percent(cntTotalN, fqFile.Count()*rLen), 3, 0, false)  << EOL;
		dout << "Reads that include 'N' relative to the total number of reads: " 
			<< sPercent(Percent(cntTotalReads, fqFile.Count()), 3, 0, false)  << EOL;
	}
	else
		dout << "No reads included 'N'\n";
}

/************************  end of class StatN ************************/

