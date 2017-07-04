#include "groupedseqs_a.h"

inline void GroupedSeqs_A::init_seq_parsed(Seq_A& seq, string desc = "N/A", string bp = "", string mid = "NNNNNNNNNN")
{
	seq.description = desc;
	seq.basepairs = bp;

	string tempstr = "";
	seq.size = bp.size();
	tempstr = bp;
	transform(tempstr.begin(), tempstr.end(), tempstr.begin(), ::toupper);

	int i = int();
	int done = int();
	while (!done && (int)tempstr.size()>(i + 30))
	{
		if (tempstr.substr(i, 10) == mid && tempstr.substr(i + 14, 2) == "CA" && tempstr.substr(i + 20, 2) == "GT" && tempstr.substr(i + 27, 2) == "GC")
		{
			seq.xx_primerid_xx = bp.substr(i + 8, 21);
			done++;
		}
		i++;
	}
	if (!done)
	{
		seq.xx_primerid_xx = "NNNNNNNNNN"; // consider mutated
	}
}

int GroupedSeqs_A::read_file_to_seqs(const char *filename, const string mid)
{
	ifstream file;
	if (!open_infile(file, filename)) {
		cout << "unable to open file: " << filename << endl;
		return 0;
	}

	string line, line2, line3, test;
	string greaterthen = ">";
	int breakout;
	Seq_A tempseq;

	while (file.good())
	{
		if (line3[0] == greaterthen[0]) line = line3;
		else getline(file, line);
		
		if (line != "end.") {
			if (line[0] == greaterthen[0]) {
				getline(file, line2);
				if (line2 != "end.") {
					breakout = 0;
					do {
						breakout++;
						getline(file, line3);
						if (line3[0] != greaterthen[0] && line3 != "end.") line2 += line3;
					} while (line3[0] != greaterthen[0] && line3 != "end." && breakout<100);
					init_seq_parsed(tempseq, line, line2, mid); 

					if (tempseq.xx_primerid_xx == "NNNNNNNNNN") seqs_mutated_id.push_back(tempseq);
					else if (tempseq.size > 350) {
						seqs_for_groups.insert(tempseq); // Set elements are unique, removin possible dups
					}
					else seqs_too_short.push_back(tempseq); // < 350 bp
				}
			}
		}
	}
	file.close();
	cout << "How many unmutated (and long enough) sequences (to be counted in groupsizes.txt): " << seqs_for_groups.size() << endl;
	cout << "How many sequences with mutated id (discarded): " << seqs_mutated_id.size() << endl;
	cout << "How many sequences too short (discarded): " << seqs_too_short.size() << endl;
	return 1;
}

GroupedSeqs_A::GroupedSeqs_A()
{;}

GroupedSeqs_A::~GroupedSeqs_A()
{
	seqs_ungrouped.shrink_to_fit();
	seqs_too_short.shrink_to_fit();
	seqs_mutated_id.shrink_to_fit();
}
