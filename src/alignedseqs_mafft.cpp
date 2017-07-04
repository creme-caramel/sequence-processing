#include "alignedseqs.h"
using namespace std;

template <class Seq>
inline void AlignedSeqs<Seq>::init_seq_mafft_aligned(Seq& seq, string bp) {
	string temp = "  ";
	seq.description = bp.substr(0, bp.find_first_of('_', 0));
	if (bp.size()>75)
		seq.basepairs = bp.substr(bp.rfind(temp) + 2, bp.rfind(' ') - bp.rfind(temp) - 2);
	else
		seq.basepairs = bp;

	seq.xx_primerid_xx = "";
	seq.size = bp.size();
}

template <class Seq>
inline void AlignedSeqs<Seq>::read_mafft_aligned_file(ifstream& file, vector<Seq>& tempgroup)
{
	Seq tempSeq;
	string line;
	string space = " ";

	if (file.good())
		getline(file, line);

	while (file.good()) {
		getline(file, line);
		getline(file, line);
		getline(file, line);

		while (line != "end." && line[0] != space[0]) {
			init_seq_mafft_aligned(tempSeq, line);
			tempgroup.push_back(tempSeq);
			getline(file, line);
		}
		getline(file, line);
		getline(file, line);
		while (line != "end." && line.size() >= 6) {
			string temp = "  ";
			for (int i = 0; i< (int)tempgroup.size(); i++) {
				tempgroup.at(i).basepairs += line.substr(line.rfind(temp) + 2, line.rfind(' ') - line.rfind(temp) - 2);
				getline(file, line);
			}
			getline(file, line);
			getline(file, line);
		}
	}
}

template <class Seq>
inline int AlignedSeqs<Seq>::use_mafft(string filein, string fileout)
{
	string CommandString; 
	CommandString = string("mafft ") + string("--thread 4 --op 0.05 --ep 0.05 --auto --clustalout ") + filein + string(" > ") + fileout;
	//CommandSting += " --op 0.17 --ep 0.07 -maxiterate 10";
	char *command = new char[CommandString.size() + 1];
	command[CommandString.size()] = 0;
	memcpy(command, CommandString.c_str(), CommandString.size());
	return system(command);
}
