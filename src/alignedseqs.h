#ifndef ALIGNEDSEQS_H
#define ALIGNEDSEQS_H

#include <vector>
#include <set>
using namespace std;

template <class Seq>
class AlignedSeqs
{
	vector<vector<Seq> > groups; // = alignedXX.aln
	vector<Seq> consensusdna; // = Consensus.fna
	vector<Seq> ConsRef;//reference + uber consensus
	vector<vector<Seq> > ConsGroups; // = consalignedXX.aln
	set<Mutation> OutputData;//this is a giant array of finished data.

	// from io.h
	friend int open_infile(ifstream&, const char*);
	friend int open_infile(ifstream&, string&);
	friend int open_outfile(ofstream&, const char*);
	friend int open_outfile(ofstream&, string&);

	// defined in alignedseqs.cpp
	void find_group_consensus(ifstream&, int);
	void find_uber_consensus(ifstream&);
	void write_group_consensus(ofstream&);
	void write_group_with_consensus(ofstream&, int);
	void write_uber_consensus(ofstream&);
	void read_group_with_consensus(ifstream&);
	void group_mut_hunt(int);

	// defined in alignedseqs_mafft.cpp
	void init_seq_mafft_aligned(Seq&, string bp = "none_      GARBIGE ");
	void read_mafft_aligned_file(ifstream&, vector<Seq>&);
	int use_mafft(string, string);
public:
	AlignedSeqs();
	~AlignedSeqs();
	int align_each_group_consenses(const int, const string&);
	int read_reference(const char*);
	int align_uber_consenses(const int, const string&);
	int write_mutations_to_file(const string&);
};
#include "alignedseqs.cpp"
#include "alignedseqs_mafft.cpp"

#endif



