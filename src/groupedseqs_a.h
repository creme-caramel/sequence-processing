#ifndef GROUPEDSEQS_A_H
#define GROUPEDSEQS_A_H

#include "groupedseqs.h"
using namespace std;

class GroupedSeqs_A : public GroupedSeqs<Seq_A> {
protected:
	inline void init_seq_parsed(Seq_A&, string, string, string);
public:
	GroupedSeqs_A();
	~GroupedSeqs_A();
	int read_file_to_seqs(const char *, const string);
};

#endif
