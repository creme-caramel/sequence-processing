/*
 * For the grouping stages:
 *
 * This is an abstract class
 * which is common to all datatypes.
 * 
 * Datatype-specific member functions will
 * inherit this class
 * and name as 'groupedseqs_a', 'groupedseqs_b' ... etc.
 */

#ifndef GROUPEDSEQS_H
#define GROUPEDSEQS_H

#include "io.h"
#include "types.h"
#include <vector>
#include <sstream>
#include <set>
#include <list>
#include <algorithm>
#include <iostream>
using namespace std;

template <typename Seq>
class GroupedSeqs {
protected:
	set<Seq> seqs_for_groups;
	vector<Seq> seqs_too_short;
	vector<Seq> seqs_ungrouped;
	vector<Seq> seqs_mutated_id;
	list<int> sizes;
	int numgroups;

	friend int open_infile(ifstream&, const char*);
	friend int open_outfile(ofstream&, string&);
	void init_seq_parsed(Seq&, string, string, string);

	void writegroupsize(ofstream& file) {
		int intsize;
		int intcount = 1;
		sizes.sort();
	
		list<int>::iterator its;
		its = sizes.begin(); 
		intsize = *its;
	
		file << "Group size" << "," << "Number of groups" << endl;
		for (int i = 0; i< (int)sizes.size(); i++) {
			if (intsize == *its)
				intcount++;
			else {
				file << intsize << "," << intcount << endl; 
				intcount = 1;
				intsize = *its;
			}
			its++;
		}
		file << intsize << "," << intcount << endl; 
	}

	inline void write_seq_vector_to_fna(ostream& file, vector<Seq>& vec) {
		Seq tempSeq1;
		for (int i = 0; i< (int)vec.size(); i++) {
			tempSeq1 = vec.at(i);
			file << tempSeq1.description << endl << tempSeq1.basepairs << endl;
		}
	}

public:
	GroupedSeqs() : numgroups(0) {}
	~GroupedSeqs() {}
	void read_file_to_seqs(const char *, const string);

	/*
	 * also erases defined seqs from the seqs_for_groups
	 */

	inline int define_and_write_each_group(const string& path) {
		vector<Seq> seqs_for_groups_copy;
		Seq tempseq1, tempseq2;
		typename set<Seq>::iterator it;
		it = seqs_for_groups.begin();
		tempseq1 = *it;
		tempseq2 = *it;
	
		do {
			seqs_for_groups_copy.push_back(tempseq2);
			seqs_for_groups.erase(it);
			if (!seqs_for_groups.empty()) {
				it = seqs_for_groups.begin();
				tempseq2 = *it;
			}
		} while (tempseq1.xx_primerid_xx == tempseq2.xx_primerid_xx && (!seqs_for_groups.empty()));
	
		int numb;
		numb = seqs_for_groups_copy.size();
		sizes.push_back(numb);

		ofstream file;
		if (seqs_for_groups_copy.size() > 1) {
			numgroups++;
			string filename = path + "group" + to_string(numgroups) + ".fna";
			if (!open_outfile(file, filename)) {
				cout << "unable to open file: " << filename << endl;
				return 0;
			}			
			write_seq_vector_to_fna(file, seqs_for_groups_copy);
			file.close();
		}
		else {
			seqs_ungrouped.push_back(seqs_for_groups_copy.at(0));
		}
		return 1;
	}
	
	int write_other_files(const string& path) {
		ofstream file;
		string filename = path + "groupsizes.txt";
		if (!open_outfile(file, filename)) {
			cout << "unable to open file: " << filename << endl;
			return 0;
		}
		writegroupsize(file);
		file.close();

		filename = path + "ungrouped.fna";
		if (!open_outfile(file, filename)) {
			cout << "unable to open file: " << filename << endl;
			return 0;
		}
		write_seq_vector_to_fna(file, seqs_ungrouped);
		file.close();

		filename = path + "too_short.fna";
		if (!open_outfile(file, filename)) {
			cout << "unable to open file: " << filename << endl;
			return 0;
		}
		write_seq_vector_to_fna(file, seqs_too_short);
		file.close();

		filename = path + "mutated_id.fna";
		if (!open_outfile(file, filename)) {
			cout << "unable to open file: " << filename << endl;
			return 0;
		}
		write_seq_vector_to_fna(file, seqs_mutated_id);
		file.close();

		seqs_ungrouped.clear();
		seqs_too_short.clear();
		seqs_mutated_id.clear();
		return 1;
	}
	
	int num_groups() {
		return numgroups;
	}
	
	bool is_empty() {
		return seqs_for_groups.empty();
	}
};

#endif
