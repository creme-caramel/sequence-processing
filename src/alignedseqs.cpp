#include "alignedseqs.h"
#include <iostream>
#include <stdlib.h>
#include <assert.h>

/*
 * group consensus
 */

template <class Seq>
inline void AlignedSeqs<Seq>::find_group_consensus(ifstream& file, int group)
{
	vector<Seq> tempgroup;
	read_mafft_aligned_file(file, tempgroup);

	string tempConStr = "";
	int minforcons = (int)(tempgroup.size() * .5 + 0.5);

	for (int j = 0; j < (int)tempgroup.at(0).basepairs.size(); j++) {
		int dash = 0;
		int A = 0;
		int T = 0;
		int G = 0;
		int C = 0;
		for (int k = 0; k< (int)tempgroup.size(); k++) {
			char tempchar = tempgroup.at(k).basepairs.at(j);
			if (tempchar == 'A' || tempchar == 'a')
				A++;
			else if (tempchar == 'T' || tempchar == 't')
				T++;
			else if (tempchar == 'G' || tempchar == 'g')
				G++;
			else if (tempchar == 'C' || tempchar == 'c')
				C++;
			else
				dash++;
		}
		if (A >= minforcons)
			tempConStr += "A";
		else if (T >= minforcons)
			tempConStr += "T";
		else if (G >= minforcons)
			tempConStr += "G";
		else if (C >= minforcons)
			tempConStr += "C";
		else
			tempConStr += "-";
	}
	Seq tempCons;
	tempCons.basepairs = tempConStr;
	tempCons.description = to_string(group);
	groups.push_back(tempgroup);
	consensusdna.push_back(tempCons);
}

template <class Seq>
void AlignedSeqs<Seq>::write_group_consensus(ofstream& file)
{
	Seq tempSeq;
	for (int i = 0; i < (int)consensusdna.size(); i++) {
		tempSeq = consensusdna.at(i);
		file << ">" << tempSeq.description << "\n" << tempSeq.basepairs << "\n";
	}
}

/*
 * uber consensus
 */

template <class Seq>
void AlignedSeqs<Seq>::find_uber_consensus(ifstream& file)
{
	vector<Seq> tempgroup;
	read_mafft_aligned_file(file, tempgroup);

	string tempConStr = "";
	int minforcons = (int)(tempgroup.size()*.49 + 0.5);

	for (int j = 0; j< (int)tempgroup.at(0).basepairs.size(); j++) {
		int dash = 0;
		int A = 0;
		int T = 0;
		int G = 0;
		int C = 0;
		char tempchar;
		for (int k = 0; k<(int)tempgroup.size(); k++) {
			if ((int)tempgroup.at(k).basepairs.size()>j)
				tempchar = tempgroup.at(k).basepairs.at(j);
			else
				tempchar = '-';

			if (tempchar == 'A' || tempchar == 'a')
				A++;
			else if (tempchar == 'T' || tempchar == 't')
				T++;
			else if (tempchar == 'G' || tempchar == 'g')
				G++;
			else if (tempchar == 'C' || tempchar == 'c')
				C++;
			else
				dash++;
		}
		if (A >= minforcons)
			tempConStr += "A";
		if (T >= minforcons)
			tempConStr += "T";
		if (G >= minforcons)
			tempConStr += "G";
		if (C >= minforcons)
			tempConStr += "C";
		else
			tempConStr += "-";
	}
	Seq tempCons;
	tempCons.basepairs = tempConStr;
	tempCons.description = "UberConsensus";
	ConsRef.push_back(tempCons);
}

template <class Seq>
void AlignedSeqs<Seq>::write_group_with_consensus(ofstream& file, int groupnum)
{
	Seq tempSeq;
	// assert((int)groups.size() >= group);

	for (int i = 0; i < (int)ConsRef.size(); i++) {
		tempSeq = ConsRef.at(i);
		file << ">" << tempSeq.description << "\n" << tempSeq.basepairs << "\n";
	}

	file << ">" << consensusdna.at(groupnum).description << "\n" << consensusdna.at(groupnum).basepairs << "\n";
	
	for (int i = 0; i < (int)groups.at(groupnum).size(); i++) {
		tempSeq = groups.at(groupnum).at(i);
		file << ">" << tempSeq.description << "\n" << tempSeq.basepairs << "\n";
	}
}

template <class Seq>
void AlignedSeqs<Seq>::write_uber_consensus(ofstream& file)
{
	Seq tempSeq;
	for (int i = 0; i< (int)ConsRef.size(); i++) {
		tempSeq = ConsRef.at(i);
		file << ">" << tempSeq.description << "\n" << tempSeq.basepairs << "\n";
	}
}

template <class Seq>
void AlignedSeqs<Seq>::read_group_with_consensus(ifstream& file)
{
	vector<Seq> tempgroup;
	read_mafft_aligned_file(file, tempgroup);
	ConsGroups.push_back(tempgroup);
}

template <class Seq>
void AlignedSeqs<Seq>::group_mut_hunt(int group)
{
	int truepos = 16534;
	int truedes = 0;
	int tempcount = 0;
	string tempmutation;
	string tempgroup;
	int counting = 1;

	for (int i = 0; i< (int)ConsGroups.at(group).at(0).basepairs.size(); i++)
	{
		counting += 1;
		if (ConsGroups.at(group).at(0).basepairs[i] != '-')
		{
			truepos++;
			if (truepos>16569)
			{
				truepos = 1;
			}
			truedes = 0;
		}
		else if (ConsGroups.at(group).at(1).basepairs[i] != '-' || ConsGroups.at(group).at(2).basepairs[i] != '-')
		{
			truedes++;
		}
		if (ConsGroups.at(group).at(1).basepairs[i] != ConsGroups.at(group).at(2).basepairs[i] && truepos != 16534)
		{
			if (ConsGroups.at(group).at(1).basepairs[i] != '-')
			{
				if (ConsGroups.at(group).at(2).basepairs[i] != '-')
				{
					for (int mutit = 3; mutit < (int)ConsGroups.at(group).size(); mutit++)
					{
						if (ConsGroups.at(group).at(2).basepairs[i] == ConsGroups.at(group).at(mutit).basepairs[i])
						{
							tempcount++;
						}
					}
					tempmutation = ConsGroups.at(group).at(1).basepairs[i];
					tempmutation += " to ";
					tempmutation += ConsGroups.at(group).at(2).basepairs[i];
					tempgroup = ConsGroups.at(group).at(2).basepairs.substr(35, 16);
					Mutation tempdata(truepos, truedes, tempmutation, tempgroup, (group + 1), tempcount, (ConsGroups.at(group).size() - 3));

cout << tempdata.posision << "\t" << tempdata.subposition << "\t" << tempdata.mutation << "\t" << tempdata.groupnum << "\t" << tempdata.FrequencyOfMutation << "\t" << tempdata.SizeOfGroup << endl;
					OutputData.insert(tempdata);
					tempcount = 0;
				}
				else
				{
					//deletion
					for (int mutit = 3; mutit < (int)ConsGroups.at(group).size(); mutit++)
					{
						if (ConsGroups.at(group).at(2).basepairs[i] == ConsGroups.at(group).at(mutit).basepairs[i])
						{
							tempcount++;
						}
					}
					tempmutation = ConsGroups.at(group).at(1).basepairs[i];
					tempmutation += " to :";
					tempgroup = ConsGroups.at(group).at(2).basepairs.substr(35, 16);
					Mutation tempdata(truepos, truedes, tempmutation, tempgroup, (group + 1), tempcount, ConsGroups.at(group).size() - 3);

cout << tempdata.posision << "\t" << tempdata.subposition << "\t" << tempdata.mutation << "\t" << tempdata.groupnum << "\t" << tempdata.FrequencyOfMutation << "\t" << tempdata.SizeOfGroup << endl;
					OutputData.insert(tempdata);
					tempcount = 0;
				}
			}
			else
			{
				//insertion
				for (int mutit = 3; mutit< (int)ConsGroups.at(group).size(); mutit++)
				{
					if (ConsGroups.at(group).at(2).basepairs[i] == ConsGroups.at(group).at(mutit).basepairs[i])
					{
						tempcount++;
					}
				}
				tempmutation = ": to ";
				tempmutation += ConsGroups.at(group).at(2).basepairs[i];
				tempgroup = ConsGroups.at(group).at(2).basepairs.substr(35, 16);
				Mutation tempdata(truepos, truedes, tempmutation, tempgroup, (group + 1), tempcount, ConsGroups.at(group).size() - 3);

cout << tempdata.posision << "\t" << tempdata.subposition << "\t" << tempdata.mutation << "\t" << tempdata.groupnum << "\t" << tempdata.FrequencyOfMutation << "\t" << tempdata.SizeOfGroup << endl;
				OutputData.insert(tempdata);
				tempcount = 0;
			}
		}
	}
}

template <class Seq>
int AlignedSeqs<Seq>::align_each_group_consenses(const int numgroups, const string& path)
{
	string fna, aln;
	ifstream infile;
	ofstream outfile;

	for (int i = 1; i <= numgroups; i++) {
		fna = path + "group" + to_string(i) + ".fna";
		aln = path + "aligned" + to_string(i) + ".aln";
		use_mafft(fna, aln);
		if (!open_infile(infile, aln)) {
			cout << "unable to open file: " << aln << endl;
			return 0;
		}
		find_group_consensus(infile, i);
		infile.close();
	}
	fna = path + "Consensus.fna";
	if (!open_outfile(outfile, fna)) {
		cout << "unable to open file: " << fna << endl;
		return 0;
	}
	write_group_consensus(outfile); // --> Consensus.fna
	outfile.close();
	return 1;
}

template <class Seq>
int AlignedSeqs<Seq>::read_reference(const char *filename)
{
	ifstream file;
	if (!open_infile(file, filename)) {
		cout << "unable to open file: " << filename << endl;
		return 0;
	}

	string line, line2, line3;
	int breakout;
	string greaterthen = ">";
	Seq tempSeq;

	while (file.good()) {
		if (line3[0] == greaterthen[0])
			line = line3;
		else
			getline(file, line);

		if (line != "end.") {
			if (line[0] == greaterthen[0]) {
				getline(file, line2);
				if (line2 != "end.") {
					breakout = 0;
					do {
						breakout++;
						getline(file, line3);

						if (line3[0] != greaterthen[0] && line3 != "end.") {
							line2 += line3;
						}
					} while (line3[0] != greaterthen[0] && line3 != "end." && breakout<100);
					tempSeq.description = line;
					tempSeq.basepairs = line2;
					ConsRef.push_back(tempSeq);
				}
			}
		}
	}
	file.close();
	return 1;
}

template <class Seq>
int AlignedSeqs<Seq>::align_uber_consenses(const int numgroups, const string& path)
{
	string fna, aln;
	ifstream infile;
	ofstream outfile;

	fna = path + "Consensus.fna";
	aln= path + "Consensus.aln";

	use_mafft(fna, aln);

	if (!open_infile(infile, aln)) {
		cout << "unable to open file: " << aln << endl;
		return 0;
	}
	find_uber_consensus(infile);
	infile.close();

	for (int j = 1; j <= numgroups; j++) {
		fna= path + "consgroup" + to_string(j) + ".fna";
		aln = path + "consaligned" + to_string(j) + ".aln";

		if (!open_outfile(outfile, fna)) {
			cout << "unable to open file: " << fna << endl;
			return 0;
		}
		write_group_with_consensus(outfile, j - 1); // --> consgroupXX.fna
		outfile.close();

		use_mafft(fna, aln); // --> consalignedXX.aln

		if (!open_infile(infile, aln)) {
			cout << "unable to open file: " << aln << endl;
			return 0;
		}
		read_group_with_consensus(infile);
		infile.close();

		group_mut_hunt(j - 1);
	}

	fna = path + "UberConsensus.fna";
	if (!open_outfile(outfile, fna)) {
		cout << "unable to open file: " << fna << endl;
		return 0;
	}
	write_uber_consensus(outfile);
	outfile.close();
	return 1;
}

template <class Seq>
int AlignedSeqs<Seq>::write_mutations_to_file(const string& path)
{
	ofstream file;
	string filename = path + "output.txt";
	if (!open_outfile(file, filename)) {
		cout << "unable to open file: " << filename << endl;
		return 0;
	}

	Mutation temp1, temp2, temp3;
		file << "Mitomap consensusdna position,Mutated to,Template number,Mutation frequence in template consensus seq,# seq that make up template consensus seq,Template number,Mutation frequence in template consensus seq,# seq that make up template consensus seq,Template number,Mutation frequence in template consensus seq,# seq that make up template consensus seq\n";
		temp1 = *OutputData.begin();

		if (temp1.subposition<10) {
			file << "\n" << to_string(temp1.posision) << "." << "00" << to_string(temp1.subposition) << "," << temp1.mutation << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
		}
		else if (temp1.subposition<100) {
			file << "\n" << to_string(temp1.posision) << "." << "0" << to_string(temp1.subposition) << "," << temp1.mutation << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
		}
		else {
			file << "\n" << to_string(temp1.posision) << "." << to_string(temp1.subposition) << "," << temp1.mutation << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
		}

		for (set<Mutation>::iterator it = (++OutputData.begin()); it != OutputData.end(); ++it) {
			temp3 = *it;
			if ((double)temp3.FrequencyOfMutation / (double)temp3.SizeOfGroup > 0.5)
			{
				temp2 = temp1;
				temp1 = temp3;
				if (temp1.posision == temp2.posision && temp1.subposition == temp2.subposition && temp1.mutation == temp2.mutation) {
					file << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
				}
				else {
					if (temp1.subposition<10) {
						file << "\n" << to_string(temp1.posision) << "." << "00" << to_string(temp1.subposition) << "," << temp1.mutation << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
					}
					else if (temp1.subposition<100) {
						file << "\n" << to_string(temp1.posision) << "." << "0" << to_string(temp1.subposition) << "," << temp1.mutation << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
					}
					else {
						file << "\n" << to_string(temp1.posision) << "." << to_string(temp1.subposition) << "," << temp1.mutation << "," << to_string(temp1.groupnum) << "," << to_string(temp1.FrequencyOfMutation) << "," << to_string(temp1.SizeOfGroup);
					}
				}
			}
		}
	cout << "done" << endl;
	file.close();
	return 1;
}

template <class Seq>
AlignedSeqs<Seq>::AlignedSeqs() {}

template <class Seq>
AlignedSeqs<Seq>::~AlignedSeqs()
{
	groups.clear();
	consensusdna.clear();
	ConsRef.clear();
	OutputData.clear();
}
