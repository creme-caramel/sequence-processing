#ifndef TYPES_H
#define TYPES_H

#include <string>
using namespace std;

// Roche 454
typedef struct Seq_A {
	string description;
	string basepairs;

	string xx_primerid_xx; // form of **[primerID]**
	long size;

	bool operator< (const Seq_A &other) const
	{
		if (xx_primerid_xx != other.xx_primerid_xx)
			return xx_primerid_xx < other.xx_primerid_xx;
		else
			return description < other.description;
	}
} Seq_A;

typedef struct Mutation {
	int posision;
	int subposition;
	string mutation;
	string group;
	int groupnum;
	int FrequencyOfMutation;
	unsigned long SizeOfGroup;

	Mutation(int pos = 0, int dec = 0, string mut = "N/A", string disc = "N/A", int num = 0, int freq = 0, int size = 0)
	{
		posision = pos;
		subposition = dec;
		mutation = mut;
		group = disc;
		groupnum = num;
		FrequencyOfMutation = freq;
		SizeOfGroup = size;
	}

	bool operator< (const Mutation &other) const
	{
		if (posision != other.posision)
			return posision < other.posision;
		else if (subposition != other.subposition)
			return subposition < other.subposition;
		else if (mutation != other.mutation)
			return mutation < other.mutation;
		else
			return groupnum < other.groupnum;
	}	
} Mutation;

#endif
