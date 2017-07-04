#include "io.h"
#include "types.h"
#include "groupedseqs_a.h"
#include "alignedseqs.h"
#include <assert.h>
using namespace std;

const char * reference = "Dloopref.fna";
const string mid_to_string(const char *num_str);

int main(int argc, const char* argv[])
{
	assert(argc == 3);

	typedef struct {
		string filepath;
		string filename;
	} InFiles;

	InFiles f;
	string in(argv[1]);

	size_t pos = in.find_last_of("\\/");
	if (string::npos == pos) {
		f.filepath = "";
		f.filename = in;
	}
	else {
		f.filepath = in.substr(0, pos + 1);
		f.filename = in.erase(0, pos + 1);
	}

	GroupedSeqs_A gs;
	AlignedSeqs<Seq_A> as;

	assert(gs.read_file_to_seqs(argv[1], mid_to_string(argv[2]))); // eg. ("/path/to/data.fna", "15")
	while (!gs.is_empty())
		assert(gs.define_and_write_each_group(f.filepath));
	assert(gs.write_other_files(f.filepath));
	assert(as.align_each_group_consenses(gs.num_groups(), f.filepath));
	assert(as.read_reference(reference)); //  --> early half of consref
	assert(as.align_uber_consenses(gs.num_groups(), f.filepath));
	assert(as.write_mutations_to_file(f.filepath)); // --> output.txt
	return 1;
}

const string mid_to_string(const char *num_str)
{
	int num = atoi(num_str);
	string mid;
	switch (num) {
	case 1:
		mid = "ACGAGTGCGT"; break;
	case 3:
		mid = "AGACGCACTC"; break;
	case 4:
		mid = "AGCACTGTAG"; break;
	case 5:
		mid = "ATCAGACACG"; break;
	case 6:
		mid = "ATATCGCGAG"; break;
	case 7:
		mid = "CGTGTCTCTA"; break;
	case 10:
		mid = "TCTCTATGCG"; break;
	case 13:
		mid = "CATAGTAGTG"; break;
	case 14:
		mid = "CGAGAGATAC"; break;
	case 15:
		mid = "ATACGACGTA"; break;
	case 16:
		mid = "TCACGTACTA"; break;
	case 17:
		mid = "CGTCTAGTAC"; break;
	case 19:
		mid = "TGTACTACTC"; break;
	case 20:
		mid = "ACGACTACAG"; break;
	case 21:
		mid = "CGTAGACTAG"; break;
	case 22:
		mid = "TACGAGTATG"; break;
	case 24:
		mid = "TAGAGACGAG"; break;
	case 30:
		mid = "AGACTATACT"; break;
	case 33:
		mid = "ATAGAGTACT"; break;
	case 36:
		mid = "CGACGTGACT"; break;
	case 39:
		mid = "TACAGATCGT"; break;
	case 40:
		mid = "TACGCTGTCT"; break;
	default:
		mid = "NNNNNNNNNN"; break;
	}
	return mid;
}
