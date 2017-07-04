#ifndef IO_H
#define IO_H

#include <fstream>
#include <string>
#include <cstring>
using namespace std;

inline int open_infile(ifstream& file, const char* to_read) {
	file.open(to_read);
	return file.is_open();
}

inline int open_infile(ifstream& file, string& to_read) {
	char *a = new char[to_read.size() + 1];
	a[to_read.size()] = 0;
	memcpy(a, to_read.c_str(), to_read.size());
	file.open(a);
	delete [] a;
	return file.is_open();
}

inline int open_outfile(ofstream& file, const char* to_write) {
	file.open(to_write);
	return file.is_open();
}

inline int open_outfile(ofstream& file, string& to_write) {
	char *a = new char[to_write.size() + 1];
	a[to_write.size()] = 0;
	memcpy(a, to_write.c_str(), to_write.size());
	file.open(a);
	delete [] a;
	return file.is_open();
}

#endif
