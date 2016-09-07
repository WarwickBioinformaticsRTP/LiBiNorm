#pragma once

#include <vector>
#include "api/BamReader.h"
#include "libCommon.h"
#include "stringEx.h"

using namespace std;
using namespace BamTools;

class bamRead
{
public:
	bamRead(){};
	bamRead(const BamAlignment & ba);
	bamRead & operator = (const BamAlignment & ba);
	void output(FILE * f);

	string readSeq,qualData;
	stringEx name; 
};

class bamReadCache : public vector<bamRead>
{
public:
	bamReadCache(size_t size):vector<bamRead>(size){};
	void clear();
	bamRead & operator[](size_t i){ return vector<bamRead>::operator[](i);};
};

class MakeFastq
{
public:
	MakeFastq(void);

	int main(int argc, char **argv);
};

