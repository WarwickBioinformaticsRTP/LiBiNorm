#pragma once

#include "api/BamReader.h"
#include "libCommon.h"
#include "stringEx.h"

using namespace std;
using namespace BamTools;

class bamRead
{
public:
	bamRead(const BamAlignment & ba);
	void output(FILE * f);

	string readSeq,qualData,name; 
};

class MakeFastq
{
public:
	MakeFastq(void);

	int main(int argc, char **argv);
};

