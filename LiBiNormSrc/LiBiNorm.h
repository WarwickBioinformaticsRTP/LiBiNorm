#ifndef LIBINORM_H
#define LIBINORM_H

#include "stringEx.h"
#include "transcriptData.h"

class LiBiNorm
{
	stringEx consFileName;
	transcriptDataMap transData;

	dataType consData;

public:
	int main(int argc, char **argv);
	int loadData();

	
};

#endif


