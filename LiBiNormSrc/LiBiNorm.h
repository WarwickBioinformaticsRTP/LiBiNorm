#ifndef LIBINORM_H
#define LIBINORM_H

#include "LiBiNormCore.h"

class LiBiNorm : protected LiBiNormCore
{
public:
	LiBiNorm()  {};



	int main(int argc, char **argv);

private:
	stringEx landscapeFilename;

};

#endif


