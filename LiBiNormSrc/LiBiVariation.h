#ifndef LIBIVARIATION_H
#define LIBIVARIATION_H

#include "LiBiNorm.h"

class LiBiVariation : protected LiBiNormCore
{
public:
	int main(int argc, char **argv);

private:
	stringEx parameterFilename;

};

#endif
