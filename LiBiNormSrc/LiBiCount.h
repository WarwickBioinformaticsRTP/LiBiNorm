#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "gtfFile.h"

using namespace std;

class region
{
public:
	size_t finish;
	vector<region *> overlaps;
	string name;
	region(	size_t finish,const string name):finish(finish),name(name){};

};


class chromosomeData : public multimap<size_t, region> 
{

};


class gtfFileEx : public gtfFile
{
	map<string,chromosomeData> chromData; 

public: 
	void index(const string & feature,const string & attribute);



};


class LiBiCount
{

public:
	int main(int argc, char **argv);


};

#endif
