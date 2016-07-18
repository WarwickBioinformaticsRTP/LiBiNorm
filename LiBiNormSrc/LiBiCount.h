#ifndef LIBICOUNT_H
#define LIBICOUNT_H

#include "gtfFile.h"

using namespace std;

class region
{
public:
	size_t finish;
	string name;
	char strand;
	vector<region *> overlaps;
	region(	size_t finish,const string & name,char strand):finish(finish),name(name),strand(strand){};
};


class chromosomeData : public multimap<size_t, region> 
{

};


class gtfFileEx : public gtfFile
{
	map<string,chromosomeData> chromData; 
public: 
	void index();
	void outputChromData(const string & filename);

};


class LiBiCount
{

public:
	int main(int argc, char **argv);


};

#endif
