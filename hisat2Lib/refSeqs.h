#pragma once

#include <string>
#include <map>
#include "FeatureFileEx.h"

class refSeqs
{
public:
	bool get(const std::string& ebwtFileBase, std::map<std::string, geneData> & genes);

private:
	std::vector<std::string> refnames;
	std::map<std::string, size_t> nameToId;

};
