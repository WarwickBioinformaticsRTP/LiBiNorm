
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "parser.h"
#include "transcriptData.h"


void transcriptDataMap::remove_invalid_values()
{
	for (auto & i : *this)
		i.second.remove_invalid_values();
}

int transcriptDataMap::loadData(const string filename)
{
	std::ifstream f;
	f.open(filename);
	if (!f.is_open())
		exitFail("Unable to open file ", filename);

	string buffer,gene,direction;
	int a = 1;
	while (!f.eof())
	{
		//	Need to check for gene and length consistency and that the gene has not appeared before
		transcriptData tempData;

		getline(f,buffer);
		parseTsv(buffer,gene,tempData.length," ",direction,tempData.counts[0].values());
		getline(f,buffer);
		parseTsv(buffer,gene,tempData.length," ",direction,tempData.counts[1].values());
		emplace(gene,move(tempData));
		a++;
	}

	return EXIT_SUCCESS;
}

void transcriptDataMap::transferTo(dataType & mcmcData,size_t maxLength)
//void dataType::consolidateWith(transcriptDataMap & transData,size_t maxLength)
{
	for (auto & gene : *this)
	{
		for (auto & counts : gene.second.counts)
		{
			counts.selectAtMost(maxLength);

			size_t insertSize = max<size_t>(maxLength,counts.size());
			mcmcData[0].append(counts);
			mcmcData[1].append(counts.size(),gene.second.length);
			//Still need data[2]
		}

	}

}

