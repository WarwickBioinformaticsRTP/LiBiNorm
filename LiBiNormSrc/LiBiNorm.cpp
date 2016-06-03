#include "libCommon.h"
#include "LiBiNorm.h"

using namespace std;

int main(int argc, char **argv)
{
	LiBiNorm lbn;
	return lbn.main(argc,argv);
}


int LiBiNorm::main(int argc, char **argv)
{


	if(argc < 1)
	{
		printf("Error: parameter wrong!\n");
		return EXIT_FAILURE;
	}
	else if(argc == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     LiBiNorm:    RNA-seq library bias normalisation   \n");
		return EXIT_SUCCESS;
	}

	int ni = 1;
	while(ni < argc)
	{
		if(strcmp(argv[ni], "-c") == 0)
		{
			consFileName = argv[++ni];
		}
		else
		{
			cout << "Invalid parameter";
		}
		ni++;
	}

	transData.loadData(consFileName);
	transData.remove_invalid_values();
	consData.consolidateWith(transData,100);



	cout << "Data loaded";
	string x;
	cin >> x;

	return EXIT_SUCCESS;
}



