/**
* This version:
* 1) Taking class Object and ObjectWithoutCoordinate as individual .h and .cpp
* 2) Using the concept of "Lizhen Wang. Efficiently Mining High Utility Co-location Patterns from Spatial Data Sets with Instance-Specific Utilities"
* 2) Utility is normalized
* 3) w1 and w2 are calculated by w1 = inter / (inter + intra), w2 = intra / (inter + intra)
* 4) Use clique-hash mining framework
*
*  Note: This version uses unordered_set to save SN and coHM
*	For example:
*   std::unordered_map<ObjWithoutCoord, std::unordered_set<ObjWithoutCoord, myHashFunc>, myHashFunc> SN;
*
* Writer: vanha tran
* Date: 31.12.2021
*/


#include <iostream>
#include <ctime>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>
#include <queue>

// For tree
#include "tree.hh"

// For collecting memory usage
#include <Windows.h>
#include <stdio.h>
#include <Psapi.h>
#pragma comment(lib, "psapi.lib")

// Import my function
#include "Object.h"
#include "printFunctions.h"
#include "readerCSV.h"
#include "toolFunctions.h"


using namespace std;


int main()
{
	// Set time
	clock_t startTime, endTime;
	clock_t timeMaterializeNeighbor;
	clock_t timeBuildTree;
	clock_t timeCalcUitility;
	clock_t timeFilHighPat;

	double time_taken;
	double totalTimeMatStarNei = 0.0f;
	double totalTimeBuildTree = 0.0f;
	double totalTimeCalcUtility = 0.0f;
	double totalTimeFilHighPat = 0.0f;

	// Begin counting time
	startTime = clock();

	// 1. Set parameters
	double dist_thres = 3200.0f;  // set the distance threshold
	double prev_thres = 0.2f; // set the prevalence
	std::vector<int> topk{ 2, 4, 6, 8, 10 }; // set the number of top-k patterns
	int numtopk = 10; // get the Qc for each size pattern in top of numtopk
	int topkpersize = 2; // get QC of top-k of each size

	// 2. Experiment number of instances
	string file_name = "./Data/quality/gongshan_plant_alphabet_sample_50_checkin.csv";
	freopen("./Data/quality/gongshan_plant_alphabet_sample_50_checkin_distance_3200_PI_02.txt", "w", stdout);

	// 2.1 Load the dataset
	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File	
	dataList.erase(dataList.begin()); // Delete the first line that is a title

	// 2.2 Get max and min of utility for nomorlizing	
	std::unordered_map<char, std::vector<float>> minmaxU;
	getMinMaxUtility(dataList, minmaxU);

	// 3. Make grid and find neighbors
	// 3.1 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid;
	makeGrid(dataList, dist_thres, grid);
	//printGrid(grid);
	dataList.clear();
	// 3.2 Find neighbors, generate A1: <B.1, C.1>
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc> SN;
	std::map<ObjWithoutCoord, int> V;
	// 3.3 Store the total utitlity of each feature 
	std::unordered_map<char, float> totalUtilityFeats;
	genStarNeighborhoods(grid, dist_thres, minmaxU, SN, V, totalUtilityFeats);
	// 3.4 Count execution time
	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);
	//printSN(SN);	
	grid.clear();
	minmaxU.clear();

	// 4. Construct I-tree and get H-cliques	
	generateCliques(V, SN);

	timeBuildTree = clock();
	totalTimeBuildTree = double(timeBuildTree - timeMaterializeNeighbor);

	std::cout << "Finish cliques:" << timeBuildTree / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Number of keys: " << coHM.size() << endl;
	
	// 5. Calculate the utility of each pattern		
	std::unordered_map<std::string, float> Pk; // saves all patterns and their utilities
	calculateUtility(totalUtilityFeats, Pk);
	
	timeCalcUitility = clock();
	totalTimeCalcUtility = totalTimeCalcUtility + double(timeCalcUitility - timeBuildTree);

	// 6. Filter high utility patterns
	std::unordered_map<std::string, float> HUPk; // saves all high utility patterns
	filterHighUtilityPatterns(Pk, HUPk, prev_thres);

	std::cout << "The number of patterns: " << HUPk.size() << endl;
	//printPatterns(HUPk);

	timeFilHighPat = clock();
	totalTimeFilHighPat = totalTimeFilHighPat + double(timeFilHighPat - timeCalcUitility);

	// 7. Print the total execution time
	endTime = clock();
	time_taken = double(endTime - startTime);

	std::cout << "Time for meterializing neighbors: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for building instance tree: " << totalTimeBuildTree / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for constructing co-location hashmap: " << totalTimeBuildHash / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for calculating utility: " << totalTimeCalcUtility / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for filtering high utility patterns: " << totalTimeFilHighPat / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time (no materialize neighbors) is : " << (totalTimeBuildTree + totalTimeBuildHash + totalTimeCalcUtility + totalTimeFilHighPat) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 8. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << " (kB)" << std::endl;
	std::cout << "Peack memory usage: " << physPeackMemUsedByMe / 1024 << " (kB)" << std::endl;

	// 9. Count the number of each size
	std::vector<int> maxsize;
	int onek;
	std::map<int, int> sizePats;
	std::map<int, int>::iterator itsizePats;
	std::unordered_map<std::string, float>::iterator itAllPats = HUPk.begin();
	while (itAllPats != HUPk.end())
	{
		onek = itAllPats->first.size();
		maxsize.push_back(onek);

		itsizePats = sizePats.find(onek);
		if (itsizePats != sizePats.end())
		{
			itsizePats->second += 1;
		}
		else
		{
			sizePats.insert({ onek, 1 });
		}
		// Terminate
		++itAllPats;
	}
	int maxSize = *std::max_element(maxsize.begin(), maxsize.end());
	std::cout << "The maximal size of patterns is: " << maxSize << endl;
	std::cout << "Patters by sizes: \n";
	printPattbySize(sizePats);


	// The next code in here



	return 0;
}








