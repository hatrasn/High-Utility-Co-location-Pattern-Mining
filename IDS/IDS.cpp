
/* This code implements the algorithm IDS that is published in 2019
* Bao, X. and Wang, L., 2019. A clique-based approach for co-location pattern mining. Information Sciences, 490, pp.244-264.
* 
* Writer: vanha tran
* Date: 1.9.2020
*/


#include <iostream>
#include <ctime>
#include <map>
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
	// 1. Set time
	clock_t startTime, endTime; // start and end time
	clock_t timeMaterializeNeighbor; // time for materializing star neighborhoods
	clock_t timeGenerateCliques;	
	clock_t timeFilPCPs;

	double time_taken; // count the total execution time
	double totalTimeMatStarNei = 0.0f; // total time for materializing star neighborhoods
	double totalTimeGenerateCliques = 0.0f;	
	double totalTimeFilPCPs = 0.0f;

	startTime = clock(); // Begin count time

	// 2. Set parameters
	float dist_thres = 30.0f;  // set the distance threshold
	float prev_thres = 0.2f; // set the prevalence

	// 3. Load data sets
	string file_name = "./Data/5k_15f_50k_ex_checkin_3.csv";
	freopen("./Data/5k_15f_50k_ex_checkin_3_distance_30_PI_02_new.txt", "w", stdout);	
	
	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File	
	dataList.erase(dataList.begin()); // Delete the first line

	// 4. Count the instance number of each feature
	countNumberInstance(dataList);	

	// 5. Make grid and find neighbors
	// 5.1 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid = makeGrid(dataList, dist_thres);	
	// 5.2 Find neighbors, generate A1: <B.1, C.1>
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN = genStarNeighborhoods(grid, dist_thres);

	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);		

	// 6. Construct I-tree, get H-cliques and construct the chash structure 
	generateCliques2(SN);
	timeGenerateCliques = clock();
	totalTimeGenerateCliques = double(timeGenerateCliques - timeMaterializeNeighbor);

	// 7. Calculate PIs and filter prevalent patterns
	timeFilPCPs = clock();
	std::map<std::string, float> PkAll; // store all prevalent patterns
	PkAll = CalPIandFilterPatterns2(prev_thres);
	totalTimeFilPCPs = double(clock() - timeFilPCPs) - totalTimeCalPIs;
	std::cout << "The number of prevalent pattern is: " << PkAll.size() << endl;
	
	// 8. Print the running time
	std::cout << "Execution time for materizaling starneighbors: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating cliques: " << (totalTimeGenerateCliques - totalTimeCreateChash) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for creating chash structure: " << totalTimeCreateChash / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for calculating PIs: " << totalTimeCalPIs / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering PCPs: " << totalTimeFilPCPs / CLOCKS_PER_SEC << "s." << endl;

	endTime = clock();
	time_taken = double(endTime - startTime);
	std::cout << "Time taken by the program is (no neighboring) : " << (time_taken - totalTimeMatStarNei) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 9. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << "(kB)" << std::endl;
	std::cout << "Peack memory usage: " << physPeackMemUsedByMe / 1024 << "(kB)" << std::endl;

	// 10. Get the maximal size of patterns
	std::vector<int> sizePats;
	std::map<std::string, float>::iterator itAllPats = PkAll.begin();
	while (itAllPats != PkAll.end())
	{
		sizePats.push_back(itAllPats->first.size());
		++itAllPats;
	}
	int maxSize = *std::max_element(sizePats.begin(), sizePats.end());
	std::cout << "The maximal size of patterns is: " << maxSize << endl;


	// The next code in here


	return 0;
}








