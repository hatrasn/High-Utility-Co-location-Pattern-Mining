
/* This version implements the paper
* Li, Y., Wang, L., Yang, P. and Li, J., 2021, September. EHUCM: An Efficient Algorithm for Mining High Utility Co-location Patterns from Spatial Datasets with 
* Feature-specific Utilities. In International Conference on Database and Expert Systems Applications (pp. 185-191). Springer, Cham.
* 
*
* Since we address utilities of instances, thus utility of a feature is equal to the average of utility of its instances
*
* Writer: Hatran
* Date: 05/12/2022
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

// For collect memory usage
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
	clock_t timeConstructFONT;

	double time_taken; // count total execution time
	double totalTimeConstructFONT = 0.0;
	double totalTimeMatStarNei = 0.0f; // total time for materializing star neighborhoods
	
	startTime = clock(); // Begin count time

	// 2. Set parameters
	float dist_thres = 200.0f;  // set the distance threshold
	float minutil = 0.2f; // set the prevalence

	// 3. Load data sets
	string file_name = "./Data/LasVegas_x_y_alphabet_version_03_2_sample_50.csv";
	freopen("./Data/LasVegas_x_y_alphabet_version_03_2_sample_50_distance_200_PI_02_new.txt", "w", stdout);

	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File	
	dataList.erase(dataList.begin()); // Delete the first line	 

	// 4. Make grid and find neighbors
	// 4.1 Get maximal utility
	std::unordered_map<char, double> maxU;
	getMaxUtility(dataList, maxU);
	// 4.2 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid = makeGrid(dataList, dist_thres);
	// 4.3 Find neighbors, generate A1: <B.1, C.1>
	std::map<char, double> tuf; // the total utility of each feature tuf = number of instances of f * utility of f
	genStarNeighborhoods(grid, dist_thres, maxU, tuf);
	
	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);

	// 5. Calculate utility of the data set
	calculateUS(tuf);	

	// 6. Find feature-object neighbor sets, FONS
	buildFONS(SN);
	
	// 7. Create FONT
	buildFONT();
	timeConstructFONT = clock();
	totalTimeConstructFONT = double(timeConstructFONT - timeMaterializeNeighbor);

	// 8. Sort features by their utility	
	sortFeatureByUtility(tuf);

	// 9. Initialize
	std::vector<std::vector<char>> c; // a candidate
	int k = 1; // the index of feature in sortF
	for (auto f : sortF)
	{
		std::vector<char> one{ f };
		c.push_back(one);
	}
	// 10. Search size k patterns
	searchHUCP(c, k, minutil);	
	
	std::cout << "The number of HUCPs: " << PkAll.size() << endl;
	//printPrevalentPattern();
	
	//11. Print the running time
	std::cout << "Execution time for materizaling starneighbors is: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for building FONT is: " << totalTimeConstructFONT / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating candidate is: " << totalGenCand / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for generating participating instances is: " << (totalGenerateParticipatingObjects - totalTimeCalUPIAndFilPrevCoLoPat) / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for calculating UPI and filtering patterns: " << totalTimeCalUPIAndFilPrevCoLoPat / CLOCKS_PER_SEC << "s." << endl;
	
	endTime = clock();
	time_taken = double(endTime - startTime);
	std::cout << "Time taken by the program is (No neighboring) : " << (time_taken - totalTimeMatStarNei) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 12. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << "(kB)" << std::endl;
	std::cout << "Peack memory usage: " << physPeackMemUsedByMe / 1024 << "(kB)" << std::endl;

	// 13. Get the maximal size of patterns
	std::vector<int> sizePats;
	std::map<std::string, double>::iterator itAllPats = PkAll.begin();
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








