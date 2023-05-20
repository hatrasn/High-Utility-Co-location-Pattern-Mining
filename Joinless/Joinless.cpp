
/**
* This code implements the paper "Jin Soung Yoo and Shashi Shekhar. 2006. A joinless approach for mining spatial colocation patterns. 18, 10 (2006), 1323–1337. https://doi.org/10.1109/TKDE.2006.150"
* Including:
* 1) size k table instances are validated by size (k-1) table instances
* 2) to reduce memory, do not product all coarse row instances, each generate one coarse row instances, use size (k-1) to check the row instance is or not a real clique
*
* Writer: vanha tran
* Date: 1.8.2020
*/


#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

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
	clock_t timeGenerateCandidate; // time for generating candidates
	clock_t timeFilStarInstance;
	clock_t timeFilCoarseColoPat;
	clock_t timeFilCliqueInstance; // time for collecting co-location instances
	clock_t timeFilPrevCoLoPat; // time for filtering co-location instances
	clock_t timeWhile;

	double time_taken; // count total execution time
	double totalTimeMatStarNei = 0.0f; // total time for materializing star neighborhoods	
	double totalTimeGenCand = 0.0f; // total time for generating candidates
	double totalTimeFilStarInst = 0.0f;
	double totalTimeFilCoarseCoLoPat = 0.0f;
	double totalTimeFilCliqueInst = 0.0f; // total time for collecting candidates
	double totalTimeFilPrevCoLoPat = 0.0f; // total time for filgering candidates

	startTime = clock(); // Begin count time

	// 2. Set parameters
	double dist_thres = 3200.0f;  // set the distance threshold
	double prev_thres = 0.1f; // set the prevalence
	std::vector<int> topk{ 1, 3, 4, 5, 6 }; // set the number of top-k patterns
	int numtopk = 1; // get the Qc for each size pattern in top of numtopk
	int topkpersize = 2; // get QC of top-k of each size

	// 3. Load data sets
	string file_name = "./Data/quality/gongshan_plant_alphabet_sample_50_checkin.csv";
	freopen("./Data/quality/gongshan_plant_alphabet_sample_50_checkin_distance_3200_PI_01_check_check_n.txt", "w", stdout);
																  
	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File
	// Delete the first line
	dataList.erase(dataList.begin());

	// 4. Count the instance number of each feature
	std::map<char, int> totalInstEachFeat = countNumberInstance(dataList);
	//printInstanceNumber(totalInstEachFeat);

	// 5. Make grids and find neighbors
	// 5.1 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid = makeGrid(dataList, dist_thres);	

	// 5.2 Find neighbors, generate A1: <B.1, C.1>
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN = genStarNeighborhoods(grid, dist_thres);	
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupSNByFeat = groupStarNeighByFeatures(SN);

	// 5.3 Count execution time
	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);
	// Free memory
	SN.clear();

	// 6. Find co-location patterns
	// 6.1 Set Pk, k, Ck
	std::map<vector<char>, float> Pk, PkAll; // Store all prevalent patterns
	std::vector<char> pattern;
	std::vector<vector<char>> Ckplus; // Store all candidate patterns	
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; // Store table isntances of size k patterns
	int k = 2;

	// 6.2 Get size 1 pattern
	for (auto const& pat : totalInstEachFeat)
	{
		pattern.push_back(pat.first);
		Pk.insert({ pattern, 1.0f });
		pattern.clear();
	}

	timeFilPrevCoLoPat = clock();
	totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeFilPrevCoLoPat - timeMaterializeNeighbor);

	int numcan = 0;
	// 6.3 Find size k patterns
	while (Pk.size())
	{
		timeWhile = clock();
		// Generate candidate patterns
		Ckplus = genCandidatePatterns(Pk, k);		
		numcan += Ckplus.size();
		timeGenerateCandidate = clock();
		totalTimeGenCand = totalTimeGenCand + double(timeGenerateCandidate - timeWhile);

		// Select coarse prevalent patterns
		if (k <= 2) // directly get table instances of candidate patterns and calculate their PIs
		{						
			CIk.clear();
			CIk = filterStarInstancesSize2(Ckplus, groupSNByFeat, k);
			timeFilCliqueInstance = clock();
			totalTimeFilCliqueInst = totalTimeFilCliqueInst + double(timeFilCliqueInstance - timeGenerateCandidate);						
			Pk = selectPrevalentPatterns(CIk, totalInstEachFeat, prev_thres, k);						
			PkAll.insert(Pk.begin(), Pk.end());
			timeFilPrevCoLoPat = clock();
			totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeFilPrevCoLoPat - timeFilCliqueInstance);
		}
		else
		{
			// Filter star instances : to store table instances of size (k+1) patterns
			std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> SIkplus
				= filterStarInstances(Ckplus, groupSNByFeat, k);
			timeFilStarInstance = clock();
			totalTimeFilStarInst = totalTimeFilStarInst + double(timeFilStarInstance - timeGenerateCandidate);
			// Use coarse prun			
			selectCoarsePrevalentPatterns(SIkplus, totalInstEachFeat, prev_thres, k);
			timeFilCoarseColoPat = clock();
			totalTimeFilCoarseCoLoPat = totalTimeFilCoarseCoLoPat + double(timeFilCoarseColoPat - timeFilStarInstance);
			// Filter real row instances			
			CIk = filterCliqueInstances2(SIkplus, CIk);
			timeFilCliqueInstance = clock();
			totalTimeFilCliqueInst = totalTimeFilCliqueInst + double(timeFilCliqueInstance - timeFilCoarseColoPat);
			// Select prevalent patterns			
			Pk = selectPrevalentPatterns(CIk, totalInstEachFeat, prev_thres, k);
			// Add Pk to the result PkAll
			PkAll.insert(Pk.begin(), Pk.end());
			timeFilPrevCoLoPat = clock();
			totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeFilPrevCoLoPat - timeFilCliqueInstance);
		}
		// Go to the next iterator
		k += 1;
		Ckplus.clear();
	}

	// 7. Print prevalent patterns	
	std::cout << "The number of canidates: " << numcan << endl;
	std::cout << "The number of prevalent patterns is: " << PkAll.size() << endl;
	// 8. Print running time
	std::cout << "Execution time for materizaling starneighbors is: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating candidates: " << totalTimeGenCand / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering star instances: " << totalTimeFilStarInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering coarse co-locations: " << totalTimeFilCoarseCoLoPat / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering clique instances: " << totalTimeFilCliqueInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering prevalent co-locations: " << totalTimeFilPrevCoLoPat / CLOCKS_PER_SEC << "s." << endl;

	// 9. Calculate time
	endTime = clock();
	time_taken = double(endTime - startTime);
	std::cout << "Time taken by the program is (No neighboring) : " <<
		(totalTimeGenCand + totalTimeFilStarInst + totalTimeFilCoarseCoLoPat + totalTimeFilCliqueInst + totalTimeFilPrevCoLoPat) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 10. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << "(kB)" << std::endl;
	std::cout << "Peak memory usage: " << physPeackMemUsedByMe / 1024 << "(kB)" << std::endl;

	// 11. Count the number patterns of each size
	std::vector<int> maxsize;
	int onek;
	std::map<int, int> sizePats;
	std::map<int, int>::iterator itsizePats;
	std::map<vector<char>, float>::iterator itAllPats = PkAll.begin();
	while (itAllPats != PkAll.end())
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
	std::cout << "Patters sorted by their sizes: \n";
	printPattbySize(sizePats);
	printPrevalentPattern(PkAll);



	// The next code in here



	return 0;
}








