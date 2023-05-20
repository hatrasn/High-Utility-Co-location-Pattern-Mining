
/* This code implements the algorithm that has been published in
* 
* LizhenWang,Wanguo Jiang, and Hongmei Chen. 2017. Efficiently mining high utility co-location patterns from spatial
* data sets with instance-specific utilities. In Proc. Int. Conf. Database Syst. Adv. Appl. 458–474. https://doi.org/10.1007/978-3-319-55699-4_28
* 
* Objective: Discovering all high utility co-location patterns
* 
* Writer: vanha tran
* Date: 1.12.2020
*/


#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

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
	clock_t timeGenerateCandidate; // time for generating candidates
	clock_t timeFilStarInstance;
	clock_t timeFilCliqueInstance; // time for collecting co-location instances
	clock_t timeComputeUPI;
	clock_t timeSelectUPI; // time for filtering co-location instances	
	clock_t timeWhile;

	double time_taken; // count total execution time
	double totalTimeMatStarNei = 0.0f; // total time for materializing star neighborhoods	
	double totalTimeComputeUPI = 0.0f;
	double totalTimeGenCand = 0.0f; // total time for generating candidates
	double totalTimeFilStarInst = 0.0f;
	double totalTimeCollectTableInstance = 0.0f; // total time for collecting candidates
	double totalTimeSelectUPI = 0.0f; // total time for filgering candidates

	startTime = clock(); // Begin count time

	// 2. Set parameters	
	double w1 = 0.6f;
	double w2 = 0.4f;
	double dist_thres = 200.0;  // set the distance threshold
	double prev_thres = 0.2; // set the prevalence
	std::vector<int> topk{ 2, 4, 6, 8, 10 }; // set the number of top-k patterns
	int numtopk = 10; // get the Qc for each size pattern in top of numtopk
	int topkpersize = 2; // get QC of top-k of each size

	// 3. Load data sets
	string file_name = "./Data/instance_real/Beijing_POI_type_0_Alpha_2_sample_35_exponent_distr_sample_20_exponent_distr.csv";
	freopen("./Data/instance_real/Beijing_POI_type_0_Alpha_2_sample_35_exponent_distr_sample_20_exponent_dist_distance_200_PI_02_log_prun.txt", "w", stdout);

	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File	
	dataList.erase(dataList.begin()); // Delete the first line

	// 4. Count the instance number of each feature
	std::map<char, float> totalUtilityEachFeat = computeUtility(dataList);

	// 5. Make grid and find neighbors
	// 5.1 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid = makeGrid(dataList, dist_thres);
	// 5.2 Find neighbors, generate A1: <B.1, C.1>
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN = genStarNeighborhoods(grid, dist_thres);
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupSNByFeat = groupStarNeighByFeatures(SN);	
	SN.clear(); // free memory
	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);
	
	// 6. Find patterns	
	std::map<std::vector<char>, float> Pk, PkAll; // size k prevalent patterns	
	std::vector<std::pair<std::vector<char>, std::vector<char>>> NonPk; // save all non-prevalent patterns
	std::vector<char> pattern;	
	std::vector<std::vector<char>> Ckplus; // size k+1 candidates	
	std::map<std::vector<char>, std::vector<float>> CTIk; // store the PRs of each candidate	
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; // store table isntances of size k patterns
	int k = 2;
	// 6.1 Get size 1 pattern
	for (auto const& pat : totalUtilityEachFeat)
	{
		pattern.push_back(pat.first);
		Pk.insert({ pattern, 1.0f });
		pattern.clear();
	}
	timeSelectUPI = clock();
	totalTimeSelectUPI = totalTimeSelectUPI + double(timeSelectUPI - timeMaterializeNeighbor);

	// 6.2 Find size k patterns
	while (k <= totalUtilityEachFeat.size())
	{		
		timeWhile = clock();
		// Generate candidate patterns		
		Ckplus = genCandidatePatterns(totalUtilityEachFeat, k);
		Pk.clear();
		// Prun candidates		
		candidatePrun(Ckplus, NonPk, k);
		timeGenerateCandidate = clock();
		totalTimeGenCand = totalTimeGenCand + double(timeGenerateCandidate - timeWhile);

		if (Ckplus.size())
		{
			// Filter star instances			
			std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> SIkplus;
			SIkplus = filterStarInstances(Ckplus, groupSNByFeat, k);			
			timeFilStarInstance = clock();
			totalTimeFilStarInst = totalTimeFilStarInst + double(timeFilStarInstance - timeGenerateCandidate);
			// Select coarse prevalent patterns
			if (k <= 2) // directly get table instances of candidate patterns and calculate their PIs
			{								
				timeFilCliqueInstance = clock();				
				CTIk = computeUPI(SIkplus, totalUtilityEachFeat, w1, w2);
				timeComputeUPI = clock();
				totalTimeComputeUPI = totalTimeComputeUPI + double(timeComputeUPI - timeFilCliqueInstance);
				CIk = SIkplus;
				SIkplus.clear();
				// Filter utility patterns add Pk to the result PkAll				
				selectUtilityPatterns(Pk, NonPk, CTIk, k, prev_thres);
				PkAll.insert(Pk.begin(), Pk.end());
				//Pk.clear();
				CTIk.clear();
				timeSelectUPI = clock();
				totalTimeSelectUPI = totalTimeSelectUPI + double(timeSelectUPI - timeComputeUPI);
			}
			else
			{
				// Filter clique instances			
				SIkplus = filterCliqueInstances(SIkplus, CIk);		
				// Update the size k+1 table instnaces for the next iterator				
				timeFilCliqueInstance = clock();
				totalTimeCollectTableInstance = totalTimeCollectTableInstance
					+ double(timeFilCliqueInstance - timeFilStarInstance);
				// Compute UPI
				CTIk = computeUPI(SIkplus, totalUtilityEachFeat, w1, w2);
				CIk = SIkplus;
				SIkplus.clear();
				timeComputeUPI = clock();
				totalTimeComputeUPI = totalTimeComputeUPI + double(timeComputeUPI - timeFilCliqueInstance);
				// Filter utility patterns add Pk to the result PkAll				
				selectUtilityPatterns(Pk, NonPk, CTIk, k, prev_thres);
				PkAll.insert(Pk.begin(), Pk.end());
				timeSelectUPI = clock();
				totalTimeSelectUPI = totalTimeSelectUPI + double(timeSelectUPI - timeComputeUPI);
			}
			SIkplus.clear();
		}
		// Go to the next iterator
		k += 1;
		Ckplus.clear();
		CTIk.clear();
	}
	// 6.3 Print prevalent patterns	
	std::cout << "The number of prevalent pattern is: " << PkAll.size() << endl;

	// 7. Print running time	
	std::cout << "Execution time for materizaling starneighbors is: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating candidates: " << totalTimeGenCand / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for generating star instances: " << totalTimeFilStarInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for collecting table instances: " << totalTimeCollectTableInstance / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for computing UPI: " << totalTimeComputeUPI / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering UPI patterns: " << totalTimeSelectUPI / CLOCKS_PER_SEC << "s." << endl;
	
	endTime = clock();
	time_taken = double(endTime - startTime);
	std::cout << "Time taken by the program is (No neighboring) : " << (time_taken - totalTimeMatStarNei) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 8. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << "(kB)" << std::endl;
	std::cout << "Peack memory usage: " << physPeackMemUsedByMe / 1024 << "(kB)" << std::endl;

	// 9. Get the maximal size of patterns
	std::vector<int> sizePats;
	std::map<vector<char>, float>::iterator itAllPats = PkAll.begin();
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








