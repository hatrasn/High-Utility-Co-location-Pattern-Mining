
/* This code implements the EPA algorithm that has been published in 
* 
* Shisheng Yang, LizhenWang, Xuguang Bao, and Junli Lu. 2016. A framework for mining spatial high utility co-location
* patterns. In Proc. 12th Int. Conf. Fuzzy Syst. Knowl. Discov. 595–601. https://doi.org/10.1109/FSKD.2015.7382010
**
* Writer: vanha tran
* Date: 11.8.2020
*/


#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <string>


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
#include "pattern.h"


using namespace std;


int main()
{
	// 1. Set time
	clock_t startTime, endTime; // start and end time
	clock_t timeMaterializeNeighbor; // time for materializing star neighborhoods	
	clock_t timeComputeUPI;
	clock_t timeGenerateCandidate; // time for generating candidates
	clock_t timeFilStarInstance;
	clock_t timeSelectUPI; // time for filtering co-location instances
	clock_t timeWhile;

	double time_taken; // count total execution time
	double totalTimeMatStarNei = 0.0; // total time for materializing star neighborhoods	
	double totalTimeComputeUPI = 0.0;
	double totalTimeGenCand = 0.0; // total time for generating candidates
	double totalTimeFilStarInst = 0.0;
	double totalTimeSelectUPI = 0.0; // total time for filgering candidates

	startTime = clock(); // Begin count time

	// 2. Set parameters
	double dist_thres = 3200.0f;  // set the distance threshold
	double prev_thres = 0.1f; // set the prevalence
	int s = 4; // the number of top-s subsets
	std::vector<int> topk{ 1, 2, 3, 4, 5 }; // set the number of top-k patterns
	int numtopk = 2; // get the Qc for each size pattern in top of numtopk
	int topkpersize = 2; // get QC of top-k of each size

	// 3. Load data sets
	string file_name = "./Data/quality/gongshan_plant_alphabet_sample_50_checkin.csv";
	freopen("./Data/quality/gongshan_plant_alphabet_sample_50_checkin_distance_3200_PI_01_new_n.txt", "w", stdout);

	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File	
	dataList.erase(dataList.begin()); // Delete the first line

	// 4. Get max and min of utilities	
	std::unordered_map<char, double> maxU;
	std::map<char, int> numInstFeat;
	std::set<char> allF;
	getMaxUtility(dataList, maxU, numInstFeat, allF);
	
	// 5. Make grids and find neighbors
	// 5.1 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid = makeGrid(dataList, dist_thres);
	std::map<char, double> totalUtilityEachFeat; //count the instance number of each feature
	// 5.2 Find neighbors, generate A1: <B.1, C.1>
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN = genStarNeighborhoods(grid, dist_thres, maxU, totalUtilityEachFeat);
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupSNByFeat = groupStarNeighByFeatures(SN);	
	SN.clear(); // Free memory
	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);

	// 6. Calculate the utiltity of all data set and each features
	double totalUti = getUtilityOFDataset(totalUtilityEachFeat);	
	// The utility of features is the average utilities of utility of instances	
	getUtilityOfFeat(totalUtilityEachFeat, numInstFeat);
	
	// 7. Find high utility co-location patterns	
	std::vector<Pattern> Pk, PkAll; // Size k prevalent patterns	
	std::map<std::string, double> NonPk; // Save all non-size-k-prevalent patterns
	std::vector<Pattern> CTIk; // Store the PI of each candidate	
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> CIk; // Store table isntances of size k patterns	
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>> TAllk; // Store table isntances of all size k patterns, // {A, B, C} : <<instances of A>, <instances of B>, <instances of C>>
	std::unordered_map<std::string, int> pruneCands; // Store all prune candidates

	int k = 2;
	// 7.1 Get size 1 pattern
	for (auto const& pat : totalUtilityEachFeat)
	{
		double uf = totalUtilityEachFeat.find(pat.first)->second;
		std::string patsize1{ pat.first };

		Pattern onesize1{ patsize1, totalUtilityEachFeat.find(pat.first)->second, 0.0 };
		Pk.push_back(onesize1);
	}
	timeSelectUPI = clock();
	totalTimeSelectUPI = totalTimeSelectUPI + double(timeSelectUPI - timeMaterializeNeighbor);

	// 7.2. Find size k patterns	
	while (k <= totalUtilityEachFeat.size())		
	{
		timeWhile = clock();		
		std::unordered_map<std::string, int> Ckplus; // Generate candidate patterns
		genCandidatePatterns2(Ckplus, totalUtilityEachFeat, k);
		Pk.clear();			
		candidatePrun(Ckplus, pruneCands); // Prun candidates	
		timeGenerateCandidate = clock();
		totalTimeGenCand = totalTimeGenCand + double(timeGenerateCandidate - timeWhile);

		if (Ckplus.size())
		{			
			if (k <= 2) // directly get table instances of candidate patterns and calculate their PIs
			{
				// Filter star instances			
				std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>> SIkplus;
				filterStarInstancesSize2(SIkplus, Ckplus, groupSNByFeat, k);
				timeFilStarInstance = clock();
				totalTimeFilStarInst = totalTimeFilStarInst + double(timeFilStarInstance - timeGenerateCandidate);				
				CTIk = computeUPI(SIkplus, totalUtilityEachFeat, totalUti, numInstFeat, TAllk);
				timeComputeUPI = clock();
				totalTimeComputeUPI = totalTimeComputeUPI + double(timeComputeUPI - timeFilStarInstance);
				CIk = SIkplus;
				SIkplus.clear();
				selectUtilityPatterns(Pk, NonPk, CTIk, k, prev_thres); // Filter utility patterns add Pk to the result PkAll	
				PkAll.insert(PkAll.end(), Pk.begin(), Pk.end());

				// Find size >k non-prvalent candidate
				findPruneCandidate(NonPk, allF, TAllk, s, totalUtilityEachFeat, totalUti, prev_thres, pruneCands);

				// Free memory
				CTIk.clear();
				timeSelectUPI = clock();
				totalTimeSelectUPI = totalTimeSelectUPI + double(timeSelectUPI - timeComputeUPI);
			}
			else
			{				
				std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>> SIkplus; // In order to save memory, this step generate row instances and validate real cliques 
				filterStarInstancesSizek(SIkplus, Ckplus, groupSNByFeat, k, CIk);
				CIk.clear();				
				timeFilStarInstance = clock();
				totalTimeFilStarInst = totalTimeFilStarInst + double(timeFilStarInstance - timeGenerateCandidate);				
				if (SIkplus.size())
				{
					CTIk = computeUPI(SIkplus, totalUtilityEachFeat, totalUti, numInstFeat, TAllk);
					CIk = SIkplus;
					SIkplus.clear();
				}
				timeComputeUPI = clock();
				totalTimeComputeUPI = totalTimeComputeUPI + double(timeComputeUPI - timeFilStarInstance);
				selectUtilityPatterns(Pk, NonPk, CTIk, k, prev_thres); // Filter utility patterns add Pk to the result PkAll	
				PkAll.insert(PkAll.end(), Pk.begin(), Pk.end());				
				findPruneCandidate(NonPk, allF, TAllk, s, totalUtilityEachFeat, totalUti, prev_thres, pruneCands);
				timeSelectUPI = clock();
				totalTimeSelectUPI = totalTimeSelectUPI + double(timeSelectUPI - timeComputeUPI);
			}
		}
		// Go to the next iterator
		k += 1;	
		CTIk.clear();		
	}

	// 7.3 Print prevalent patterns	
	std::cout << "The number of prevalent pattern is: " << PkAll.size() << endl;
	
	//8. Print the running time
	std::cout << "Execution time for materizaling starneighbors is: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating candidates: " << totalTimeGenCand / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for generating star instances: " << totalTimeFilStarInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for computing UPI: " << totalTimeComputeUPI / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering UPI patterns: " << totalTimeSelectUPI / CLOCKS_PER_SEC << "s." << endl;

	// 9. Calculate the total execution time
	endTime = clock();
	time_taken = double(endTime - startTime);
	std::cout << "Time taken by the program is (No neighboring) : " <<
		(totalTimeGenCand + totalTimeFilStarInst + totalTimeComputeUPI + totalTimeSelectUPI) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 10. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << "(kB)" << std::endl;
	std::cout << "Peak memory usage: " << physPeackMemUsedByMe / 1024 << "(kB)" << std::endl;

	// 11. Get the maximal size of patterns
	int maxSizePatt;
	std::map<int, int> sizePats;
	classifyPatternAndGetMaxSizeofPatterns(PkAll, sizePats, maxSizePatt);
	std::cout << "The maximal size of patterns is: " << maxSizePatt << endl;
	
		

	// The next code in here




	return 0;
}








