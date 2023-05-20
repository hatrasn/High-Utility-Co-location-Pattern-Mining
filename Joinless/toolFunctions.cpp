#include <iostream>
#include <map>
#include <vector>
#include "toolFunctions.h"
#include "Object.h"
#include <string>
#include <algorithm>
#include <math.h>
#include <cstdlib>
#include <numeric>
#include <boost/algorithm/string.hpp>


using namespace std;

/**
* @brief: This function counts the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
std::map<char, int> countNumberInstance(std::vector<std::vector<std::string>> dataList)
{
	std::map<char, int> totalInstNumEachFeat;
	char feature;
	int instance;
	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		feature = vec[0][0]; // get the feature	
		if (totalInstNumEachFeat.empty())
		{
			totalInstNumEachFeat.insert({ feature, 1 });
		}
		else if (totalInstNumEachFeat.find(feature) == totalInstNumEachFeat.end())
		{
			totalInstNumEachFeat.insert({ feature, 1 });
		}
		else
		{
			instance = totalInstNumEachFeat.find(feature)->second;
			instance += 1;
			totalInstNumEachFeat[feature] = instance;
		}
	}

	return totalInstNumEachFeat;
}


/**
* @brief This function makes a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(
	std::vector<std::vector<std::string>> & dataList, 
	float dist_thres)
{
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid;
	int cell_x, cell_y;
	std::pair<int, int> cell_key;
	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		cell_x = ceil(std::stof(vec[2]) / dist_thres); // x coordinate
		cell_y = ceil(std::stof(vec[3]) / dist_thres); // y coordinate		
													   // make keys
		cell_key = make_pair(cell_x, cell_y);		
		std::vector<ObjWithCoord> value;
		ObjWithCoord instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[2]), std::stof(vec[3]) };
		value.push_back(instance);
		// check if this key is exsiting?
		if (grid.empty()) //The grid is empty
		{
			grid[cell_key] = value;
		}
		else if (grid.find(cell_key) == grid.end()) // if the key is not exist
		{
			grid[cell_key] = value;
		}
		else // the key has already existed
		{
			std::vector<ObjWithCoord> old_value = grid.find(cell_key)->second;
			old_value.push_back(instance);
			grid[cell_key] = old_value;
		}
	}

	return grid;
}


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst)
{
	return  sqrt((currentInst.x - checkInst.x)*(currentInst.x - checkInst.x)
		+ (currentInst.y - checkInst.y)*(currentInst.y - checkInst.y));
}



/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> & grid, 
	float dist_thres)
{
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN;
	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	float dist;
	std::vector<ObjWithCoord> fiveCellInst; // save all instance in the five cells and two neighboring instances.
	std::vector<ObjWithCoord> neighborPair;
	std::set<ObjWithoutCoord> value; // save all neighbors of the center instance
	// Loop each cell in grid and computation neighbors in five cells
	for (std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator it = grid.begin(); it != grid.end(); ++it)
	{		
		i = it->first.first;
		j = it->first.second;
		fiveCells.push_back(it->first);
		fiveCells.push_back(make_pair(i, j + 1));
		fiveCells.push_back(make_pair(i + 1, j + 1));
		fiveCells.push_back(make_pair(i + 1, j));
		fiveCells.push_back(make_pair(i + 1, j - 1));
		for (auto const& cell : fiveCells)
		{
			if (grid.count(cell))
			{
				fiveCellInst.insert(fiveCellInst.end(), grid.find(cell)->second.begin(), grid.find(cell)->second.end());
			}
		}
		// Sort all instances in the five cells
		std::sort(fiveCellInst.begin(), fiveCellInst.end());
		// Iterator each instance in the currentCellInst and check its neighbors		
		for (auto const& currentInst : it->second)
		{
			for (auto const& checkInst : fiveCellInst) // check with each instance in the five cells
			{
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{					
					dist = calculateDistanceTwoInstances(currentInst, checkInst);
					if (dist <= dist_thres) // the two instances have neighbor relationship
					{						
						neighborPair.push_back(currentInst);
						neighborPair.push_back(checkInst);
						std::sort(neighborPair.begin(), neighborPair.end());						
						ObjWithoutCoord key = { neighborPair[0].feature, neighborPair[0].instance };
						ObjWithoutCoord starN = { neighborPair[1].feature, neighborPair[1].instance };
						if (SN.find(key) != SN.end()) // this instance has already existed
						{
							// update value				
							SN.find(key)->second.insert(starN);
						}
						else // This feature has not existed in SN, directly put into SN	
						{
							// Put into SN							
							SN.insert({ key, std::set<ObjWithoutCoord> {starN} });
						}
						// Clear all values for the next iterator
						neighborPair.clear();
					}
				}
			}
		}
		// Clear all elements for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}

	return SN;
}


/**
* @brief: This function groups star neighbor instances by the feature type of the key.
* @param: SN: a hashmap of star neighborhoods*
* @retval: groupSNByFeat: a hash map that stores as <feature, <instance: neighborhoods>>.
*/
std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupStarNeighByFeatures(
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>& SN)
{
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupSNByFeat;
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> outValue;
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>::iterator itGroupNei;
	for (std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>::iterator it = SN.begin(); it != SN.end(); ++it)
	{		
		itGroupNei = groupSNByFeat.find(it->first.feature);
		if (itGroupNei != groupSNByFeat.end())
		{
			// update			
			itGroupNei->second.insert({ it->first, it->second });
		}
		else
		{
			// put as new item
			outValue.insert({ it->first, it->second });
			groupSNByFeat.insert({ it->first.feature, outValue });

		}
		outValue.clear();
	}
	return groupSNByFeat;
}



/**
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*         k:
* @retval: Cnk
*/
template<typename T>
std::vector<std::vector<T>> combinekOfn(std::vector<T>& c, int k)
{
	std::vector<std::vector<T>> Cnk;
	int n = c.size();
	std::vector<bool> v(n);
	std::fill(v.begin(), v.begin() + k, true);
	do {
		std::vector<T> onec;
		for (int i = 0; i < n; ++i) {
			if (v[i])
			{
				onec.push_back(c[i]);
			}
		}
		Cnk.push_back(onec);
	} while (std::prev_permutation(v.begin(), v.end()));
	// Return
	return Cnk;
}



/**
* @brief: This function generates candidates
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::vector<std::vector<char>> genCandidatePatterns(
	std::map<std::vector<char>, float> Pk, int k)
{
	std::vector<vector<char>> Ck;
	if (k <= 2)
	{
		// Directly combination		
		std::vector<char> features;
		for (auto const& ele : Pk)
		{
			features.push_back(ele.first[0]);
		}
		// Combination to generate size 2 patterns	
		Ck = combinekOfn(features, k);
	}
	else // k>2
	{		
		std::vector<vector<char>> keys; // Save size k patterns
		std::vector<vector<char>> subPatterns; // Save size k patterns of a size k+1 pattern
		std::vector<char> sizeKPlusPattern;	// Save size k+1 pattern
		int flag = 0; // flag the subsets of size (k+1) are all in the size k pattern
		// Retrieve all keys form Pk
		for (auto const& ele : Pk)
		{
			keys.push_back(ele.first);
		}		
		std::vector<vector<char>> tempKeys = keys;
		for (auto const& pattern : keys)
		{
			// Delete this pattern from tempKeys
			tempKeys.erase(std::remove(tempKeys.begin(), tempKeys.end(), pattern), tempKeys.end());						
			// Check (k-1) is the same
			for (auto const& tempPattern : tempKeys)
			{
				bool checkFirstK = std::equal(pattern.begin(), pattern.begin() + (k - 2), tempPattern.begin());
				// check (k-1) features are the same	
				if (checkFirstK) 
				{
					// compose to size k pattern					
					sizeKPlusPattern = pattern;
					sizeKPlusPattern.push_back(tempPattern.back());
					std::sort(sizeKPlusPattern.begin(), sizeKPlusPattern.end());
					// Check sub-sets are prevalent					
					subPatterns = combinekOfn(sizeKPlusPattern, k - 1);
					for (auto const& subPat : subPatterns)
					{
						if (std::find(keys.begin(), keys.end(), subPat) != keys.end())
						{
							flag += 1;
						}
						else
						{
							break;
						}
					}
					if (flag == subPatterns.size()) // if all sub patterns of size k+1 is prevalent, size k+1 is a candidate
					{
						Ck.push_back(sizeKPlusPattern);
					}
					// Clear all for the next iterator
					sizeKPlusPattern.clear();
					subPatterns.clear();
					flag = 0;
				}
			}
		}
	}	
	return Ck;
}



/**
* @brief: This function groups a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeature3(
	std::set<ObjWithoutCoord> instSet,
	std::vector<char>& remainCandidate)
{
	std::map<char, std::vector<ObjWithoutCoord>> gp;
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgp;
	// Loop each neihbor instance
	for (auto const& inst : instSet)
	{
		// Only put the instance has feature type is in remainCandidate
		if (std::find(remainCandidate.begin(), remainCandidate.end(), inst.feature) != remainCandidate.end())
		{
			itgp = gp.find(inst.feature);
			if (itgp != gp.end()) // the fearture has already existed
			{
				itgp->second.push_back(inst);				
			}
			else
			{				
				std::vector<ObjWithoutCoord> innverValueVec{ inst };
				gp.insert({ inst.feature, innverValueVec });
			}
		}
	}

	return gp;
}



/**
* @brief: This function generates cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::vector<std::vector<ObjWithoutCoord>> v)
{
	std::vector<std::vector<ObjWithoutCoord>> resultCartesian;
	auto product = [](long long a, std::vector<ObjWithoutCoord>& b)
	{
		return a*b.size();
	};

	const long long N = accumulate(v.begin(), v.end(), 1LL, product);
	std::vector<ObjWithoutCoord> u(v.size());
	for (long long n = 0; n<N; ++n)
	{
		lldiv_t q{ n, 0 };
		for (long long i = v.size() - 1; 0 <= i; --i) {
			q = div(q.quot, v[i].size());
			u[i] = v[i][q.rem];
		}		
		resultCartesian.push_back(u);
	}
	return resultCartesian;
}



/**
* @brief: This function filters star instances of a candidate
* @param: totalInstEachFeat: the instance number of each feature
* @retval: Nope
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterStarInstancesSize2(
	std::vector<std::vector<char>> Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>& groupSNByFeat,
	int k)
{
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> SIk; // Save table instances of size k patterns
	char firstFeature; // The first feature of the candidate pattern
	std::vector<std::vector<ObjWithoutCoord>> innerValue;
	std::map<char, std::vector<ObjWithoutCoord>> groupInst; // group instances by feature types;
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgroupInst;
	std::vector<char> remainCandidate; // Save execpt first element of the candidate
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>::iterator itSN;
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>>::iterator itSIk;
	// Check each candidate and take back its table instance
	for (auto const& candidate : Ckplus)
	{
		remainCandidate.clear();
		// Retrieve the first feature		
		firstFeature = *(candidate.begin());
		// Get the remain feature in the candidate
		remainCandidate.insert(remainCandidate.end(), candidate.begin() + 1, candidate.end());
		// Get the star neighbor of the first feature		
		itSN = groupSNByFeat.find(firstFeature);
		if (itSN != groupSNByFeat.end())// if the first feature has star neighborhoods
		{
			// Package to vector<vector<ObjWithoutCoordinate>>
			for (auto const& star : itSN->second) // iterator each star neighborhoos
			{
				// Put instances of the other feature into innverValue				
				groupInst = groupInstanceByFeature3(star.second, remainCandidate);
				if (groupInst.size() == remainCandidate.size())
				{
					// Put the instance of the first feature into the innerValue					
					std::vector<std::vector<ObjWithoutCoord>> tempInnerValue(k); // Save innver value in SIk structure
					tempInnerValue[0].push_back(star.first);
					// Update SIk		
					for (size_t t = 1; t < k; ++t)
					{
						itgroupInst = groupInst.find(candidate[t]);
						tempInnerValue[t].insert(tempInnerValue[t].end(),
							itgroupInst->second.begin(),
							itgroupInst->second.end());
					}
					// Cartesian Product tempInnerValue to generate row instances
					innerValue = cartesianProduct(tempInnerValue);
					tempInnerValue.clear();
					// Put into SIk
					if (innerValue.size())
					{
						itSIk = SIk.find(candidate);
						if (itSIk != SIk.end())
						{
							// Update SIk fast version
							itSIk->second.insert(itSIk->second.end(), innerValue.begin(), innerValue.end());
						}
						else  // if SIk is empty
						{
							// Put as a new value
							SIk.insert({ candidate, innerValue });
						}
					}
				}
				// Clear temporatory varibles to next iterator				
				groupInst.clear();
				innerValue.clear();
			}
		}
	}
	return SIk;
}




/**
* @brief: This function filters star instances of a candidate
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> filterStarInstances(
	std::vector<vector<char>> Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> & groupSNByFeat,
	int k)
{	
	// Save table instances of size k patterns
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> SIk;
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>::iterator itSIk;
	// The first feature of the candidate pattern
	char firstFeature;
	std::vector<std::vector<ObjWithoutCoord>> innerValue;
	// group instances by feature types;
	std::map<char, std::vector<ObjWithoutCoord>> groupInst; 
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgroupInst;
	// Save execpt first element of the candidate
	std::vector<char> remainCandidate; 
	// Check each candidate and take back its table instance
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>::iterator itSN;	
	for (auto const& candidate : Ckplus)
	{
		remainCandidate.clear();
		// Retrieve the first feature		
		firstFeature = *(candidate.begin());
		// Get the remain feature in the candidate
		remainCandidate.insert(remainCandidate.end(), candidate.begin() + 1, candidate.end());
		// Get the star neighbor of the first feature		
		itSN = groupSNByFeat.find(firstFeature);
		if (itSN != groupSNByFeat.end())// if the first feature has star neighborhoods
		{
			// Package to vector<vector<ObjWithoutCoordinate>>
			for (auto const& star : itSN->second) // iterator each star neighborhoos
			{
				// Put instances of the other feature into innverValue				
				groupInst = groupInstanceByFeature3(star.second, remainCandidate);
				if (groupInst.size() == remainCandidate.size())
				{
					// Put the instance of the first feature into the innerValue					
					std::vector<std::vector<ObjWithoutCoord>> tempInnerValue(k); // Save innver value in SIk structure
					tempInnerValue[0].push_back(star.first);
					// Update SIk		
					for (size_t t = 1; t < k; ++t)
					{
						itgroupInst = groupInst.find(candidate[t]);
						tempInnerValue[t].insert(tempInnerValue[t].end(),
							itgroupInst->second.begin(),
							itgroupInst->second.end());
					}
					// Put this instance rows into the result
					itSIk = SIk.find(candidate);
					if (itSIk != SIk.end())
					{
						itSIk->second.push_back(tempInnerValue);
					}
					else
					{
						SIk.insert({ candidate, 
							std::vector<std::vector<std::vector<ObjWithoutCoord>>> {tempInnerValue} });
					}
				}
			}
		}
	}
	return SIk;
}



/**
* @brief: This function calculates PIs of candidate patterns and filters prevalent patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         totalInstEachFeat: the number of instances of each features.
*         prev_thres: a minimum prevalence threshold
*         k: the size of patterns
* @retval: Pk: a hash map of prevalent patterns.
*/
std::map<std::vector<char>, float> selectPrevalentPatterns(
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> & CIkplus,
	std::map<char, int> & totalInstEachFeat,
	float prev_thres,
	int k)
{
	std::map<vector<char>, float> Pk;
	std::vector<float> PRs;
	float PI;
	std::vector<char> pattern;	
	std::vector<std::set<ObjWithoutCoord>> uniqueTableInstance;
	int numberInst, totalInst;
	for (std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>>::iterator iter = CIkplus.begin(); 
		iter != CIkplus.end(); ++iter)
	{
		// Retrieve patterns and table instances
		pattern = iter->first;
		// Loop each row instance in tableInstance
		for (auto const& rowInst : iter->second)
		{
			if (uniqueTableInstance.empty())
			{
				for (int i = 0; i < k; ++i)
				{
					std::set<ObjWithoutCoord> temp;
					temp.insert(rowInst[i]);
					uniqueTableInstance.push_back(temp);
				}
			}
			else
			{
				for (int i = 0; i < k; ++i)
				{
					uniqueTableInstance[i].insert(rowInst[i]);
				}
			}
		}
		// Loop each element
		for (int i = 0; i < k; ++i)
		{
			totalInst = totalInstEachFeat.find(pattern[i])->second;			
			numberInst = uniqueTableInstance[i].size();
			PRs.push_back((float)numberInst / (float)totalInst);
		}
		// Find the minimum elemnet
		PI = *std::min_element(std::begin(PRs), std::end(PRs));
		// Check prevalent patterns
		if (PI >= prev_thres)
		{
			Pk.insert({ pattern, PI });
		}
		// Clear all temporary varibles
		uniqueTableInstance.clear();
		PRs.clear();
		pattern.clear();		
	}
	return Pk;
}


/**
* @brief: This function selects coarse prevalent patterns.
* @param: SIK: a hash map of table instances
* @retval: CIk: a hash map of table instances which satify coarse min_prev
*/
void selectCoarsePrevalentPatterns(
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIkplus,
	std::map<char, int>& totalInstEachFeat,
	float prev_thres,
	int k)
{
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>::iterator iter = SIkplus.begin();
	while(iter != SIkplus.end())
	{
		// Loop each row instance to get tableInstance
		std::vector<std::set<ObjWithoutCoord>> uniqueTableInstance(k);
		for (auto const& rows : iter->second)
		{
			for (size_t t = 0; t < k; ++t)
			{
				uniqueTableInstance[t].insert(rows[t].begin(), rows[t].end());
			}
		}		
		// Calculate PRs
		std::vector<float> PRs;
		for (size_t t = 0; t < k; ++t)
		{
			int totalInst = totalInstEachFeat.find(iter->first[t])->second;			
			int numberInst = uniqueTableInstance[t].size();			
			PRs.push_back((float)numberInst / (float)totalInst);
		}
		// Find the minimum elemnet
		float PI = *std::min_element(std::begin(PRs), std::end(PRs));		
		// Check prevalent patterns
		if (PI < prev_thres)
		{
			iter = SIkplus.erase(iter);
		}
		else
		{
			++iter;
		}				
	}
}



/**
* @brief: This function finds real row instances of candidate patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         CIk: a hash map of table instances of size k candiate patterns.
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterCliqueInstances2(
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> & SIkplus,
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk)
{
	// Save the size (k+1) table instances
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> tabInstkplus;
	// Save the sub pattern of pattern of size (k+1). Example (A, B, C) -> patternkSub = (B, C)
	std::vector<char> patternkSub; 
	std::vector<std::vector<ObjWithoutCoord>> innerValue, realTableInstkplus;	
	std::vector<ObjWithoutCoord> checkRowkPlus;
	// Loop each element of CIkplus
	for (auto const& item : SIkplus)
	{
		// Get sub pattern (execpt the first element)
		patternkSub.insert(patternkSub.end(), item.first.begin() + 1, item.first.end());		
		// Check each row instance in coarseTabInstkPlus
		if (CIk.find(patternkSub) != CIk.end())
		{	
			for (auto const& rows : item.second)
			{
				// Product coarse row instances			
				innerValue = cartesianProduct(rows);
				// Check real row instances				
				for (auto const& rowInstKPlus : innerValue)
				{					
					checkRowkPlus.insert(checkRowkPlus.end(), rowInstKPlus.begin() + 1, rowInstKPlus.end());
					if (std::find(CIk.find(patternkSub)->second.begin(), 
						CIk.find(patternkSub)->second.end(), checkRowkPlus) 
						!= CIk.find(patternkSub)->second.end())
					{
						// The current row instance is a clique
						realTableInstkplus.push_back(rowInstKPlus);
					}
					checkRowkPlus.clear();
				}
			}			
			// Save the table instance of the current candidate
			if (realTableInstkplus.size() > 0)
			{
				// Put into the final result
				tabInstkplus.insert({ item.first, realTableInstkplus });
				realTableInstkplus.clear();
			}			
		}
		// Clear all
		patternkSub.clear();
	}
	return tabInstkplus;
}



