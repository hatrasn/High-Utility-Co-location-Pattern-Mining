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
* @brief: This function computes the utility of each feature type
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its utility.
*/
std::map<char, float> computeUtility(
	std::vector<std::vector<std::string>>& dataList)
{
	std::map<char, float> totalUtilityEachFeat;
	std::map<char, float>::iterator it;
	char feature;
	for (std::vector<std::string> vec : dataList)
	{		
		feature = vec[0][0];
		if (totalUtilityEachFeat.empty())
		{
			totalUtilityEachFeat.insert({ feature, stof(vec[4]) });
		}
		else
		{
			it = totalUtilityEachFeat.find(vec[0][0]);
			if (it != totalUtilityEachFeat.end())
			{
				it->second = it->second + stof(vec[4]);
			}
			else
			{
				totalUtilityEachFeat.insert({ feature, stof(vec[4]) });
			}
		}
	}
	return totalUtilityEachFeat;
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
		cell_key = make_pair(cell_x, cell_y);		
		std::vector<ObjWithCoord> value;
		ObjWithCoord instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[2]), std::stof(vec[3]), std::stof(vec[4]) };
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
	for (std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator it = grid.begin(); it != grid.end(); ++it)
	{
		// Create the five cells
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
					// Calculate the distance
					dist = calculateDistanceTwoInstances(currentInst, checkInst);
					if (dist <= dist_thres) // the two instances have neighbor relationship
					{						
						// Put them into SN
						neighborPair.push_back(currentInst);
						neighborPair.push_back(checkInst);
						std::sort(neighborPair.begin(), neighborPair.end());
						// convert Objects to ObjWithoutCoordinate
						ObjWithoutCoord key = { neighborPair[0].feature, neighborPair[0].instance, neighborPair[0].u };
						ObjWithoutCoord starN = { neighborPair[1].feature, neighborPair[1].instance, neighborPair[1].u };
						if (SN.find(key) != SN.end()) // this instance has already existed
						{
							// update value				
							SN.find(key)->second.insert(starN);
						}
						else // This feature has not existed in SN, directly put into SN	
						{												
							SN.insert({ key, std::set<ObjWithoutCoord> {starN} });
						}
						// Clear all values for the next iterator
						neighborPair.clear();
					}
				}
			}
		}
		// clear all element for the next iterator
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
			// put as a new item
			outValue.insert({ it->first, it->second });
			groupSNByFeat.insert({ it->first.feature, outValue });

		}
		outValue.clear();
	}
	return groupSNByFeat;
}


/*
*@brief: This function generates all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::vector<char> getCombinationChar(const T& c, int combo)
{
	std::vector<char> result;
	int n = c.size();
	for (int i = 0; i < n; ++i) {
		if ((combo >> i) & 1)
			result.push_back(c[i]);
	}
	return result;
}


/**
* @brief: This function generates a subset of a vector.
* @param: c, k
* @r
*/
template<typename T>
std::vector<std::vector<char>> comboChar(const T& c, int k)
{
	std::vector<std::vector<char>> combination;

	int n = c.size();
	int combo = (1 << k) - 1;       // k bit sets
	while (combo < 1 << n)
	{
		combination.push_back(getCombinationChar(c, combo));
		int x = combo & -combo;
		int y = combo + x;
		int z = (combo & ~y);
		combo = z / x;
		combo >>= 1;
		combo |= y;
	}

	return combination;
}


/**
* @brief: This function generates size k candidates
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::vector<std::vector<char>> genCandidatePatterns(
	std::map<char, float> & totalUtilityEachFeat,
	int k)
{		
	std::vector<char> allFeats; // Get all feature types
	for (auto const& item : totalUtilityEachFeat)
	{
		allFeats.push_back(item.first);
	}
	// Combo to get size k candidates
	std::vector<std::vector<char>> Ck;
	Ck = comboChar(allFeats, k);

	// Return
	return Ck;
}



/**
* @brief: This function group a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeature(
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
				// sort
				//std::sort(itgp->second.begin(), itgp->second.end());
			}
			else
			{
				// Directly put into gp
				std::vector<ObjWithoutCoord> innverValueVec{ inst };				
				gp.insert({ inst.feature, innverValueVec });
			}			
		}
	}

	return gp;

}

/**
* @brief: This function generate cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::vector<std::vector<ObjWithoutCoord>>& v)
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
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterStarInstances(
	std::vector<vector<char>> Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> & groupSNByFeat,
	int k)
{
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> SIk; // Save table instances of size k patterns
	char firstFeature; // The first feature of the candidate pattern	
	std::vector<std::vector<ObjWithoutCoord>> innerValue;		
	std::map<char, std::vector<ObjWithoutCoord>> groupInst;
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgroupInst;	
	std::vector<char> remainCandidate; // Save execpt first element of the candidate
	std::vector<char> groupInstKeys; // Save all feature types of the star neighborhoods 
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
			for (auto const& star : itSN->second) // iterator each star neighborhoos
			{				
				// Group a vector of instances by their features, return a hashmap <Feature: Instance>, then put they into innverValue
				groupInst = groupInstanceByFeature(star.second, remainCandidate);
				if(groupInst.size() == remainCandidate.size())				
				{										
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
				groupInstKeys.clear();				
				innerValue.clear();
			}
		}
	}
	return SIk;
}



/**
* @brief: This function calculates PIs of candidate patterns and filter prevalent patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         totalInstEachFeat: the number of instances of each features.
*         prev_thres: a minimum prevalence threshold
*         k: the size of patterns
* @retval: Pk: a hash map of prevalent patterns.
*/
std::map<std::vector<char>, std::vector<float>> computeUPI(
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>>& CIkplus,
	std::map<char, float>& totalUtilityEachFeat,
	float w1,
	float w2)
{
	std::map<std::vector<char>, std::vector<float>> CTIk;
	for (auto const& item : CIkplus)
	{
		// Get instances of each feature
		std::vector<std::set<ObjWithoutCoord>> Tc(item.first.size());
		for (auto const& row : item.second)
		{
			for (size_t t = 0; t < row.size(); ++t)
			{
				Tc[t].insert(row[t]);
			}
		}
		// Calculate utility
		std::vector<float> PRs;
		for (auto const& item : Tc)
		{
			float u = 0.0f;
			for (auto const& inst : item)
			{
				u += inst.u;
			}
			PRs.push_back(u);
		}
		// Calculate InterUR
		std::vector<float> InterUR;
		for (size_t i = 0; i < item.first.size(); ++i)
		{
			float suminterTc = 0.0f;
			float suminterAll = 0.0f;
			if (i == 0)
			{
				for (size_t j = 1; j < item.first.size(); ++j)
				{
					suminterTc += PRs[j];
					suminterAll += totalUtilityEachFeat.find(item.first[j])->second;
				}
				InterUR.push_back(suminterTc / suminterAll);
			}
			else
			{
				for (size_t j = 0; (j != i) && (j < item.first.size()); ++j)
				{
					suminterTc += PRs[j];
					suminterAll += totalUtilityEachFeat.find(item.first[j])->second;
				}
				InterUR.push_back(suminterTc / suminterAll);
			}
		}
		// Calculate IntraUR
		std::vector<float>::iterator itPR = PRs.begin();
		int t = 0;
		while (itPR != PRs.end())
		{
			(*itPR) = (*itPR) / totalUtilityEachFeat.find(item.first[t])->second;
			++t;

			++itPR;
		}
		// Calculate UPR
		itPR = PRs.begin();
		t = 0;
		while (itPR != PRs.end())
		{
			(*itPR) = w1 * (*itPR) + w2 * InterUR[t];
			++t;
			++itPR;
		}
		// Save this candidate
		CTIk.insert({ item.first, PRs });
	}
	// return
	return CTIk;
}



/**
* @brief: This function selects hight utility patterns.
* @param: Pk: a hash map of hight utility patterns
*			NonPK: a hash map of non hight utility patterns
* NonPk = [not high utility pattern: [feature that UPR < threshold], ....]
*			CTIk:
*			prev_thres
* @retval: PK and NonPk
*/
void selectUtilityPatterns(
	std::map<std::vector<char>, float>& Pk,
	std::vector<std::pair<std::vector<char>, std::vector<char>>>& NonPk,
	std::map<std::vector<char>, std::vector<float>>& CTIk,
	float k,
	float prev_thres)
{
	std::map<std::vector<char>, std::vector<float>>::iterator it = CTIk.begin();
	while (it != CTIk.end())
	{
		// Calculate the UPI
		float PI = *std::min_element(it->second.begin(), it->second.end());
		if (PI >= prev_thres)
		{
			Pk.insert({ it->first, PI });
		}
		else
		{			
			// Get non-hight-utility patterns for prunning
			std::vector<char> oneNonHU;
			for (size_t t = 0; t < it->first.size(); ++t)
			{
				if (it->second[t] <= prev_thres)
				{
					oneNonHU.push_back(it->first[t]);
				}
			}			
			std::pair<std::vector<char>, std::vector<char>> nonuti = std::make_pair(it->first, oneNonHU);
			NonPk.push_back(nonuti);
		}		
		// Terminate
		++it;
	}
}



/**
* @brief: This function finds real row instances of candidate patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         CIk: a hash map of table instances of size k candiate patterns.
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterCliqueInstances(
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> & CIkplus,
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> & CIk)

{
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> tabInstkplus; // Save the size (k+1) table instances
	std::vector<char> patternkSub; // Save the sub pattern of pattern of size (k+1). Example (A, B, C) -> patternkSub = (B, C)
	std::vector<std::vector<ObjWithoutCoord>> realTableInstkplus;	
	std::vector<ObjWithoutCoord> checkRowkPlus;
	// Loop each element of CIkplus
	for (auto const& item : CIkplus)
	{
		// Get sub pattern (execpt the first element)
		patternkSub.insert(patternkSub.end(), item.first.begin() + 1, item.first.end());		
		// Check each row instance in coarseTabInstkPlus
		if (CIk.count(patternkSub)) // If the sub pattern is prevalent
		{			
			// Loop each row instance in coarseTabInstkPlus to check if it is a clique
			for (auto const& rowInstKPlus : item.second)
			{				
				checkRowkPlus.insert(checkRowkPlus.end(), rowInstKPlus.begin() + 1, rowInstKPlus.end());				
				if (std::find(CIk.find(patternkSub)->second.begin(), CIk.find(patternkSub)->second.end(), checkRowkPlus) != CIk.find(patternkSub)->second.end())
				{
					// The current row instance is a clique
					realTableInstkplus.push_back(rowInstKPlus);
				}
				checkRowkPlus.clear();
			}
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



/**
* @brief: This function pruns candidates
* @param:
* @retval: PK and NonPk
*/
void candidatePrun(
	std::vector<std::vector<char>>& Ckplus, 
	std::vector<std::pair<std::vector<char>, std::vector<char>>>& NonPk,
	int k)
{
	if (k > 2)
	{		
		for (size_t i = 0; i < NonPk.size(); ++i)
		{
			for (size_t j = i + 1; j < NonPk.size(); ++j)
			{
				// Get common feature
				std::vector<char> commonF;
				std::set_intersection(NonPk[i].first.begin(), NonPk[i].first.end(),
					NonPk[j].first.begin(), NonPk[j].first.end(),
					std::back_inserter(commonF));
				// Only have one common feature
				if (commonF.size() == 1)
				{
					// Check this commomn feature is or not a non high utility feature
					if (std::find(NonPk[j].second.begin(), NonPk[j].second.end(), commonF[0])
						!= NonPk[j].second.end())
					{						
						// generate c = c1 and c2	
						std::vector<char> prunPat;
						std::set_union(NonPk[i].first.begin(), NonPk[i].first.end(),
							NonPk[j].first.begin(), NonPk[j].first.end(),
							std::back_inserter(prunPat));
						std::sort(prunPat.begin(), prunPat.end());
						// Only combine to a size-k candidate
						if (prunPat.size() == k)
						{
							// Find this candidate from candidate set and delete it
							auto itfind = std::find(Ckplus.begin(), Ckplus.end(), prunPat);
							if (itfind != Ckplus.end())
							{
								Ckplus.erase(itfind);
							}						
						}
					}					
				}
			}
		}
	}
}



