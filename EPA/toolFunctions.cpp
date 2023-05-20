#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <cstdlib>
#include <numeric>
#include <unordered_map>
#include <sstream>
#include <iomanip>


#include "toolFunctions.h"
#include "Object.h"

#include <boost/algorithm/string.hpp>


using namespace std;



/**
* @brief: This function computes the utility of each feature type
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its utility.
*/
std::map<char, double> computeUtility(
	std::vector<std::vector<std::string>>& dataList)
{
	std::map<char, double> totalUtilityEachFeat;
	std::map<char, double>::iterator it;
	char feature;
	int numofinst = 0;
	std::map<char, int> numInstEachFeat;
	std::map<char, int>::iterator itnumInstEachFeat;
	for (std::vector<std::string> vec : dataList)
	{
		// get the feature	
		feature = vec[0][0];
		// Check this feature for utility
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
		// Check this feature for number of instances
		if (numInstEachFeat.empty())
		{
			numInstEachFeat.insert({ feature, 1 });
		}
		else
		{
			itnumInstEachFeat = numInstEachFeat.find(feature);
			if (itnumInstEachFeat != numInstEachFeat.end())
			{
				itnumInstEachFeat->second += 1;
			}
			else
			{
				numInstEachFeat.insert({ feature, 1 });
			}
		}

	}
	// The utility of each feature is the average utility of sum its instances
	it = totalUtilityEachFeat.begin();
	while (it != totalUtilityEachFeat.end())
	{
		it->second = it->second / (double)numInstEachFeat.find(it->first)->second;
		++it;
	}
	// Return
	return totalUtilityEachFeat;
}



/**
* @brief: This function gets the maximal and minimal values of the utility of each feature type
* @param: dataList: an input dataset;
* maxU: the maximal utility
* numInstFeat: the total number instanes of each feature
* allF: a set of feature types
* @retval: maxminU: <A, [minU, maxU]>, <B, [minU, maxU]>, ...
*/
void getMaxUtility(
	std::vector<std::vector<std::string>>& dataList,
	std::unordered_map<char, double>& maxU,
	std::map<char, int> & numInstFeat,
	std::set<char>& allF)
{
	//1. Loop each instance and get its utility
	std::unordered_map<char, std::vector<double>> allU;
	std::unordered_map<char, std::vector<double>>::iterator itallU;
	std::map<char, int>::iterator itnumInstFeat;
	for (auto const& row : dataList)
	{
		// Get feature
		allF.insert(row[0][0]);
		// Get utility
		itallU = allU.find(row[0][0]);
		if (itallU != allU.end())
		{
			itallU->second.push_back(std::stof(row[4]));
		}
		else
		{
			allU.insert({ row[0][0], std::vector<double>{std::stof(row[4])} });
		}
		// Get number of instances of features
		itnumInstFeat = numInstFeat.find(row[0][0]);
		if (itnumInstFeat != numInstFeat.end())
		{
			itnumInstFeat->second += 1;
		}
		else
		{
			numInstFeat.insert({ row[0][0], 1});
		}
	}
	//2. Get min and max
	for (auto const& item : allU)
	{		
		double mU = *std::max_element(item.second.begin(), item.second.end());

		maxU.insert({ item.first, mU });
	}
}



/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(
	std::vector<std::vector<std::string>> & dataList, 
	double dist_thres)
{
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid;
	int cell_x, cell_y;
	std::pair<int, int> cell_key;
	for (std::vector<std::string> vec : dataList)
	{		
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
double calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst)
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
	double dist_thres,
	std::unordered_map<char, double>& maxU,
	std::map<char, double>& totalUtilityEachFeat)
{
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN;
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>::iterator itSN;
	std::unordered_map<char, double>::iterator itmaxU;
	std::map<char, double>::iterator ittotalU;
	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	double dist;
	// save all instance in the five cells and two neighboring instances.
	std::vector<ObjWithCoord> fiveCellInst;	
	for (std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator it = grid.begin();
		it != grid.end();
		++it)
	{
		// Create the five cells
		i = it->first.first;
		j = it->first.second;
		fiveCells.push_back(it->first);
		fiveCells.push_back(make_pair(i, j + 1));
		fiveCells.push_back(make_pair(i + 1, j + 1));
		fiveCells.push_back(make_pair(i + 1, j));
		fiveCells.push_back(make_pair(i + 1, j - 1));
		// Get all instances in the current cell and five cells
		std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator itgrid;
		for (auto const& cell : fiveCells)
		{
			itgrid = grid.find(cell);
			if (itgrid != grid.end())
			{
				fiveCellInst.insert(fiveCellInst.end(),
					itgrid->second.begin(), itgrid->second.end());
			}
		}
		// Sort all instances in the five cells
		std::sort(fiveCellInst.begin(), fiveCellInst.end());
		// Iterator each instance in the currentCellInst and check its neighbors		
		for (auto const& currentInst : it->second)
		{		
			// Calculate nomorlized utiltiy of this instances
			itmaxU = maxU.find(currentInst.feature);			
			double uc = log10(currentInst.u) / log10(itmaxU->second);			
			// Save this uc to corresponding feature type
			ittotalU = totalUtilityEachFeat.find(currentInst.feature);
			if (ittotalU != totalUtilityEachFeat.end())
			{
				ittotalU->second += uc;
			}
			else
			{
				totalUtilityEachFeat.insert({ currentInst.feature, uc });
			}
			// Save this instance to SN
			ObjWithoutCoord key = {
				currentInst.feature, // feature
				currentInst.instance, // instance
				uc};			

			// Create SN for this instance
			itSN = SN.find(key);
			if (itSN == SN.end())
			{
				SN.insert({ key, std::set<ObjWithoutCoord> {} });
			}
			// Check neighboring with other instances in the current cell
			for (auto const& checkInst : fiveCellInst)
			{
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{
					// Calucalate distance between the two instances
					dist = calculateDistanceTwoInstances(currentInst, checkInst);
					if (dist <= dist_thres) // the two instances have neighbor relationship
					{
						// Put them into SN
						// Calculate nomorlized utiltiy of this instances
						itmaxU = maxU.find(checkInst.feature);
						double uc = log10(checkInst.u) / log10(itmaxU->second);						
						// Save this uc to corresponding feature type
						ObjWithoutCoord starN = {
							checkInst.feature,
							checkInst.instance,
							uc};
						if (key < starN)
						{
							itSN = SN.find(key);
							if (itSN != SN.end())
							{
								itSN->second.insert(starN);
							}
							else
							{
								SN.insert({ key, std::set<ObjWithoutCoord> {starN} });
							}
						}
						else
						{
							itSN = SN.find(starN);
							if (itSN != SN.end())
							{
								itSN->second.insert(key);
							}
							else
							{
								SN.insert({ starN, std::set<ObjWithoutCoord> {key} });
							}
						}
					}
				}
			}
		}
		// Clear all element for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}
	return SN;
}


/**
* @bref This function calculates the utility of each feature
* utility of feature = all utility of instances of the feature / number of instances.
* @param utility of each feature
* @retval
*/
void getUtilityOfFeat(std::map<char, double>& totalUtilityEachFeat,
	std::map<char, int>& numInstFeat)
{
	std::map<char, double>::iterator ituti = totalUtilityEachFeat.begin();
	while (ituti != totalUtilityEachFeat.end())
	{
		ituti->second = ituti->second / (double)numInstFeat.find(ituti->first)->second;
		
		++ituti;
	}
}



/**
* @brief: This function calculates the total utility of the input data set.
* @param: total utilities of features
* @retval:
*/
double getUtilityOFDataset(std::map<char, double> & totalUtilityEachFeat)
{
	double totalUti = 0.0;
	for (auto const& item : totalUtilityEachFeat)
	{
		totalUti = totalUti + item.second;
	}

	return totalUti;
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



/**
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
template<typename T>
std::list<std::vector<T>> combinekOfn(std::vector<T>& c, int k)
{
	std::list<std::vector<T>> Cnk;
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
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
template<typename T>
std::vector<std::vector<T>> combinekOfnVec(std::vector<T>& c, int k)
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
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
std::vector<std::string> combinekOfnStr(std::string & c, int k)
{
	std::vector<std::string> Cnk;
	int n = c.size();
	std::vector<bool> v(n);

	std::fill(v.begin(), v.begin() + k, true);
	do {
		std::string onec;
		for (int i = 0; i < n; ++i) {
			if (v[i])
			{
				onec.push_back(c[i]);
				//onec += c[i];
			}
		}
		Cnk.push_back(onec);
	} while (std::prev_permutation(v.begin(), v.end()));
	// Return
	return Cnk;
}



/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombinationChar(const T& c, int combo)
{
	std::string result;
	int n = c.size();
	for (int i = 0; i < n; ++i) {
		if ((combo >> i) & 1)
			result.push_back(c[i]);
	}
	return result;
}


/**
* @brief: This function generate a subset of a vector.
* @param: c, k
* @r
*/
template<typename T>
std::vector<std::string> comboChar(const T& c, int k)
{
	std::vector<std::string> combination;

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
* @brief: This function print the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void genCandidatePatterns2(
	std::unordered_map<std::string, int> & Ckplus,
	std::map<char, double> & totalUtilityEachFeat,
	int k)
{	
	// Get all feature types
	std::string allFeats;
	for (auto const& item : totalUtilityEachFeat)
	{
		allFeats.push_back(item.first);
	}
	// Combo to get size k candidates
	int n = allFeats.size();
	std::vector<bool> v(n);
	std::fill(v.begin(), v.begin() + k, true);
	do {
		std::string onec;
		for (int i = 0; i < n; ++i) 
		{
			if (v[i])
			{
				onec.push_back(allFeats[i]);
				//onec += c[i];
			}
		}
		//std::cout << "onec: " << onec << endl;
		Ckplus.insert({ onec, 1 });

	} while (std::prev_permutation(v.begin(), v.end()));
}




/**
* @brief: This function groups a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<int>> groupInstanceByFeature(std::set<ObjWithoutCoord> instSet)
{
	std::map<char, std::vector<int>> gp;
	std::set<int> innverValueSet;
	std::vector<int> innverValueVec;
	for (auto const& inst : instSet)
	{		
		if (gp.find((inst.feature)) != gp.end()) // the fearture has already existed
		{
			// Retrieve the old value and update it
			innverValueVec = gp.find(inst.feature)->second;
			innverValueVec.push_back(inst.instance);

			// Delete duplicate elements invector			
			innverValueSet.insert(innverValueVec.begin(), innverValueVec.end());
			innverValueVec.clear();
			innverValueVec.insert(innverValueVec.end(), innverValueSet.begin(), innverValueSet.end());
			innverValueSet.clear();
			// Update
			gp[inst.feature] = innverValueVec;
		}
		else
		{
			// Directly put into gp
			innverValueVec.push_back(inst.instance);
			gp.insert({ inst.feature, innverValueVec });
		}
		// Clear all temporatory for next loop
		innverValueVec.clear();
	}
	return gp;
}


/**
* @brief: This function groups a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeature2(
	std::set<ObjWithoutCoord> instSet,
	std::string & remainCandidate)
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
				// Directly put into gp
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
void filterStarInstancesSize2(
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>& SIk,
	std::unordered_map<std::string, int>& Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>& groupSNByFeat,
	int k)
{
	char firstFeature; // The first feature of the candidate pattern
	std::vector<std::vector<ObjWithoutCoord>> innerValue;
	std::map<char, std::vector<ObjWithoutCoord>> groupInst; // group instances by feature types;
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgroupInst;
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>::iterator itSN;
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>::iterator itSIk;
	// Check each candidate and take back its table instance
	for (auto const& candidate : Ckplus)
	{
		// Retrieve the first feature		
		firstFeature = candidate.first[0];
		// Get the remain feature in the candidate
		std::string remainCandidate{ candidate.first.begin() + 1, candidate.first.end() };
		// Get the star neighbor of the first feature		
		itSN = groupSNByFeat.find(firstFeature);
		if (itSN != groupSNByFeat.end())// if the first feature has star neighborhoods
		{			
			for (auto const& star : itSN->second) // iterator each star neighborhoos
			{
				// Put instances of the other feature into innverValue
				if (star.second.size())
				{
					// Group a vector of instances by their features, return a hashmap <Feature: Instance>, then put they into innverValue
					groupInst = groupInstanceByFeature2(star.second, remainCandidate);
					if (groupInst.size() == remainCandidate.size())
					{
						// Put the instance of the first feature into the innerValue					
						std::vector<std::vector<ObjWithoutCoord>> tempInnerValue(k); // Save innver value in SIk structure
						tempInnerValue[0].push_back(star.first);
						// Update SIk		
						for (size_t t = 1; t < k; ++t)
						{
							itgroupInst = groupInst.find(candidate.first[t]);
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
							itSIk = SIk.find(candidate.first);
							if (itSIk != SIk.end())
							{
								// Update SIk fast version
								itSIk->second.insert(itSIk->second.end(), innerValue.begin(), innerValue.end());
							}
							else  // if SIk is empty
							{
								// Put as a new value
								SIk.insert({ candidate.first, innerValue });
							}
						}
					}
				}
				// Clear temporatory varibles to next iterator				
				groupInst.clear();				
				innerValue.clear();
			}
		}
	}
}



/**
* @brief: This function group a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeaturek(
	std::set<ObjWithoutCoord> instSet,
	std::string& remainCandidate)
{
	std::map<char, std::vector<ObjWithoutCoord>> gp;
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgp;

	// Loop each neihbor instance
	if (instSet.size() >= remainCandidate.size())
	{
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
	}
	return gp;

}



/**
* @brief: This function filter star instances of a candidate
* @param: totalInstEachFeat: instance number of each feature
* CIk: size k table instances
* 
* @retval: Nope
*/
void filterStarInstancesSizek(
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>& SIkplus,
	std::unordered_map<std::string, int>& Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> & groupSNByFeat,
	int k,
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>>& CIk)
{
	// Save real table instance of one candidate
	std::vector<std::vector<ObjWithoutCoord>> realTableInst;
	char firstFeature; // The first feature of the candidate pattern	
	std::vector<std::vector<ObjWithoutCoord>> innerValue;
	std::map<char, std::vector<ObjWithoutCoord>> groupInst; // group instances by feature types;
	std::map<char, std::vector<ObjWithoutCoord>>::iterator itgroupInst;
	// Save execpt first element of the candidate	
	std::vector<ObjWithoutCoord> checkRowkPlus;
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>::iterator itSN;
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>::iterator itSIk;
	// Check each candidate and take back its table instance
	for (auto const& candidate : Ckplus)
	{			
		// Retrieve the first feature		
		firstFeature = candidate.first[0];
		// Get the remain feature in the candidate
		std::string remainCandidate{ candidate.first.begin() + 1, candidate.first.end() };
		// Get the star neighbor of the first feature		
		itSN = groupSNByFeat.find(firstFeature);
		if (itSN != groupSNByFeat.end())// if the first feature has star neighborhoods
		{			
			for (auto const& star : itSN->second) // iterator each star neighborhoos
			{				
				// Group a vector of instances by their features, return a hashmap <Feature: Instance>, then put they into innverValue				
				if (star.second.size() > 1)
				{
					groupInst = groupInstanceByFeaturek(star.second, remainCandidate);
					if (groupInst.size() == remainCandidate.size())
					{
						// Put the instance of the first feature into the innerValue					
						std::vector<std::vector<ObjWithoutCoord>> tempInnerValue(k);
						tempInnerValue[0].push_back(star.first);
						// Update SIk		
						for (size_t t = 1; t < k; ++t)
						{
							itgroupInst = groupInst.find(candidate.first[t]);
							tempInnerValue[t].insert(tempInnerValue[t].end(),
								itgroupInst->second.begin(),
								itgroupInst->second.end());
						}
						// Cartesian Product tempInnerValue to generate row instances
						innerValue = cartesianProduct(tempInnerValue);
						tempInnerValue.clear();
						// Loop each row innerValue to check real row instances						
						if (CIk.find(remainCandidate) != CIk.end()) // If the sub pattern is prevalent
						{
							// Loop each row instance in coarseTabInstkPlus to check if it is a clique
							for (auto const& rowInstKPlus : innerValue)
							{								
								checkRowkPlus.insert(checkRowkPlus.end(), rowInstKPlus.begin() + 1, rowInstKPlus.end());								
								if (std::find(CIk.find(remainCandidate)->second.begin(),
									CIk.find(remainCandidate)->second.end(), checkRowkPlus)
									!= CIk.find(remainCandidate)->second.end())
								{
									// The current row instance is a clique
									realTableInst.push_back(rowInstKPlus);
								}
								checkRowkPlus.clear();
							}
						}
					}
				}
				// Clear temporatory varibles to next iterator				
				groupInst.clear();							
				innerValue.clear();
			}
		}
		// Put the real table instance
		if (realTableInst.size() > 0)
		{
			// Put into the final result
			SIkplus.insert({ candidate.first, realTableInst });
			realTableInst.clear();
		}
	}	
}



/*
*@brief: This function saves bits decimal places
*@param:
*@retval:
*/
double roundNum(double number, unsigned int bits)
{
	stringstream ss;
	ss << fixed << setprecision(bits) << number;
	ss >> number;
	return number;
}



/**
* @brief: This function calculate PIs of candidate patterns and filter prevalent patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         totalInstEachFeat: the number of instances of each features.
*         prev_thres: a minimum prevalence threshold
*         k: the size of patterns
* @retval: Pk: a hash map of prevalent patterns.
*/
std::vector<Pattern> computeUPI(
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>& CIkplus,
	std::map<char, double>& totalUtilityEachFeat,
	double & totalUti,
	std::map<char, int>& numInstFeat,
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>> & TAllk)
{
	std::vector<Pattern> CTIk;
	double sumuinst; // sum utility of all instances in the table instance
	double sumufeat; // sum utility of all features in the candidate

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
		// Save this table instances to all
		TAllk.insert({item.first, Tc});
		// Calculate utility ratio
		std::vector<double> PRs;
		sumuinst = 0.0;
		sumufeat = 0.0;
		int m = 0;
		for (auto const& instofFeat : Tc)
		{
			// Calculate utility of the candidate
			double u = 0.0;
			for (auto const& inst : instofFeat)
			{
				u += inst.u;				
				sumuinst += inst.u;
			}
			// Calculate pattern utility
			double uf = instofFeat.size() * totalUtilityEachFeat.find(item.first[m])->second;
			PRs.push_back(uf);
			m += 1;
		}
		// Calculate QC
		for (size_t t = 0; t < item.first.size(); ++t)
		{			
			sumufeat = sumufeat + totalUtilityEachFeat.find(item.first[t])->second * (double)numInstFeat.find(item.first[t])->second;
		}
		double QC = sumuinst / sumufeat;				
		// Calculate utility of candidate
		double utiOfc = 0.0f;
		for (auto const& pr : PRs)
		{
			utiOfc += pr;
		}			
		// Calculate pattern utility ratio
		double pur = utiOfc / totalUti;			
		// Save 2 decimal places		
		pur = roundNum(pur, 2);
		// Package one pattern
		Pattern patt{item.first, pur, QC};
		// Save this candidate
		CTIk.push_back(patt);
	}
	// return
	return CTIk;
}



/**
* @brief: This function select hight utility patterns.
* @param: Pk: a hash map of hight utility patterns
*			NonPK: a hash map of non hight utility patterns
*			CTIk:
*			prev_thres
* @retval: PK and NonPk
*/
void selectUtilityPatterns(
	std::map<std::vector<char>, double>& Pk,
	std::vector<std::vector<char>>& NonPk,
	std::map<std::vector<char>, std::vector<double>>& CTIk,
	double k,
	double prev_thres)
{
	if (k <= 2)
	{
		std::map<std::vector<char>, std::vector<double>>::iterator it = CTIk.begin();
		while (it != CTIk.end())
		{
			// Calculate the UPI
			double PI = *std::min_element(it->second.begin(), it->second.end());
			if (PI >= prev_thres)
			{
				Pk.insert({ it->first, PI });
			}			
			// Terminate
			++it;
		}
	}
	else
	{		
		std::map<std::vector<char>, std::vector<double>>::iterator it = CTIk.begin();
		while (it != CTIk.end())
		{			
			// Calculate the UPI
			double PI = *std::min_element(it->second.begin(), it->second.end());
			if (PI >= prev_thres)
			{
				Pk.insert({ it->first, PI });
			}
			// Get the non-high utility feature set (Definition 7)
			std::vector<char> oneNonUHU;
			// Get maximal UPR, need to have one UPR > min_prev
			double maxUPR = *std::max_element(it->second.begin(), it->second.end());		
			if (maxUPR >= prev_thres)
			{
				// Get the feature < min_prev
				for (size_t t = 0; t < k; ++t)
				{
					if (it->second[t] < prev_thres)
					{
						oneNonUHU.push_back(it->first[t]);
					}
				}
				// Only get size larger than two
				if (oneNonUHU.size() >= 2)
				{					
					NonPk.push_back(oneNonUHU);
				}
			}
			// Terminate
			++it;
		}		
	}
}


/**
* @brief: This function calculates EPUR.
* @param:
* @retval:
*/
double calculateEUPR(
	std::vector<std::string>& cfi,
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>>& TAllk,
	std::map<char, double>& totalUtilityEachFeat,
	double& totalUti)
{
	double sumu = 0.0;
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>>::iterator it;
	for (auto const& c : cfi)
	{
		it = TAllk.find(c);
		if (it != TAllk.end())
		{
			std::vector<double> PRs;
			int m = 0;
			for (auto const& instofFeat : it->second)
			{
				// Calculate pattern utility
				double uf = instofFeat.size() * totalUtilityEachFeat.find(c[m])->second;
				PRs.push_back(uf);
				m += 1;
			}			
			// Calculate utility of candidate
			double utiOfc = 0.0f;
			for (auto const& pr : PRs)
			{
				utiOfc += pr;
			}
			// Calculate pattern utility ratio
			double pur = utiOfc / totalUti;
			sumu += pur;
		}
	}
}


/**
* @brief: This function finds the non-prevalent candidates from non-size-k prevalent patterns.
* @param:allF: a set of all feature types
*			TAllk: the participating instances of size k patterns
*			s: the top-s subset candidates
* @retval: NonPk: a set of non-prevalent candidates
*/
void findPruneCandidate(
	std::map<std::string, double>& NonPk,
	std::set<char>& allF,
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>>& TAllk,
	int s,
	std::map<char, double>& totalUtilityEachFeat, 
	double & totalUti,
	double & prev_thres,
	std::unordered_map<std::string, int>& pruneCands)
{
	for (auto const& nonc : NonPk)
	{
		// 1. Get the different between two sets		
		std::string Fsubc;
		std::set_difference(allF.begin(), allF.end(),
			nonc.first.begin(), nonc.first.end(),
			std::back_inserter(Fsubc));		
		// 2. Get vss of c and fi
		double epurc = 0.0;
		std::vector<std::string> prunc;

		for (auto const& fi : Fsubc)
		{
			std::vector<std::string> cfi;
			// Get top-s subset that include fi
			std::string c{ nonc.first.begin(), nonc.first.end() };
			int k = nonc.first.size();
			for (int i = k - 1; i >= 2; --i)
			{
				// Get sub of combok
				std::vector<std::string> combok = combinekOfnStr(c, i);
				// Add fi in Fsubc
				for (auto j = 0; j < combok.size(); ++j)
				{
					combok[j].push_back(fi);
					std::sort(combok[j].begin(), combok[j].end());
				}
				// Sort combok by size and alphabet
				std::sort(combok.begin(), combok.end(), 
					[](const std::string& lhs, const std::string& rhs) 
					{
						return lhs.size() >= rhs.size();
					});
				// check combok containing top-s
				if (cfi.size() == 0)
				{
					if (combok.size() >= s)
					{
						cfi.insert(cfi.end(), combok.begin(), combok.begin() + s);
						break;
					}
					else
					{
						cfi.insert(cfi.end(), combok.begin(), combok.end());
					}
				}
				else 
				{
					int n = 0;
					while (cfi.size() < s)
					{
						cfi.push_back(combok[n]);
						++n;
					}
				}
			}
			// Calculate EUPR
			double epur = calculateEUPR(cfi, TAllk, totalUtilityEachFeat, totalUti);			
			double fepur = nonc.second + epur;
			epurc += fepur;
			// Get super set of c
			c.push_back(fi);
			std::sort(c.begin(), c.end());
			prunc.push_back(c);
		}
		// Check for super set
		if (epurc < prev_thres)
		{
			for (auto const& pc : prunc)
			{
				pruneCands.insert({ pc, 1 });
			}			
		}
	}
}



/**
* @brief: This function selects hight utility patterns.
* @param: Pk: a hash map of hight utility patterns
*			NonPK: a hash map of non hight utility patterns
*			CTIk:
*			prev_thres
* @retval: PK and NonPk
*/
void selectUtilityPatterns(
	std::vector<Pattern>& Pk,
	std::map<std::string, double>& NonPk,
	std::vector<Pattern>& CTIk,
	double k,
	double prev_thres)
{
	std::vector<Pattern>::iterator it = CTIk.begin();	
	while (it != CTIk.end())
	{
		// Get high utility
		if ((*it).upr >= prev_thres)
		{
			Pk.push_back(*it);
		}
		else
		{
			// Get non-high utility patterns
			NonPk.insert({(*it).c, (*it).upr});
		}
		// Terminate
		++it;
	}
}



/**
* @brief: This function pruns candidates
* @param:
* @retval: PK and NonPk
*/
void candidatePrun(
	std::unordered_map<std::string, int>& Ckplus,
	std::unordered_map<std::string, int>& pruneCands)
{
	// check can in Ckplus is in prenneCands
	//std::cout << "numof prune:" << pruneCands.size() << endl;
	if (pruneCands.size())
	{
		//std::cout << "prune.." << endl;
		std::unordered_map<std::string, int>::iterator itprune;

		std::unordered_map<std::string, int>::iterator it = Ckplus.begin();
		while (it != Ckplus.end())
		{
			itprune = pruneCands.find(it->first);
			if (itprune != pruneCands.end())
			{
				it = Ckplus.erase(it);
			}
			else
			{
				++it;
			}
		}
	}
}



/**
* @brief: This function gets the maximal size patterns in the result.
* @param: PkAll all prevalent patterns
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
int getMaxSizeofPatterns(std::vector<Pattern>& PkAll)
{
	std::vector<int> sizePats;
	std::vector<Pattern>::iterator itAllPats = PkAll.begin();
	while (itAllPats != PkAll.end())
	{
		sizePats.push_back((*itAllPats).c.size());
		++itAllPats;
	}
	int maxSize = *std::max_element(sizePats.begin(), sizePats.end());

	return maxSize;
}


/**
* @brief: This function calculate the Qc of top k patterns.
* @param: Qctopk
*		topk
* @retval: nope
*/
void calculateQcofTopk(
	std::vector<Pattern>& PkAll,
	std::vector<double>& Qctopk,
	std::vector<int>& topk)
{
	int i = 0;
	int numPatts = PkAll.size();

	for (auto const& tk : topk) // get top-k number
	{
		// Get utility of the tk-th element
		double uptitk = PkAll[tk].upr;
		// Get Qc of top-k patterns
		double qc = 0.0;
		for (size_t j = 0; j < numPatts; ++j)
		{
			// Sum QC
			if (PkAll[j].upr >= uptitk)
			{
				qc += PkAll[j].Qc;
			}
		}
		// Save
		Qctopk[i] = qc;
		++i;
	}
}



/*
*@brief: This function calculates the Qc of top-k patterns of each size
*@param: HUPk: all prevalent patterns
*  topkpersize : top-k
*@retval: QcTopkbySize
*/
void calculateQCofTopkofEachSize(
	std::vector<Pattern>& PkAll,
	int topkpersize,
	std::map<int, double>& QcTopkbySize)
{
	//  1 Classify patterns by sizes
	std::map<int, std::vector<Pattern>> allPatbySize;
	std::map<int, std::vector<Pattern>>::iterator itallPatbySize;
	for (auto const& p : PkAll)
	{
		itallPatbySize = allPatbySize.find(p.c.size());
		if (itallPatbySize != allPatbySize.end())
		{
			itallPatbySize->second.push_back(p);
		}
		else
		{
			allPatbySize.insert({ p.c.size(), std::vector<Pattern> {p} });
		}
	}
	//  2 Sort patterns by uti	
	while (itallPatbySize != allPatbySize.end())
	{
		std::sort(itallPatbySize->second.begin(), itallPatbySize->second.end());
		// Terminate
		++itallPatbySize;
	}
	//  3 Get top-k of each size
	for (auto const& item : allPatbySize)
	{
		double avgsumQc = 0.0;

		if (item.second.size() >= topkpersize)
		{
			// Get the uti of the k-th element
			double pikth = item.second[topkpersize].upr;
			int nk = 0;
			for (auto const& pat : item.second)
			{
				if (pat.upr >= pikth)
				{
					avgsumQc += pat.Qc;
					nk += 1;
				}
			}
			// Calculate average
			avgsumQc = avgsumQc / nk;
			// Save 2 decimal places
			avgsumQc = roundNum(avgsumQc, 2);
			// Save
			QcTopkbySize.insert({ item.first, avgsumQc });
		}
		else // The number of size k is small than topk, direct calculate
		{
			for (auto const& pat : item.second)
			{
				avgsumQc += pat.Qc;
			}
			// Calculate average
			avgsumQc = avgsumQc / item.second.size();
			// Save 2 decimal places
			avgsumQc = roundNum(avgsumQc, 2);
			// Save
			QcTopkbySize.insert({ item.first, avgsumQc });
		}
	}
}



/**
* @brief: This function calculate the Qc of top k patterns by size.
* @param: PkAll
*		PkAll
* @retval: nope
*/
void calculateQCofTopkClassifybySize(
	std::vector<Pattern>& PkAll,
	int numtopk,
	std::map<int, double>& QcbySize)
{
	int numPatts = PkAll.size();
	std::map<int, std::vector<Pattern>> PtopkbySize;
	std::map<int, std::vector<Pattern>>::iterator itPtopkbySize;
	if (numtopk < PkAll.size())
	{
		// Get the uti of the numtopk-th pattern
		double upinumtopk = PkAll[numtopk].upr;
		// Check top k patterns
		for (size_t t = 0; t < numPatts; ++t)
		{
			if (PkAll[t].upr >= upinumtopk)
			{
				itPtopkbySize = PtopkbySize.find(PkAll[t].c.size());
				if (itPtopkbySize != PtopkbySize.end())
				{
					itPtopkbySize->second.push_back(PkAll[t]);
				}
				else
				{
					PtopkbySize.insert({ PkAll[t].c.size(), std::vector<Pattern>{ PkAll[t]} });
				}
			}
		}
	}
	// Save average QC of patterns in a certain size
	for (auto const& item : PtopkbySize)
	{
		// Calculate total utility of pattern
		double uofsize = 0.0;
		for (auto const& pt : item.second)
		{
			uofsize += pt.Qc;
		}
		// Calculate average QC
		uofsize = uofsize / (double)item.second.size();
		// Save 2 decimal places
		uofsize = roundNum(uofsize, 2);
		// Save the result
		QcbySize.insert({ item.first, uofsize });
	}
}



/**
* @brief: This function classifies patterns by sizes and get the maximal size of patterns
* @param: PkAll
*		PkAll
* @retval: nope
*/
void classifyPatternAndGetMaxSizeofPatterns(
	std::vector<Pattern> & PkAll,
	std::map<int, int> & sizePats,
	int& maxSizePatt)
{
	std::vector<int> maxsize;
	int onek;
	std::map<int, int>::iterator itsizePats;

	std::vector<Pattern>::iterator itAllPats = PkAll.begin();
	while (itAllPats != PkAll.end())
	{
		// Get the size of the current pattern
		onek = (*itAllPats).c.size();

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

	maxSizePatt = *std::max_element(maxsize.begin(), maxsize.end());
}


