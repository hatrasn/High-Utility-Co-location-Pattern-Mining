#include <iostream>
#include <map>
#include <vector>
#include "Object.h"
#include <set>


#ifndef TOOLFUNCTIONS_H
#define TOOLFUNCTIONS_H


/**
* @brief: This function counts the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
std::map<char, int> countNumberInstance(std::vector<std::vector<std::string>> dataList);

/**
* @brief This function makes a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(
	std::vector<std::vector<std::string> >& dataList, 
	float dist_thres);


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>>& grid,
	float dist_thres);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst);


/**
* @brief: This function groups star neighbor instances by the feature type of the key.
* @param: SN: a hashmap of star neighborhoods*
* @retval: groupSNByFeat: a hash map that stores as <feature, <instance: neighborhoods>>.
*/
std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupStarNeighByFeatures(
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> & SN);



/**
* @brief: This function generates candidates
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::vector<std::vector<char>> genCandidatePatterns(std::map<std::vector<char>, float> Pk, int k);


/**
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
template<typename T>
std::vector<std::vector<T>> combinekOfn(std::vector<T>& c, int k);



/**
* @brief: This function filters star instances of a candidate
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> filterStarInstances(
	std::vector<std::vector<char>> Ck,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> & groupSNByFeat,
	int k);


/**
* @brief: This function filters star instances of a candidate
* @param: totalInstEachFeat: the instance number of each feature
* @retval: Nope
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterStarInstancesSize2(
	std::vector<std::vector<char>> Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>& groupSNByFeat,
	int k);


/**
* @brief: This function group a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeature3(
	std::set<ObjWithoutCoord> instSet,
	std::vector<char>& remainCandidate);



/**
* @brief: This function generates cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::vector<std::vector<ObjWithoutCoord>> v);


/**
* @brief: This function calculates PIs of candidate patterns and filters prevalent patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         totalInstEachFeat: the number of instances of each features.
*         prev_thres: a minimum prevalence threshold
*         k: the size of patterns
* @retval: Pk: a hash map of prevalent patterns.
*/
std::map<std::vector<char>, float> selectPrevalentPatterns(
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>>& CIkplus,
	std::map<char, int>& totalInstEachFeat,
	float prev_thres,
	int k);


/**
* @brief: This function selects coarse prevalent patterns.
* @param: SIK: a hash map of table instances
* @retval: CIk: a hash map of table instances which satify coarse min_prev
*/
void selectCoarsePrevalentPatterns(
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIkplus,
	std::map<char, int> & totalInstEachFeat,
	float prev_thres,
	int k);


/**
* @brief: This function finds real row instances of candidate patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         CIk: a hash map of table instances of size k candiate patterns.
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterCliqueInstances2(
	std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> & SIkplus,
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk);



#endif