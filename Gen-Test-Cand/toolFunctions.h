#include <iostream>
#include <map>
#include <vector>
#include "Object.h"
#include <set>


#ifndef TOOLFUNCTIONS_H
#define TOOLFUNCTIONS_H




/**
* @brief: This function computes the utility of each feature type
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its utility.
*/
std::map<char, float> computeUtility(
	std::vector<std::vector<std::string>>& dataList);


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
* @brief: This function generates  size k candidates
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::vector<std::vector<char>> genCandidatePatterns(
	std::map<char, float>& totalUtilityEachFeat,
	int k);


/*
*@brief: This function generates all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::vector<char> getCombinationChar(const T& c, int combo);

/**
* @brief: This function generates a subset of a vector.
* @param: c, k
* @retval: a vector
*/
template<typename T>
std::vector<std::vector<char>> comboChar(const T& c, int k);


/**
* @brief: This function filters star instances of a candidate
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterStarInstances(
	std::vector<std::vector<char>> Ck,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> & groupSNByFeat,
	int k);


/**
* @brief: This function groups a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeature(
	std::set<ObjWithoutCoord> instSet,
	std::vector<char>& remainCandidate);


/**
* @brief: This function generate cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::vector<std::vector<ObjWithoutCoord>>& v);


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
	float w2);


/**
* @brief: This function selects hight utility patterns.
* @param: Pk: a hash map of hight utility patterns
*			NonPK: a hash map of non hight utility patterns
*			CTIk:
*			prev_thres
* @retval: PK and NonPk
*/
void selectUtilityPatterns(
	std::map<std::vector<char>, float>& Pk,
	std::vector<std::pair<std::vector<char>, std::vector<char>>> & NonPk,
	std::map<std::vector<char>, std::vector<float>>& CTIk,
	float k,
	float prev_thres);



/**
* @brief: This function finds real row instances of candidate patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         CIk: a hash map of table instances of size k candiate patterns.
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> filterCliqueInstances(
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> & CIkplus,
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> & CIk);



/**
* @brief: This function pruns candidates
* @param:
* @retval: PK and NonPk
*/
void candidatePrun(
	std::vector<std::vector<char>>& Ckplus,
	std::vector<std::pair<std::vector<char>, std::vector<char>>>& NonPk,
	int k);



#endif