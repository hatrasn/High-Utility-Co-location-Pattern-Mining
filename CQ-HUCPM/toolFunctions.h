#include <iostream>
#include <map>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include<string>


#include "Object.h"

#include"tree.hh"

#include <boost/algorithm/string.hpp>


#ifndef TOOLFUNCTIONS_H
#define TOOLFUNCTIONS_H


// Definition global variable
extern std::unordered_map<std::string, std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>>> coHM;
extern clock_t timeBuildHash;
extern double totalTimeBuildHash;


/**
* @brief: This function gets the maximal and minimal values of the utility of each feature type
* @param: dataList: an input dataset;
* @retval: maxminU: <A, [minU, maxU]>, <B, [minU, maxU]>, ...
*/
void getMinMaxUtility(
	std::vector<std::vector<std::string>>& dataList,
	std::unordered_map<char, std::vector<float>>& maxminU);


/**
* @brief This function makes a grid on the input dataset
* @param dataList: an input dataset
*		dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it
*/
void makeGrid(
	std::vector<std::vector<std::string>> & dataList, 
	float dist_thres,
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> & grid);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst);


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>>& grid,
	float dist_thres,
	std::unordered_map<char, std::vector<float>>& minmaxU,
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc>& SN,
	std::map<ObjWithoutCoord, int>& V,
	std::unordered_map<char, float>& totalUtilityFeats);


/**
* @brief: This function constructs co-location Hash map
* @param: onecl: one row instance
* @retval: coloHM: a co-location hashmap storing all row instances
*/
void createCoLoHM(std::vector<ObjWithoutCoord> const& onecl);


/**
* @brief: This function generates all cliques
* @param: SN: star neighbors
* @retval: CLs: a vector of all cliques
*/
void generateCliques(
	std::map<ObjWithoutCoord, int>& V,
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc>& SN);



/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombination(const T& c, int combo);

template<typename T>
std::vector<std::string> combo(const T& c, int k);


/*
*@brief: This function calculates the utility of the current pattern
*@param: pattern: the current pattern which need to calculate its utility.
*       coHM: a co-location pattern hash map
*@retval: UT: the utility of the current pattern
*/
float calcUtilityWithSuperPatterns2Phases(
	std::string pattern,
	std::unordered_map<char, float>& totalUtilityFeats,
	int sizeofPatttern);


/**
* @brief: This function calculates the utility of each pattern
* @param: V and SN
* @retval: coloHM: a co-location hashmap storing all row instances
*/
void calculateUtility(std::unordered_map<char, float>& totalUtilityFeats,
	std::unordered_map<std::string, float>& Pk);


/*
*@brief: This function filters high utility patterns
*@param: Pk: all patterns with their utility values.
*       high_thres: threshold given by users
*@retval: HUPk: the high utility patterns
*/
void filterHighUtilityPatterns(
	std::unordered_map<std::string, float>& Pk,
	std::unordered_map<std::string, float>& HUPk,
	float high_thres);



#endif