#include <iostream>
#include <map>
#include <ctime>
#include<unordered_map>
#include <vector>
#include "Object.h"
#include <set>
#include<string>

#include"tree.hh"

#include <boost/algorithm/string.hpp>

#ifndef TOOLFUNCTIONS_H
#define TOOLFUNCTIONS_H


// Definition global variable
extern std::map<char, int> totalInstEachFeat;
extern std::unordered_map<std::string, std::unordered_map<char, std::set<ObjWithoutCoord>>> chash;
extern float PI;

extern clock_t timeCreateChash; // time for generating candidates
extern clock_t timeTimeCalPIs; // time for calculating PIS

extern double totalTimeCreateChash;
extern double totalTimeCalPIs;


/**
* @brief: This function count the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
void countNumberInstance(std::vector<std::vector<std::string>> dataList);

/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(std::vector<std::vector<std::string> > dataList, float dist_thres);


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> genStarNeighborhoods(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid, float dist_thres);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst);


/**
* @brief: This function generates all cliques.
* @param: SN: star neighbors
* @retval: CLs: a vector of all cliques
*/
void generateCliques(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN);



/**
* @brief: This function generates all cliques.
* @param: SN: star neighbors
* @retval: CLs: a vector of all cliques
*/
void generateCliques2(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN);


/**
* @brief: This function creates the chash to store all cliques.
* @param: Cls: all cliques
* @retval: chash: a hash map
*/
void createChash(std::set<ObjWithoutCoord> Cls);


/**
* @brief: This function creates the chash to store all cliques.
* @param: Cls: all cliques
* @retval: chash: a hash map
*/
void createChash2(std::vector<std::set<ObjWithoutCoord>> Cls);


/**
* @brief: This function calculates PIs and filter prevalent patterns.
* @param: chash: a hash map of all cliques
* @retval: PkAll: all prevalent patterns
*/
std::map<std::string, float> CalPIandFilterPatterns(float prev_thres);


// Function to sort the map according 
// to value in a (key-value) pairs 
std::unordered_map<std::string, int> sortMap(std::unordered_map<std::string, int>& M);


/**
* @brief: This function calculates PIs and filter prevalent patterns.
* @param: chash: a hash map of all cliques
* @retval: PkAll: all prevalent patterns
*/
std::map<std::string, float> CalPIandFilterPatterns2(float prev_thres);

/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombination(const T& c, int combo);

template<typename T>
std::vector<std::string> combo(const T& c, int k);


/**
* @brief: This function calculates PIs.
* @param: currCandidate: one candidate, chash: a hash map
* @retval: PI: the pi of the candidate
*/
void calculatePI(std::string currCandidate);


#endif