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
extern std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN;
extern tree<std::vector<char>> CandTree; // generate candidates
extern tree<ObjWithoutCoord> FONT;
extern tree<ObjWithoutCoord>::iterator rootFONT;
extern std::vector<char> sortF;  // the set of features that are sorted by their utilities
extern std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; // table instances of size k patterns
extern std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIkplus; // table instances of size k patterns
extern std::map<char, double> vuf; // the utility of each feature
extern std::map<std::string, double> PkAll;
extern double US; // utility of all data set
extern std::map<ObjWithoutCoord, std::map<char, std::set<ObjWithoutCoord>>> FONS;
extern std::vector<std::vector<ObjWithoutCoord>> realTableInst;

//extern clock_t startFilHUCMP;
extern double totalGenerateParticipatingObjects;
extern double totalGenCand;
extern double totalTimeCalUPIAndFilPrevCoLoPat;


/**
* @brief This function makes a grid on an input dataset
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
void genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid, 
	float dist_thres, 
	std::unordered_map<char, double> maxU,
	std::map<char, double> & tuf);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst);


/**
* @bref This function calculates the utility of the input data set
* @param tuf
* @retval US
*/
void calculateUS(std::map<char, double> tuf);


/**
* @brief: This function gets the maximal and minimal values of the utility of each feature type
* @param: dataList: an input dataset;
* @retval: maxminU: <A, [minU, maxU]>, <B, [minU, maxU]>, ...
*/
void getMaxUtility(
	std::vector<std::vector<std::string>>& dataList,
	std::unordered_map<char, double>& maxU);


/**
* @bref This function sorts features by their utilities.
* @param tuf: utilities of features
* @retval std::unordered_map<char, double>& sortF
*/
void sortFeatureByUtility(std::map<char, double> tuf);



/**
* @bref This function builds the candidate tree
* @param sortF: a set of sorted features
* @retval tree<std::vector<char>> CandTree;
*/
std::vector<std::vector<char>> genCands(std::vector<std::vector<char>> oldCands);


/**
* @bref This function finds HUCP
* @param c: a candidates
* @retval 
*/
void searchHUCP(std::vector<std::vector<char>> c, int k, double minutil);


/**
* @bref This function computes feature-object neighbor sets.
* @param SN: star neighbor
* @retval FONS: neighbors that are classified by features and instances
*/
void buildFONS(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN);


/**
* @bref This function builds feature-object neighbor tree (Figure 3-2 in Reference[2])
* @param std::map<ObjWithoutCoord, std::map<char, std::set<ObjWithoutCoord>>> FONS: in global level
* @retval FONT: in globel level
*/
void buildFONT();


/**
* @bref This function executes searching to get high utility patterns
* @param c, k
* @retval
*/
bool searchPatterns(std::vector<char> c, float minutil);


/**
* @bref This function calculates the utility of candidate c
* @param c: candidate
*		CSHBS: the participating instances of features in c
* @retval: uc: the utility of c
*/
double calculateUtilityOfc(std::vector<char> c, std::unordered_map<char, std::set<ObjWithoutCoord>> CSHBS);


/**
* @bref This function calculates the pattern utility ratio of c
* @param uc: the utility of c
* @retval: lambdac
*/
double calculateLambda(double uc);


/**
* @bref This function calculates luc to prune
* @param uc: the utility of c
* @retval: luc
*/
double calculateLuc(std::vector<char> c, double uc);


/**
* @bref This function calculates luc to prune
* @param uc: the utility of c
* @retval: luc
*/
double calculateUbc(std::vector<char> c);


/**
* @bref This function calculates euc
* @param ubc and lambdac
* @retval: euc
*/
double calculateEuc(double ubc, double lambdac);


/**
* @bref This function generates the participating instances of c
* @param FONT, c
* @retval: CSHBS
*/
void generateParticipatingInstance(std::vector<char> c, std::unordered_map<char, std::set<ObjWithoutCoord>>& CSHBS);



/**
* @bref This function checks wwhether a vector is a subset of another
* @param A, B: two vectors
* @retval: bool
*/
template <typename T>
bool IsSubset(std::vector<T> A, std::vector<T> B);


/**
* @brief: This function generates cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::unordered_map<char, std::vector<ObjWithoutCoord>> tempInner);


/**
* @brief: This function validates real row instances
* @param: innervalue: candidate row instances
* std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; // table instances of size k patterns
* @retval: realTableInst: real row instances, global level, std::vector<std::vector<ObjWithoutCoord>> rowInstSet;
*/
void validateRealRowInst(std::vector<char> c, std::vector<std::vector<ObjWithoutCoord>> innervalue);



#endif