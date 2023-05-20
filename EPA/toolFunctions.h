
#ifndef TOOLFUNCTIONS_H
#define TOOLFUNCTIONS_H


#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>


#include <boost/algorithm/string.hpp>


#include "Object.h"
#include "pattern.h"



/**
* @brief: This function computes the utility of each feature type
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its utility.
*/
std::map<char, double> computeUtility(std::vector<std::vector<std::string>>& dataList);


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
	std::map<char, int>& numInstFeat,
	std::set<char>& allF);



/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(
	std::vector<std::vector<std::string> >& dataList, 
	double dist_thres);


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>>& grid,
	double dist_thres,
	std::unordered_map<char, double>& maxU,
	std::map<char, double>& totalUtilityEachFeat);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
double calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst);


/**
* @brief: This function calculates the total utility of the input data set.
* @param: total utilities of features
* @retval:
*/
double getUtilityOFDataset(std::map<char, double>& totalUtilityEachFeat);


/**
* @bref This function calculates the utility of each feature
* utility of feature = all utility of instances of the feature / number of instances.
* @param utility of each feature
* @retval 
*/
void getUtilityOfFeat(std::map<char, double>& totalUtilityEachFeat,
	std::map<char, int>& numInstFeat);


/**
* @brief: This function groups star neighbor instances by the feature type of the key.
* @param: SN: a hashmap of star neighborhoods*
* @retval: groupSNByFeat: a hash map that stores as <feature, <instance: neighborhoods>>.
*/
std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupStarNeighByFeatures(
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> & SN);



/**
* @brief: This function print the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void genCandidatePatterns2(
	std::unordered_map<std::string, int>& Ckplus,
	std::map<char, double>& totalUtilityEachFeat,
	int k);



/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombinationChar(const T& c, int combo);


/**
* @brief: This function generate a subset of a vector.
* @param: c, k
* @r
*/
template<typename T>
std::vector<std::string> comboChar(const T& c, int k);


/**
* @brief: This function filter star instances of a candidate
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void filterStarInstancesSizek(
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>& SIkplus,
	std::unordered_map<std::string, int>& Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>& groupSNByFeat,
	int k,
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>>& CIk);


/**
* @brief: This function filters star instances of a candidate
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void filterStarInstancesSize2(
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>>& SIk,
	std::unordered_map<std::string, int>& Ckplus,
	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>>& groupSNByFeat,
	int k);



/**
* @brief: This function groups a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeature2(
	std::set<ObjWithoutCoord> instSet,
	std::string & remainCandidate);


/**
* @brief: This function group a set of instances by their feature types.
* @param: instSet: a set of instances.
* @retval: gp: a hash map that stores instances are grouped by feature type.
*/
std::map<char, std::vector<ObjWithoutCoord>> groupInstanceByFeaturek(
	std::set<ObjWithoutCoord> instSet,
	std::string & remainCandidate);



/**
* @brief: This function generate cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::vector<std::vector<ObjWithoutCoord>>& v);


/*
*@brief: This function saves bits decimal places
*@param:
*@retval:
*/
double roundNum(double number, unsigned int bits);



/**
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
template<typename T>
std::list<std::vector<T>> combinekOfn(std::vector<T>& c, int k);


/**
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
template<typename T>
std::vector<std::vector<T>> combinekOfnVec(std::vector<T>& c, int k);


/**
* @brief: This function creates all possible k combinations of n items
* @param: c: a vector
*           k:
* @retval: Cnk
*/
std::vector<std::string> combinekOfnStr(std::string & c, int k);


/**
* @brief: This function calculates EPUR.
* @param:
* @retval:
*/
double calculateEUPR(
	std::vector<std::string>& cfi,
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>>& TAllk,
	std::map<char, double>& totalUtilityEachFeat,
	double& totalUti);


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
	double& totalUti,
	double& prev_thres,
	std::unordered_map<std::string, int>& pruneCands);



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
	double& totalUti,
	std::map<char, int>& numInstFeat,
	std::map<std::string, std::vector<std::set<ObjWithoutCoord>>>& TAllk);


/**
* @brief: This function select hight utility patterns.
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
	double prev_thres);


/**
* @brief: This function pruns candidates
* @param:
* @retval: PK and NonPk
*/
void candidatePrun(
	std::unordered_map<std::string, int>& Ckplus,
	std::unordered_map<std::string, int>& pruneCands);



/**
* @brief: This function gets the maximal size patterns in the result.
* @param: PkAll all prevalent patterns
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
int getMaxSizeofPatterns(std::vector<Pattern>& PkAll);


/**
* @brief: This function calculate the Qc of top k patterns.
* @param: Qctopk
*		topk
* @retval: nope
*/
void calculateQcofTopk(
	std::vector<Pattern>& PkAll,
	std::vector<double>& Qctopk,
	std::vector<int>& topk);



/*
*@brief: This function calculates the Qc of top-k patterns of each size
*@param: HUPk: all prevalent patterns
*  topkpersize : top-k
*@retval: QcTopkbySize
*/
void calculateQCofTopkofEachSize(
	std::vector<Pattern>& PkAll,
	int topkpersize,
	std::map<int, double>& QcTopkbySize);


/**
* @brief: This function calculate the Qc of top k patterns by size.
* @param: PkAll
*		PkAll
* @retval: nope
*/
void calculateQCofTopkClassifybySize(
	std::vector<Pattern>& PkAll,
	int numtopk,
	std::map<int, double>& QcbySize);



/**
* @brief: This function classifies patterns by sizes and get the maximal size of patterns
* @param: PkAll
*		PkAll
* @retval: nope
*/
void classifyPatternAndGetMaxSizeofPatterns(
	std::vector<Pattern> & PkAll,
	std::map<int, int>& sizePats,
	int& maxSizePatt);



#endif