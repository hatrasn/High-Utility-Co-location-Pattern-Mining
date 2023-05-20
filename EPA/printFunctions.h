#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

#include <iostream>
#include <map>>
#include <vector>
#include <set>
#include <unordered_map>

#include "Object.h"
#include "pattern.h"



using namespace std;



/**
* @brief: This function prints the utility of each feature type.
* @param: totalUtilityEachFeat: utility of each feature
* @retval: Nope
*/
void printUtilityFeature(std::map<char, double>& totalUtilityEachFeat);



/**
* @brief: This function prints the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int>& totalInstEachFeat);

/**
* @brief: This function prints the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printGrid(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid);


/**
* @brief: This function prints the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printSN(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN);


/**
* @brief: This function prints the group star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> SN);

/**
* @brief: This function prints all candidate patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printCandidates(std::unordered_map<std::string, int> Ck);


/**
* @brief: This function prints candidate patterns and their table instances
* @param: SIk: a hash map of candidate pateterns and their table instasnces
* @retval: Nope
*/
void printStarInstance(std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> SIk);


/**
* @brief: This function prints all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern(std::vector<Pattern> Pk);



/**
* @brief: This function prints Qc of topk patterns
* @param: topk: number of top-k
*			Qc: the quality of each top-k pattern set
* @retval: Nope
*/
void printQc(std::vector<int> topk, std::vector<double> Qc);



/**
* @brief: This function prints average Qc of topk patterns by sizes
* @param: QcbySize
* @retval: Nope
*/
void printAvgQCbySize(std::map<int, double> QcbySize);



/**
* @brief: This function prints QC of top-k patterns by each size
* @param: QcbySize
* @retval: Nope
*/
void printQCofTopkofSize(std::map<int, double>& QcTopkbySize);



/**
* @brief: This function prints patterns by sizes
* @param: sizePats
* @retval: Nope
*/
void printPattbySize(std::map<int, int>& sizePats);





#endif
