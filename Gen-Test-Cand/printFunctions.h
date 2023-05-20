#include <iostream>
#include <map>>
#include <vector>
#include "Object.h"
#include <set>

#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

using namespace std;

/**
* @brief: This function prints the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int> totalInstEachFeat);

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
void printCandidates(std::vector<std::vector<char>> Ck);


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
void printPrevalentPattern(std::map<vector<char>, float> Pk);





#endif
