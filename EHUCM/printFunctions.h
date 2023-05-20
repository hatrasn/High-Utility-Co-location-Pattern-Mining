#include <iostream>
#include <map>>
#include <vector>
#include<unordered_map>
#include <set>

#include "Object.h"

#include"tree.hh"

#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

using namespace std;

/**
* @brief: This function print the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int> totalInstEachFeat);

/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printGrid(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid);


/**
* @brief: This function prints the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printSN();


/**
* @brief: This function prints sorted tuf.
* @param: sortF
* @retval: Nope
*/
void printSortTUF();


/**
* @brief: This function prints a set of instances.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printInstances(std::vector<ObjWithoutCoord > SN);


/**
* @brief: This function print the group star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> SN);


/**
* @brief: This function prints feature-object neighbor sets
* @param: std::map<ObjWithoutCoord, std::map<char, std::set<ObjWithoutCoord>>> FONS
* @retval: Nope
*/
void printFONS();


/**
* @brief: This function prints FONT
* @param: FONT
* @retval: Nope
*/
void printFONT();



/**
* @brief: This function prints candidate tree
* @param:tree < std::vector<char>>  CandTree
* @retval: Nope
*/
void printCandTree();



/**
* @brief: This function prints candidates
* @param: c: a candiate
* @retval: Nope
*/
void printCandidate(std::vector<char> c);


/**
* @brief: This function prints all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern();


/**
* @brief: This function prints all participating instances of c
* @param:  CSHBS
* @retval: Nope
*/
void printCSHBS(std::unordered_map<char, std::set<ObjWithoutCoord>> CSHBS);


/**
* @brief: This function prints all row instances of c
* @param:  rowInstSet
* @retval: Nope
*/
void printRowInstSet(std::vector<std::vector<ObjWithoutCoord>> rowInstSet);




#endif
