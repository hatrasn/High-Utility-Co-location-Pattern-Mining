#include <iostream>
#include <map>>
#include <vector>
#include "Object.h"
#include <set>
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
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printSN(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN);


/**
* @brief: This function print the group star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> SN);

/**
* @brief: This function print the row instance tree of an instance. This is a pre_order_iterator tree
* @param: tree<string>: a row instance tree.
* @retval: Nope
*/
void printITree(tree<ObjWithoutCoord> ITree);


/**
* @brief: This function print all obtained cliques
* @param: Cls: a vector of all cliques
* @retval: Nope
*/
void printCliques(std::vector<std::set<ObjWithoutCoord>> Cls);


/**
* @brief: This function print all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern(std::map<std::string, float> Pk);


/**
* @brief: This function prints the hash map.
* @param: chash: a hash map of all cliques
* @retval: Nope
*/
void printChash();




#endif
