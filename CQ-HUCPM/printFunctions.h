#include <iostream>
#include <map>>
#include <vector>
#include <set>
#include <unordered_map>


#include "Object.h"

#include"tree.hh"


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
* @brief: This function prints the row instance tree of an instance. This is a pre_order_iterator tree
* @param: tree<string>: a row instance tree.
* @retval: Nope
*/
void printITree(tree<ObjWithoutCoord> ITree);


/**
* @brief: This function prints patterns by sizes
* @param: sizePats
* @retval: Nope
*/
void printPattbySize(std::map<int, int>& sizePats);


/**
* @brief: This function prints the high utility patterns
* @param: HUPk
* @retval: Nope
*/
void printPatterns(std::unordered_map<std::string, float>& HUPk);


#endif
