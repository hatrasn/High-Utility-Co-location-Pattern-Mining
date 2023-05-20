#include <iostream>
#include <map>
#include <vector>
#include<string>
#include <set>

#include "printFunctions.h"
#include "Object.h"
#include"toolFunctions.h"
#include"tree.hh"


using namespace std;

/**
* @brief: This function print the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int> totalInstEachFeat)
{
	std::cout << "The instance number of each feature is: " << std::endl;
	for (auto const& inst : totalInstEachFeat)
	{
		std::cout << inst.first << ":" << inst.second << ",";
	}
	std::cout << endl;
}


/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printGrid(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid)
{
	for (auto const& cell : grid)
	{
		std::cout << "The cell is: " << "(" << cell.first.first << cell.first.second << "): ";
		std::cout << "Instances: ";
		for (auto const& obj : cell.second)
		{
			std::cout << obj.feature << "." << obj.instance << ", ";
		}
		std::cout << endl;
	}

}



/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printSN(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN)
{
	std::cout << "The star neighbor: " << endl;
	for (auto const& inst : SN)
	{
		std::cout << inst.first.feature << "." << inst.first.instance << " : ";
		for (auto const& nei : inst.second)
		{
			std::cout << nei.feature << "." << nei.instance << ", ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}

/**
* @brief: This function print the group star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> SN)
{
	std::cout << "The star neighborhoods are: " << endl;
	for (auto const& feat : SN)
	{
		std::cout << "Feature: " << feat.first << endl;
		for (auto const& inst : feat.second)
		{
			std::cout << "Instance: " << inst.first.feature << "." << inst.first.instance << ": ";
			for (auto const& neigh : inst.second)
			{
				std::cout << neigh.feature << "." << neigh.instance << ", ";
			}
			std::cout << endl;
		}
		std::cout << endl;
	}
}


/**
* @brief: This function print the row instance tree of an instance. This is a pre_order_iterator tree
* @param: tree<string>: a row instance tree.
* @retval: Nope
*/
void printITree(tree<ObjWithoutCoord> ITree)
{
	//std::cout << "The instance tree: " << endl;
	tree<ObjWithoutCoord>::pre_order_iterator it = ITree.begin();
	while (it != ITree.end())
	{
		for (int i = 0; i < ITree.depth(it); ++i)
			cout << "   ";
		std::cout << (*it).feature<<"."<<(*it).instance << endl;
		++it;
	}
	std::cout << endl << endl;
}


/**
* @brief: This function print all obtained cliques
* @param: Cls: a vector of all cliques
* @retval: Nope
*/
void printCliques(std::vector<std::set<ObjWithoutCoord>> Cls)
{
	std::cout << "All cliques: " << endl;
	for (auto const& cl : Cls)
	{
		for (auto const& inst : cl)
		{
			std::cout << inst.feature << "." << inst.instance << ", ";
		}
		std::cout << endl;
	}

}


/**
* @brief: This function print all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern(std::map<std::string, float> Pk)
{
	std::cout << "The prevalent patterns are: " << endl;
	for (auto const& pattern : Pk)
	{
		std::cout << pattern.first << " : " << pattern.second << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints the hash map.
* @param: chash: a hash map of all cliques
* @retval: Nope
*/
void printChash()
{
	for (auto const& item : chash)
	{
		std::cout << item.first << endl;
		for (auto const& elemt : item.second)
		{
			std::cout << elemt.first << " : ";
			for (auto const& inst : elemt.second)
			{
				std::cout << inst.feature << "." << inst.instance << ", ";
			}
			std::cout << endl;
		}
		std::cout << endl;
	}
	std::cout << endl;
}

