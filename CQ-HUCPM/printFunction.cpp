#include <iostream>
#include <map>
#include <vector>
#include<string>
#include <set>
#include <unordered_map>

#include "printFunctions.h"
#include "Object.h"
#include"toolFunctions.h"
#include"tree.hh"


using namespace std;

/**
* @brief: This function prints the instance number of each feature.
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
* @brief: This function prints the star neighborhoods.
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
* @brief: This function prints the star neighborhoods.
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
* @brief: This function prints the row instance tree of an instance. This is a pre_order_iterator tree
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
		std::cout << (*it).feature<<"."<<(*it).instance<<" : "<< (*it).fl << endl;
		++it;
	}
	std::cout << endl << endl;
}


/**
* @brief: This function prints patterns by sizes
* @param: sizePats
* @retval: Nope
*/
void printPattbySize(std::map<int, int>& sizePats)
{
	std::cout << "Size	# patterns" << endl;
	for (auto const& item : sizePats)
	{
		std::cout << item.first << " : " << item.second << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints the high utility patterns
* @param: HUPk
* @retval: Nope
*/
void printPatterns(std::unordered_map<std::string, float>& HUPk)
{
	for (auto const& item : HUPk)
	{
		std::cout << item.first << " : " << item.second << endl;
	}
	std::cout << endl;
}

