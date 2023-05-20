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
void printSN()
{
	//std::cout << "The star neighbor: " << endl;
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
* @brief: This function prints a set of instances.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printInstances(std::vector<ObjWithoutCoord > SN)
{
	for (auto const& inst : SN)
	{
		std::cout << inst.feature << "." << inst.instance << ", ";
	}
	std::cout << endl;
}


/**
* @brief: This function prints sorted tuf.
* @param: sortF
* @retval: Nope
*/
void printSortTUF()
{
	std::cout << "Sorted utility of features: \n";
	for (auto const& item : sortF)
	{
		std::cout << item <<", ";
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
* @brief: This function prints FONT
* @param: FONT: tree<ObjWithoutCoord> FONT;
* @retval: Nope
*/
void printFONT()
{
	std::cout << "The FONT tree: " << endl;
	tree<ObjWithoutCoord>::pre_order_iterator it = FONT.begin();
	while (it != FONT.end())
	{
		for (int i = 0; i < FONT.depth(it); ++i)
			cout << "   ";
		std::cout << (*it).feature << "." << (*it).instance << endl;
		++it;
	}
	std::cout << endl << endl;
}


/**
* @brief: This function prints candidate tree
* @param:tree < std::vector<char>>  CandTree
* @retval: Nope
*/
void printCandTree()
{
	std::cout << "The candidate tree: " << endl;
	tree < std::vector<char>> ::pre_order_iterator it = CandTree.begin();
	while (it != CandTree.end())
	{
		for (int i = 0; i < CandTree.depth(it); ++i)
			cout << "   ";
		for (auto f : *it)
		{
			std::cout << f << ", ";
		}
		std::cout << endl;
		++it;
	}
	std::cout << endl << endl;
}


/**
* @brief: This function prints feature-object neighbor sets
* @param: std::map<ObjWithoutCoord, std::map<char, std::set<ObjWithoutCoord>>> FONS
* @retval: Nope
*/
void printFONS()
{
	for (auto const& item : FONS)
	{
		std::cout << item.first.feature << "." << item.first.instance << ": \n";
		for (auto const& inner : item.second)
		{
			std::cout << inner.first << ": ";
			for (auto const& nei : inner.second)
			{
				std::cout << nei.feature << "." << nei.instance << ", ";
			}
			std::cout << endl;
		}
		std::cout << endl;
	}
	std::cout << endl;
}




/**
* @brief: This function prints candidates
* @param: c: a candiate
* @retval: Nope
*/
void printCandidate(std::vector<char> c)
{
	//std::cout << "Candidate is: ";
	for (auto const& f : c)
	{
		std::cout << f << ", ";
	}
	std::cout << endl;
}


/**
* @brief: This function print all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern()
{
	std::cout << "The HUCP are: " << endl;
	for (auto const& pattern : PkAll)
	{
		std::cout << pattern.first << " : " << pattern.second << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints all participating instances of c
* @param:  CSHBS
* @retval: Nope
*/
void printCSHBS(std::unordered_map<char, std::set<ObjWithoutCoord>> CSHBS)
{
	std::cout << "The participating instances: \n";
	for (auto row : CSHBS)
	{
		std::cout << row.first << " : ";
		for (auto pin : row.second)
		{
			std::cout << pin.feature << "." << pin.instance << ", ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}



/**
* @brief: This function prints all row instances of c
* @param:  rowInstSet
* @retval: Nope
*/
void printRowInstSet(std::vector<std::vector<ObjWithoutCoord>> rowInstSet)
{
	std::cout << "Row instances are: \n";
	for (auto const row : rowInstSet)
	{
		for (auto inst : row)
		{
			std::cout << inst.feature << "." << inst.instance << ", ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}