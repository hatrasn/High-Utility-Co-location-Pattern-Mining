#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>


#include "printFunctions.h"
#include "Object.h"
#include "pattern.h"


using namespace std;



/**
* @brief: This function prints the utility of each feature type.
* @param: totalUtilityEachFeat: utility of each feature
* @retval: Nope
*/
void printUtilityFeature(std::map<char, double>& totalUtilityEachFeat)
{
	std::cout << "The utility of features: \n";
	for (auto const& f : totalUtilityEachFeat)
	{
		std::cout << f.first << " : " << f.second << endl;
	}
	std::cout << endl;
}



/**
* @brief: This function prints the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int>& totalInstEachFeat)
{
	std::cout << "The instance number of each feature is: " << std::endl;
	for (auto const& inst : totalInstEachFeat)
	{
		std::cout << inst.first << ":" << inst.second << endl;
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
* @brief: This function prints the group star neighborhoods.
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
* @brief: This function prints all candidate patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printCandidates(std::unordered_map<std::string, int> Ck)
{
	std::cout << "The candidates are: " << endl;
	for (auto const& pattern : Ck)
	{
		std::cout << pattern.first << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints candidate patterns and their table instances
* @param: SIk: a hash map of candidate pateterns and their table instasnces
* @retval: Nope
*/
void printStarInstance(std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> SIk) 
{
	std::cout << "The star instances of candidate patterns: " << endl;
	for (auto const& si : SIk)
	{
		std::cout << "Pattern is: ";
		for (auto const& patt : si.first)
		{
			std::cout << patt;
		}
		std::cout << endl << "Its table instance is: " << endl;
		for (auto const& tbist : si.second)
		{
			for (auto const& inst : tbist)
			{
				std::cout << inst.feature<<"."<< inst.instance << ",";
			}
			std::cout << endl;
		}
	}
}


/**
* @brief: This function prints all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern(std::vector<Pattern> Pk)
{
	std::cout << endl << "The prevalent patterns are (sorted by PI): " << endl;
	std::cout << "Pattern		UPI		Utility" << endl;
	for (auto const& pattern : Pk)
	{		
		for (auto const& feature : pattern.c)
		{
			std::cout << feature;
		}		
		std::cout <<"   " << pattern.upr << "   " << pattern.Qc << endl;
	}
}



/**
* @brief: This function prints Qc of topk patterns
* @param: topk: number of top-k
*			Qc: the quality of each top-k pattern set
* @retval: Nope
*/
void printQc(std::vector<int> topk, std::vector<double> Qc)
{
	for (size_t t = 0; t < topk.size(); ++t)
	{
		std::cout << "top: " << topk[t] << " : " << Qc[t] << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints average Qc of topk patterns by sizes
* @param: QcbySize
* @retval: Nope
*/
void printAvgQCbySize(std::map<int, double> QcbySize)
{
	for (auto const& item : QcbySize)
	{
		std::cout << item.first << " : " << item.second << endl;
	}
	std::cout << endl;
}



/**
* @brief: This function prints QC of top-k patterns by each size
* @param: QcbySize
* @retval: Nope
*/
void printQCofTopkofSize(std::map<int, double>& QcTopkbySize)
{
	std::cout << "top-k	Qc" << endl;
	for (auto const& item : QcTopkbySize)
	{
		std::cout << " " << item.first << "	" << item.second << endl;
	}
	std::cout << endl;
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









