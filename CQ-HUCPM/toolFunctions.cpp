
#include <string>
#include <ctime>
#include <algorithm>
#include <math.h>
#include <cstdlib>
#include <numeric>
#include <iostream>
#include <map>
#include <vector>
#include<unordered_set>

#include "toolFunctions.h"
#include "Object.h"
#include "printFunctions.h"

#include"tree.hh"

#include <boost/algorithm/string.hpp>


using namespace std;

// Definition global variables
std::unordered_map<std::string, std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>>> coHM;
clock_t timeBuildHash;
double totalTimeBuildHash = 0.0f;


// ------ For compare of set_intersection ------//
struct compa_set_intersec
{
	template<typename T>
	bool operator()(T const& obT, ObjWithoutCoord const& obU) const {
		if ((*obT).feature < obU.feature) return true;
		else if ((*obT).feature == obU.feature && (*obT).instance < obU.instance) return true;
		else return false;
	}
	template<typename T>
	bool operator()(ObjWithoutCoord const& obU, T const& obT) const {
		if (obU.feature < (*obT).feature) return true;
		else if (obU.feature == (*obT).feature && obU.instance < (*obT).instance) return true;
		else return false;
	}
};


/**
* @brief: This function gets the maximal and minimal values of the utility of each feature type
* @param: dataList: an input dataset;
* @retval: maxminU: <A, [minU, maxU]>, <B, [minU, maxU]>, ...
*/
void getMinMaxUtility(
	std::vector<std::vector<std::string>>& dataList,
	std::unordered_map<char, std::vector<float>>& maxminU)
{
	//1. Loop each instance and get its utility
	std::unordered_map<char, std::vector<float>> allU;
	std::unordered_map<char, std::vector<float>>::iterator itallU;
	for (auto const& row : dataList)
	{
		itallU = allU.find(row[0][0]);
		if (itallU != allU.end())
		{
			itallU->second.push_back(std::stof(row[4]));
		}
		else
		{
			allU.insert({ row[0][0], std::vector<float>{std::stof(row[4])} });
		}
	}
	//2. Get min and max
	for (auto const& item : allU)
	{
		std::vector<float> maxmin(2);
		maxmin[0] = *std::min_element(item.second.begin(), item.second.end());
		maxmin[1] = *std::max_element(item.second.begin(), item.second.end());

		maxminU.insert({ item.first, maxmin });
	}
}



/**
* @brief This function makes a grid on the input dataset
* @param dataList: an input dataset
*		dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it
*/
void makeGrid(
	std::vector<std::vector<std::string>>& dataList, 
	float dist_thres,
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> & grid)
{
	int cell_x, cell_y;
	std::pair<int, int> cell_key;
	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		cell_x = ceil(std::stof(vec[2]) / dist_thres); // x coordinate
		cell_y = ceil(std::stof(vec[3]) / dist_thres); // y coordinate		
													   // make keys
		cell_key = make_pair(cell_x, cell_y);		
		std::vector<ObjWithCoord> value;
		ObjWithCoord instance = {
			vec[0][0], // feature
			std::stoi(vec[1]), // instance
			std::stof(vec[2]), // x
			std::stof(vec[3]), // y
			std::stof(vec[4]) // utility 
		};
		value.push_back(instance);
		// check if this key is exsiting?
		if (grid.empty()) //The grid is empty
		{
			grid[cell_key] = value;
		}
		else if (grid.find(cell_key) == grid.end()) // if the key is not exist
		{
			grid[cell_key] = value;
		}
		else // the key has already existed
		{
			std::vector<ObjWithCoord> old_value = grid.find(cell_key)->second;
			old_value.push_back(instance);
			grid[cell_key] = old_value;
		}
	}	
}


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst)
{
	return  sqrt((currentInst.x - checkInst.x)*(currentInst.x - checkInst.x)
		+ (currentInst.y - checkInst.y)*(currentInst.y - checkInst.y));
}


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>>& grid,
	float dist_thres,
	std::unordered_map<char, std::vector<float>>& minmaxU,
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc>& SN,
	std::map<ObjWithoutCoord, int>& V,
	std::unordered_map<char, float>& totalUtilityFeats)
{	
	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	float dist;
	// save all instance in the five cells and two neighboring instances.
	std::vector<ObjWithCoord> fiveCellInst;
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc>::iterator itSN;
	std::unordered_map<char, std::vector<float>>::iterator itminmaxU;
	std::unordered_map<char, float>::iterator ittotalU;
	for (std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator it = grid.begin();
		it != grid.end();
		++it)
	{
		// Create the five cells
		i = it->first.first;
		j = it->first.second;

		fiveCells.push_back(it->first);
		fiveCells.push_back(make_pair(i, j + 1));
		fiveCells.push_back(make_pair(i + 1, j + 1));
		fiveCells.push_back(make_pair(i + 1, j));
		fiveCells.push_back(make_pair(i + 1, j - 1));
		// Get all instances in the current cell and five cells
		std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator itgrid;
		for (auto const& cell : fiveCells)
		{
			itgrid = grid.find(cell);
			if (itgrid != grid.end())
			{
				fiveCellInst.insert(fiveCellInst.end(),
					itgrid->second.begin(), itgrid->second.end());
			}
		}
		// Sort all instances in the five cells
		std::sort(fiveCellInst.begin(), fiveCellInst.end());
		// Iterator each instance in the currentCellInst and check its neighbors		
		for (auto const& currentInst : it->second)
		{			
			// Calculate nomorlized utiltiy of this instances
			itminmaxU = minmaxU.find(currentInst.feature);
			float uc = (currentInst.u - itminmaxU->second[0]) / (itminmaxU->second[1] - itminmaxU->second[0]);
			// Save this uc to corresponding feature type
			ittotalU = totalUtilityFeats.find(currentInst.feature);
			if (ittotalU != totalUtilityFeats.end())
			{
				ittotalU->second += uc;
			}
			else
			{
				totalUtilityFeats.insert({ currentInst.feature, uc });
			}
			// Save this instance to SN
			ObjWithoutCoord key = {
				currentInst.feature, // feature
				currentInst.instance, // instance
				uc,  // utility
				0 }; // flag
			// Save to V
			V.insert({ key, 1 });
			// Create SN for this instance
			itSN = SN.find(key);
			if (itSN == SN.end())
			{
				SN.insert({ key, std::set<ObjWithoutCoord> {} });
			}
			// Check neighboring with other instances in the current cell
			for (auto const& checkInst : fiveCellInst)
			{
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{
					// Calucalate distance between the two instances
					dist = calculateDistanceTwoInstances(currentInst, checkInst);
					if (dist <= dist_thres) // the two instances have neighbor relationship
					{						
						// Calculate nomorlized utiltiy of this instances
						itminmaxU = minmaxU.find(checkInst.feature);
						uc = (checkInst.u - itminmaxU->second[0]) / (itminmaxU->second[1] - itminmaxU->second[0]);
						// Save this uc to corresponding feature type
						ObjWithoutCoord starN = {
							checkInst.feature,
							checkInst.instance,
							uc,
							0 };
						if (key < starN)
						{
							itSN = SN.find(key);
							if (itSN != SN.end())
							{
								itSN->second.insert(starN);
							}
							else
							{
								SN.insert({ key, std::set<ObjWithoutCoord> {starN} });
							}
						}
						else
						{
							itSN = SN.find(starN);
							if (itSN != SN.end())
							{
								itSN->second.insert(key);
							}
							else
							{
								SN.insert({ starN, std::set<ObjWithoutCoord> {key} });
							}
						}
					}
				}
			}
		}
		// clear all elements for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}
}


/**
* @brief: This function constructs co-location Hash map
* @param: onecl: one row instance
* @retval: coloHM: a co-location hashmap storing all row instances
*/
void createCoLoHM(std::vector<ObjWithoutCoord> const& onecl)
{
	std::string pattern;
	// Get all feature of instances in cl to build key
	for (auto const& inst : onecl)
	{
		pattern += inst.feature;
	}
	// Check this key in coloHM
	std::unordered_map<std::string, std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>>>::iterator itcoloHM = coHM.find(pattern);
	if (itcoloHM != coHM.end()) // This key is in the coloHM
	{
		// update
		for (size_t t = 0; t < onecl.size(); ++t)
		{
			itcoloHM->second[t].insert(onecl[t]);
		}
	}
	else // this key is not in the chash
	{
		std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>> newItem;
		for (auto const& item : onecl)
		{
			newItem.push_back(std::unordered_set<ObjWithoutCoord, myHashFunc>{item});
		}
		coHM.insert({ pattern, newItem });
	}
}



/**
* @brief: This function generates all cliques
* @param: SN: star neighbors
* @retval: CLs: a vector of all cliques
*/
void generateCliques(
	std::map<ObjWithoutCoord, int>& V,
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc>& SN)
{
	// 1. Create a vector to store all cliques	
	std::vector<ObjWithoutCoord> onecl; // ----- use vector and sort ------ //
	// 2. Create an root
	tree<ObjWithoutCoord> ITree;
	tree<ObjWithoutCoord>::iterator root, headNode, nextsib, ancestor, delNode, addNode;
	// 3. Initial tree
	root = ITree.begin();
	std::unordered_map<ObjWithoutCoord, std::set<ObjWithoutCoord>, myHashFunc>::iterator itSN;
	// 4. Loop each instance to get Hs-cliques
	while (V.size())
	{
		// 4.1 Get one instances	
		std::vector<ObjWithoutCoord> vertices{ V.begin()->first };
		// 4.2 Delete the instance from V
		V.erase(V.begin());
		// 4.3 And add its star neighboring intances
		itSN = SN.find(vertices[0]);
		if (itSN != SN.end())
		{
			vertices.insert(vertices.end(), itSN->second.begin(), itSN->second.end());
		}
		// 4.5 Check each star neigboring instance to delete from V
		std::sort(vertices.begin() + 1, vertices.end());
		int n = vertices.size();
		if (n > 2)
		{
			// Get the vertices
			for (size_t t = 1; t < n; ++t)
			{
				itSN = SN.find(vertices[t]);
				// If one instance, its star neighboring instances are all in vertices,  it can be deleted from V
				if (itSN->second.size())
				{
					if (std::includes(vertices.begin(), vertices.end(),
						itSN->second.begin(), itSN->second.end()))
					{
						V.erase(vertices[t]);
					}
				}
				else // This instance has no neighbor instances, it can be deleted directly
				{
					V.erase(vertices[t]);
				}
			}			
			headNode = ITree.insert(root, vertices[0]);			
			// Add neighboring instances of the center instance as children
			for (size_t i = 1; i < n; ++i)
			{
				ITree.append_child(headNode, vertices[i]);
			}			
			tree<ObjWithoutCoord>::iterator postit = ITree.begin();
			++postit; // skip the headNode
			// Loop queue
			while (ITree.size())
			{				
				std::vector<tree<ObjWithoutCoord>::iterator> childrenNodes;
				// get all sibling of currNode
				std::vector<tree<ObjWithoutCoord>::iterator> sibVec;
				nextsib = postit;
				while (nextsib.node->next_sibling != nullptr)
				{
					sibVec.push_back(nextsib.node->next_sibling);
					nextsib = nextsib.node->next_sibling;
				}
				// get intersection of sibling and starnei
				if (sibVec.size())
				{
					std::set_intersection(sibVec.begin(), sibVec.end(),
						SN.find(*postit)->second.begin(), SN.find(*postit)->second.end(),
						std::back_inserter(childrenNodes),
						compa_set_intersec{}
					);
				}
				// check these children
				if (childrenNodes.empty()) // The currNode is a leaf node, it has no children
				{
					if (postit->fl == 0) // only combine a clique it the current node .fl = 0
					{
						// build cliques by get all ancestors of currNode and add to the results
						onecl.push_back(*postit);
						ancestor = postit.node->parent;
						while (ancestor != nullptr)
						{
							onecl.push_back(*ancestor);
							ancestor = ancestor.node->parent; // get the next pattern
						}
						if (onecl.size() > 1)
						{
							// Sort cliques
							std::sort(onecl.begin(), onecl.end());
							// Create chash			
							timeBuildHash = clock();
							createCoLoHM(onecl);
							totalTimeBuildHash = totalTimeBuildHash + double(clock() - timeBuildHash);
						}
						onecl.clear();
					}
					// Remove the ancestor if it has only one child
					ancestor = postit.node->parent;
					// Delete currNode
					postit = ITree.erase(postit);
					while (ancestor != nullptr)
					{
						if (ancestor.number_of_children() == 0) // this ancestor has only one child is currNode, it is deleted
						{
							delNode = ancestor.node->parent;
							ITree.erase(ancestor);
							ancestor = delNode;
						}
						else
						{
							break;
						}
					}	
				}
				else
				{
					// The currNode has children, add the children into the tree and put them into queue
					for (auto const& onechild : childrenNodes)
					{
						addNode = ITree.append_child(postit, *onechild);
						addNode->fl = 0;					
						onechild->fl = 1;
					}
					childrenNodes.clear();
					++postit;
				}
			}			
			// Clear tree for the next instances
			ITree.clear();
		}
	}
}


/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombination(const T& c, int combo)
{
	std::string result;
	int n = c.size();
	for (int i = 0; i < n; ++i) {
		if ((combo >> i) & 1)
			result.push_back(c[i]);
	}
	return result;
}


template<typename T>
std::vector<std::string> combo(const T& c, int k)
{
	std::vector<std::string> combination;

	int n = c.size();
	int combo = (1 << k) - 1;       // k bit sets

	while (combo < 1 << n)
	{
		combination.insert(combination.end(), getCombination(c, combo));
		int x = combo & -combo;
		int y = combo + x;
		int z = (combo & ~y);
		combo = z / x;
		combo >>= 1;
		combo |= y;
	}
	return combination;
}



/*
*@brief: This function calculates the utility of the current pattern
*@param: pattern: the current pattern which need to calculate its utility.
*       coHM: a co-location pattern hash map
*@retval: UT: the utility of the current pattern
*/
float calcUtilityWithSuperPatterns2Phases(
	std::string pattern,
	std::unordered_map<char, float>& totalUtilityFeats,
	int sizeofPatttern)
{
	std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>> UTAll(sizeofPatttern);
	bool isInHashMap;
	int indexFeat;
	// 1. Get instances of the super patterns	
	for (
		std::unordered_map<std::string, std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>>>::iterator it = coHM.begin();
		it != coHM.end();
		++it)
	{
		// Check if the key in cpHashMap is super set of the current pattern
		if (it->first.size() >= sizeofPatttern)
		{
			isInHashMap = std::includes(it->first.begin(), it->first.end(), pattern.begin(), pattern.end());
			if (isInHashMap)					
			{
				// Get Tc
				for (size_t t = 0; t < sizeofPatttern; t++)
				{
					// First: find the index of the feature in the key
					indexFeat = lower_bound(it->first.begin(), it->first.end(), pattern[t]) - it->first.begin();
					// Get the utility of the current feature in the candidate
					UTAll[t].insert(it->second[indexFeat].begin(), it->second[indexFeat].end());
				}
			}
		}
	}

	// 2. Calculate intra utility
	std::vector<float> intraU(sizeofPatttern);
	for (size_t t = 0; t < sizeofPatttern; ++t)
	{
		float utf = 0.0f;
		for (auto const& inst : UTAll[t])
		{
			utf += inst.u;
		}
		intraU[t] = utf / totalUtilityFeats.find(pattern[t])->second;
	}

	// 3. Calculate inter utility	
	std::vector<float> interU(sizeofPatttern);
	for (size_t t = 0; t < sizeofPatttern; ++t)
	{
		// 3.1 Get utility of the remain features
		float u1 = 0.0f, u2 = 0.0f;
		for (size_t j = 0; j < sizeofPatttern; ++j)
		{
			if (j != t)
			{
				for (auto const& inst : UTAll[j])
				{
					u1 += inst.u;
				}
				u2 += totalUtilityFeats.find(pattern[j])->second;
			}
		}
		// 3.2 Obtain interU of one feature
		interU[t] = u1 / u2;
	}

	// 4. Calcualte the weight w1 and w2
	std::vector<float> UPRs(sizeofPatttern);
	for (size_t t = 0; t < sizeofPatttern; ++t)
	{
		float w1 = intraU[t] / (intraU[t] + interU[t]);
		float w2 = interU[t] / (intraU[t] + interU[t]);
		UPRs[t] = w1 * intraU[t] + w2 * interU[t];
	}

	// 5. Calculate UPI
	float UPI = *std::min_element(std::begin(UPRs), std::end(UPRs));

	// Return
	return UPI;
}




/**
* @brief: This function calculates the utility of each pattern
* @param: V and SN
* @retval: coloHM: a co-location hashmap storing all row instances
*/
void calculateUtility(std::unordered_map<char, float>& totalUtilityFeats,
	std::unordered_map<std::string, float>& Pk)
{
	// 1. Get all keys from coHM
	std::unordered_map<std::string, int> keys;
	for (auto const& item : coHM)
	{
		keys.insert({item.first, 1});		
	}
	// 2. Devide the keys in colocHM into a feature of group
	
	std::unordered_map<std::string, int>::iterator itkeys = keys.begin(), itcheksubc;

	// 3. Take one key and calculate its utility
	std::unordered_map<std::string, float>::iterator itPk;
	while (itkeys != keys.end())
	{
		// 3.1 Check this is has already calcualed
		itPk = Pk.find(itkeys->first);
		int k = itkeys->first.size();

		if (itPk == Pk.end()) // This cand need to calculate its utitlity
		{
			float upi = calcUtilityWithSuperPatterns2Phases(itkeys->first, totalUtilityFeats, k);
			Pk.insert({ itkeys->first, upi });
		}
		// 3.2 Calculate its all subsets		
		if (k > 2)
		{
			for (size_t t = 2; t < k; ++t)
			{
				std::vector<std::string> subCand = combo(itkeys->first, t);
				for (auto const& subc : subCand)
				{
					itPk = Pk.find(subc);
					if (itPk == Pk.end())
					{
						// Check in keys
						itcheksubc = keys.find(subc);
						if (itcheksubc == keys.end())
						{
							float upi = calcUtilityWithSuperPatterns2Phases(subc, totalUtilityFeats, subc.size());
							Pk.insert({ subc, upi });
						}
					}
				}
			}
		}
		// Terminate
		++itkeys;
	}
}


/*
*@brief: This function filters high utility patterns
*@param: Pk: all patterns with their utility values.
*       high_thres: threshold given by users
*@retval: HUPk: the high utility patterns
*/
void filterHighUtilityPatterns(
	std::unordered_map<std::string, float>& Pk,
	std::unordered_map<std::string, float>& HUPk,
	float high_thres)
{
	for (auto const& item : Pk)
	{
		if (item.second >= high_thres)
		{
			HUPk.insert({ item.first, item.second });
		}
	}
}


