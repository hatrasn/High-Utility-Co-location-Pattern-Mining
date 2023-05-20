
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
std::map<char, int> totalInstEachFeat;
std::unordered_map<std::string, std::unordered_map<char, std::set<ObjWithoutCoord>>> chash;
float PI;

clock_t timeCreateChash; // time for generating candidates
clock_t timeTimeCalPIs; // time for calculating PIS

double totalTimeCreateChash = 0.0f;
double totalTimeCalPIs = 0.0f;



/**
* @brief: This function count the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
void countNumberInstance(std::vector<std::vector<std::string>> dataList)
{
	char feature;
	int instance;

	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		feature = vec[0][0]; // get the feature	

		if (totalInstEachFeat.empty())
		{
			totalInstEachFeat.insert({ feature, 1 });
		}
		else if (totalInstEachFeat.find(feature) == totalInstEachFeat.end())
		{
			totalInstEachFeat.insert({ feature, 1 });
		}
		else
		{
			instance = totalInstEachFeat.find(feature)->second;
			instance += 1;
			totalInstEachFeat[feature] = instance;
		}
	}	
}


/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(std::vector<std::vector<std::string> > dataList, float dist_thres)
{
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid;

	int cell_x, cell_y;
	std::pair<int, int> cell_key;

	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		cell_x = ceil(std::stof(vec[2]) / dist_thres); // x coordinate
		cell_y = ceil(std::stof(vec[3]) / dist_thres); // y coordinate		
													   // make keys
		cell_key = make_pair(cell_x, cell_y);
		//cout << "One cell key: " << cell_key.first << ":" << cell_key.second << endl;

		// package the value
		std::vector<ObjWithCoord> value;
		ObjWithCoord instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[2]), std::stof(vec[3]) };
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

	return grid;
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
std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> genStarNeighborhoods(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid, float dist_thres)
{
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN;

	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	float dist;
	std::vector<ObjWithCoord> fiveCellInst; // save all instance in the five cells and two neighboring instances.
	std::vector<ObjWithCoord> neighborPair;
	std::set<ObjWithoutCoord> value; // save all neighbors of the center instance		

										  // Loop each cell in grid and computation neighbors in five cells
	for (std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator it = grid.begin(); it != grid.end(); ++it)
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

		for (auto const& cell : fiveCells)
		{
			if (grid.count(cell))
			{
				fiveCellInst.insert(fiveCellInst.end(), grid.find(cell)->second.begin(), grid.find(cell)->second.end());
			}
		}

		// Sort all instances in the five cells
		std::sort(fiveCellInst.begin(), fiveCellInst.end());

		// Iterator each instance in the currentCellInst and check its neighbors		
		for (auto const& currentInst : it->second)
		{
			// Create this key
			ObjWithoutCoord key = { currentInst.feature, currentInst.instance };
			if (SN.find(key) == SN.end()) // this instance has not already existed
			{
				// update value				
				SN.insert({ key, std::set<ObjWithoutCoord>{} });
			}
			// Check neighboring instances
			for (auto const& checkInst : fiveCellInst) // check with each instance in the five cells
			{
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{
					//std::cout << "Check is neigborhoods";

					dist = calculateDistanceTwoInstances(currentInst, checkInst);

					if (dist <= dist_thres) // the two instances have neighbor relationship
					{				
						// convert Objects to ObjWithoutCoordinate						
						ObjWithoutCoord starN = { checkInst.feature, checkInst.instance };
						if (key < starN)
						{
							if (SN.find(key) != SN.end()) // this instance has already existed
							{
								// update value				
								SN.find(key)->second.insert(starN);
							}
							else // This feature has not existed in SN, directly put into SN	
							{
								// Put into SN							
								SN.insert({ key, std::set<ObjWithoutCoord> {starN} });
							}
						}
						else
						{
							if (SN.find(starN) != SN.end()) // this instance has already existed
							{
								// update value				
								SN.find(starN)->second.insert(key);
							}
							else // This feature has not existed in SN, directly put into SN	
							{
								// Put into SN							
								SN.insert({ starN, std::set<ObjWithoutCoord> {key} });
							}
						}
					}
				}
			}
		}
		// clear all element for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}

	return SN;

}


/**
* @brief: This function generates all cliques.
* @param: SN: star neighbors         
* @retval: CLs: a vector of all cliques
*/
void generateCliques(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN)
{
	// 1. Create a vector to store all cliques	
	std::set<ObjWithoutCoord> onecl;
	// 2. Create an root
	tree<ObjWithoutCoord> ITree;
	tree<ObjWithoutCoord>::iterator root, headNode, nextsib, ancestor, delNode, addNode;
	// Initial tree
	root = ITree.begin();

	// 5.2 Loop each instance to get Hs-cliques
	for (auto const& s : SN)
	{
		// Create a queue to store nodes
		queue<tree<ObjWithoutCoord>::iterator> Q;
		// Add headNode
		headNode = ITree.insert(root, s.first);
		// Put headNode to the end of queue
		Q.push(headNode);
		// Loop queue
		while (!Q.empty())
		{
			// Get the first element in queue
			tree<ObjWithoutCoord>::iterator currNode;
			currNode = Q.front(); // get the first element	
			Q.pop(); // delete the first element
			// Get all children of this node (Lemma 3)
			std::vector<ObjWithoutCoord> childrenNodes;
			if (ITree.depth(currNode) == 0) // this currNode is a head-cliques, only get the star neighbor of it as childrenNodes
			{
				childrenNodes.insert(childrenNodes.end(),
					SN.find(*currNode)->second.begin(), SN.find(*currNode)->second.end());
			}
			else // this node is a normal node
			{
				// get intersection of sibling and star neighbor
				// 1. get all sibling of currNode
				std::vector<ObjWithoutCoord> sibVec;
				nextsib = currNode;
				while (nextsib.node->next_sibling != nullptr)
				{
					sibVec.push_back(nextsib.node->next_sibling->data);
					nextsib = nextsib.node->next_sibling;
				}
				// 2. get intersection of sibling and starnei
				if (sibVec.size())
				{
					std::set_intersection(sibVec.begin(), sibVec.end(),
						SN.find(*currNode)->second.begin(), SN.find(*currNode)->second.end(),
						back_inserter(childrenNodes));
				}
				else
				{
					childrenNodes.clear();
				}
			}
			// Check these children
			if (childrenNodes.empty()) // The currNode is a leaf node, it has no children
			{
				// build cliques by get all ancestors of currNode and add to the results
				onecl.insert(*currNode);
				ancestor = currNode.node->parent;
				while (ancestor != nullptr)
				{
					onecl.insert(*ancestor);
					ancestor = ancestor.node->parent; // get the next pattern
				}
				if (onecl.size() > 1)
				{
					// Create chash
					timeCreateChash = clock();
					createChash(onecl);
					totalTimeCreateChash = totalTimeCreateChash + double(clock() - timeCreateChash);
				}
				onecl.clear();

				// Remove the ancestor if it has only one child
				ancestor = currNode.node->parent;
				// Delete currNode
				ITree.erase(currNode);
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
					addNode = ITree.append_child(currNode, onechild);
					Q.push(addNode);
				}
				childrenNodes.clear();
			}
		}
		// Clear tree for the next instances
		ITree.clear();
	}	
}



/**
* @brief: This function generates all cliques.
* @param: SN: star neighbors
* @retval: CLs: a vector of all cliques
*/
void generateCliques2(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN)
{
	// 1. Create a vector to store all cliques	
	std::set<ObjWithoutCoord> onecl;
	// 2. Create an root
	tree<ObjWithoutCoord> ITree;
	tree<ObjWithoutCoord>::iterator root, headNode, nextsib, ancestor, delNode, addNode;
	// Initial tree
	root = ITree.begin();

	// 5.2 Loop each instance to get Hs-cliques
	for (auto const& s : SN)
	{
		if (s.second.size()) // Only process the instance which has neighbor instances
		{
			//std::cout << "current node: " << s.first.feature << "." << s.first.instance << endl;
			// 5.2.1 Create a tree
			// Add headNode
			headNode = ITree.insert(root, s.first);
			// Add all its instances
			for (auto const& nei : s.second)
			{
				ITree.append_child(headNode, nei);
			}
			//std::cout << "Inital tree \n";
			//printITree(ITree);

			tree<ObjWithoutCoord>::iterator postit = ITree.begin_breadth_first();
			++postit;

			// 5.2.3 Build tree
			while (ITree.size())
			{
				// Get all children of this node (Lemma 3)
				// 1. get all sibling of currNode
				std::vector<ObjWithoutCoord> sibVec;
				std::vector<ObjWithoutCoord> childrenNodes;

				nextsib = postit;
				while (nextsib.node->next_sibling != nullptr)
				{
					sibVec.push_back(nextsib.node->next_sibling->data);
					nextsib = nextsib.node->next_sibling;
				}
				// 2. get intersection of sibling and starnei
				if (sibVec.size())
				{
					std::set_intersection(sibVec.begin(), sibVec.end(),
						SN.find(*postit)->second.begin(), SN.find(*postit)->second.end(),
						back_inserter(childrenNodes));
				}

				// Check these children
				if (childrenNodes.empty()) // The currNode is a leaf node, it has no children
				{
					// build cliques by get all ancestors of currNode and add to the results
					onecl.insert(*postit);
					ancestor = postit.node->parent;
					while (ancestor != nullptr)
					{
						onecl.insert(*ancestor);
						ancestor = ancestor.node->parent; // get the next pattern
					}
					if (onecl.size() > 1)
					{
						// Create chash
						timeCreateChash = clock();
						createChash(onecl);
						totalTimeCreateChash = totalTimeCreateChash + double(clock() - timeCreateChash);
					}
					onecl.clear();

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
						ITree.append_child(postit, onechild);
					}
					childrenNodes.clear();

					//std::cout << "Add tree \n";
					//printITree(ITree);

					++postit;
				}
			}
			// Clear tree for the next instances
			ITree.clear();
		}
	}
}



/**
* @brief: This function creates the chash to store all cliques.
* @param: Cls: all cliques
* @retval: chash: a hash map
*/
void createChash(std::set<ObjWithoutCoord> onecl)
{	
	std::string newKey;
	// Get all feature of instances in cl to build key
	for (auto const& inst : onecl)
	{
		newKey += inst.feature;
	}
	// Check this key in chash
	std::unordered_map<std::string, std::unordered_map<char, std::set<ObjWithoutCoord>>>::iterator itchash = chash.find(newKey);
	if (itchash != chash.end()) // This key is in the chash
	{
		// update
		for (auto const& elemt : onecl)
		{
			itchash->second.find(elemt.feature)->second.insert(elemt);
		}
	}
	else // this key is not in the chash
	{
		// create a new item
		std::unordered_map<char, std::set<ObjWithoutCoord>> newItem;
		for (auto const& item : onecl)
		{
			newItem.insert({ item.feature, std::set<ObjWithoutCoord>{item} });
		}
		chash.insert({newKey, newItem});
	}		
}


/**
* @brief: This function creates the chash to store all cliques.
* @param: Cls: all cliques
* @retval: chash: a hash map
*/
void createChash2(std::vector<std::set<ObjWithoutCoord>> Cls)
{
	std::unordered_map<std::string, std::unordered_map<char, std::set<ObjWithoutCoord>>>::iterator itchash;
	for (auto const& cl : Cls)
	{
		std::string newKey;
		// Get all feature of instances in cl to build key
		for (auto const& inst : cl)
		{
			newKey += inst.feature;
		}

		// Check this key in chash
		itchash = chash.find(newKey);
		if (itchash == chash.end()) // this key is not in chash, create a new key
		{
			std::unordered_map<char, std::set<ObjWithoutCoord>> newValue;
			chash.insert({newKey, newValue});
		}

		// Add value into key
		for (auto const& elemt: cl)
		{
			chash.find(newKey)->second.find(elemt.feature)->second.insert(elemt);
		}		
	}
}


/**
* @brief: This function calculates PIs.
* @param: currCandidate: one candidate, chash: a hash map
* @retval: PI: the pi of the candidate
*/
void calculatePI(std::string currCandidate)
{
	// Get all supersets of the current candidate
	std::vector<std::string> supersets;
	std::unordered_map<std::string, std::unordered_map<char, std::set<ObjWithoutCoord>>>::iterator itchash = chash.begin();
	while (itchash != chash.end())
	{
		// check currCandidate is not or subset of keys in chash
		if (std::includes(itchash->first.begin(), itchash->first.end(),
			currCandidate.begin(), currCandidate.end()))
		{
			supersets.push_back(itchash->first);
		}
		++itchash;
	}
	// Loop each cp in supersets
	std::vector<std::unordered_set<ObjWithoutCoord, myHashFunc>> insts(currCandidate.size()); // save instances
	for (auto const& cp : supersets)
	{
		for (size_t t = 0; t < currCandidate.size(); ++t)
		{
			insts[t].insert(chash[cp].find(currCandidate[t])->second.begin(),
				chash[cp].find(currCandidate[t])->second.end());
		}
	}
	
	// Calculate PI
	std::vector<float> PRs;
	for (size_t t = 0; t< currCandidate.size(); ++t)
	{
		PRs.push_back((float)insts[t].size() / (float)totalInstEachFeat.find(currCandidate[t])->second);		
	}

	// find the minimum elemnet
	PI = *std::min_element(std::begin(PRs), std::end(PRs));	
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


/**
* @brief: This function calculates PIs and filter prevalent patterns.
* @param: chash: a hash map of all cliques
* @retval: PkAll: all prevalent patterns
*/
std::map<std::string, float> CalPIandFilterPatterns(float prev_thres)
{
	std::map<std::string, float> PkAll;
	//1. Get all candidate
	std::vector<std::string> candidates, checked;
	for (auto const& item : chash)
	{
		candidates.push_back(item.first);
	}
	// Sort by size of element
	std::sort(candidates.begin(), candidates.end(), [](const std::string& pat1, const std::string& pat2)
		{
			return pat1.size() > pat2.size();
		});	
	while (!candidates.empty())
	{
		// get one pattern
		std::string currCandidate = candidates.front();		
		candidates.erase(candidates.begin());
		// Marked as checked patterns
		checked.push_back(currCandidate);

		// cal PI
		calculatePI(currCandidate);		
		// Check prevelence
		if (PI >= prev_thres)
		{
			// this currCandidate and its subsets are prevalent patterns
			PkAll.insert({ currCandidate, PI });
			for (size_t t = 2; t < currCandidate.size(); ++t)
			{
				// generate subsets
				std::vector<std::string> subsets = combo(currCandidate, t);
				// calcualte PIs and add to the result
				for (auto const& onesub : subsets)
				{
					
					calculatePI(onesub);
					PkAll.insert({onesub, PI});
					
					checked.push_back(onesub);
					// Remove from candidate
					for (std::vector<std::string>::iterator iter = candidates.begin(); iter != candidates.end(); ++iter)
					{
						if (*iter == onesub)
						{
							candidates.erase(iter);
							break;
						}
					}
				}
			}
		}
		else
		{
			// Generate the direcsubset of the currcandidate
			std::vector<std::string> directsubsets = combo(currCandidate, currCandidate.size()-1);
			// Add the direcsubset as new candidate to the candidate
			for (auto const& newCand : directsubsets)
			{
				if (std::find(checked.begin(), checked.end(), newCand) == checked.end())
				{
					candidates.push_back(newCand);
				}
			}
			//candidates.insert(candidates.end(), directsubsets.begin(), directsubsets.end());
			// Sort
			std::sort(candidates.begin(), candidates.end(), [](const std::string& pat1, const std::string& pat2)
				{
					return pat1.size() > pat2.size();
				});
		}
	}

	return PkAll;
}

/**
* @brief: This function sorts elements in map by the size of keys.
* @param: candidate: a hash map of all cliques
* @retval: PkAll: all prevalent patterns
*/
struct cmpByStringLength {
	bool operator()(const std::string& a, const std::string& b) const {
		return a.size() >= b.size();
	}
};


// Comparator function to sort pairs 
// according to second value 
bool cmp(pair<string, int>& a,
	pair<string, int>& b)
{
	return a.first.size() >= b.first.size();
}


// Function to sort the map according 
// to value in a (key-value) pairs 
std::unordered_map<std::string, int> sortMap(std::unordered_map<std::string, int>& M)
{
	// Declare vector of pairs 
	std::vector<std::pair<std::string, int> > A;

	// Copy key-value pair from Map 
	// to vector of pairs 
	for (auto& it : M) {
		A.push_back(it);
	}

	// Sort using comparator function 
	sort(A.begin(), A.end(), cmp);

	std::unordered_map<std::string, int> sortedMap;	
	// Print the sorted value 
	for (auto& it : A) 
	{
		sortedMap.insert(it);
	}

	return sortedMap;
}


/**
* @brief: This function compares the size of two patterns
* @param: cand1 and cand2: two candidates
* @retval: bool
*/
struct compareBySize
{
	bool operator()(std::string const& cand1, std::string const& cand2)
	{
		return cand1.size() < cand2.size();
	}
};



/**
* @brief: This function calculates PIs and filter prevalent patterns.
* @param: chash: a hash map of all cliques
* @retval: PkAll: all prevalent patterns
*/
std::map<std::string, float> CalPIandFilterPatterns2(float prev_thres)
{
	std::map<std::string, float> PkAll;
	//1. Get all candidate and they are sorted by sizes of keys	
	std::vector<std::string> candidates, checked;

	// Sort the candiates by sizes
	for (auto const& item : chash)
	{
		candidates.push_back(item.first);
	}
	// Sort candidates by sizes of keys
	std::sort(candidates.begin(), candidates.end(), compareBySize());

	// Loop each candidate
	while (!candidates.empty())
	{
		// 4.1 Get the first element
		std::string currCandidate = candidates.back();
		
		// ------------------Do not use this to keep slower------------------ //
		checked.push_back(currCandidate);

		// 4.2 Delete the element from candidates
		candidates.pop_back(); 

		// 4.3 Cal PI
		timeTimeCalPIs = clock();		
		calculatePI(currCandidate);
		totalTimeCalPIs += double(clock() - timeTimeCalPIs);
		// Check prevelence
		if (PI >= prev_thres)
		{
			// this currCandidate and its subsets are prevalent patterns
			PkAll.insert({ currCandidate, PI });
			for (size_t t = 2; t < currCandidate.size(); ++t)
			{
				// generate subsets
				std::vector<std::string> subsets = combo(currCandidate, t);
				// calcualte PIs and add to the result
				for (auto const& onesub : subsets)
				{
					// Calcualte the PI of the candiate
					timeTimeCalPIs = clock();
					calculatePI(onesub);
					totalTimeCalPIs += double(clock() - timeTimeCalPIs);
					PkAll.insert({ onesub, PI });
					checked.push_back({ onesub, 1 });					
					// Remove from candidate
					candidates.erase(std::remove(candidates.begin(), candidates.end(), onesub), candidates.end());
				}
			}
		}
		else
		{
			// Generate the direcsubset of the currcandidate
			if (currCandidate.size() > 2)
			{
				std::vector<std::string> directsubsets = combo(currCandidate, currCandidate.size() - 1);
				// Add the direcsubset as new candidate to the candidate
				for (auto const& newCand : directsubsets)
				{
					// check if this pattern already calculated, do not put it in the candidate
					if (std::find(candidates.begin(), candidates.end(), newCand) == candidates.end()
						&& std::find(checked.begin(), checked.end(), newCand) == checked.end())
					{						
						candidates.push_back(newCand);
					}
				}
				// sort candidates
				std::sort(candidates.begin(), candidates.end(), compareBySize());
			}			
		}
	}
	return PkAll;
}