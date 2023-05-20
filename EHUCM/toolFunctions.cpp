
#include <string>
#include <algorithm>
#include <ctime>
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
std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN;
tree<std::vector<char>> CandTree;
tree<ObjWithoutCoord> FONT;
tree<ObjWithoutCoord>::iterator rootFONT;

std::vector<char> sortF;  // the set of features that are sorted by their utilities
std::map<char, double> vuf; // the utility of each feature

std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; // table instances of size k patterns
std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIkplus; // table instances of size k patterns

std::map<std::string, double> PkAll; // store all HUCP
double US; // utility of all data set
std::map<ObjWithoutCoord, std::map<char, std::set<ObjWithoutCoord>>> FONS; 

std::vector<std::vector<ObjWithoutCoord>> realTableInst;

//clock_t startFilHUCMP;
double totalGenerateParticipatingObjects = 0.0f;
double totalGenCand = 0.0f;
double totalTimeCalUPIAndFilPrevCoLoPat = 0.0f;


/**
* @brief This function makes a grid on an input dataset
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
		cell_key = make_pair(cell_x, cell_y);		
		std::vector<ObjWithCoord> value;
		ObjWithCoord instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[2]), std::stof(vec[3]), std::stof(vec[4]) };
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
	return  sqrt((currentInst.x - checkInst.x) * (currentInst.x - checkInst.x)
		+ (currentInst.y - checkInst.y) * (currentInst.y - checkInst.y));
}



/**
* @brief: This function gets the maximal and minimal values of the utility of each feature type
* @param: dataList: an input dataset;
* @retval: maxminU: <A, [minU, maxU]>, <B, [minU, maxU]>, ...
*/
void getMaxUtility(
	std::vector<std::vector<std::string>>& dataList,
	std::unordered_map<char, double>& maxU)
{
	//1. Loop each instance and get its utility
	std::unordered_map<char, std::vector<double>> allU;
	std::unordered_map<char, std::vector<double>>::iterator itallU;
	for (auto const& row : dataList)
	{
		itallU = allU.find(row[0][0]);
		if (itallU != allU.end())
		{
			itallU->second.push_back(std::stof(row[4]));
		}
		else
		{
			allU.insert({ row[0][0], std::vector<double>{std::stof(row[4])} });
		}
	}
	//2. Get min and max
	for (auto const& item : allU)
	{
		double mU = *std::max_element(item.second.begin(), item.second.end());

		maxU.insert({ item.first, mU });
	}
}



/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid, 
	float dist_thres, 
	std::unordered_map<char, double> maxU,
	std::map<char, double>& tuf)
{	
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>::iterator itSN;
	std::map<char, int> totalInsPerFeat;
	std::map<char, int>::iterator ittotalInsPerFeat;
	std::unordered_map<char, double>::iterator itmaxU;	
	std::map<char, double>::iterator ittuf;	
	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	double dist;
	// save all instance in the five cells and two neighboring instances.
	std::vector<ObjWithCoord> fiveCellInst;
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
			// calcualte normal utility of currentInst
			itmaxU = maxU.find(currentInst.feature);
			double uc = log10(currentInst.u) / log10(itmaxU->second);
			// Save this instance to SN
			ObjWithoutCoord key = {
				currentInst.feature, // feature
				currentInst.instance, // instance
				uc}; // utility 
			// Create SN for this instance
			itSN = SN.find(key);
			if (itSN == SN.end())
			{
				SN.insert({ key, std::set<ObjWithoutCoord> {} });
			}
			// Get utility
			ittuf = tuf.find(key.feature);
			if (ittuf != tuf.end())
			{
				ittuf->second += key.u;
			}
			else
			{
				tuf.insert({ key.feature, key.u });
			}
			// Get instance number
			ittotalInsPerFeat = totalInsPerFeat.find(key.feature);
			if (ittotalInsPerFeat != totalInsPerFeat.end())
			{
				ittotalInsPerFeat->second += 1;
			}
			else
			{
				totalInsPerFeat.insert({ key.feature, 1 });
			}
			// Check neighboring with other instances in the current cell
			for (auto const& checkInst : fiveCellInst)
			{
				// Calculate utility
				itmaxU = maxU.find(checkInst.feature);
				double uc2 = log10(checkInst.u) / log10(itmaxU->second);
				// Create to SN
				ObjWithoutCoord check = {
				checkInst.feature, // feature
				checkInst.instance, // instance
				uc2 }; // utility 				
				itSN = SN.find(check);
				if (itSN == SN.end())
				{
					SN.insert({ check, std::set<ObjWithoutCoord> {} });
				}
				// Check neighbor relationship
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{
					// Calucalate distance between the two instances
					dist = calculateDistanceTwoInstances(currentInst, checkInst);
					if (dist <= dist_thres) // the two instances have neighbor relationship
					{
						// Put them into SN
						itSN = SN.find(key);
						if (itSN != SN.end())
						{
							itSN->second.insert(check);
						}
						itSN = SN.find(check);
						if (itSN != SN.end())
						{
							itSN->second.insert(key);
						}
					}
				}
			}
		}
		// clear all element for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}
	// Calculate vuf: the utility of each feature that is the average utility
	for (auto const& item : tuf)
	{
		ittotalInsPerFeat = totalInsPerFeat.find(item.first);
		if (ittotalInsPerFeat != totalInsPerFeat.end())
		{
			vuf.insert({ item.first, item.second / double(ittotalInsPerFeat->second) });
		}		
	}
}


/**
* @bref This function calculates the utility of the input data set
* @param tuf
* @retval US
*/
void calculateUS(std::map<char, double> tuf) 
{
	US = 0.0;
	for (auto const& item: tuf)
	{
		US += item.second;
	}	
}


/**
* @bref This function compares maps by values
* @param 
* @retval
*/
bool cmp(pair<char, double>& a,
	pair<char, double>& b)
{
	return a.second > b.second;
}


/**
* @bref This function sorts features by their utilities.
* @param tuf: utilities of features
* @retval  std::unordered_map<char, double>& sortF
*/
void sortFeatureByUtility(std::map<char, double> tuf)
{
	// Declare vector of pairs
	std::vector<pair<char, double> > A;	
	for (auto& it : tuf) {
		A.push_back(it);
	}
	// Sort using comparator function
	sort(A.begin(), A.end(), cmp);
	// Put  the result the sorted value
	for (auto& it : A) {
		sortF.push_back(it.first);
	}
}


/**
* @bref This function finds HUCP
* @param c: a set of candidates
* @retval
*/
void searchHUCP(std::vector<std::vector<char>> c, int k, double minutil)
{
	clock_t timeStart, timeGenCand, timeGenPartInst;	
	bool flag; // flag=true -> the candidate can be extended, flag=false -> the candidate can not be extended
	while (c.size())
	{
		timeStart = clock();
		std::vector<std::vector<char>> extendedC;		
		// check utility of each candidate in size k candidate c
		for (int t = 0; t < c.size(); t++)
		{				
			flag = searchPatterns(c[t], minutil);		
			if (flag)
			{
				extendedC.push_back(c[t]);
			}			
		}		
		timeGenPartInst = clock();
		totalGenerateParticipatingObjects += double(timeGenPartInst - timeStart);

		// finish k size candidates
		c.clear();
		CIk.clear();
		CIk = CIkplus;
		CIkplus.clear();
		// Generate new candidates and put into c
		c = genCands(extendedC);
		timeGenCand = clock();
		totalGenCand += double(timeGenCand - timeGenPartInst);
	}
}



/**
* @bref This function builds the candidate tree
* @param sortF: a set of sorted features
* @retval tree<std::vector<char>> CandTree;
*/
std::vector<std::vector<char>> genCands(std::vector<std::vector<char>> oldCands)
{
	std::vector<std::vector<char>> newCands;
	// find the index of the last element in the oldCand
	int index;	
	for (auto c : oldCands)
	{
		auto it = find(sortF.begin(), sortF.end(), c.back());
		// If element was found
		if (it != sortF.end())
		{
			// calculating the index		
			index = it - sortF.begin();			
		}
		// Generate new cands
		for (auto t = index + 1; t < sortF.size(); t++)
		{
			c.push_back(sortF[t]);
			newCands.push_back(c);
			c.pop_back();
		}
	}
	// Return
	return newCands;
}



/**
* @bref This function computes feature-object neighbor sets.
* @param SN: star neighbor
* @retval FONS: neighbors that are classified by features and instances
*/
void buildFONS(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN)
{	
	for (auto const& item : SN)
	{
		std::map<char, std::set<ObjWithoutCoord>> inner;
		std::map<char, std::set<ObjWithoutCoord>>::iterator itinner;
		for (auto const& nei : item.second)
		{
			itinner = inner.find(nei.feature);
			if(itinner != inner.end()) // find out
			{
				itinner->second.insert(nei);
			}
			else
			{
				inner.insert({ nei.feature, std::set<ObjWithoutCoord>{nei} });
			}
		}
		FONS.insert({ item.first, inner });
	}
}



/**
* @bref This function builds feature-object neighbor tree (Figure 3-2 in Reference[2])
* @param std::map<ObjWithoutCoord, std::map<char, std::set<ObjWithoutCoord>>> FONS: in global level
* @retval FONT: in globel level
*/
void buildFONT()
{
	tree<ObjWithoutCoord>::iterator begintree, headNode, childeNode, leafNode;
	// Initial tree
	begintree = FONT.begin();
	// Create a root
	ObjWithoutCoord rootNode = { 'R', 0, 0.0f };
	rootFONT = FONT.insert(begintree, rootNode);
	// Add all features
	for (auto const& f : vuf) // vuf saves utility of each feature
	{
		// create node
		ObjWithoutCoord nodef = { f.first, 0, 0.0f };
		headNode = FONT.append_child(rootFONT, nodef);
		// add instances in the FONS of the current feature node
		for (auto const& item : FONS)
		{
			if (item.first.feature == f.first && item.second.size()>0) // get instances of the current feature
			{
				// add neighbor feature
				childeNode = FONT.append_child(headNode, item.first); // e.g., child D1
				for (auto const& nei : item.second)
				{
					ObjWithoutCoord child = { nei.first, 0, 0.0f };
					leafNode = FONT.append_child(childeNode, child);					
					for (auto const& fnei : nei.second)
					{
						FONT.append_child(leafNode, fnei);
					}
				}				
			}
		}

	}	
}


/**
* @bref This function executes searching to get high utility patterns
* @param c, k
* @retval
*/
bool searchPatterns(std::vector<char> c, float minutil)
{
	bool flag;	
	// 1. If this is size 1 candidate, return true
	if (c.size() <= 1)
	{		
		flag = true;
		return flag;
	}
	else // 2. The size of candidate is larger than 1
	{		
		std::unordered_map<char, std::set<ObjWithoutCoord>> CSHBS;		
		generateParticipatingInstance(c, CSHBS);		
		
		clock_t timecalUtility = clock();
		if (CSHBS.size() > 0) // has real row instances
		{
			// calculate utiltiy of c
			//startFilHUCMP = clock();
			double uc = calculateUtilityOfc(c, CSHBS);			
			double lambdac = calculateLambda(uc);
			if (lambdac >= minutil)
			{
				std::string pattern;
				for (auto const& f : c)
				{
					pattern += f;
				}
				PkAll.insert({ pattern, lambdac });
				// go to the next candidate
				flag = true;
				return flag;				
			}
			else
			{				
				// calculating lu(c) needs vuf, uc, US
				double luc = calculateLuc(c, uc);				
				if (luc <= (1 - minutil))
				{
					// Calculate ubc, enc
					double ubc = calculateUbc(c);									
					double euc = calculateEuc(ubc, lambdac);					
					if (euc >= minutil)
					{
						flag = true;
						return flag;
					}
					else
					{
						flag = true;
						return flag;
					}
				}
				else
				{
					flag = false;
					return flag;
				}				
			}			
		}
		totalTimeCalUPIAndFilPrevCoLoPat += double(clock() - timecalUtility);
	}	
}



/**
* @bref This function calculates the utility of candidate c
* @param c: candidate
*		CSHBS: the participating instances of features in c
* @retval: the utility of c
*/
double calculateUtilityOfc(std::vector<char> c, std::unordered_map<char, std::set<ObjWithoutCoord>> CSHBS)
{
	double uc=0.0;
	// Get utility of each feature in c
	std::vector<float> vf(c.size());
	int t = 0; // index of vf
	for (auto const& item : CSHBS)
	{
		vf[t] = vuf.find(item.first)->second * float(item.second.size());
		t += 1;
	}
	// Calculate uc
	for (auto const& u : vf)
	{
		uc += u;
	}

	return uc;
}



/**
* @bref This function calculates the pattern utility ratio of c
* @param uc: the utility of c
* @retval: lambdac
*/
double calculateLambda(double uc)
{
	double lambdac = uc / US;
	return lambdac;
}


/**
* @bref This function calculates luc to prune
* @param uc: the utility of c
* @retval: luc
*/
double calculateLuc(std::vector<char> c, double uc)
{
	std::map<char, double>::iterator itvuf; // vuf: utility of each feature
	double para1 = 0.0;
	for (auto const& f : c)
	{
		itvuf = vuf.find(f);
		if (itvuf != vuf.end())
		{
			para1 += itvuf->second;
		}
	}
	double molecule = para1 - uc;
	double luc = molecule / US;

	return luc;
}


/**
* @bref This function calculates luc to prune
* @param uc: the utility of c
* @retval: luc
*/
double calculateUbc(std::vector<char> c)
{
	double ubc=0.0;
	// Find Ec: a set of extended features of c
	std::vector<char> Ec;
	// only needs to check the last feature in c in F
	for (int t = 0; t < sortF.size(); t++)
	{
		if (sortF[t] == c.back())
		{
			// get the remain features
			for (int j = t+1; j < sortF.size(); j++)
			{
				Ec.push_back(sortF[j]);
			}
			break;
		}
	}
	// Find nei(f) f in Ec
	std::map<char, std::set<ObjWithoutCoord>> neiEC;
	std::map<char, std::set<ObjWithoutCoord>>::iterator itneiEC;
	std::map<char, std::set<ObjWithoutCoord>>::iterator itfc;
	for (auto f : Ec)
	{
		for (auto item : FONS)
		{
			if (item.first.feature == f)
			{
				// check
				bool flag = true;
				for (auto fc : c)
				{
					itfc = item.second.find(fc);
					if (itfc == item.second.end() || itfc->second.size() == 0)
					{
						flag = false;
						break;
					}					
				}				
				if (flag)
				{
					itneiEC = neiEC.find(f);
					if (itneiEC != neiEC.end())
					{
						itneiEC->second.insert(item.first);
					}
					else
					{						
						neiEC.insert({ f, std::set<ObjWithoutCoord> {item.first} });
					}
				}
			}
		}
	}	
	// calculate ubc 	
	for (auto item : neiEC)
	{
		ubc = ubc + (vuf.find(item.first)->second * item.second.size());
	}
	ubc = ubc / US;

	return ubc;
}



/**
* @bref This function calculates euc
* @param ubc and lambdac
* @retval: euc
*/
double calculateEuc(double ubc, double lambdac)
{
	double euc = lambdac + ubc;
	return euc;
}


/**
* @bref This function generates the participating instances of c
* @param FONT, c
* @retval: CSHBS: save participating instances of each feature
*/
void generateParticipatingInstance(std::vector<char> c, std::unordered_map<char, std::set<ObjWithoutCoord>>& CSHBS)
{
	std::vector<std::vector<ObjWithoutCoord>> innervalue; // to save candidate row instances	
	tree<ObjWithoutCoord>::sibling_iterator ittree = FONT.begin(rootFONT);	
	while (ittree != FONT.end(rootFONT))
	{		
		// Only get the instances that matchs with c[0] the first feature in c		
		if (ittree.node->data.feature == c[0])
		{				
			tree<ObjWithoutCoord>::sibling_iterator itfirstF = FONT.begin(ittree);
			while (itfirstF != FONT.end(ittree))
			{
				// Get all children of this node				
				std::vector<char> childF;
				childF.push_back(itfirstF.node->data.feature);
				tree<ObjWithoutCoord>::sibling_iterator itsib = FONT.begin(itfirstF);
				while (itsib != FONT.end(itfirstF))
				{					
					childF.push_back(itsib.node->data.feature);
					++itsib;
				}
				// check all features of the current node is in c?			
				if (IsSubset(childF, c))
				{
					// get candidate row instances					
					std::unordered_map<char, std::vector<ObjWithoutCoord>> tempInnerValue;
					std::unordered_map<char, std::vector<ObjWithoutCoord>>::iterator ittempInnerValue;
					tempInnerValue.insert({ c[0], std::vector<ObjWithoutCoord>{itfirstF.node->data} }); // put the firt					
					// create an empty for the other features
					for (int j = 1; j < c.size(); j++)
					{
						tempInnerValue.insert({ c[j], std::vector<ObjWithoutCoord> {} });
					}									
					// put the child of the other feature in c
					itsib = FONT.begin(itfirstF);
					while (itsib != FONT.end(itfirstF))
					{
						if (std::find(c.begin(), c.end(), itsib.node->data.feature) != c.end()) // find out the instance that its feature is in c, then get all children					
						{
							// get all children
							std::vector<ObjWithoutCoord> instonef;
							tree<ObjWithoutCoord>::sibling_iterator itnei = FONT.begin(itsib);
							while (itnei != FONT.end(itsib))
							{
								instonef.push_back(itnei.node->data);
								++itnei;
							}
							// put into tempvalue
							ittempInnerValue = tempInnerValue.find(itsib.node->data.feature);
							if (ittempInnerValue != tempInnerValue.end())
							{
								ittempInnerValue->second.insert(ittempInnerValue->second.end(), instonef.begin(), instonef.end());
							}
						}
						++itsib;
					}
					// run cartesianProduct to get row instances
					innervalue = cartesianProduct(tempInnerValue);
					tempInnerValue.clear();					
					// check real row instances and update the table instance					
					validateRealRowInst(c, innervalue);	
					innervalue.clear();					
				}	
				// Terminal
				++itfirstF;
			}			
		}
		++ittree;
	}	
	// Put in CSHBS	
	if (realTableInst.size())
	{
		std::vector<std::set<ObjWithoutCoord>> Tb(c.size());		
		for (auto row : realTableInst)
		{			
			for (int t = 0; t < c.size(); t++)
			{
				Tb[t].insert(row[t]);
			}
		}		
		// put table instance of the current candidate into CIKplus
		CIkplus.insert({ c, realTableInst });
		// clear
		realTableInst.clear();
		// put into CSHBS		
		for (auto i = 0; i < c.size(); i++)
		{
			CSHBS.insert({ c[i], Tb[i]});
		}		
	}	
}



/**
* @bref This function checks wwhether a vector is a subset of another
* @param A, B: two vectors
* @retval: bool
*/
template <typename T>
bool IsSubset(std::vector<T> A, std::vector<T> B)
{
	std::sort(A.begin(), A.end());
	std::sort(B.begin(), B.end());
	return std::includes(A.begin(), A.end(), B.begin(), B.end());
}



/**
* @brief: This function generates cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(std::unordered_map<char, std::vector<ObjWithoutCoord>> tempInner)
{
	std::vector<std::vector<ObjWithoutCoord>> v;
	for (auto const& item : tempInner)
	{
		v.push_back(item.second);
	}
	std::vector<std::vector<ObjWithoutCoord>> resultCartesian;

	auto product = [](long long a, std::vector<ObjWithoutCoord>& b)
	{
		return a * b.size();
	};
	const long long N = accumulate(v.begin(), v.end(), 1LL, product);
	std::vector<ObjWithoutCoord> u(v.size());
	for (long long n = 0; n < N; ++n)
	{
		lldiv_t q{ n, 0 };
		for (long long i = v.size() - 1; 0 <= i; --i) {
			q = div(q.quot, v[i].size());
			u[i] = v[i][q.rem];
		}
		resultCartesian.push_back(u);
	}

	return resultCartesian;
}


/**
* @brief: This function validates real row instances
* @param: innervalue: candidate row instances
* std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; // table instances of size k patterns
* @retval: realTableInst: real row instances, global level, std::vector<std::vector<ObjWithoutCoord>> rowInstSet;
*/
void validateRealRowInst(std::vector<char> c, std::vector<std::vector<ObjWithoutCoord>> innervalue)
{
	int sizek = innervalue[0].size();
	if (sizek == 2)
	{
		CIk.insert({ c, innervalue });
		realTableInst.insert(realTableInst.end(), innervalue.begin(), innervalue.end());
	}
	else // size k>2 patterns, need to check relation ship
	{
		bool flag = true;
		// get the table instances of the sub candiate c[1-end]
		std::vector<char> remainc;
		remainc.insert(remainc.end(), c.begin() + 1, c.end());

		std::vector<ObjWithoutCoord> checkRowkPlus;

		std::map<std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>>::iterator itCIk;
		itCIk = CIk.find(remainc);
		if (itCIk != CIk.end()) // if c is a pattern, realTableInst has its table instances,if not, realTableInst has not table instances
		{
			for (auto row : innervalue)
			{
				checkRowkPlus.insert(checkRowkPlus.end(), row.begin() + 1, row.end());
				if (std::find(itCIk->second.begin(), itCIk->second.end(), checkRowkPlus) != itCIk->second.end())
				{
					// The current row instance is a clique					
					realTableInst.push_back(row);
				}
				checkRowkPlus.clear();
			}
		}
		else // realTableInst has not its table instance, we need to check
		{
			bool flag2 = true;
			for (auto row : innervalue)
			{
				for (int t = 0; t < sizek - 1; t++)
				{
					flag2 = std::includes(SN.find(row[t])->second.begin(), SN.find(row[t])->second.end(), row.begin() + t + 1, row.end());
					if (flag2 == false)
					{
						break;
					}
				}
				// check
				if (flag2) // real row instance
				{
					realTableInst.push_back(row);
				}
			}
		}			
	}
}


