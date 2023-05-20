#pragma once

#ifndef PATTERN_H
#define PATTERN_H

#include <vector>
#include <string>

// With coordinate
struct Pattern
{
	std::string c; // pattern name	
	double upr; // the pattern utility ratio
	double Qc;  // the quality of a pattern
	
	// Compare two objects
	bool operator==(const Pattern& obj) const;
	bool operator<(const Pattern& obj) const;

};



#endif