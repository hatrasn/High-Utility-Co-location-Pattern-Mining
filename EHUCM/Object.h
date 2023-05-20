
#ifndef OBJECT_H
#define OBJECT_H

// With coordinate
struct ObjWithCoord
{
	char feature;
	int instance;	
	float x;
	float y;
	double u;

	// Compare two objects
	bool operator==(const ObjWithCoord& obj) const;
	bool operator<(const ObjWithCoord& obj) const;

};

// No coordinate
struct ObjWithoutCoord
{
	char feature;
	int instance;
	double u;

	bool operator==(const ObjWithoutCoord & obj) const;
	bool operator<(const ObjWithoutCoord & obj) const;

};

// Hash function for unordered_map
struct myHashFuncWithoutCoord
{
	std::size_t operator()(const ObjWithoutCoord& inst) const;
};


#endif