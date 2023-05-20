
#ifndef OBJECT_H
#define OBJECT_H

// With coordinate
struct ObjWithCoord
{
	char feature;
	int instance;
	double x;
	double y;
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


#endif