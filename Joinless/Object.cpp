#include <iostream>
#include "Object.h"


// For ObjWithCoord
bool ObjWithCoord::operator==(const ObjWithCoord& obj) const
{
	return feature == obj.feature && instance == obj.instance;
}


bool ObjWithCoord::operator<(const ObjWithCoord& obj) const
{
	if (feature < obj.feature) return true;
	else if (feature == obj.feature && instance < obj.instance) return true;
	else return false;
}


// For ObjWithoutCoord
bool ObjWithoutCoord::operator==(const ObjWithoutCoord& obj) const
{
	return this->feature == obj.feature && this->instance == obj.instance;
}


bool ObjWithoutCoord::operator<(const ObjWithoutCoord& obj) const
{
	if (this->feature < obj.feature) return true;
	else if (this->feature == obj.feature && this->instance < obj.instance) return true;
	else return false;
}

