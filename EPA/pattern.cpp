#include <iostream>

#include "pattern.h"


bool Pattern::operator==(const Pattern& p) const
{
	return this->c == p.c;
}


bool Pattern::operator<(const Pattern& p) const
{
	if (this->upr > p.upr) return true;	
	else return false;
}

