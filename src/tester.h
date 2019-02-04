/*
 * tester.h
 *
 *  Created on: 26/05/2014
 *      Author: mac
 */

#ifndef TESTER_H_
#define TESTER_H_

//#include "utility/alignment.h"

namespace flatbush {

typedef float real;

class tester {
private:
	real value;
//	Alignment* A;
public:
	tester(real r = 0.0) : value(r) {}
	virtual ~tester();
	void testCountPatterns();
};

}	// end namespace
#endif /* TESTER_H_ */
