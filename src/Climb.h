/*
 * Climb.h
 *
 *  Created on: 16 Dec. 2018
 *      Author: mc7
 */

#ifndef SRC_CLIMB_H_
#define SRC_CLIMB_H_

#include <ostream>
#include "project.h"

namespace flatbush {

class Climb {
public:
	leafset start;
	leafset finish;
	unsigned int walkLength;
	double startScore;
	double finalScore;
	Climb();
	Climb(const leafset& s, double zs, const leafset& f, unsigned int w, double zf);
	virtual ~Climb();
};

std::ostream& operator<<(std::ostream& os, const Climb& cl);

} /* namespace flatbush */

#endif /* SRC_CLIMB_H_ */
