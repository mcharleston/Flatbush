/*
 * Climb.cpp
 *
 *  Created on: 16 Dec. 2018
 *      Author: mc7
 */

#include "Climb.h"

namespace flatbush {

Climb::Climb(const leafset& s, double zs, const leafset& f, unsigned int w, double zf) :
		start(s), startScore(zs), finish(f), walkLength(w), finalScore(zf) {

}

Climb::Climb() {
	// TODO Auto-generated constructor stub

}

Climb::~Climb() {
	// TODO Auto-generated destructor stub
}

std::ostream& operator<<(std::ostream& os, const Climb& cl) {
	os << strmLeafset(cl.start) << ',' << cl.startScore << ',' << strmLeafset(cl.finish) << ',' << cl.finalScore << ',' << cl.walkLength;
	return os;
}

} /* namespace flatbush */
