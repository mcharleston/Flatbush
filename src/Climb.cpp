/*
 * Climb.cpp
 *
 *  Created on: 16 Dec. 2018
 *      Author: mc7
 */

#include "Climb.h"

namespace flatbush {

leafset start;
leafset finish;
unsigned int walkLength;
double startScore;
double finalScore;

Climb::Climb(const leafset& s, double zs, const leafset& f, unsigned int w, double zf) :
		start(s), finish(f), walkLength(w), startScore(zs), finalScore(zf) {
}

Climb::Climb() :
		walkLength(0.0), startScore(0.0), finalScore(0.0) {
	// TODO Auto-generated constructor stub

}

Climb::~Climb() {
	// TODO Auto-generated destructor stub
}

std::ostream& operator<<(std::ostream& os, const Climb& cl) {
	os << strmLeafset(cl.start) << ',' << cl.startScore << ',' << strmLeafset(cl.finish) << ',' << cl.finalScore
			<< ',' << card(cl.start) << ',' << cl.walkLength;
	return os;
}

} /* namespace flatbush */
