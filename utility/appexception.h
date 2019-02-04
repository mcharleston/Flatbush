/*
 * appexception.h
 *
 *  Created on: 15Aug.,2016
 *      Author: mac
 */

#ifndef UTILITY_APPEXCEPTION_H_
#define UTILITY_APPEXCEPTION_H_

#include <iostream>
#include <stdexcept>

namespace flatbush {

class app_exception: public std::runtime_error {
public:
	virtual ~app_exception() {}
	explicit app_exception(const char* msg) : runtime_error(msg) {
		std::cout << msg << std::endl;
	}
	explicit app_exception(const std::string& msg) : runtime_error(msg) {
		std::cout << msg.c_str() << std::endl;
	}
};

} /* namespace flatbush */

#endif /* UTILITY_APPEXCEPTION_H_ */
