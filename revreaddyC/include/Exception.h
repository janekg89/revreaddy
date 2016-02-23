/* Exception.h - derived from std::exception, throwing a string. */

#ifndef __EXCEPTION_H_INCLUDED__
#define __EXCEPTION_H_INCLUDED__
#include <exception>
#include <string>

struct Exception : public std::exception {
	Exception(const std::string& inMsg) : msg(inMsg) {}
	std::string msg;
	const char* what() const throw() { return msg.c_str(); }
};

#endif //__EXCEPTION_H_INCLUDED__