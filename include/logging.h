/* logging.h - using compile-time flags to turn on/off boost.logging
 * since the simulation is all about performance. 
 * Severity levels are those of boost: trace, debug, info, warning, error, fatal 
 * There is no filtering involved. The severity is set upon compilation by setting the flags
 * __TRACE__, __DEBUG__, __INFO__, __WARNING__, __ERROR__, __FATAL__ */
 
#ifndef __LOGGING_H_INCLUDED__
#define __LOGGING_H_INCLUDED__

#include <boost/log/trivial.hpp>

#ifndef __TRACE__
#define LOG_TRACE(X)
#endif
#ifdef __TRACE__
#define LOG_TRACE(X) BOOST_LOG_TRIVIAL(trace) << X ;
#define __DEBUG__
#define __INFO__
#define __WARNING__
#define __ERROR__
#define __FATAL__
#endif

#ifndef __DEBUG__
#define LOG_DEBUG(X)
#endif
#ifdef __DEBUG__
#define LOG_DEBUG(X) BOOST_LOG_TRIVIAL(debug) << X ;
#define __INFO__
#define __WARNING__
#define __ERROR__
#define __FATAL__
#endif

#ifndef __INFO__
#define LOG_INFO(X)
#endif
#ifdef __INFO__
#define LOG_INFO(X) BOOST_LOG_TRIVIAL(info) << X ;
#define __WARNING__
#define __ERROR__
#define __FATAL__
#endif

#ifndef __WARNING__
#define LOG_WARNING(X)
#endif
#ifdef __WARNING__
#define LOG_WARNING(X) BOOST_LOG_TRIVIAL(warning) << X ;
#define __ERROR__
#define __FATAL__
#endif

#ifndef __ERROR__
#define LOG_ERROR(X)
#endif
#ifdef __ERROR__
#define LOG_ERROR(X) BOOST_LOG_TRIVIAL(error) << X ;
#define __FATAL__
#endif

#ifndef __FATAL__
#define LOG_FATAL(X)
#endif
#ifdef __FATAL__
#define LOG_FATAL(X) BOOST_LOG_TRIVIAL(fatal) << X ;
#endif

#endif // __LOGGING_H_INCLUDED__