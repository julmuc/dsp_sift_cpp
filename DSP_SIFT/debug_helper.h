/** @internal
 ** @file   debug_helper.h
 ** @brief   debug_helper
 ** @author Julian heuser
 **/
#ifndef DEBUGHELPER_H
#define DEBUGHELPER_H

#include <iostream>


#ifdef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...) std::cout << "Func.: " << __FUNCTION__ << "  L.: " << __LINE__ << "  Value: " << M << std::endl;
#endif



#endif // DEBUGHELPER_H