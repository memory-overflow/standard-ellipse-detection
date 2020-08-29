/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-27
 *Author:zhengaohong@minieye.cc
 *Version 1.0.1
*/


#ifndef _INCLUDE_CVCANNYAPI_H_
#define _INCLUDE_CVCANNYAPI_H_
#include <stdint.h>

namespace MINI {

bool calculateGradient3(const uint8_t *image, int row, int col, double * angles);

} // namespace MINI


#endif // _INCLUDE_CVCANNYAPI_H_
