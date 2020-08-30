/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/


#ifndef _INCLUDE_DEFINES_H_
#define _INCLUDE_DEFINES_H_

#include <cmath>

namespace zgh {

#ifndef M_LN10
#define M_LN10 2.30258509299404568402    // ln10
#endif /* !M_LN10 */


#define ANGLE_NOT_DEF -1024.0

#define NOTUSED 0

#define USED 1

#define RELATIVE_ERROR_FACTOR 100.0

#define NONE_POL 0

#define SAME_POL 1

#define OPP_POL -1

#define TABSIZE 100000

#define log_gamma(x) ((x) > 15.0? log_gamma_windschitl(x) : log_gamma_lanczos(x))

#define PI 3.14159265358979323846

#define PI_2 1.57079632679489661923

#define PI_4 0.78539816339744830962

#define PI_8 0.392699081

const double DEPS = 1e-8;

const double EPS = DEPS;

const float FEPS = 1e-4;

const double ANGLE_TH = 22.5;

const double GRAD_THRESHOLD = 2.0 / sin(PI * ANGLE_TH / 180.0);

const double MIN_ELLIPSE_THRESHOLD_LENGTH = 2.0;

const double REGION_LIMITATION_DIS_TOLERACE = -3.0 * MIN_ELLIPSE_THRESHOLD_LENGTH;

} // namespace zgh

#endif // _INCLUDE_DEFINES_H_
