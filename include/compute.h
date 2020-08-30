/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/


#ifndef _INCLUDE_COMPUTE_H_
#define _INCLUDE_COMPUTE_H_

#include "types.hpp"

namespace zgh {

bool calculateGradient(const uint8_t *image, int row, int col, double *angles);

std::shared_ptr<Ellipse> calcElliseParam(const std::shared_ptr<Arc> &arc1, 
                                        const std::shared_ptr<Arc> &arc2,
                                        const double *angles,
                                        int row,
                                        int col);

bool regionLimitation(const std::shared_ptr<Arc> &arc1, const std::shared_ptr<Arc> &arc2);

bool gaussianSampler(const uint8_t *ori_data, int ori_row, int ori_col,
                     double *data, int row, int col,
                     double scale, double sigma_scale);

std::shared_ptr<Ellipse> fitEllipse(const std::vector<Pixel> &points);

} // namespace zgh

#endif // _INCLUDE_COMPUTE_H_
