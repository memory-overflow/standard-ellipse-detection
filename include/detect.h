/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-22
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#ifndef _INCLUDE_DETECT_H_
#define _INCLUDE_DETECT_H_

#include "types.hpp"


namespace zgh {

bool lineSegmentDetection(const double *image, int row, int col,
                          std::vector<std::shared_ptr<Lined> > &lines);

bool lsdgroups(const double *image, int row, int col, double scale,
               std::vector<std::shared_ptr<Arc> >& arcs);

bool getValidInitialEllipseSet(const uint8_t *image,
                               const double *angles,
                               int row, int col, 
                               std::vector<std::shared_ptr<Ellipse> > &ells,
                               int polarity = 0);
                              

bool generateEllipseCandidates(const uint8_t *image,
                               const double *angles,
                               int row, int col, 
                               std::vector<std::shared_ptr<Ellipse> > &ells, int polarity);


bool detectEllipse(const uint8_t *image, int row, int col,
                   std::vector<std::shared_ptr<Ellipse> > &ells,
                   int polarity = 0, double width = 2.0);

}
//namespace zgh


#endif // _INCLUDE_DETECT_H_
