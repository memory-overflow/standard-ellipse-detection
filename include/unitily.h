/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/


#ifndef _INCLUDE_UNITILT_H_
#define _INCLUDE_UNITILT_H_

namespace zgh {

bool equal(int x, int y);

bool euqal(float x, float y);

bool equal(double x, double y);

double angle2rad(double angle);

double rad2angle(double rad);

double angle_diff(double a, double b);

double angle_diff_signed(double a, double b);

bool inRect(int h, int w, int x, int y);

double rotateAngle(double start_angle, double end_angle, int polarity);

double sqr(double a);

int sqr(int a);

double log_gamma_windschitl(double x);

double log_gamma_lanczos(double x);

double nfa(int n, int k, double p, double logNT);

} // namespace zgh


#endif // _INCLUDE_UNITILT_H_
