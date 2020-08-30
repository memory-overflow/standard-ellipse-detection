/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/


#ifndef _INCLUDE_TYPES_H_
#define _INCLUDE_TYPES_H_

#include <memory>
#include <vector>
#include <stack>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "defines.h"
#include "unitily.h"

namespace zgh {


template <typename T>
class Point_ {
 public:
  T x, y;
  Point_();
  Point_(T _x, T _y);
  ~Point_();
  Point_(const Point_<T>& pt);
  T dot(const Point_<T>& pt) const;
  T cross(const Point_<T>& pt) const;
  double length() const;
  double angle(const Point_<T>& pt) const;
  Point_<T> rotate(double rad) const;
  bool operator < (const Point_<T>& pt) const;
  bool operator == (const Point_<T>& pt) const;
  Point_<T> operator + (const Point_<T>& pt) const;
  Point_<T> operator - (const Point_<T>& pt) const;
  Point_<T> operator * (T k) const;
  Point_<T> operator / (T k) const;
  void rotation(double rad);
  void normal();
  Point_<T>& operator = (const Point_<T>& pt);
  Point_<T>& operator += (const Point_<T>& pt);
  Point_<T>& operator -= (const Point_<T>& pt);
  Point_<T>& operator *= (T k);
  Point_<T>& operator /= (T k);
};

template<typename T>
using Vector_ = Point_<T>;

using Pointi = Point_<int>;
using Pointf = Point_<float>;
using Pointd = Point_<double>;
using Vectori = Vector_<int>;
using Vectorf = Vector_<float>;
using Vectord = Vector_<double>;
using Pixel = Point_<int>;



template <typename T>
std::istream& operator >> (std::istream& in, Point_<T>& pt);

template <typename T>
std::istream& operator >> (std::istream& in, Point_<T>&& pt);

template <typename T>
std::ostream& operator << (std::ostream& out, const Point_<T>& pt);

template <typename T>
std::ostream& operator << (std::ostream& out, const Point_<T>&& pt);



template <typename T>
class Line_ {
 public:
  Point_<T> sp, ep;
  Point_<T> center;
  double length;
  int polarity;
  double width;
  Vector_<T> dir;
  Line_();
  Line_(const Point_<T> _sp, const Point_<T> _ep, int _pol = 0);
  Line_(T _x1, T _y1, T _x2, T _y2, int _pol);
  Line_(const Line_<T>& lt);
  Line_(Line_<T>&& lt);
  ~Line_();
  double angle(const Line_<T>& lt) const;
  const std::vector<Pixel>& getregs() const;
  void addpixel(const Pixel& pix);
  void setregs(const std::vector<Pixel>& pv);
  void setregs(std::vector<Pixel>&& pv);
 private:
  std::vector<Pixel> regions;
};

using Linei = Line_<int>;
using Linef = Line_<float>;
using Lined = Line_<double>;


class Arc {
 public:
  double coverages;
  int polarity;
  double fitMat[6][6];
  std::vector<std::shared_ptr<Lined> > lines;
  Arc(std::shared_ptr<Lined>& lptr);
  Arc(std::shared_ptr<Lined>&& lptr);
  Arc(const Arc& arc);
  ~Arc();
  const std::vector<Pixel>& getregs() const;
  void setregs(const std::vector<Pixel>& pv);
  void setregs(std::vector<Pixel>&& pv);
  void merge(std::shared_ptr<Lined>& lptr);
  void merge(std::shared_ptr<Lined>&& lptr);
  void updateMatrix(const std::shared_ptr<Lined>& lptr);
 private:
  std::vector<Pixel> regions;
};

class Ellipse {
 public:
  Pointd o; // center point
  double a, b; // short, long axis length
  double phi;
  int polarity;
  double coefficients[6];
  double goodness;
  double coverangle;
  std::vector<Pixel> inliers;
  Ellipse();
  Ellipse(const Ellipse &e);
  Ellipse(Pointd center, double a, double b, double phi);
  Ellipse(double a, double b, double c, double d, double e, double f);
  bool isCircle();
  bool isLegal();
  bool equal(const std::shared_ptr<Ellipse> &eptr,
             double centers_distabce__threshold = MIN_ELLIPSE_THRESHOLD_LENGTH, 
             double semimajor_errorratio = 0.1,
             double semiminor_errorratio = 0.1,
             double angle_errorratio = 0.1,
             double iscircle_ratio = 0.9);
  double distopoint(Pointd p);
  Vectord getTangent(Pointd p);
  Vectord getTangent(Pointi p);
 private:
  
};

class RectIter {
 public:
  RectIter(const std::shared_ptr<Lined> &line);
  Pointd polys[4];
  Pixel np;
  RectIter operator ++ ();
  RectIter operator ++ (int);
  bool isEnd();
 private:
  int xmin, xmax;
  int ymin, ymax;
  double ang;
  void doinc();
  void calcYAxisRange(int x);
};

class EllipseIter {
 public:
  EllipseIter(const std::shared_ptr<Ellipse> &ell_t, double distance_tolerance);
  Pixel np;
  EllipseIter operator ++ ();
  EllipseIter operator ++ (int);
  bool isEnd();
 private:
  double a, b, c, d, e, f;
  const std::shared_ptr<Ellipse> &ell;
  double distance_tolerance;
  int dir;
  int xmin, xmax;
  std::stack<Pixel> points;
  bool isend;
  void doinc();
  void calcYAxis(int x);
};


template<typename T1, typename T2>
class TemplateSameType {
 public:
  operator bool () {
    return false;
  }
};

template<typename T1>
class TemplateSameType<T1, T1> {
 public:
  operator bool () {
    return true;
  }
};

template <typename return_type>
class FuncTimerDecorator {
 public:
  static double times;
  FuncTimerDecorator(std::string funcname) : funcname(funcname) {

  }
  template <typename T, typename... Args>
  return_type operator () (T func, Args&&... args){
    clock_t start = clock();
    return_type res = func(args...);
    clock_t finish = clock();
    double millis = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    std::cerr << "Using " << millis << " ms to call function " << funcname << std::endl;
    return res;
  }
 private:
  std::string funcname;
};


} // namespace zgh


namespace zgh {

/*--------------------------------------------------------------*/

template <typename T>
T min(const T& x, const T& y) {
  return x < y ? x : y;
}

template <typename T>
T max(const T& x, const T& y) {
  return x > y ? x : y;
}

template <typename T>
Point_<T>::Point_(): x(0), y(0) {

}

template <typename T>
Point_<T>::Point_(T _x, T _y): x(_x), y(_y) {

}

template <typename T>
Point_<T>::~Point_<T>() {

}

template <typename T>
Point_<T>::Point_(const Point_<T>& pt): x(pt.x), y(pt.y) {

}

template <typename T>
T Point_<T>::dot(const Point_<T>& pt) const {
  return x * pt.x + y * pt.y;
}

template <typename T>
T Point_<T>::cross(const Point_<T>& pt) const {
  return x * pt.y - y * pt.x;
}

template <typename T>
double Point_<T>::length() const {
  return sqrt(max(decltype(dot(*this))(0), dot(*this)));
}


template <typename T>
double Point_<T>::angle(const Point_<T>& pt) const {
  double temp = 1.0 * dot(pt) / length() / pt.length();
  return acos(min(max(0.0, temp), 1.0));
}

template <typename T>
Point_<T> Point_<T>::rotate(double rad) const {
  return Point_<T>(x * cos(rad) - y * sin(rad), x * sin(rad) + y * cos(rad));
}

template <typename T>
bool Point_<T>::operator < (const Point_<T>& pt) const {
  if (equal(x, pt.x)) {
    return y < pt.y;
  }
  return x < pt.x;
}

template <typename T>
bool Point_<T>::operator == (const Point_<T>& pt) const {
  return equal(x, pt.x) && equal(y, pt.y);
}

template <typename T>
Point_<T> Point_<T>::operator + (const Point_<T>& pt) const {
  return Point_<T>(x + pt.x, y + pt.y);
}

template <typename T>
Point_<T> Point_<T>::operator - (const Point_<T>& pt) const {
  return Point_<T>(x - pt.x, y - pt.y);
}

template <typename T>
Point_<T> Point_<T>::operator * (T k) const {
  return Point_<T>(x * k, y * k);
}

template <typename T>
Point_<T> Point_<T>::operator / (T k) const {
  if (equal(k, 0 * k)) {
    return Point_<T>(0, 0);
  }
  return Point_<T>(x / k, y / k);
}
template <typename T>
void Point_<T>::rotation(double rad) {
  T xx = x, yy = y;
  x = xx * cos(rad) - yy * sin(rad);
  y = xx * sin(rad) + yy * cos(rad);
}

template <typename T>
void Point_<T>::normal() {
  double len = length();
  if (equal(len, 0.0)) {
    x = y = 0;
  }
  x /= len;
  y /= len;
}

template <typename T>
Point_<T>& Point_<T>::operator = (const Point_<T> &pt) {
  x = pt.x;
  y = pt.y;
  return *this;
}

template <typename T>
Point_<T>& Point_<T>::operator += (const Point_<T> &pt) {
  x += pt.x;
  y += pt.y;
  return *this;
}

template <typename T>
Point_<T>& Point_<T>::operator -= (const Point_<T> &pt) {
  x -= pt.x;
  y -= pt.y;
  return *this;
}

template <typename T>
Point_<T>& Point_<T>::operator *= (T k) {
  x *= k;
  y *= k;
  return *this;
}

template <typename T>
Point_<T>& Point_<T>::operator /= (T k) {
  if (equal(k, 0 * k)) {
    x = y = 0;
    return *this;
  }
  x /= k;
  y /= k;
  return *this;
}


/*--------------------------------------------------------------*/

template <typename T>
Line_<T>::Line_() {
  width = 2.0;
}

template <typename T>
Line_<T>::Line_(const Point_<T> _sp, const Point_<T> _ep, int _pol): sp(_sp), ep(_ep), polarity(_pol) {
  dir = ep - sp;
  length = dir.length();
  if (!equal(length, 0.0)) {
    dir.normal();
  }
  width = 2.0;
}

template <typename T>
Line_<T>::Line_(T _x1, T _y1, T _x2, T _y2, int _pol): sp(_x1, _y1), ep(_x2, _y2), polarity(_pol) {
  dir = ep - sp;
  length = dir.length();
  if (!equal(length, 0.0)) {
    dir.normal();
  }
  width = 2.0;
}

template <typename T>
Line_<T>::Line_(const Line_<T>& lt): sp(lt.sp), ep(lt.ep), length(lt.length), polarity(lt.polarity),
  dir(lt.dir), regions(lt.regions) {
  width = 2.0;
}

template <typename T>
Line_<T>::Line_(Line_<T>&& lt): sp(lt.sp), ep(lt.ep), length(lt.length), polarity(lt.polarity),
  dir(lt.dir), regions(std::forward<Line_<T> >(lt.regions)) {
  width = 2.0;
}

template <typename T>
Line_<T>::~Line_() {
  std::vector<Pixel>().swap(regions);
}

template <typename T>
double Line_<T>::angle(const Line_<T>& lt) const {
  if (!equal(length, 0.0) && !equal(lt.length, 0.0)) {
    return dir.angle(lt.dir);
  } else {
    return NAN;
  }
}

template <typename T>
const std::vector<Pixel>& Line_<T>::getregs() const {
  return regions;
}

template <typename T>
void Line_<T>::addpixel(const Pixel& pix) {
  regions.push_back(pix);
}

template <typename T>
void Line_<T>::setregs(const std::vector<Pixel>& pv) {
  std::vector<Pixel>().swap(regions);
  regions = pv;
}

template <typename T>
void Line_<T>::setregs(std::vector<Pixel>&& pv) {
  std::vector<Pixel>().swap(regions);
  regions = std::forward<std::vector<Pixel> >(pv);
}



/*--------------------------------------------------------------*/

inline Arc::Arc(std::shared_ptr<Lined>& lptr): coverages(0), polarity(lptr->polarity) {
  memset(fitMat, 0, sizeof(double) * 6 * 6);
  updateMatrix(lptr);
  lines.push_back(lptr);
}

inline Arc::Arc(std::shared_ptr<Lined>&& lptr): coverages(0), polarity(lptr->polarity) {
  memset(fitMat, 0, sizeof(double) * 6 * 6);
  updateMatrix(lptr);
  lines.push_back(std::forward<std::shared_ptr<Lined> >(lptr));
}

inline Arc::Arc(const Arc& arc): coverages(arc.coverages), polarity(arc.polarity), lines(arc.lines) {
  memcpy(fitMat, arc.fitMat, sizeof(double) * 6 * 6);
}

inline Arc::~Arc() {
  
}

inline void Arc::updateMatrix(const std::shared_ptr<Lined>& lptr) {

  double D[2][6];
  // start point
  D[0][0] = sqr(lptr->sp.x);
  D[0][1] = lptr->sp.x * lptr->sp.y;
  D[0][2] = sqr(lptr->sp.y);
  D[0][3] = lptr->sp.x;
  D[0][4] = lptr->sp.y;
  D[0][5] = 1.0;


  // end point
  D[1][0] = sqr(lptr->ep.x);
  D[1][1] = lptr->ep.x * lptr->ep.y;
  D[1][2] = sqr(lptr->ep.y);
  D[1][3] = lptr->ep.x;
  D[1][4] = lptr->ep.y;
  D[1][5] = 1.0;

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      fitMat[i][j] += D[0][i] * D[0][j];
      fitMat[i][j] += D[1][i] * D[1][j];
    }
  }

}
inline void Arc::merge(std::shared_ptr<Lined>& lptr) {
  if (lptr->polarity != polarity) {
    std::cerr << "Different polarity segment cannot be merge!" << std::endl;
    return;
  }
  updateMatrix(lptr);
  lines.push_back(lptr);
  double start_angle = atan2(lines[lines.size() - 2]->dir.y, lines[lines.size() - 2]->dir.x);
  double end_angle = atan2(lptr->dir.y, lptr->dir.x);
  coverages += rotateAngle(start_angle, end_angle, polarity);
}

inline void Arc::merge(std::shared_ptr<Lined>&& lptr) {
  if (lptr->polarity != polarity) {
    std::cerr << "Different polarity segment cannot be merge!" << std::endl;
    return;
  }
  updateMatrix(lptr);
  double start_angle = atan2(lines[lines.size() - 1]->dir.y, lines[lines.size() - 1]->dir.x);
  double end_angle = atan2(lptr->dir.y, lptr->dir.x);
  coverages += rotateAngle(start_angle, end_angle, polarity);
  lines.push_back(std::forward<std::shared_ptr<Lined> >(lptr));
}



/*--------------------------------------------------------------*/

inline Ellipse::Ellipse() {
  a = b = phi = 0;
}

inline Ellipse::Ellipse(Pointd center, double a, double b, double phi) {
  this->o = center;
  this->a = a;
  this->b = b;
  this->phi = phi;
  coefficients[0] = sqr(b * cos(phi)) + sqr(a * sin(phi));
  coefficients[1] = 2.0 * (b * b - a * a) * sin(phi) * cos(phi);
  coefficients[2] = sqr(b * sin(phi)) + sqr(a * cos(phi));

  coefficients[3] = -2.0 * (sqr(b * cos(phi)) + sqr(a * sin(phi))) * o.x +
                     2.0 * (a * a - b * b) * sin(phi) * cos(phi) * o.y;

  coefficients[4] = -2.0 * (sqr(b * sin(phi)) + sqr(a * cos(phi))) * o.y +
                     2.0 * (a * a - b * b) * sin(phi) * cos(phi) * o.x;

  coefficients[5] = (sqr(b * cos(phi)) + sqr(a * sin(phi))) * sqr(o.x) + 
                    (sqr(b * sin(phi)) + sqr(a * cos(phi))) * sqr(o.y) + 
                    2.0 * (b * b - a * a) * o.x * o.y * sin(phi) * cos(phi) -
                    sqr(a * b);

}

inline Ellipse::Ellipse(double a, double b, double c, double d, double e, double f) {

  // ax^2 + bxy + cy^2 + dx + ey + f = 0;
  coefficients[0] = a;
  coefficients[1] = b;
  coefficients[2] = c;
  coefficients[3] = d;
  coefficients[4] = e;
  coefficients[5] = f;
  
  double thetarad = 0.5 * atan2(b, a - c);
  double cost = cos(thetarad);
  double sint = sin(thetarad);
  double sin_squared = sqr(sint);
  double cos_squared = sqr(cost);
  double cos_sin = sint * cost;
  double Ao = f;
  double Au = d * cost + e * sint;
  double Av = -d * sint + e * cost;
  double Auu = a * cos_squared + c * sin_squared + b * cos_sin;
  double Avv = a * sin_squared + c * cos_squared - b * cos_sin;
  if (zgh::equal(0.0, Auu) || zgh::equal(0.0, Avv)) {
    this->a = this->b = this->phi = 0;
  } else {
    double tuCentre = -Au / (2.0 * Auu);
    double tvCentre = -Av / (2.0 * Avv);
    double wCentre = Ao - Auu * tuCentre * tuCentre - Avv * tvCentre * tvCentre;
    double uCentre = tuCentre * cost - tvCentre * sint;
    double vCentre = tuCentre * sint + tvCentre * cost;
    double Ru = -wCentre / Auu;
    double Rv = -wCentre / Avv;
    if (Ru <= 0 || Rv <= 0) {
      this->a = this->b = this->phi = 0;
      return;
    }
    Ru = sqrt(Ru);
    Rv = sqrt(Rv);
    this->o = Pointd(uCentre, vCentre);
    this->a = Ru;
    this->b = Rv;
    this->phi = thetarad;
    if (Ru < Rv) {
      std::swap(this->a, this->b);
      if (this->phi < 0) {
        this->phi += PI_2;
      } else {
        this->phi -= PI_2;
      }
      if (this->phi < -PI_2) {
        this->phi += PI;
      }
      if (this->phi > PI_2) {
        this->phi -= PI;
      }
    }
  }
}

inline Ellipse::Ellipse(const Ellipse &e): o(e.o), a(e.a), b(e.b), polarity(e.polarity) {
  memcpy(coefficients, e.coefficients, sizeof(double) * 6);
}

inline bool Ellipse::isLegal() {
  if (min(a, b) <= 3 * MIN_ELLIPSE_THRESHOLD_LENGTH) {
    return false;
  }
  return true;
}
inline bool Ellipse::isCircle() {
  return isLegal() && std::fabs(a - b) <= FEPS;
}

inline bool Ellipse::equal(const std::shared_ptr<Ellipse> &eptr,
                          double centers_distabce_threshold, 
                          double semimajor_errorratio,
                          double semiminor_errorratio,
                          double angle_errorratio,
                          double iscircle_ratio) {
  
  if (!eptr) {
    return false;
  }
  bool con1 = std::fabs(this->o.x - eptr->o.x) < centers_distabce_threshold
           && std::fabs(this->o.y - eptr->o.y) < centers_distabce_threshold
           && std::fabs(this->a - eptr->a) / max(this->a, eptr->a) < semimajor_errorratio
           && std::fabs(this->b - eptr->b) / min(this->b, eptr->b) < semiminor_errorratio;

  double temp = angle_diff(this->phi, eptr->phi);
  double phi_diff = min(PI - temp, temp);
  bool con2 = (this->b / this->a) >= iscircle_ratio;
  bool con3 = (eptr->b / eptr->a) >= iscircle_ratio;
  bool con4 = ((con2 && con3) || (!con2 && !con3 &&  phi_diff <= angle_errorratio * PI));

  return con1 && con4;
}

inline double Ellipse::distopoint(Pointd p) {

  Vectord pt = p - o;
  pt.rotation(-phi);
  double ae2 = a * a;
  double be2 = b * b;
  double fe2 = a * a - b * b;
  if (isCircle()) {
    return std::fabs((p - o).length() - a);
  }

  double X = sqr(pt.x);
  double Y = sqr(pt.y);

  double delta = sqr(X + Y + fe2) - 4 * fe2 * X;
  double A = (X + Y + fe2 - sqrt(delta)) / 2.0;
  double ah = sqrt(A);
  double bh2 = fe2 - A;
  double term = A * be2 + ae2 * bh2;
  double xi = ah * sqrt(ae2 * (be2 + bh2) / term);
  double yi = b * sqrt(bh2 * (ae2 - A) / term);
  double d = sqr(pt.x - xi) + sqr(pt.y - yi);
  d = min(d, sqr(pt.x + xi) + sqr(pt.y - yi));
  d = min(d, sqr(pt.x - xi) + sqr(pt.y + yi));
  d = min(d, sqr(pt.x + xi) + sqr(pt.y + yi));
  return sqrt(d);
}


inline Vectord Ellipse::getTangent(Pointd p) {
  double dx = coefficients[1] * p.x + 2.0 * p.y * coefficients[2] + coefficients[4];
  double dy = -1.0 * (coefficients[3] + 2.0 * coefficients[0] * p.x + coefficients[1] * p.y);
  Vectord tan(dx, dy);
  tan.rotation(-PI_2);
  tan.normal();
  return tan;
}

inline Vectord Ellipse::getTangent(Pointi p) {
  return getTangent(Pointd(1.0 * p.x, 1.0 * p.y));
}



/*--------------------------------------------------------------*/

inline RectIter::RectIter(const std::shared_ptr<Lined> &line) {
  Pointd sp = line->sp, ep = line->ep;
  Vectord dir = line->ep - line->sp;
  if (dir.y < 0) {
    dir.x = -dir.x;
    dir.y = -dir.y;
    std::swap(sp, ep);
  }
  dir.normal();
  ang = atan2(dir.x, dir.y);
  if (ang < 0) {
    /* 
        --------------------- -> y
        |                   |
        |         /         |
        |        /          |
        |       /           |
        |      /            |
        |                   |
        ---------------------
        |
        v
        x
     */

    // Clockwise
    polys[0] = sp + dir.rotate(-PI / 2) * line->width / 2.0;
    polys[1] = sp + dir.rotate(PI / 2) * line->width / 2.0;
    polys[2] = ep + dir.rotate(PI / 2) * line->width / 2.0;
    polys[3] = ep + dir.rotate(-PI / 2) * line->width / 2.0;
  } else {
    /* 
        --------------------- -> y
        |                   |
        |     \             |
        |      \            |
        |       \           |
        |        \          |
        |                   |
        ---------------------
        |
        v
        x
     */

    // Counterclockwise
    polys[0] = ep + dir.rotate(-PI / 2) * line->width / 2.0;
    polys[1] = ep + dir.rotate(PI / 2) * line->width / 2.0;
    polys[2] = sp + dir.rotate(PI / 2) * line->width / 2.0;
    polys[3] = sp + dir.rotate(-PI / 2) * line->width / 2.0;
  }
  xmin = (int)ceil(polys[2].x);
  xmax = (int)floor(polys[0].x);
  np.x = xmin;
  while (np.x <= xmax) {
    calcYAxisRange(np.x);
    if (ymin <= ymax) {
      break;
    }
    ++np.x;
  }
  np.y = ymin;
}


inline bool RectIter::isEnd() {
  return np.x > xmax;
}


inline RectIter RectIter::operator ++ () {
  doinc();
  return *this;
}

inline RectIter RectIter::operator ++ (int) {
  auto res = *this;
  doinc();
  return res;
}

inline void RectIter::doinc() {
  ++np.y;
  if (np.y > (double) ymax) {
    ++np.x;
    calcYAxisRange(np.x);
    np.y = ymin;
  }
}

inline void RectIter::calcYAxisRange(int x) {
  if (std::fabs(ang) <= FEPS || std::fabs(std::fabs(ang) - PI_2) <= FEPS) {
    ymin = min(min(polys[0].y, polys[1].y), polys[2].y);
    ymax = max(max(polys[0].y, polys[1].y), polys[2].y);
    return;
  }
  if (ang < 0) {
    // calculate ymin 
    if ((double)x <= polys[1].x) {
      double k = ((double)x - polys[2].x) / (polys[1].x - polys[2].x);
      Pointd p = polys[2] + (polys[1] - polys[2]) * k;
      ymin = (int)ceil(p.y);
    } else {
      double k = ((double)x - polys[1].x) / (polys[0].x - polys[1].x);
      Pointd p = polys[1] + (polys[0] - polys[1]) * k;
      ymin = (int)ceil(p.y);
    }
    // calculate ymax
    if ((double)x <= polys[3].x) {
      double k = ((double)x - polys[2].x) / (polys[3].x - polys[2].x);
      Pointd p = polys[2] + (polys[3] - polys[2]) * k;
      ymax = (int)floor(p.y);
    } else {
      double k = ((double)x - polys[3].x) / (polys[0].x - polys[3].x);
      Pointd p = polys[3] + (polys[0] - polys[3]) * k;
      ymax = (int)floor(p.y);
    }
  } else {
    // calculate ymin 
    if ((double)x <= polys[3].x) {
      double k = ((double)x - polys[2].x) / (polys[3].x - polys[2].x);
      Pointd p = polys[2] + (polys[3] - polys[2]) * k;
      ymin = (int)ceil(p.y);
    } else {
      double k = ((double)x - polys[3].x) / (polys[0].x - polys[3].x);
      Pointd p = polys[3] + (polys[0] - polys[3]) * k;
      ymin = (int)ceil(p.y);
    }
    // calculate ymax
    if ((double)x <= polys[1].x) {
      double k = ((double)x - polys[2].x) / (polys[1].x - polys[2].x);
      Pointd p = polys[2] + (polys[1] - polys[2]) * k;
      ymax = (int)floor(p.y);
    } else {
      double k = ((double)x - polys[1].x) / (polys[0].x - polys[1].x);
      Pointd p = polys[1] + (polys[0] - polys[1]) * k;
      ymax = (int)floor(p.y);
    }
  }
}

/*--------------------------------------------------------------*/

inline EllipseIter::EllipseIter(const std::shared_ptr<Ellipse> &ell_t,
                         double _distance_tolerance): ell(ell_t) {
  distance_tolerance = _distance_tolerance;
  a = ell->coefficients[0];
  b = ell->coefficients[1];
  c = ell->coefficients[2];
  d = ell->coefficients[3];
  e = ell->coefficients[4];
  f = ell->coefficients[5];

  if (b * b - 4 * a * c >= 0) {
    std::cerr << "error" << std::endl;
    isend = true;
    return;
  } 
  double B = (2 * b * e - 4 * d * c) / (b * b - 4 * a * c);
  double C = (e * e - 4 * f * c) / (b * b - 4 * a * c);
  double delta = B * B - 4 * C;
  if (delta < 0) {
    isend = true;
    return;
  }
  xmin = (int)ceil((-B - sqrt(delta)) / 2.0);
  xmax = (int)floor((-B + sqrt(delta)) / 2.0);
  np.x = xmin;
  dir = 1;
  isend = false;
  calcYAxis(np.x);
  while (np.x < xmax && points.empty()) {
    calcYAxis(++np.x);
  }
  if (!points.empty()) {
    np = points.top();
    points.pop();
  }
}

inline bool EllipseIter::isEnd() {
  return isend;
}

inline EllipseIter EllipseIter::operator ++ () {
  doinc();
  return *this;
}

inline EllipseIter EllipseIter::operator ++ (int) {
  auto res = *this;
  doinc();
  return res;
}

inline void EllipseIter::doinc() {
  if (points.empty()) {
    if (dir == 1) {
      while (np.x < xmax && points.empty()) {
        calcYAxis(++np.x);
      }
      if (points.empty()) {
        np.x = xmax + 1;
        dir = -1;
      }
    }

    if (dir == -1) {
      while (np.x > xmin && points.empty()) {
        calcYAxis(--np.x);
      }
      if (points.empty()) {
        isend = true;
      }
    }
  }

  if (!isend) {
    np = points.top();
    points.pop();
  }
}

inline void EllipseIter::calcYAxis(int x) {
  double A = c;
  double B = (b * x + e);
  double C = f + d * x + a * x * x;
  double delta = B * B - 4 * A * C;
  if (delta < 0) {
    return;
  }
  double ya = (-B + sqrt(delta)) / (2.0 * A);
  double yb = (-B - sqrt(delta)) / (2.0 * A);
  double y;
  if (dir == 1) {
    y = max(ya, yb);
  } else {
    y = min(ya, yb);
  }
  int ymid = (int)floor((ya + yb) / 2.0);
  int yy = (int)y;
  for (int idy = yy + 1; dir == 1 || idy < ymid; ++idy) {
    if (ell->distopoint(Pointd(1.0 * x, 1.0 * idy)) <= distance_tolerance) {
      points.push(Pointi(x, idy));
    } else {
      break;
    }
  }
  for (int idy = yy; dir == -1 || idy >= ymid; --idy) {
    if (ell->distopoint(Pointd(1.0 * x, 1.0 * idy)) <= distance_tolerance) {
      points.push(Pointi(x, idy));
    } else {
      break;
    }
  }
  if ((dir == 1 && x == xmax) || (dir == -1 && x == xmin)) {
    for (int id = 1; id <= (int)ceil(distance_tolerance); ++id) {
      for (int idy = yy + 1; ; ++idy) {
        if (ell->distopoint(Pointd(1.0 * (x + id * dir), 1.0 * idy)) <= distance_tolerance) {
          points.push(Pointi(x + id * dir, idy));
        } else {
          break;
        }
      }

      for (int idy = yy; ; --idy) {
        if (ell->distopoint(Pointd(1.0 * (x + id * dir), 1.0 * idy)) <= distance_tolerance) {
          points.push(Pointi(x + id * dir, idy));
        } else {
          break;
        }
      }
    }
  }

}

/*--------------------------------------------------------------*/

template <typename T>
std::istream& operator >> (std::istream& in, Point_<T>& pt) {
  in >> pt.x >> pt.y;
  return in;
}

template <typename T>
std::ostream& operator << (std::ostream& out, const Point_<T>& pt) {
  out << std::fixed << std::setprecision(5) << "[" << pt.x << ", " << pt.y << "]";
  return out;
}

template <typename T>
std::istream& operator >> (std::istream& in, Point_<T>&& pt) {
  in >> pt.x >> pt.y;
  return in;
}

template <typename T>
std::ostream& operator << (std::ostream& out, const Point_<T>&& pt) {
  out << std::fixed << std::setprecision(5) << "[" << pt.x << ", " << pt.y << "]";
  return out;
}

} // namespace zgh



#endif // _INCLUDE_TYPES_H_
