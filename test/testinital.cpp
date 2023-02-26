/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-6-24
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <iostream>
#include <vector>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "types.hpp"
#include "detect.h"
#include "compute.h"
#include "cvcannyapi.h"

using namespace std;
using namespace cv;
using namespace zgh;


int main(int argc, char* argv[]) {
  Mat board = imread(argv[1]);
  Mat image = imread(argv[1], IMREAD_GRAYSCALE);
  imshow("car", image);
  // calc the gradient
  int row = image.rows;
  int col = image.cols;
  double* angles = new double[row * col];
  calculateGradient3(image.data, row, col, angles);
  vector<shared_ptr<Ellipse> > ells;
  getValidInitialEllipseSet(image.data, angles, row, col, ells);
  cout << "Find " << ells.size() << " initial ellipse" << endl;
  for (int i = 0; i < (int)ells.size(); ++i) {
    auto ell = ells[i];
    // std::cout << ell->o << " " << ell->a << " " << ell->b << " " << ell->phi << " " << endl;
    ellipse(board,
      Point(ell->o.y, ell->o.x),
      Size(ell->a, ell->b),
      rad2angle(PI_2 - ell->phi),
      0,
      360,
      Scalar(0, 255, 0),
      1,
      8,
      0);
  }
  imshow("initial ellpise", board);
  waitKey(0);

  return 0;
}

