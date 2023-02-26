/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "types.hpp"
#include "detect.h"
#include "compute.h"

using namespace std;
using namespace cv;
using namespace zgh;



int main(int argc, char* argv[]) {
  Mat image = imread(argv[1], IMREAD_GRAYSCALE);
  imshow("origin", image);
  vector<shared_ptr<Lined> > lines;
  double scale = 0.8;
  double sigma_scale = 0.6;
  int row = (int)ceil(1.0 * image.rows * scale);
  int col = (int)ceil(1.0 * image.cols * scale);
  double *data = new double[row * col];
  gaussianSampler(image.data, image.rows, image.cols, data, row, col, scale, sigma_scale);
  image.release();
  
  lineSegmentDetection(data, row, col, lines);

  std::cout << "Detect " << lines.size() << " lines" << endl;
  
  Mat boarda(row, col, CV_8UC3, Scalar(255, 255, 255));
  
  Mat boardb(row, col, CV_8UC3, Scalar(255, 255, 255));
  
  for (auto& li : lines) {
    Scalar color(0, rand() % 256, rand() % 256);
    line(boarda, Point2f(li->sp.y, li->sp.x), Point2f(li->ep.y, li->ep.x), color, 2, 8, 0);
    for (auto& pix : li->getregs()) {
      boardb.at<Vec3b>(pix.x, pix.y)[0] = color[0];
      boardb.at<Vec3b>(pix.x, pix.y)[1] = color[1];
      boardb.at<Vec3b>(pix.x, pix.y)[2] = color[2];
    }
  }
  imshow("line", boarda);
  imshow("region", boardb);
  waitKey(0);
  delete [] data; 
  return 0;
}

