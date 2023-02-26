/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-6-31
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <iostream>
#include <string>
#include <vector>

#include "detect.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "types.hpp"

using namespace std;
using namespace zgh;
using namespace cv;

int main(int argc, char *argv[])
{
  // if (argc <= 1) {
  //   std::cout << "[Usage]: testdetect [image_dir1] [image_dir2] [image_dir3] ..." << std::endl;
  // }
#if 0
  VideoCapture cap(0);
  Mat board;
  // int cnt = 0;
  while (1)
  {
    cap >> board;
    // resize(board, board, Size(board.cols / 3, board.rows / 3));
    if (board.empty())
    {
      break;
    }
    Mat image;
    cvtColor(board, image, COLOR_BGR2GRAY);
    // cv::imwrite("test" + to_string(cnt++) + ".jpg", image);
    vector<shared_ptr<Ellipse>> ells;
    int row = image.rows;
    int col = image.cols;
    double width = 2.0;
    FuncTimerDecorator<int>("detectEllipse")(detectEllipse, image.data, row, col, ells, NONE_POL, width);
    cout << "Find " << ells.size() << " ellipse" << endl;
    for (int i = 0; i < (int)ells.size(); ++i)
    {
      auto ell = ells[i];
      std::cout << "coverangle : " << ell->coverangle << ",\tgoodness : " << ell->goodness << ",\tpolarity : " << ell->polarity << endl;
      ellipse(board,
              Point(ell->o.y, ell->o.x),
              Size(ell->a, ell->b),
              rad2angle(PI_2 - ell->phi),
              0,
              360,
              Scalar(0, 255, 0),
              width,
              8,
              0);
    }
    imshow("origin", image);
    imshow("detect", board);
    waitKey(0);
  }
#else
  for (int i = 1; i < argc; ++i)
  {
    cv::Mat board = cv::imread(argv[i]);
    cv::Mat image = cv::imread(argv[i], cv::IMREAD_GRAYSCALE);
    // cv::resize(image, image, cv::Size(image.cols / 2, image.rows / 2));
    // cv::resize(board, board, cv::Size(board.cols / 2, board.rows / 2));
    cv::imshow("origin image", image);
    vector<shared_ptr<Ellipse>> ells;
    int row = image.rows;
    int col = image.cols;
    double width = 2.0;
    FuncTimerDecorator<int>("detectEllipse")(detectEllipse, image.data, row, col, ells, NONE_POL, width);
    cout << "Find " << ells.size() << " ellipse(s)" << endl;
    for (int i = 0; i < (int)ells.size(); ++i)
    {
      auto ell = ells[i];
      std::cout << "coverangle : " << ell->coverangle
                << ",\tgoodness : " << ell->goodness
                << ",\tpolarity : "
                << ell->polarity
                << endl;
      cv::ellipse(board,
                  cv::Point(ell->o.y, ell->o.x),
                  cv::Size(ell->a, ell->b),
                  rad2angle(PI_2 - ell->phi),
                  0,
                  360,
                  cv::Scalar(0, 255, 0),
                  width,
                  8,
                  0);
    }
    cv::imshow("detected result", board);
    cv::waitKey(0);
  }
#endif
  return 0;
}
