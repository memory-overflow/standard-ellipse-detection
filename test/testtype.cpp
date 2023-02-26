#include <cassert>
#include "types.hpp"


#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace zgh;
using namespace cv;


int main() {
  // test unitily.h
  assert(equal(1, 3) == false);

  assert(equal(1, 1) == true);

  assert(equal(-0, 0) == true);

  assert(equal(1.3, 1.300000000000001) == true);

  assert(equal(-1.56, 5.546) == false);

  assert(equal(1.0, 1.34) == false);

  assert(equal(angle2rad(45.0), PI_4) == true);

  assert(equal(rad2angle(PI_4), 45.0) == true);


  // test types.hpp
  Pointi a(1, 2), b(3, 4);

  assert(Pointi(a) == a);

  assert(a.dot(b) == b.dot(a));

  assert(a.dot(b) == 11);

  assert(a.cross(b) == -2);

  assert(a.cross(b) == -b.cross(a));

  assert(equal(a.length(), sqrt(5.0)) == true);

  Pointd c(1.0, sqrt(3.0)), d(1.0, 0);

  assert(equal(rad2angle(c.angle(d)), 60.0) == true);

  assert(c.rotate(angle2rad(-30.0)) == Pointd(sqrt(3.0), 1.0));

  c.rotation(angle2rad(-30.0));


  assert(c == Pointd(sqrt(3.0), 1.0));

  assert(a + b == Pointi(4, 6));

  assert(a - b == Pointi(-2, -2));


  assert(a / 1 == a);

  assert(a * 2 == Pointi(2, 4));

  (a += b) *= 2;

  assert(a == Pointi(8, 12));

  auto e = Pointd(-4.0, 4.0);

  e.normal();

  assert(e == Pointd(-sqrt(0.5), sqrt(0.5)));
  
  assert(fabs(nfa(193, 145, 0.125, 13.6128) - 74.2627) <= 1e-3);

  cout << "All test passed" << endl;


  // test RectIter

  Mat boardb(720, 1280, CV_8UC3, Scalar(255, 255, 255));
  shared_ptr<Lined> line = make_shared<Lined>(200, 200, 200, 500, 1);
  line->width = 4;
  for (RectIter iter(line); !iter.isEnd(); ++iter) {
    Scalar color(0, 0, 255);
    auto pix = iter.np;
    boardb.at<Vec3b>(pix.x, pix.y)[0] = color[0];
    boardb.at<Vec3b>(pix.x, pix.y)[1] = color[1];
    boardb.at<Vec3b>(pix.x, pix.y)[2] = color[2];
  }
  cv::line(boardb, Point2f(line->sp.y, line->sp.x), Point2f(line->ep.y, line->ep.x), Scalar(0, 0, 0), 2, 8, 0);
  
  line = make_shared<Lined>(200, 200, 500, 200, 1);
  line->width = 4;
  for (RectIter iter(line); !iter.isEnd(); ++iter) {
    Scalar color(0, 255, 0);
    auto pix = iter.np;
    boardb.at<Vec3b>(pix.x, pix.y)[0] = color[0];
    boardb.at<Vec3b>(pix.x, pix.y)[1] = color[1];
    boardb.at<Vec3b>(pix.x, pix.y)[2] = color[2];
  }
  cv::line(boardb, Point2f(line->sp.y, line->sp.x), Point2f(line->ep.y, line->ep.x), Scalar(0, 0, 0), 2, 8, 0);
  
  
  imshow("test", boardb);
  waitKey();

  std::shared_ptr<Ellipse> ell = make_shared<Ellipse>(Pointd(200, 200), 100, 50, 0.344);
  ellipse(boardb,
          Point(ell->o.y, ell->o.x),
          Size(ell->a, ell->b),
          rad2angle(PI_2 - ell->phi),
          0,
          360,
          Scalar(0, 255, 0),
          1,
          8,
          0);
  Pointd p(ell->o);
  assert(equal(ell->distopoint(Pointd(ell->o.x, ell->o.y)), ell->b));
  
  // test EllipseIter
  for (EllipseIter iter(ell, 1.0); !iter.isEnd(); ++iter) {
    Scalar color(255, 0, 0);
    auto pix = iter.np;
    if (inRect(720, 1280, pix.x, pix.y)) {
      boardb.at<Vec3b>(pix.x, pix.y)[0] = color[0];
      boardb.at<Vec3b>(pix.x, pix.y)[1] = color[1];
      boardb.at<Vec3b>(pix.x, pix.y)[2] = color[2];
    }
  }

  for (EllipseIter iter(ell, 1.0); !iter.isEnd(); ++iter) {
    Scalar color(255, 0, 0);
    auto pix = iter.np;
    if (inRect(720, 1280, pix.x, pix.y)) {
      cv::Mat temp = boardb.clone();
      Vectord tan = ell->getTangent(pix);
      Point2f st(pix.y, pix.x);
      Point2f en(pix.y + tan.y * 50, pix.x + tan.x * 50);
      cv::line(temp, st, en, color, 2, 8, 0);
      
      Vectord tan1 = tan.rotate(PI_2);
      Point2f st1(pix.y - tan1.y * 50, pix.x - tan1.x * 50);
      Point2f en1(pix.y + tan1.y * 50, pix.x + tan1.x * 50);
      cv::line(temp, st1, en1, color, 2, 8, 0);

      imshow("test", temp);
      waitKey(30);
    }
  }



  return 0;
}
