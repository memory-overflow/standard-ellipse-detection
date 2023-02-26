/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-6-27
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <cstring>
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "cvcannyapi.h"
#include "defines.h"

namespace zgh
{

    using namespace cv;
#include "stdio.h"
    // #include <iostream>

    bool calculateGradient3(const uint8_t* data, int row, int col, double* angles)
    {

        int low_thresh, high_thresh; //canny thresh

        cv::Mat DX, DY;
        cv::Mat gray(row, col, CV_8UC1, (uint8_t*)data);
        cv::Mat Angle(row, col, CV_64FC1, angles);
        cv::Mat Angle_masked = cv::Mat::zeros(row, col, CV_64FC1);
        cv::Mat Edge = cv::Mat::zeros(row, col, CV_8UC1);
        cv::Size size = gray.size();

        Sobel(gray, DX, CV_64FC1, 1, 0, 3);
        Sobel(gray, DY, CV_64FC1, 0, 1, 3);

        double maxVal;
        cv::Mat magnitude;
        cv::magnitude(DX, DY, magnitude);

        // 自动设置canny阈值
        {
            double minVal;
            cv::Point minLoc;
            cv::Point maxLoc;

            cv::minMaxLoc(magnitude, &minVal, &maxVal, &minLoc, &maxLoc);
            cv::normalize(magnitude, magnitude, 0.0, 1.0, cv::NORM_MINMAX);

        }
        {
            // set magic numbers
            const int NUM_BINS = 64;
            const double percent_of_pixels_not_edges = 0.65;
            const double threshold_ratio = 0.3;
            // compute histogram
            int bin_size = cvFloor(maxVal / float(NUM_BINS) + 0.5f) + 1;
            if (bin_size < 1) bin_size = 1;
            int bins[NUM_BINS] = { 0 };
            magnitude.forEach<double>([&bins, bin_size](double& pixel, const int* position) {
                int hgf = int(pixel);
                (bins[hgf / bin_size])++;
                });
            // / % Select the thresholds
            float total(0.f);
            float target =
                float(size.height * size.width * percent_of_pixels_not_edges);
            int low_thresh, high_thresh(0);

            while (total < target) {
                total += bins[high_thresh];
                high_thresh++;
            }
            high_thresh *= bin_size;
            low_thresh = cvFloor(threshold_ratio * float(high_thresh));
        }
        Canny(gray, Edge, low_thresh, high_thresh);

        for (int idx = 0; idx < row; ++idx) {
            double* _dx = DX.ptr<double>(idx);
            double* _dy = DY.ptr<double>(idx);
            uchar* _e = Edge.ptr<uchar>(idx);
            for (int idy = 0; idy < col; ++idy) {
                if (_e[idy] > 0) {
                    angles[idx * col + idy] =
                        atan2((double)_dy[idy], (double)_dx[idy]);  // calculate gradient
                }
                else {
                    angles[idx * col + idy] = ANGLE_NOT_DEF;
                }
            }
        }

        {
            cv::normalize(DX, DX, 0.0, 1, cv::NORM_MINMAX);
            cv::normalize(DY, DY, 0.0, 1, cv::NORM_MINMAX);
            imshow("gray", gray);
            imshow("DX", DX);
            imshow("DY", DY);
            imshow("magnitude", magnitude);

            imshow("angle", Angle);
            imshow("canny", Edge);
        }

        return true;
    }

} // namespace zgh

