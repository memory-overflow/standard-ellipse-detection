# 介绍
由于业界没有好的椭圆检测算法，所以我写一个库用来做椭圆检测。这个库对于图像中的标准、明显、完整、大小在 100x100 像素以上的椭圆的检测效果非常好，速度也很快。
这个库的实现参考了论文 ttps://arxiv.org/pdf/1810.03243v4.pdf。

# ubuntu 的使用方法
1. 首先需要安装两个库的支持，opencv 库，这个可以搜一下网上的教程安装一下。第二个库是一个矩阵运算的库lapack，需要源码安装。
先下载[lapack源码](https://github.com/Reference-LAPACK/lapack/archive/v3.9.0.tar.gz)，这个库是gfortran写的，所以要先`sudo apt-get install gfortran`安装gfortran。
然后
```
tar -xzvf lapack-3.9.0.tar.gz && cd lapack-3.9.0
mkdir build && cd build
cmake ..
make -j7
sudo make install
sudo ldconfig
sudo cp sudo cp LAPACKE/include/*.h /usr/local/include/
```

2. 安装ellipse-detection库
```
git clone https://github.com/memory-overflow/standard-ellipse-detection.git
cd standard-ellipse-detection
mkdir build && cd build
cmake ..
make
sudo make install
```

3. 测试
提供了1个测试工具，可以查看算法效果。
```
cmake .. -DBUILD_TESTING=ON
make
./bin/testdetect [image_dir1] [image_dir2] [image_dir3] ...
```

4. 接口和使用方法
代码中引用头文件`#include "ellipse_detection/detect.h"`，接口说明
```
bool detectEllipse(const uint8_t *image, int height, int width,
                   std::vector<std::shared_ptr<Ellipse> > &ells,
                   int polarity = 0, double line_width = 2.0);
```
- 输入:
    - image 图像原始数据，按照"BRG"排列
    - height 图像高度
    - width 图像宽度
    - polarity 表示椭圆极，默认为 0。
    - line_width 椭圆线宽
- 输出
    - ells 检测到的椭圆列表

关于 Ellipse 结构的说明
```
Pointd o; // 椭圆中心点
double a, b; // 短半轴线，长半轴长度
double phi; // 偏角，单位为弧度
int polarity; // 极性
double goodness; // 椭圆评分
double coverangle; // 椭圆角度完整程度
std::vector<Pixel> inliers;
```

# 效果图
[图1](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test12_result.jpg)

[图2](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test6_result.jpg)

[图3](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test9_result.jpg)