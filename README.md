# 介绍
本库提供了一个可以在工业中稳定使用的椭圆检测方法。对于图像中的标准、明显、完整、大小在 100x100 像素以上的椭圆的检测效果非常好，速度也很快。
这个库的实现参考了论文 https://arxiv.org/pdf/1810.03243v4.pdf。

# 安装
## Ubuntu
### install opencv
opencv 需要通过源码安装，opencv 的安装可以参考博客 https://blog.csdn.net/Arthur_Holmes/article/details/100146463, 注意本项目中一定要安装 opencv3 的版本，否则可能会有兼容性问题。

### install lapack
lapack 是一个矩阵运算的库，也是需要源码安装。
先下载[lapack源码](https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.9.1.tar.gz)，lapack 是 gfortran 写的，要先`sudo apt-get install gfortran`安装gfortran 才能正常编译。

执行下面的操作完成安装
```
tar -xzvf lapack-3.9.1.tar.gz && cd lapack-3.9.1
mkdir build && cd build
cmake ..
make -j7
sudo make install
sudo ldconfig
cd ..
sudo cp LAPACKE/include/*.h /usr/local/include/
```

### install ellipse-detection
执行下面的操作完成安装
```
git clone https://github.com/memory-overflow/standard-ellipse-detection.git
cd standard-ellipse-detection
mkdir build && cd build
cmake ..
make
sudo make install
```

### 接口和使用方法
代码中引用头文件`#include "ellipse_detection/detect.h"`，然后引入namespace zgh。

接口说明
```
bool detectEllipse(const uint8_t *image, int height, int width,
                   std::vector<std::shared_ptr<Ellipse> > &ells,
                   int polarity = 0, double line_width = 2.0);
```
- 输入参数:
    - image 图像原始数据，灰度图，彩色图需要先转换成灰度图，并且转换成一维数组输入
    - height 图像高度
    - width 图像宽度
    - polarity 表示椭圆极性，-1、0、1, 默认为 0，检测所有极性。
    - line_width 椭圆线宽，单位像素，推荐使用默认值
- 输出
    - ells 检测到的椭圆列表

关于 Ellipse 结构的说明
```
Pointd o; // 椭圆中心点坐标
double a, b; // 短半轴长度，长半轴长度
double phi; // 椭圆偏角，单位为弧度
int polarity; // 椭圆极性
double goodness; // 椭圆评分
double coverangle; // 椭圆角度完整程度
std::vector<Pixel> inliers; // 构成的像素点
```

### 测试
提供了1个测试工具，可以查看算法效果。需要桌面版的操作系统才能显示图片，如果是服务器版本的操作系统，需要注释掉 imshow 部分。
```
cmake .. -DBUILD_TESTING=ON
make
./bin/testdetect [image_dir1] [image_dir2] [image_dir3] ...
```



# 效果展示
![图1](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test12_result.jpg)

![图2](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test6_result.jpg)

![图3](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test9_result.jpg)

有问题欢迎联系 zhengaohong@gmail.com, wechat: islands___
