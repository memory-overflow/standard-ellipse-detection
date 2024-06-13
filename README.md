# 介绍
本库提供了一个可以在工业中稳定使用的椭圆检测方法。对于图像中的标准、明显、完整、大小在 100x100 像素以上的椭圆的检测效果非常好，速度也很快。
这个库的实现参考了论文 https://arxiv.org/pdf/1810.03243v4.pdf。

微信小程序搜索"椭圆识别"可以体验效果（小程序首次使用需要启动服务，第一张图片可能会失败，多试几次）。

![image](https://user-images.githubusercontent.com/15645203/226512866-71bbaab5-6e43-41ef-bac7-23b61733269d.png)

也可以线上体验：http://43.154.37.202




# 效果展示


![图1](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test12_result.jpg)

![图2](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test6_result.jpg)

![图3](https://github.com/memory-overflow/standard-ellipse-detection/blob/master/images/test9_result.jpg)



# 安装
## Ubuntu
### Install opencv
opencv 需要通过源码安装，ubuntu opencv 的安装可以参考博客 https://blog.csdn.net/Arthur_Holmes/article/details/100146463, 注意本项目中一定要安装 opencv3 的版本，否则可能会有兼容性问题。服务器版本无需安装图像化相关的库。

### Install lapack
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

### Install ellipse-detection
执行下面的操作完成安装
```
git clone https://github.com/memory-overflow/standard-ellipse-detection.git
cd standard-ellipse-detection
mkdir build && cd build
cmake ..
make
sudo make install
```

## Centos
### Install opencv
opencv 需要通过源码安装，centos 的 opencv 的安装可以参考博客 https://www.jianshu.com/p/1cb1ca235eb3, 注意本项目中一定要安装 opencv3 的版本，否则可能会有兼容性问题。服务器版本无需安装图像化相关的库。


### Install lapack
lapack 是一个矩阵运算的库，也是需要源码安装。
先下载[lapack源码](https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.9.1.tar.gz)，lapack 是 gfortran 写的，要先`yum install gcc-gfortran`安装 gfortran 才能正常编译。

由于 lapack 的编译需要 cmake3, 先安装 cmake3, `yum install cmake3`

执行下面的操作完成 lapack 的安装
```
tar -xzvf lapack-3.9.1.tar.gz && cd lapack-3.9.1
mkdir build && cd build
cmake3 ..
make -j7
sudo make install
sudo ldconfig
cd ..
sudo cp LAPACKE/include/*.h /usr/local/include/
```

### Install ellipse-detection
执行下面的操作完成安装
```
git clone https://github.com/memory-overflow/standard-ellipse-detection.git
cd standard-ellipse-detection
mkdir build && cd build
cmake3 ..
make
sudo make install
```
提供一个 centos server 的镜像，可以作为基础镜像开发，里面装好的依赖和库。
```
docker pull jisuanke/standard-ellipse-detection-centos7
```


## Other system
不建议在 mac 和 windows 上使用此库。mac 的兼容性比较差，我尝试在 mac 上编译，最后因为 gfortran 链接库不兼容放弃。
windows 上的编译和使用的IDE有关，如需要使用，按照 opencv，lapack, ellipsedetect 的顺序依次编译。

其他类unix的操作系统比如安卓可以自行参考 ubuntu 和 centos 的流程编译安装。

# 接口和使用方法
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
输出的图像坐标系
![image](https://github.com/memory-overflow/standard-ellipse-detection/assets/15645203/a6f388df-3a51-4dac-8d27-a106d496362a)



# 测试
提供了1个测试工具，可以查看算法效果。需要桌面版的操作系统才能显示图片，如果是服务器版本的操作系统，需要注释掉 imshow 部分。
```
cmake3 .. -DBUILD_TESTING=ON
make
./bin/testdetect [image_dir1] [image_dir2] [image_dir3] ...
```




有问题欢迎联系 zhengaohong@gmail.com, wechat: islands___
