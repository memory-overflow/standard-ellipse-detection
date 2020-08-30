/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <lapacke.h>


#include "compute.h"
#include "defines.h"


namespace zgh {

bool calculateGradient(const uint8_t *image, int row, int col, double *angles) {
  double *mod = new double[row * col];
  double normsum = 0.0;
  for (int idx = 0; idx < row; ++idx) {
    angles[idx * col] = ANGLE_NOT_DEF;
    angles[idx * col + col - 1] = ANGLE_NOT_DEF;
    mod[idx * col] = ANGLE_NOT_DEF;
    mod[idx * col + col - 1] = ANGLE_NOT_DEF;
  }
  for (int idy = 0; idy < col; ++idy) {
    angles[idy] = ANGLE_NOT_DEF;
    angles[(row - 1) * col + idy] = ANGLE_NOT_DEF;
    mod[idy] = ANGLE_NOT_DEF;
    mod[(row - 1) * col + idy] = ANGLE_NOT_DEF;
  }
  for (int idy = 1; idy < col - 1; ++idy) {
    for (int idx = 1; idx < row - 1; ++idx) {
      int addr = idx * col + idy;
      /*
        Norm 2 computation using 3x3 pixel window:
		        A B C
		        D E F
			      G H I
		     and
		       com1 = C - G,  com2 = I - A
		     Then
		       gx = C + 2F + I - (A + 2D + G) = com1 + com2 + 2(F - D)   horizontal difference
		       gy = G + 2H + I - (A + 2B + C) = -com1 + com2 + 2(H - B)   vertical difference
		     com1 and com2 are just to avoid 2 additions.
		   */
      double com1, com2;
      double dx, dy;
      com1 = (double)image[addr - col + 1] - (double)image[addr + col - 1];
      com2 = (double)image[addr + col + 1] - (double)image[addr - col - 1];
      dx = com1 + com2 + 2.0 * ((double)image[addr + 1] - (double)image[addr - 1]);
      dy = -com1 + com2 + 2.0 * ((double)image[addr + col] - (double)image[addr - col]);
      dx /= (8.0 * 255);
      dy /= (8.0 * 255);
      double norm_square = dx * dx + dy * dy;
      normsum += norm_square;
      double norm = sqrt(norm_square);
      mod[addr] = norm;
      angles[addr] = atan2(dy, dx);
    }
  }
  
  double threshold = sqrt(normsum / (row * col)) / 3.0;  // auto threshold
  for (int idx = 1; idx < row - 1; ++idx) {
    for (int idy = 1; idy < col - 1; ++idy) {
      int addr = idx * col + idy;
      if (mod[addr] < threshold) {
        angles[addr] = ANGLE_NOT_DEF;
        continue;
      }
      double val = angles[addr];
      if ((val > -PI_8 && val <= PI_8) || val <= -7.0 * PI_8 || val > 7.0 * PI_8) {
        if (mod[addr] <= mod[addr + 1] || mod[addr] <= mod[addr - 1]) {
          angles[addr] = ANGLE_NOT_DEF;
        }
      } else if ((val > PI_8 && val <= 3.0 * PI_8) || (val > -7.0 * PI_8 && val <= -5.0 * PI_8)) {
        if (mod[addr] <= mod[addr - col - 1] || mod[addr] <= mod[addr + col + 1]) {
          angles[addr] = ANGLE_NOT_DEF;
        }
      } else if ((val > 3.0 * PI_8 && val <= 5.0 * PI_8) || (val > -5.0 * PI_8 && val <= -3.0 * PI_8)) {
        if (mod[addr] <= mod[addr - col] || mod[addr] <= mod[addr + col]) {
          angles[addr] = ANGLE_NOT_DEF; 
        }
      } else {
        if (mod[addr] <= mod[addr - col + 1] || mod[addr] < mod[addr + col - 1]) {
          angles[addr] = ANGLE_NOT_DEF;
        }
      }
    }
  }

  delete [] mod;
  return true;
}


static std::shared_ptr<Ellipse> fitEllipse(double *S) {
  double C[36];
  memset(C, 0, sizeof(double) * 36);
  C[0 * 6 + 2] = 2;
  C[1 * 6 + 1] = -1;
  C[2 * 6 + 0] = 2;


  double alphar[6], alphai[6], beta[6];
  double vl[36] = {0};
  double vr[36] = {0};
	char JOBVL = 'N';
	char JOBVR = 'V';
	int fitN = 6;
	double fitWork[64];
	int workLen = 64;
	int info;
  dggev_(&JOBVL, &JOBVR, &fitN, S, &fitN, C, &fitN, alphar, alphai,
            beta, vl, &fitN, vr, &fitN, fitWork, &workLen, &info);
  if (info == 0) {
    int index = -1;
    for (int i = 0; i < 6; ++i) {
      if (alphar[i] >= -2.2204460492503131e-014 && alphai[i] == 0 && beta[i] != 0) {
        index = i;
      }
    }
    if (index == -1) {
      double temp = -0.005;
      for (int i = 0; i < 6; ++i) {
        if (alphar[i] >= temp  && alphai[i] == 0 && beta[i] != 0) {
				  temp = alphar[i];
				  index = i; // vr[:,i], vr 第 i 列对应的特征向量则为拟合参数
			  }
      }
    }
    if (index != -1) {
      if(vr[6 * index + 0] < 0) {
        return std::make_shared<Ellipse>(-vr[6 * index + 0], -vr[6 * index + 1], -vr[6 * index + 2],
                                         -vr[6 * index + 3], -vr[6 * index + 4], -vr[6 * index + 5]);
			} else {
        return std::make_shared<Ellipse>(vr[6 * index + 0], vr[6 * index + 1], vr[6 * index + 2],
                                         vr[6 * index + 3], vr[6 * index + 4], vr[6 * index + 5]);
			}
    }
  }
  return nullptr;
}

std::shared_ptr<Ellipse> fitEllipse(const std::vector<Pixel> &points) {
  int pnum = points.size();
  double *D = new double[6 * pnum];
  for (int id = 0; id < pnum; ++id) {
    D[6 * id + 0] = sqr(points[id].x);
    D[6 * id + 1] = points[id].x * points[id].y;
    D[6 * id + 2] = sqr(points[id].y);
    D[6 * id + 3] = points[id].x;
    D[6 * id + 4] = points[id].y;
    D[6 * id + 5] = 1;
  }
  double S[36];
  for (int idx = 0; idx < 6; ++idx) {
    for (int idy = 0; idy < 6; ++idy) {
      S[idx * 6 + idy] = 0;
      for (int k = 0; k < pnum; ++k) {
        S[idx * 6 + idy] += D[k * 6 + idx] * D[k * 6 + idy];
      }
    }
  }
  delete [] D;
  return fitEllipse(S);
}

static std::shared_ptr<Ellipse> fitEllipse(const std::vector<Pixel> &points1, const std::vector<Pixel> &points2) {
  int pnum = points1.size() + points2.size();
  double *D = new double[6 * pnum];
  for (int id = 0; id < (int)points1.size(); ++id) {
    D[6 * id + 0] = sqr(points1[id].x);
    D[6 * id + 1] = points1[id].x * points1[id].y;
    D[6 * id + 2] = sqr(points1[id].y);
    D[6 * id + 3] = points1[id].x;
    D[6 * id + 4] = points1[id].y;
    D[6 * id + 5] = 1;
  }

  int offset = (int)points1.size();
  for (int id = offset; id < pnum; ++id) {
    D[6 * id + 0] = sqr(points2[id - offset].x);
    D[6 * id + 1] = points2[id - offset].x * points2[id - offset].y;
    D[6 * id + 2] = sqr(points2[id - offset].y);
    D[6 * id + 3] = points2[id - offset].x;
    D[6 * id + 4] = points2[id - offset].y;
    D[6 * id + 5] = 1;
  }

  double S[36];
  for (int idx = 0; idx < 6; ++idx) {
    for (int idy = 0; idy < 6; ++idy) {
      S[idx * 6 + idy] = 0;
      for (int k = 0; k < pnum; ++k) {
        S[idx * 6 + idy] += D[k * 6 + idx] * D[k * 6 + idy];
      }
    }
  }
  delete [] D;
  return fitEllipse(S);
}

std::shared_ptr<Ellipse> calcElliseParam(const std::shared_ptr<Arc> &arc1, const std::shared_ptr<Arc> &arc2,
                                         const double *angles, int row, int col) {
  
  double S[36];
  if (!arc2) {
    for (int idx = 0; idx < 6; ++idx) {
      for (int idy = 0; idy < 6; ++idy) {
        S[idx * 6 + idy] = arc1->fitMat[idx][idy];
      }
    }
  } else {
    for (int idx = 0; idx < 6; ++idx) {
      for (int idy = 0; idy < 6; ++idy) {
        S[idx * 6 + idy] = arc1->fitMat[idx][idy] + arc2->fitMat[idx][idy];
      }
    }
  }
  std::shared_ptr<Ellipse> ell = fitEllipse(S);
  if (!ell || !ell->isLegal()) {
    return nullptr;
  }
  
  if (!inRect(row, col, ell->o.x, ell->o.y) || max(ell->a, ell->b) > min(row, col)) {
    return nullptr;
  }
  
  int support_cnt, inlier_cnt;
  // verification first arc1
  bool significant1 = true;
  std::vector<Pixel> support_regs1;
  for (auto &line : arc1->lines) {
    support_cnt = inlier_cnt = 0;
    line->width = 3 * MIN_ELLIPSE_THRESHOLD_LENGTH;
    for (RectIter iter(line); !iter.isEnd(); ++iter) {
      auto pix = iter.np;
      if (inRect(row, col, pix.x, pix.y)) {
        double temp = angles[pix.x * col + pix.y];
        if (!equal(temp, ANGLE_NOT_DEF)) {
          double point_normalx = ell->coefficients[0] * pix.x + 
                                (ell->coefficients[1] * pix.y + ell->coefficients[3]) / 2.0;
          double point_normaly = ell->coefficients[2] * pix.y + 
                                (ell->coefficients[1] * pix.x + ell->coefficients[4]) / 2.0;
          double point_normal;
          if (line->polarity == SAME_POL) {
            point_normal = atan2(-point_normalx, -point_normaly);
          } else {
            point_normal = atan2(point_normalx, point_normaly);
          }
          ++inlier_cnt;
          if (angle_diff(point_normal, temp) <= PI / 10.0) {
            ++support_cnt;
            support_regs1.push_back(pix);
          }
        }
      }
    }
    if (support_cnt == 0
        || (1.0 * support_cnt <= 0.7 * line->length
        && 1.0 * support_cnt / inlier_cnt <= 0.6)) {
      significant1 = false;
      break;
    }
  }

  // end verification first arc1

  if (!arc2) {
    if (significant1) {
      // 再次拟合提高质量
      std::shared_ptr<Ellipse> ell_t = fitEllipse(support_regs1);
      if (ell_t && ell_t->equal(ell, 3 * MIN_ELLIPSE_THRESHOLD_LENGTH, 0.1, 0.1, 0.1, 0.9)) {
        ell = ell_t;
      }
    }
    ell->polarity = arc1->polarity;
    return ell;
  }
  
  if (!significant1) {
    return nullptr;
  }

  
  // verification second arc2
  bool significant2 = true;
  std::vector<Pixel> support_regs2;
  for (auto &line : arc2->lines) {
    support_cnt = inlier_cnt = 0;
    line->width = 3 * MIN_ELLIPSE_THRESHOLD_LENGTH;
    for (RectIter iter(line); !iter.isEnd(); ++iter) {
      auto pix = iter.np;
      if (inRect(row, col, pix.x, pix.y)) {
        double temp = angles[pix.x * col + pix.y];
        if (!equal(temp, ANGLE_NOT_DEF)) {
          double point_normalx = ell->coefficients[0] * pix.x + 
                                (ell->coefficients[1] * pix.y + ell->coefficients[3]) / 2.0;
          double point_normaly = ell->coefficients[2] * pix.y + 
                                (ell->coefficients[1] * pix.x + ell->coefficients[4]) / 2.0;
          double point_normal;
          if (line->polarity == 1) {
            point_normal = atan2(-point_normalx, -point_normaly);
          } else {
            point_normal = atan2(point_normalx, point_normaly);
          }
          ++inlier_cnt;
          if (angle_diff(point_normal, temp) <= PI / 10.0) {
            ++support_cnt;
            support_regs2.push_back(pix);
          }
        }
      }
    }
    if (support_cnt == 0
        || (1.0 * support_cnt <= 0.7 * line->length
        && 1.0 * support_cnt / inlier_cnt <= 0.6)) {
      significant2 = false;
      break;
    }
  }
  
  // end verification first arc2

  if (significant1 && significant2) {
    std::shared_ptr<Ellipse> ell_t = fitEllipse(support_regs1, support_regs2);
    if (!ell_t) {
      return nullptr;
    }
    double semimajor_errorratio, semiminor_errorratio, iscircle_ratio;
    if (ell->a <= 50) {
      semimajor_errorratio = 0.25;
    } else if (ell->a <= 100) {
      semimajor_errorratio = 0.15;
    } else {
      semimajor_errorratio = 0.1;
    }

    if (ell->b <= 50) {
      semiminor_errorratio = 0.25;
    } else if (ell->b <= 100) {
      semiminor_errorratio = 0.15;
    } else {
      semiminor_errorratio = 0.1;
    }

    if (ell->a <= 50 && ell->b <= 50) {
      iscircle_ratio = 0.75;
    } else if (50 <= ell->a && ell->a <= 100 && 50 <= ell->b && ell->b <= 100) {
      iscircle_ratio = 0.85;
    } else {
      iscircle_ratio = 0.9;
    }
    if (ell_t && ell_t->equal(ell, 3 * MIN_ELLIPSE_THRESHOLD_LENGTH, 
        semimajor_errorratio, semiminor_errorratio, 0.1, iscircle_ratio)) {
      ell = ell_t;
      ell->polarity = arc1->polarity;
      return ell;
    }
  }
  return nullptr;
}

static bool regionLimitation_reverse(const std::shared_ptr<Arc> &arc1, const std::shared_ptr<Arc> &arc2) {
  if (!arc1 || !arc2) {
    return false;
  }
  Lined line1((*arc1->lines.begin())->sp, (*arc1->lines.rbegin())->ep, 0);
  Lined line2((*arc2->lines.begin())->sp, (*arc2->lines.rbegin())->ep, 0);
  Vectord arc_dir_s1, arc_dir_e1, arc_dir_m1;
  Vectord arc_dir_s2, arc_dir_e2, arc_dir_m2; 
  if (arc1->polarity == SAME_POL) {
    arc_dir_s1 = (*arc1->lines.begin())->dir.rotate(PI_2);
    arc_dir_e1 = (*arc1->lines.rbegin())->dir.rotate(PI_2);
    arc_dir_m1 = line1.dir.rotate(PI_2);
    arc_dir_s2 = (*arc2->lines.rbegin())->dir.rotate(-PI_2);
    arc_dir_e2 = (*arc2->lines.begin())->dir.rotate(-PI_2);
    arc_dir_m2 = line2.dir.rotate(-PI_2);
  } else {
    arc_dir_s1 = (*arc1->lines.begin())->dir.rotate(-PI_2);
    arc_dir_e1 = (*arc1->lines.rbegin())->dir.rotate(-PI_2);
    arc_dir_m1 = line1.dir.rotate(-PI_2);
    arc_dir_s2 = (*arc2->lines.rbegin())->dir.rotate(PI_2);
    arc_dir_e2 = (*arc2->lines.begin())->dir.rotate(PI_2);
    arc_dir_m2 = line2.dir.rotate(PI_2);
  }
  Vectord test_vec1 = (*arc2->lines.begin())->sp - (*arc1->lines.begin())->sp;
  Vectord test_vec2 = (*arc2->lines.rbegin())->ep - (*arc1->lines.rbegin())->ep;
  Vectord test_vec3 = (test_vec1 + test_vec2) / 2.0;
  
  double t1, t2, t3, t4, t5, t6;
  t1 = arc_dir_s1.dot(test_vec1);
  t2 = arc_dir_e1.dot(test_vec2);
  t3 = arc_dir_m1.dot(test_vec3);
  t4 = -arc_dir_e2.dot(test_vec1);
  t5 = -arc_dir_s2.dot(test_vec2);
  t6 = -arc_dir_m2.dot(test_vec3);
  return t1 >= REGION_LIMITATION_DIS_TOLERACE &&
         t2 >= REGION_LIMITATION_DIS_TOLERACE &&
         t3 >= REGION_LIMITATION_DIS_TOLERACE &&
         t4 >= REGION_LIMITATION_DIS_TOLERACE &&
         t5 >= REGION_LIMITATION_DIS_TOLERACE &&
         t6 >= REGION_LIMITATION_DIS_TOLERACE;
}


bool regionLimitation(const std::shared_ptr<Arc> &arc1, const std::shared_ptr<Arc> &arc2) {
  if (!arc1 || !arc2) {
    return false;
  }
  Lined line1((*arc1->lines.begin())->sp, (*arc1->lines.rbegin())->ep, 0);
  Lined line2((*arc2->lines.begin())->sp, (*arc2->lines.rbegin())->ep, 0);
  Vectord arc_dir_s1, arc_dir_e1, arc_dir_m1;
  Vectord arc_dir_s2, arc_dir_e2, arc_dir_m2; 
  if (arc1->polarity == SAME_POL) {
    arc_dir_s1 = (*arc1->lines.begin())->dir.rotate(PI_2);
    arc_dir_e1 = (*arc1->lines.rbegin())->dir.rotate(PI_2);
    arc_dir_m1 = line1.dir.rotate(PI_2);
    arc_dir_s2 = (*arc2->lines.begin())->dir.rotate(PI_2);
    arc_dir_e2 = (*arc2->lines.rbegin())->dir.rotate(PI_2);
    arc_dir_m2 = line2.dir.rotate(PI_2);
  } else {
    arc_dir_s1 = (*arc1->lines.begin())->dir.rotate(-PI_2);
    arc_dir_e1 = (*arc1->lines.rbegin())->dir.rotate(-PI_2);
    arc_dir_m1 = line1.dir.rotate(-PI_2);
    arc_dir_s2 = (*arc2->lines.begin())->dir.rotate(-PI_2);
    arc_dir_e2 = (*arc2->lines.rbegin())->dir.rotate(-PI_2);
    arc_dir_m2 = line2.dir.rotate(-PI_2);
  }
  Vectord test_vec1 = (*arc2->lines.rbegin())->ep - (*arc1->lines.begin())->sp;
  Vectord test_vec2 = (*arc2->lines.begin())->sp - (*arc1->lines.rbegin())->ep;
  Vectord test_vec3 = (test_vec1 + test_vec2) / 2.0;
  
  double t1, t2, t3, t4, t5, t6;
  t1 = arc_dir_s1.dot(test_vec1);
  t2 = arc_dir_e1.dot(test_vec2);
  t3 = arc_dir_m1.dot(test_vec3) * arc1->polarity * arc2->polarity;
  t4 = -arc_dir_e2.dot(test_vec1);
  t5 = -arc_dir_s2.dot(test_vec2);
  t6 = -arc_dir_m2.dot(test_vec3) * arc1->polarity * arc2->polarity;
  
  return (t1 >= REGION_LIMITATION_DIS_TOLERACE &&
         t2 >= REGION_LIMITATION_DIS_TOLERACE &&
         t3 >= REGION_LIMITATION_DIS_TOLERACE &&
         t4 >= REGION_LIMITATION_DIS_TOLERACE &&
         t5 >= REGION_LIMITATION_DIS_TOLERACE &&
         t6 >= REGION_LIMITATION_DIS_TOLERACE) ||
         regionLimitation_reverse(arc1, arc2);
}

static void gaussianKernel(std::vector<double> &kernel, double sigma, double mean) {
  double sum = 0.0;
  for (int i = 0; i < (int)kernel.size(); ++i) {
    double val = ((double)i - mean) / sigma;
    kernel[i] = exp(-0.5 * val * val);
    sum += kernel[i];
  }
  if (sum >= 0.0) {
    for (int i = 0; i < (int)kernel.size(); ++i) {
      kernel[i] /= sum;
    }
  }
}

bool gaussianSampler(const uint8_t *ori_data, int ori_row, int ori_col,
                     double *data, int row, int col,
                     double scale, double sigma_scale) {
  
  double *aux = new double[ori_row * col];
  double sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
  double prec = 3.0;
  int h = (int)ceil(sigma * sqrt(2.0 * prec * log(10.0)));
  int n = 1 + 2 * h;
  std::vector<double> kernel(n);

  for (int idy = 0; idy < col; ++idy) {
    double yy = (double)idy / scale;
    int yc = (int)floor(yy + 0.5);
    gaussianKernel(kernel, sigma, (double)h + yy - (double)yc);
    for (int idx = 0; idx < ori_row; ++idx) {
      double sum = 0.0;
      for (int dim = 0; dim < n; ++dim) {
        int j = yc - h + dim;
        while (j < 0) {
          j += 2 * ori_col;
        }
        while (j >= 2 * ori_col) {
          j -= 2 * ori_col;
        }
        if (j >= ori_col) {
          j = 2 * ori_col - 1 - j;
        }
        sum += (double)ori_data[idx * ori_col + j] * kernel[dim];
      }
      aux[idx * col + idy] = sum;
    }
  }

  for (int idx = 0; idx < row; ++idx) {
    double xx = (double)idx / scale;
    int xc = (int)floor(xx + 0.5);
    gaussianKernel(kernel, sigma, (double)h + xx - (double)xc);
    for (int idy = 0; idy < col; ++idy) {
      double sum = 0.0;
      for (int dim = 0; dim < n; ++dim) {
        int j = xc - h + dim;
        while (j < 0) {
          j += 2 * ori_row;
        }
        while (j >= 2 * ori_row) {
          j -= 2 * ori_row;
        }
        if (j >= ori_row) {
          j = 2 * ori_row - 1 - j;
        }
        sum += aux[j * col + idy] * kernel[dim];
      }
      data[idx * col + idy] = sum;
    }
  }
  delete [] aux;
  return true;
}


} // zgh