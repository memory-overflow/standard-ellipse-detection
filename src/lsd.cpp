/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-24
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <vector>
#include <algorithm>
#include "types.hpp"
#include "detect.h"
#include "unitily.h"

namespace zgh {
  

static bool ll_angle(const double *data, double *angles, double *mod, std::vector<Pixel> &sort_points,
                     int nbins, int row, int col) {
  double max_grad = 0.0;
  for (int idx = 0; idx < row; ++idx) {
    angles[idx * col + col - 1] = ANGLE_NOT_DEF;
  }
  for (int idy = 0; idy < col; ++idy) {
    angles[(row - 1) * col + idy] = ANGLE_NOT_DEF;
  }
  
  for (int idy = 0; idy < col - 1; ++idy) {
    for (int idx = 0; idx < row - 1; ++idx) {
      int adr = idx * col + idy;
      /*
       * Norm 2 computation using 2X2 pixel window:
       *   A B
       *   C D
       * and
       *   com1 = D - A,  com2 = B - C.
       * Then
       *   gx = B + D - (A + C)   horizontal difference
       *   gy = C + D - (A + B)   vertical difference
       * com1 and com2 are just to avoid 2 additions
       */
      double com1, com2;
      com1 = data[adr + col + 1] - data[adr];
      com2 = data[adr + 1] - data[adr + col];
      double dx = com1 + com2, dy = com1 - com2;
      double norm = sqrt(1.0 * (dx * dx + dy * dy) / 4.0);
      mod[adr] = norm;
      if (norm <= GRAD_THRESHOLD) {
        angles[adr] = ANGLE_NOT_DEF;
      } else {
        angles[adr] = atan2(dx, -dy);
        sort_points.push_back(Pixel(idx, idy));
        if (norm > max_grad) {
          max_grad = norm;
        }
      }
    }
  }

  std::vector<int> bucket_ranks(nbins, 0);
  for (auto &pix : sort_points) {
    int adr = pix.x * col + pix.y;
    int bucket_id = (int)floor(mod[adr] * (double)nbins / max_grad);
    if (bucket_id >= nbins) {
      bucket_id = nbins - 1;
    }
    ++bucket_ranks[bucket_id];
  }
  for (int bucket_id = nbins - 2; bucket_id >= 0; --bucket_id) {
    bucket_ranks[bucket_id] += bucket_ranks[bucket_id + 1];
  }

  std::vector<Pixel> temp_points(sort_points.size());
  for (int id = (int)sort_points.size() - 1; id >= 0; --id) {
    auto &pix = sort_points[id];
    int adr = pix.x * col + pix.y;
    int bucket_id = (int)floor(mod[adr] * (double)nbins / max_grad);
    if (bucket_id >= nbins) {
      bucket_id = nbins - 1;
    }
    temp_points[--bucket_ranks[bucket_id]] = std::move(pix);
  }
  sort_points.clear();
  sort_points = std::move(temp_points);
  return true;
}

static bool isaligned(double angle, double line_angle, double prec) {
  if (equal(angle, ANGLE_NOT_DEF)) {
    return false;
  }
  line_angle -= angle;
  if (line_angle < 0) {
    line_angle = -line_angle;
  }
  if (line_angle > 3.0 * PI_2) {
    line_angle = 2.0 * PI - line_angle;
    if(line_angle < 0) {
      line_angle = -line_angle;
    }
  }
  return line_angle <= prec;
}

static std::shared_ptr<Lined> regionGrow(Pixel sp, const double *angles, bool *used,
                                         double prec, int row, int col) {
  if (!inRect(row, col, sp.x, sp.y)) {
    return nullptr;
  }
  std::shared_ptr<Lined> line = std::make_shared<Lined>();
  line->addpixel(sp);
  int addr = sp.x * col + sp.y;
  double line_angle = angles[addr];
  used[addr] = USED;
  double sumdx = sin(line_angle);
  double sumdy = cos(line_angle);

  // bfs search
  auto &regs = line->getregs();
  for (int i = 0; i < (int)regs.size(); ++i) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        int x = regs[i].x + dx;
        int y = regs[i].y + dy;
        int addr = x * col + y;
        if (inRect(row, col, x, y) && used[addr] == NOTUSED &&
            isaligned(angles[addr], line_angle, prec)) {
          used[addr] = USED;
          line->addpixel(Pixel(x, y));
          sumdx += sin(angles[addr]);
          sumdy += cos(angles[addr]);
          line_angle = atan2(sumdx, sumdy);
        }
      }
    }
  }
  line->dir = Pointd(sin(line_angle), cos(line_angle));
  return line;
}

static double getTheta(const std::vector<Pixel> &regs, double x, double y, const double *mod,
                     double line_angle, double prec, int col) {
  
  double Ixx = 0.0;
  double Iyy = 0.0;
  double Ixy = 0.0;
  for (auto &pix : regs) {
    int addr = pix.x * col + pix.y;
    Ixx += sqr((double)pix.x - x) * mod[addr];
    Iyy += sqr((double)pix.y - y) * mod[addr];
    Ixy -= ((double)pix.x - x) * ((double)pix.y - y) * mod[addr];
  }
  double lambda = 0.5 * (Ixx + Iyy - sqrt(sqr(Ixx - Iyy) + 4.0 * sqr(Ixy)));
  double theta = fabs(Ixx) > fabs(Iyy) ? atan2(lambda - Ixx, Ixy) : atan2(Ixy, lambda - Iyy);
  double temp1 = angle_diff(theta, line_angle);
  if (temp1 > prec) {
    double temp2 = angle_diff(theta + PI, line_angle);
    if (temp2 < prec) {
      theta += PI;
      if (theta > PI) {
        theta -= 2.0 * PI;
      }
    } else {
      theta = (temp2 < temp1) ? (theta + PI) : theta;
      while (theta <= -PI) {
        theta += 2.0 * PI;
      }
      while (theta > PI) {
        theta -= 2.0 * PI;
      }
    }
  }
  return theta;
}

static void region2rect(std::shared_ptr<Lined> &line, const double *mod, int col) {
  double x = 0.0, y = 0.0, sum = 0.0;
  for (auto &pix : line->getregs()) {
    int addr = pix.x * col + pix.y;
    x += (double)pix.x * mod[addr];
    y += (double)pix.y * mod[addr];
    sum += mod[addr];
  }
  x /= sum;
  y /= sum;
  double line_angle = atan2(line->dir.x, line->dir.y);
  double theta = getTheta(line->getregs(), x, y, mod, line_angle, angle2rad(ANGLE_TH), col);
  double dx = sin(theta);
  double dy = cos(theta);
  double lmin = 0.0;
  double lmax = 0.0;
  double wmin = 0.0;
  double wmax = 0.0;
  for (auto &pix : line->getregs()) {
    double l = ((double)pix.x - x) * dx + ((double)pix.y - y) * dy;
    double w = -((double)pix.x - x) * dy + ((double)pix.y - y) * dx;
    lmax = max(lmax, l);
    lmin = min(lmin, l);
    wmax = max(wmax, w);
    wmin = min(wmin, w);
  }
  line->sp.x = x + lmin * dx;
  line->sp.y = y + lmin * dy;
  line->ep.x = x + lmax * dx;
  line->ep.y = y + lmax * dy;
  line->length = (line->ep - line->sp).length();
  line->width = wmax - wmin;
  if (line->width < 1.0) {
    line->width = 1.0;
  }
  line->center.x = x;
  line->center.y = y;
  line->dir.x = dx;
  line->dir.y = dy;

}


static bool isArcSegment(std::shared_ptr<Lined> &line, const double *angles, int8_t *pol, int col) {

  int same_pol_cnt = 0, opp_pol_cnt = 0;
  for (auto &pix : line->getregs()) {
    int addr = pix.x * col + pix.y;
    if (pol[addr] == (uint8_t)SAME_POL) {
      ++same_pol_cnt;
    } else if (pol[addr] == (uint8_t)OPP_POL) {
      ++opp_pol_cnt;
    }
  }
  if ((same_pol_cnt + opp_pol_cnt) > (int)line->getregs().size() / 2) {
    if (same_pol_cnt > opp_pol_cnt) {
      line->polarity = SAME_POL;
    } else {
      line->polarity = OPP_POL;
    }
    return true;
  }

  double angle_up = 0.0, angle_down = 0.0, angle_main;
  double reg_up_sin = 0.0, reg_up_cos = 0.0;
  double reg_down_sin = 0.0, reg_down_cos = 0.0;

  for (auto &pix : line->getregs()) {
    int addr = pix.x * col + pix.y;
    Pointd p(pix.x, pix.y);
    if (line->dir.dot(p - line->center) >= 0) {
      reg_up_sin += sin(angles[addr]);
      reg_up_cos += cos(angles[addr]);
    } else {
      reg_down_sin += sin(angles[addr]);
      reg_down_cos += cos(angles[addr]);
    }
  }
  angle_up = atan2(reg_up_sin, reg_up_cos);
  angle_down = atan2(reg_down_sin, reg_down_cos);
  angle_main = atan2(reg_up_sin + reg_down_sin, reg_up_cos + reg_down_cos);

  double temp1 = angle_diff_signed(angle_up, angle_main);
  double temp2 = angle_diff_signed(angle_down, angle_main);
  // 实验结果最好的阈值
  if (temp1 >= PI_8 / 10.0 && temp2 <= -PI_8 / 10.0) {
    line->polarity = OPP_POL;
    for (auto &pix : line->getregs()) {
      pol[pix.x * col + pix.y] = (uint8_t)OPP_POL;
    }
  } else if (temp1 <= -PI_8 / 10.0 && temp2 >= PI_8 / 10.0) {
    line->polarity = SAME_POL;
    for (auto &pix : line->getregs()) {
      pol[pix.x * col + pix.y] = (uint8_t)SAME_POL;
    }
  } else {
    // not a Arc support segment
    return false;
  }
  return true;
}

static bool reduce_region_radius(std::shared_ptr<Lined> &line, const double *mod, bool *used,
                                 double density_th, int col) {

  double density = (double)line->getregs().size() / (line->length * line->width);
  if (density >= density_th) {
    return true;
  }
  Pointd p(line->getregs()[0].x, line->getregs()[0].y);
  double rad1 = (line->sp - p).length();
  double rad2 = (line->ep - p).length();
  double rad = max(rad1, rad2);
  while (density < density_th) {
    auto &regs = line->getregs();
    int cnt = regs.size();
    rad *= 0.75;
    std::vector<Pixel> temppix;
    for (int id = 0; id < cnt; ++id) {
      auto pix = regs[id];
      if ((pix - regs[0]).length() > rad) {
        used[pix.x * col + pix.y] = NOTUSED;
      } else {
        temppix.push_back(pix);
      }
    }
    if ((int)temppix.size() < 2) {
      return false;
    }
    line->setregs(std::move(temppix));
    region2rect(line, mod, col);
    density = (double)line->getregs().size() / (line->length * line->width);
  }
  return true;
}

static bool refine(std::shared_ptr<Lined> &line, const double *angles, const double *mod,
                   bool *used, double prec, double density_th, int row, int col) {
  double density = (double)line->getregs().size() / (line->length * line->width);
  if (density >= density_th) {
    return true;
  }
  auto &regs = line->getregs();
  int xc = regs[0].x;
  int yc = regs[0].y;
  double sum = 0.0;
  double s_sum = 0.0;
  int n = 0;
  for (auto &pix : regs) {
    int addr = pix.x * col + pix.y;
    used[addr] = NOTUSED;
    if ((regs[0] - pix).length() < line->width) {
      double ang_d = angle_diff_signed(angles[addr], angles[xc * col + yc]);
      sum += ang_d;
      s_sum += ang_d * ang_d;
      ++n;
    }
  }

  double mean_angle = sum / (double)n;
  //  以 2 倍标准差作为新的角度容忍度，最开始为 22.5° * pi / 180
  double tau = 2.0 * sqrt((s_sum - 2.0 * mean_angle * sum) / (double)n + mean_angle * mean_angle);
  int pol = line->polarity;
  line = regionGrow(regs[0], angles, used, tau, row, col);
  line->polarity = pol;
  if (!line || line->getregs().size() < 2) {
    return false;
  }
  region2rect(line, mod, col);
  density = (double)line->getregs().size() / (line->length * line->width);
  
  if (density < density_th) {
    if (reduce_region_radius(line, mod, used, density_th, col)) {
      return true;
    } else {
      return false;
    }
  } else {
    return true;
  }
  return true;
}

static double rect_nfa(const std::shared_ptr<Lined> &line, const double *angles, double logNT,
                       double prec, int row, int col) {
  int pts = 0;
  int alg = 0;
  double line_angle = atan2(line->dir.x, line->dir.y);
  for (RectIter iter(line); !iter.isEnd(); ++iter) {
    if (inRect(row, col, iter.np.x, iter.np.y)) {
      ++pts;
      if (isaligned(angles[iter.np.x * col + iter.np.y], line_angle, prec)) {
        ++alg;
      }
    }
  }
  return nfa(pts, alg, prec / PI, logNT);
}

static double rect_improve(std::shared_ptr<Lined> &line, const double *angles, double logNT,
                            double log_eps, int row, int col) {
  double delta = 0.5;
  double delta_2 = delta / 2.0;
  double prec = angle2rad(ANGLE_TH);
  double log_nfa = rect_nfa(line, angles, logNT, prec, row, col);

  if (log_nfa > log_eps) {
    return log_nfa;
  }

  for (int i = 0; i < 5; ++i) {
    prec /= 2.0;
    double log_nfa_new = rect_nfa(line, angles, logNT, prec, row, col);
    if (log_nfa_new > log_nfa) {
      log_nfa = log_nfa_new;
    }
  }
  if (log_nfa > log_eps) {
    return log_nfa;
  }

  double min_width = line->width;
  for (int i = 0; i < 5; ++i) {
    if (line->width - delta >= 0.5) {
      line->width -= delta;
      double log_nfa_new = rect_nfa(line, angles, logNT, prec, row, col);
      if(log_nfa_new > log_nfa) {
        log_nfa = log_nfa_new;
        min_width = line->width;
      }
    }
  }
  line->width = min_width;
  if (log_nfa > log_eps) {
    return log_nfa;
  }


  int idx = 0;
  for (int i = 1; i <= 5; ++i) {
    if (line->width - i * delta >= 0.5) {
      line->sp.x += -(double)i * line->dir.y * delta_2;
      line->sp.y +=  (double)i * line->dir.x * delta_2;
      line->ep.x += -(double)i * line->dir.y * delta_2;
      line->ep.y +=  (double)i * line->dir.x * delta_2;
      line->width -= (double)i * delta;
      double log_nfa_new = rect_nfa(line, angles, logNT, prec, row, col);
      if(log_nfa_new > log_nfa) {
        log_nfa = log_nfa_new;
        idx = i;
      }
      line->sp.x -= -(double)i * line->dir.y * delta_2;
      line->sp.y -=  (double)i * line->dir.x * delta_2;
      line->ep.x -= -(double)i * line->dir.y * delta_2;
      line->ep.y -=  (double)i * line->dir.x * delta_2;
      line->width += (double)i * delta;;
    }
  }
  line->sp.x += -(double)idx * line->dir.y * delta_2;
  line->sp.y +=  (double)idx * line->dir.x * delta_2;
  line->ep.x += -(double)idx * line->dir.y * delta_2;
  line->ep.y +=  (double)idx * line->dir.x * delta_2;
  line->width -= (double)idx * delta;
  if (log_nfa > log_eps) {
    return log_nfa;
  }


  idx = 0;
  for (int i = 1; i <= 5; ++i) {
    if (line->width - i * delta >= 0.5) {
      line->sp.x -= -(double)i * line->dir.y * delta_2;
      line->sp.y -=  (double)i * line->dir.x * delta_2;
      line->ep.x -= -(double)i * line->dir.y * delta_2;
      line->ep.y -=  (double)i * line->dir.x * delta_2;
      line->width -= (double)i * delta;
      double log_nfa_new = rect_nfa(line, angles, logNT, prec, row, col);
      if(log_nfa_new > log_nfa) {
        log_nfa = log_nfa_new;
        idx = i;
      }
      line->sp.x += -(double)i * line->dir.y * delta_2;
      line->sp.y +=  (double)i * line->dir.x * delta_2;
      line->ep.x += -(double)i * line->dir.y * delta_2;
      line->ep.y +=  (double)i * line->dir.x * delta_2;
      line->width += (double)i * delta;;
    }
  }
  line->sp.x -= -(double)idx * line->dir.y * delta_2;
  line->sp.y -=  (double)idx * line->dir.x * delta_2;
  line->ep.x -= -(double)idx * line->dir.y * delta_2;
  line->ep.y -=  (double)idx * line->dir.x * delta_2;
  line->width -= (double)idx * delta;
  if (log_nfa > log_eps) {
    return log_nfa;
  }

  for (int i = 0; i < 5; ++i) {
    prec /= 2.0;
    double log_nfa_new = rect_nfa(line, angles, logNT, prec, row, col);
    if (log_nfa_new > log_nfa) {
      log_nfa = log_nfa_new;
    }
  }
  return log_nfa;
}

bool lineSegmentDetection(const double *image, int row, int col, std::vector<std::shared_ptr<Lined> > &lines) {
  int nbins = 1024;
  double logNT = 5.0 * (log10((double)row) + log10((double)col)) / 2.0 + log10(11.0);
  double min_reg_size = (int)(-logNT / log10(ANGLE_TH / 180.0)); // 每个矩形区域内 align point 最小数量
  
  double *angles = new double[row * col];
  double *mod = new double[row * col];
  bool *used = new bool[row * col];
  int8_t *pol = new int8_t[row * col];
  std::vector<Pixel> sort_points;

  ll_angle(image, angles, mod, sort_points, nbins, row, col);
  memset(used, NOTUSED, sizeof(bool) * row * col);
  memset(pol, 0, sizeof(int8_t) * row * col);

  for (auto &pix : sort_points) {
    int addr = pix.x * col + pix.y;
    if (used[addr] == NOTUSED && !equal(angles[addr], ANGLE_NOT_DEF)) {
      std::shared_ptr<Lined> line = regionGrow(pix, angles, used, angle2rad(ANGLE_TH), row, col);

      if (!line || line->getregs().size() < min_reg_size) {
        continue;
      }
      region2rect(line, mod, col);
      
      if (!isArcSegment(line, angles, pol, col)) {
        continue;
      }

      // 提纯，通过重新生长区域来达到期望的密度阈值
      if (!refine(line, angles, mod, used, angle2rad(ANGLE_TH), 0.7, row, col)) {
        continue;
      }
      
      double log_eps = 0.0;
      double log_nfa = rect_improve(line, angles, logNT, log_eps, row, col);
      if (log_nfa <= log_eps) {
        continue;
      }

      line->sp += Pointd(0.5, 0.5);
      line->ep += Pointd(0.5, 0.5);
      lines.push_back(std::move(line));
    }
  }
  delete [] angles;
  delete [] mod;
  delete [] used;
  delete [] pol;
  return true;
}

} // namespace zgh