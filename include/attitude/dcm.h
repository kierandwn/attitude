// Copyright (c) 2020 Kieran Downie. All rights reserved.
//
// This file is part of attitude.
//
// attitude is free software : you can redistribute it and /
// or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License,
// or (at your option) any later version.
//
// attitude is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with attitude.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef ATT_DCM_H_
#define ATT_DCM_H_

#include <cmath> // TODO: Compiler instructions - option for lightwight trig implementation?

#include "attitude/matrix.h"

namespace attitude {
namespace dcm {

// ZERO (function)
// Returns rotation matrix for a zero rotation (identity).
template <typename Tp>
mn_matrix<Tp, 3, 3> ZERO() {
  return mn_matrix<Tp, 3, 3>{1., 0., 0., 
                          0., 1., 0.,
                          0., 0., 1.};
}


// R1 (function)
// Returns rotation matrix for a theta radian rotation about the first axis.
template <typename Tp>
mn_matrix<Tp, 3, 3> R1(Tp theta) {
  return mn_matrix<Tp, 3, 3>{1.,          0.,         0.,
                          0.,  cos(theta), sin(theta),
                          0., -sin(theta), cos(theta)};
}


// R2 (function)
// Returns rotation matrix for a theta radian rotation about the second axis.
template <typename Tp>
mn_matrix<Tp, 3, 3> R2(Tp theta) {
  return mn_matrix<Tp, 3, 3>{cos(theta), 0., -sin(theta), 
                                  0., 1.,          0., 
                          sin(theta), 0.,  cos(theta)};
}


// R3 (function)
// Returns rotation matrix for a theta radian rotation about the third axis.
template <typename Tp>
mn_matrix<Tp, 3, 3> R3(Tp theta) {
  return mn_matrix<Tp, 3, 3>{ cos(theta), sin(theta), 0., 
                          -sin(theta), cos(theta), 0.,
                                   0.,         0., 1.};
}


// AXIS (function)
// Returns rotation matrix for a theta radian rotation about the specified axis.
template <typename Tp>
mn_matrix<Tp, 3, 3> AXIS(int axis, Tp theta) {
  switch (axis) {
    case 1:
      return R1(theta);
    case 2:
      return R2(theta);
    case 3:
      return R3(theta);
    default:
      return ZERO<Tp>();
  }
}

} // namespace dcm
} // namespace attitude
#endif // ATT_DCM_H_