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
// -
// shepherds.h
// Contains implementation of Shepherd's Rule for extracting quaternion components
// from rotation matrix. Used by quaternion and rodriguez parameter classes.
//
// File namespace attitude::
//
#ifndef ATT_SHEPHERDS_H_
#define ATT_SHEPHERDS_H_

#include "matrix.h"

namespace attitude {

template <typename Tp>
int argmax(Tp * arr, int length) {
  int max_index = 0;

  for (int k = 1; k < length; k++)
    if (arr[k] > arr[max_index]) {
      max_index = k;
    }

  return max_index;
}

template<typename Tp> 
Tp * shepherds_rule(matrix<Tp, 3, 3> matrix, Tp quaternion[4]) {
  Tp _4q0q0 = 1 + matrix[0][0] + matrix[1][1] + matrix[2][2];
  Tp _4q1q1 = 1 + matrix[0][0] - matrix[1][1] - matrix[2][2];
  Tp _4q2q2 = 1 - matrix[0][0] + matrix[1][1] - matrix[2][2];
  Tp _4q3q3 = 1 - matrix[0][0] - matrix[1][1] + matrix[2][2];

  Tp vals[4]{abs(_4q0q0), abs(_4q1q1), abs(_4q2q2), abs(_4q3q3)};
  int i = argmax(vals, 4);

  switch (i) {
    case 0: {
      Tp _4q0q1 = matrix[1][2] - matrix[2][1];
      Tp _4q0q2 = matrix[2][0] - matrix[0][2];
      Tp _4q0q3 = matrix[0][1] - matrix[1][0];
      Tp _q0 = sqrt(0.25 * _4q0q0);

      quaternion[0] = _q0;
      quaternion[1] = 0.25 * _4q0q1 / _q0;
      quaternion[2] = 0.25 * _4q0q2 / _q0;
      quaternion[3] = 0.25 * _4q0q3 / _q0;
    } break;

    case 1: {
      Tp _4q0q1 = matrix[1][2] - matrix[2][1];
      Tp _4q1q2 = matrix[0][1] + matrix[1][0];
      Tp _4q1q3 = matrix[2][0] + matrix[0][2];
      Tp _q1 = sqrt(0.25 * _4q1q1);

      quaternion[0] = 0.25 * _4q0q1 / _q1;
      quaternion[1] = _q1;
      quaternion[2] = 0.25 * _4q1q2 / _q1;
      quaternion[3] = 0.25 * _4q1q3 / _q1;
    } break;

    case 2: {
      Tp _4q0q2 = matrix[2][0] - matrix[0][2];
      Tp _4q1q2 = matrix[0][1] + matrix[1][0];
      Tp _4q2q3 = matrix[1][2] + matrix[2][1];
      Tp _q2 = sqrt(0.25 * _4q2q2);

      quaternion[0] = 0.25 * _4q0q2 / _q2;
      quaternion[1] = 0.25 * _4q1q2 / _q2;
      quaternion[2] = _q2;
      quaternion[3] = 0.25 * _4q2q3 / _q2;
    } break;

    case 3: {
      Tp _4q0q3 = matrix[0][1] - matrix[1][0];
      Tp _4q1q3 = matrix[2][0] + matrix[0][2];
      Tp _4q2q3 = matrix[1][2] + matrix[2][1];
      Tp _q3 = sqrt(0.25 * _4q3q3);

      quaternion[0] = 0.25 * _4q0q3 / _q3;
      quaternion[1] = 0.25 * _4q1q3 / _q3;
      quaternion[2] = 0.25 * _4q2q3 / _q3;
      quaternion[3] = _q3;
    } break;
  }
  return quaternion;
}


} // namespace atttitude;
