#ifndef EULER_H
#define EULER_H

#include <cmath>

#include "dcm.h"

namespace attitude
{
    template<typename T>
    int argmax(T* arr, int length)
    {
        int max_index = 0;

        for (int k = 1; k < length; k++)
            if (arr[k] > arr[max_index]){ max_index = k; }

        return max_index;
    }

    template<typename T>
    class quaternion_
    {
    private:
        T _q[4] = { T(1.), T(0.), T(0.), T(0.) };
        dcm_<T> _R;

        T _get(int i) { return _q[i]; }
        void _set(int i, T val) { _q[i] = val; }


    public:
        quaternion_() : _q{ T(1.), T(0.), T(0.), T(0.) },
                        _R(to_dcm(*this)) {}

        quaternion_(T q0, T q1, T q2, T q3) : _q{ q0, q1, q2, q3 },
                                              _R(to_dcm(*this)) {}

        quaternion_(dcm_<T> R) : _R(R) { from_dcm(R); }

        void from_dcm(dcm_<T> R)
        {
            // Shepherd's Rule
            T _4q0q0 = 1 + R[0][0] + R[1][1] + R[2][2];
            T _4q1q1 = 1 + R[0][0] - R[1][1] - R[2][2];
            T _4q2q2 = 1 - R[0][0] + R[1][1] - R[2][2];
            T _4q3q3 = 1 - R[0][0] - R[1][1] + R[2][2];

            T vals[4]{ abs(_4q0q0), abs(_4q1q1), abs(_4q2q2), abs(_4q3q3) };
            int i = argmax(vals, 4);

            if (i == 0)
            {
                T _4q0q1 = R[1][2] - R[2][1];
                T _4q0q2 = R[2][0] - R[0][2];
                T _4q0q3 = R[0][1] - R[1][0];
                T _q0 = sqrt(0.25 * _4q0q0);

                _set(0, _q0);
                _set(1, 0.25 * _4q0q1 / _q0);
                _set(2, 0.25 * _4q0q2 / _q0);
                _set(3, 0.25 * _4q0q3 / _q0);
            }
            else if (i == 1)
            {
                T _4q0q1 = R[1][2] - R[2][1];
                T _4q1q2 = R[0][1] + R[1][0];
                T _4q1q3 = R[2][0] + R[0][2];
                T _q1 = sqrt(0.25 * _4q1q1);

                _set(0, 0.25 * _4q0q1 / _q1);
                _set(1, _q1);
                _set(2, 0.25 * _4q1q2 / _q1);
                _set(3, 0.25 * _4q1q3 / _q1);
            }
            else if (i == 2)
            {
                T _4q0q2 = R[2][0] - R[0][2];
                T _4q1q2 = R[0][1] + R[1][0];
                T _4q2q3 = R[1][2] + R[2][1];
                T _q2 = sqrt(0.25 * _4q2q2);

                _set(0, 0.25 * _4q0q2 / _q2);
                _set(1, 0.25 * _4q1q2 / _q2);
                _set(2, _q2);
                _set(3, 0.25 * _4q2q3 / _q2);
            }
            else if (i == 3)
            {
                T _4q0q3 = R[0][1] - R[1][0];
                T _4q1q3 = R[2][0] + R[0][2];
                T _4q2q3 = R[1][2] + R[2][1];
                T _q3 = sqrt(0.25 * _4q3q3);

                _set(0, 0.25 * _4q0q3 / _q3);
                _set(1, 0.25 * _4q1q3 / _q3);
                _set(2, 0.25 * _4q2q3 / _q3);
                _set(3, _q3);
            }
        }

        dcm_<T> dcm() { return _R; }
        quaternion_<T> reverse() { return quaternion_<T>( dcm().reverse() ); }
        
        quaternion_<T> operator+(quaternion_<T> q)
        {
            return quaternion_<T>( 
                q[0] * _get(0) + -1. * q[1] * _get(1) + -1. * q[2] * _get(2) + -1. * q[3] * _get(3),
                q[1] * _get(0) +       q[0] * _get(1) +       q[3] * _get(2) +       q[2] * _get(3),
                q[2] * _get(0) + -1. * q[3] * _get(1) +       q[0] * _get(2) +       q[1] * _get(3),
                q[3] * _get(0) +       q[2] * _get(1) + -1. * q[1] * _get(2) +       q[0] * _get(3) );
        }

        quaternion_<T> operator-(quaternion_<T> q)
        {
            return quaternion_<T>( 
                      q[0] * _get(0) +       q[1] * _get(1) +       q[2] * _get(2) +       q[3] * _get(3),
                -1. * q[1] * _get(0) +       q[0] * _get(1) + -1. * q[3] * _get(2) +       q[2] * _get(3),
                -1. * q[2] * _get(0) +       q[3] * _get(1) +       q[0] * _get(2) + -1. * q[1] * _get(3),
                -1. * q[3] * _get(0) + -1. * q[2] * _get(1) +       q[1] * _get(2) +       q[0] * _get(3) );
        }

        quaternion_<T>* operator+=(quaternion_<T> q)
        {
            quaternion_( 
                q[0] * _get(0) + -1. * q[1] * _get(1) + -1. * q[2] * _get(2) + -1. * q[3] * _get(3),
                q[1] * _get(0) +       q[0] * _get(1) +       q[3] * _get(2) +       q[2] * _get(3),
                q[2] * _get(0) + -1. * q[3] * _get(1) +       q[0] * _get(2) +       q[1] * _get(3),
                q[3] * _get(0) +       q[2] * _get(1) + -1. * q[1] * _get(2) +       q[0] * _get(3) );
            return this;
        }

        quaternion_<T>* operator-=(quaternion_<T> q)
        {
            quaternion_( 
                      q[0] * _get(0) +       q[1] * _get(1) +       q[2] * _get(2) +       q[3] * _get(3),
                -1. * q[1] * _get(0) +       q[0] * _get(1) + -1. * q[3] * _get(2) +       q[2] * _get(3),
                -1. * q[2] * _get(0) +       q[3] * _get(1) +       q[0] * _get(2) + -1. * q[1] * _get(3),
                -1. * q[3] * _get(0) + -1. * q[2] * _get(1) +       q[1] * _get(2) +       q[0] * _get(3) );
            return this;
        }

        T operator[](int i) { return _get(i); } 

        // Comparison
        bool operator==(quaternion_<T> q) { return dcm() == q.dcm(); }
        bool operator==(dcm_<T> R) { return dcm() == R; }
    };

    // Conversions
    template<typename T>
    dcm_<T> to_dcm(quaternion_<T> q)
    {
        return dcm_<T>(
            pow(q[0], 2.) + pow(q[1], 2.) - pow(q[2], 2.) - pow(q[3], 2.), 2. * (q[1] * q[2] + q[0] * q[3]), 2. * (q[1] * q[3] - q[0] * q[2]),
            2. * (q[1] * q[2] - q[0] * q[3]), pow(q[0], 2.) - pow(q[1], 2.) + pow(q[2], 2.) - pow(q[3], 2.), 2. * (q[2] * q[3] + q[0] * q[1]),
            2. * (q[1] * q[3] + q[0] * q[2]), 2. * (q[2] * q[3] - q[0] * q[1]), pow(q[0], 2.) - pow(q[1], 2.) - pow(q[2], 2.) + pow(q[3], 2.)
        );
    }
}

template<typename T>
void display(attitude::quaternion_<T> q)
{
    std::cout
        << "[ "
        << q[0] << ","
        << q[1] << ","
        << q[2] << ","
        << q[3]
        << " ]"
        << std::endl;
}

#endif // EULER_H