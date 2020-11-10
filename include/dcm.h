#ifndef DCM_H
#define DCM_H

#include "common.h"

// Compiler instructions - option for lightwight trig implementation?
#include <cmath>

namespace attitude
{
    template<typename T>
    class dcm_
    {
    private:
        math::matrix9<T> _R;

        dcm_() : _R{ 1, 0, 0,
                     0, 1, 0,
                     0, 0, 1 }
        {}

    public:
        dcm_(math::matrix9<T> R) : _R(R)
        {}

        dcm_(   T _1, T _2, T _3,
                T _4, T _5, T _6,
                T _7, T _8, T _9 ) : _R{ _1, _2, _3, _4, _5, _6, _7, _8, _9 }
        {}

        dcm_<T> get_reverse()
        {
            return _R.transpose();
        }

        math::matrix9<T> get_matrix()
        {
            return _R;
        }

        dcm_<T> operator+(dcm_<T> R)
        {
            return dcm_<T>(_R * R.get_matrix());
        }

        dcm_<T> operator-(dcm_<T> R)
        {
            dcm_<T> R3 = operator+(R.get_reverse());
            return R3;
        }

        bool operator==(dcm_<T> R)
        {
            math::matrix9<T> diff = _R - R.get_matrix();
            return diff == T(0);
        }

        bool operator==(math::matrix9<T> R)
        {
            return _R == R;
        }

        std::array<T, 3> operator[](int i)
        {
            return _R[i];
        }
    };

    template<typename T>
    dcm_<T> ZERO()
    {
        return dcm_<T>(T(1.), T(0.), T(0.),
                       T(0.), T(1.), T(0.),
                       T(0.), T(0.), T(1.) );
    }

    template<typename T>
    dcm_<T> R1(T theta)
    {
        return dcm_<T>( T(1.), T(0.), T(0.),
                        T(0.),  cos(theta), sin(theta),
                        T(0.), -sin(theta), cos(theta) );
    }

    template<typename T>
    dcm_<T> R2(T theta)
    {
        return dcm_<T>( cos(theta), T(0.), -sin(theta),
                        T(0.), T(1.), T(0.),
                        sin(theta), T(0.),  cos(theta) );
    }

    template<typename T>
    dcm_<T> R3(T theta)
    {
        return dcm_<T>(  cos(theta), sin(theta), T(0.),
                        -sin(theta), cos(theta), T(0.),
                         T(0.), T(0.), T(1.) );
    }
}

template<typename T>
void display(attitude::dcm_<T> R)
{
    attitude::math::matrix9<T> _R = R.get_matrix();

    for (int i = 0; i < 3; i++)
    {
        std::array<T, 3> row = _R[i];

        std::cout
            << "[ "
            << row[0] << ","
            << row[1] << ","
            << row[2]
            << " ]"
            << std::endl;
    }
}

#endif // DCM_H