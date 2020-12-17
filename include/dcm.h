#ifndef DCM_H
#define DCM_H

// Compiler instructions - option for lightwight trig implementation?
#include <cmath>

#include "matrix.h"

namespace attitude
{
    template<typename T>
    class dcm_
    {
    private:
        linalg::square<T> _R = linalg::square<T>(3);

        dcm_() : _R{ 3, 1, 0, 0,
                        0, 1, 0,
                        0, 0, 1 } {}

        T* _get(int i) { return _R[i]; }

    public:
        dcm_(linalg::square<T> R) : _R(R) {}

        dcm_(T _1, T _2, T _3,
             T _4, T _5, T _6,
             T _7, T _8, T _9 ) : _R{ 3, _1, _2, _3, _4, _5, _6, _7, _8, _9 }
        {}

        dcm_<T> reverse() { return _R.transpose(); }

        linalg::square<T> matrix() { return _R; }

        dcm_<T> operator+(dcm_<T> R)
        {
            return dcm_<T>(_R * R.matrix());
        }

        dcm_<T> operator-(dcm_<T> R)
        {
            dcm_<T> R3 = operator+(R.reverse());
            return R3;
        }

        dcm_<T>* operator+=(dcm_<T> R)
        {
            _R *= R.matrix();
            return this;
        }

        dcm_<T>* operator-=(dcm_<T> R)
        {
            _R += R.matrix().transpose();
            return this;
        }

        bool operator==(dcm_<T> R)
        {
            linalg::square<T> diff = _R - R.matrix();
            return diff == T(0);
        }

        bool operator==(linalg::matrix<T> R) { return _R == R; } // TODO: remove, and add explicit flags to avoid implicit conversion from square to dcm

        T* operator[](int i){ return _get(i); }
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

    template<typename T>
    dcm_<T> AXIS(int axis, T theta)
    {
        switch (axis) {
            case 1:
                return R1(theta);
            case 2:
                return R2(theta);
            case 3:
                return R3(theta);
            default:
                return ZERO<T>();
        }
    }
}

template<typename T>
void display(attitude::dcm_<T> R)
{
    attitude::linalg::square<T> _R = R.matrix();

    for (int i = 0; i < 3; i++)
    {
        T* row = _R[i];

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