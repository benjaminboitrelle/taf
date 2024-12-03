#pragma once

#include <iostream>
#include <variant>
#include <type_traits>

template <typename T>
struct Limit
{
    T min;
    T max;

    Limit(T min_val, T max_val)
    {
        if (min_val > max_val)
        {
            std::swap(min_val, max_val);
        }
        min = min_val;
        max = max_val;
    }
};

template <typename T>
struct CartesianCoord2D
{
    T x;
    T y;

    CartesianCoord2D(T x_val, T y_val) : x{x_val}, y{y_val} {}
};

template <typename T>
struct CartesianCoord3D
{
    T x;
    T y;
    T z;

    CartesianCoord3D(T x_val, T y_val, T z_val) : x{x_val}, y{y_val}, z{z_val} {}
};

template <typename T>
struct SensorCoord2D
{
    T u;
    T v;

    SensorCoord2D(T u_val, T v_val) : u(u_val), v(v_val) {}
};

template <typename T>
struct SensorCoord3D
{
    T u;
    T v;
    T w;

    SensorCoord3D(T u_val, T v_val, T w_val) : u(u_val), v(v_val), w(w_val) {}
};

template <typename T>
struct CartesianCoordLimits
{
    Limit<T> x;
    Limit<T> y;
    std::variant<CartesianCoord2D<T>, CartesianCoord3D<T>> coord;

    // Constructor for 2D
    CartesianCoordLimits(T xmin, T xmax, T ymin, T ymax, T x_val, T y_val)
        : x(xmin, xmax), y(ymin, ymax), coord(CartesianCoord2D<T>(x_val, y_val)) {}

    // Constructor for 3D
    CartesianCoordLimits(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax, T x_val, T y_val, T z_val)
        : x(xmin, xmax), y(ymin, ymax), coord(CartesianCoord3D<T>(x_val, y_val, z_val)) {}
};