#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include "parser.h"

template<class T>
T multiply(T first, float scalar) {
    T result;
    result.x = first.x * scalar;
    result.y = first.y * scalar;
    result.z = first.z * scalar;
    return result;
}

template<class T>
T multiplyTwo(T first, T second) {
    T result;
    result.x = first.x * second.x;
    result.y = first.y * second.y;
    result.z = first.z * second.z;
    return result;
}

template<class T>
T add(T first, T second) {
    T result;
    result.x = first.x + second.x;
    result.y = first.y + second.y;
    result.z = first.z + second.z;
    return result;
}

template<class T, class N>
T add(T first, N num) {
    T result;
    result.x = first.x + num;
    result.y = first.y + num;
    result.z = first.z + num;
    return result;
}

template<class T>
T subtract(T first, T second) {
    T result;
    result.x = first.x - second.x;
    result.y = first.y - second.y;
    result.z = first.z - second.z;
    return result;
}

template<class T>
float dot(T first, T second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

template<class T>
T cross(T vec1, T vec2) {
    T result;
    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.z * vec2.x - vec1.x * vec2.z;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return result;
}

template<class T, class N>
T divide(T vec, N number) {
    T result;
    result.x = vec.x / number;
    result.y = vec.y / number;
    result.z = vec.z / number;
    return result;
}

template<class T>
float length(T vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

template<class T>
T normalize(T vector) {
    float dotProduct = dot(vector, vector);
    float k = 1.0 / sqrt(dotProduct);
    return multiply(vector, k);
}

parser::Vec3f convert(parser::Vec3i i) {
    parser::Vec3f result;
    result.x = i.x;
    result.y = i.y;
    result.z = i.z;
    return result;
}

template<class T>
void printVec(T vec) {
    std::cout << vec.x << " " << vec.y << " " << vec.z << std::endl;
}

template<class T>
float determinant(T i00, T i01, T i02, T i10, T i11, T i12, T i20, T i21, T i22) {
    return i00 * (i11 * i22 - i12* i21) + i10 * (i02 * i21 - i01 * i22) + i20 * (i01 * i12 - i11 * i02);
}

template<class T>
T scale(T vec) {
    return add(multiply(vec, 255), 0.5);
}

template<class T>
T clamp(T a) {
    if (a > 255) {
        return 255;
    } else if (a < 0) {
        return 0;
    }
    return a;
}

template<class T>
T clampVec(T vec) {
    T result;
    result.x = clamp(vec.x);
    result.y = clamp(vec.y);
    result.z = clamp(vec.z);
    return result;
}

bool parser::Vec3f::operator!=(Vec3f a) {
    return x != a.x || y != a.y || z == a.z;
}

#endif