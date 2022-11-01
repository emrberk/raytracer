#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>

template<class T>
T multiply(T first, float scalar) {
    T result;
    result.x = first.x * scalar;
    result.y = first.y * scalar;
    result.z = first.z * scalar;
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

template<class T>
float dot(T first, T second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

template<class T>
T cross(T vec1, T vec2) {
    T result;
    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.x * vec2.z - vec1.z * vec2.x;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return result;
}

template<class T>
float length(T vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

template<class T>
T normalize(T vector) {
    float dot = dot(vector, vector);
    float k = 1.0 / sqrt(dot);
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

#endif