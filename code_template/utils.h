#ifndef _UTILS_H
#define _UTILS_H

#include <math.h>

template<class T>
T multiply(T first, float scalar);

template<class T>
T add(T first, T second);

template<class T>
float dot(T first, T second);

template<class T>
T cross(T vec1, T vec2);

template<class T>
float length(T vec);

template<class T>
T normalize(T vector);

#endif