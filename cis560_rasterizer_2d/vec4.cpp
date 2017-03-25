//written by Jiawei Wang
//University of Pennsylvania
#include "vec4.h"
#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

vec4::vec4(){
    for(int i=0;i<4;i++){
        data[i]=0.0;
    }
}

vec4::vec4(float x, float y, float z, float w){
    data[0] = x;
    data[1] = y;
    data[2] = z;
    data[3] = w;
}

// copy constructor
vec4::vec4(const vec4 &v2) {
    data[0] = v2.data[0];
    data[1] = v2.data[1];
    data[2] = v2.data[2];
    data[3] = v2.data[3];
}

/// Returns the value at index
float vec4::operator[](unsigned int index) const{
    assert(index >= 0);
    assert(index <= 3);
    return data[index];
}

/// Returns a reference to the value at index
float &vec4::operator[](unsigned int index){
    assert(index >= 0);
    assert(index <= 3);
    return data[index];
}

vec4 &vec4::operator =(const vec4 &v2){
    data[0] = v2.data[0];
    data[1] = v2.data[1];
    data[2] = v2.data[2];
    data[3] = v2.data[3];
    return *this;
}

/// Test for equality
bool vec4::operator==(const vec4 &v2) const{
    float epsilon = 0.00001;
    for(int i = 0; i < 4; i++){
        if(fabs(data[i] - v2.data[i]) > epsilon)
            return false;
    }
    return true;
}

/// Test for inequality
bool vec4::operator!=(const vec4 &v2) const{
    float epsilon = 0.00001;
    for(int i = 0; i < 4; i++){
        if(fabs(data[i] - v2.data[i]) > epsilon)
            return true;
    }
    return false;
}

/// Arithmetic:
/// e.g. += adds v2 to this and return this (like regular +=)
///      +  returns a new vector that is sum of this and v2
vec4 &vec4::operator+=(const vec4 &v2){
    data[0] += v2.data[0];
    data[1] += v2.data[1];
    data[2] += v2.data[2];
    data[3] += v2.data[3];
    return *this;
}

vec4 &vec4::operator-=(const vec4 &v2){
    data[0] -= v2.data[0];
    data[1] -= v2.data[1];
    data[2] -= v2.data[2];
    data[3] -= v2.data[3];
    return *this;
}

// multiplication by a scalar
vec4 &vec4::operator*=(float c){
    data[0] *= c;
    data[1] *= c;
    data[2] *= c;
    data[3] *= c;
    return *this;
}

// division by a scalar
vec4 &vec4::operator/=(float c){
    data[0] /= c;
    data[1] /= c;
    data[2] /= c;
    data[3] /= c;
    return *this;
}

vec4 vec4::operator+(const vec4 &v2) const{
    vec4 v;
    v.data[0] = data[0] + v2.data[0];
    v.data[1] = data[1] + v2.data[1];
    v.data[2] = data[2] + v2.data[2];
    v.data[3] = data[3] + v2.data[3];
    return v;
}
vec4 vec4::operator-(const vec4 &v2) const{
    vec4 v;
    v.data[0] = data[0] - v2.data[0];
    v.data[1] = data[1] - v2.data[1];
    v.data[2] = data[2] - v2.data[2];
    v.data[3] = data[3] - v2.data[3];
    return v;

}

// multiplication by a scalar
vec4 vec4::operator*(float c) const{
    vec4 v;
    v.data[0] = data[0] * c;
    v.data[1] = data[1] * c;
    v.data[2] = data[2] * c;
    v.data[3] = data[3] * c;
    return v;
}

// division by a scalar
vec4 vec4::operator/(float c) const{
    vec4 v;
    v.data[0] = data[0] / c;
    v.data[1] = data[1] / c;
    v.data[2] = data[2] / c;
    v.data[3] = data[3] / c;
    return v;
}

/// Dot Product
float dot(const vec4 &v1, const vec4 &v2){
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]);
}

//Compute the result of v1 x v2 using only their X, Y, and Z elements.
//The fourth element of the resultant vector should be 0.
/// Cross Product
vec4 cross(const vec4 &v1, const vec4 &v2){
    vec4 v;
    v[0] = v1[1] * v2[2] - v2[1] * v1[2];
    v[1] = v1[2] * v2[0] - v2[2] * v1[0];
    v[2] = v1[0] * v2[1] - v2[0] * v1[1];
    return v;
}

/// Returns the geometric length of the input vector
float length(const vec4 &v){
    double lenlen;
    lenlen = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
    return sqrt(lenlen);
}

/// Scalar Multiplication (c * v)
vec4 operator*(float c, const vec4 &v){
    return v * c;
}

vec4 normalize(const vec4& v){
    float len = length(v);
    return v / len;
}

/// Prints the vector to a stream in a nice format
std::ostream &operator<<(std::ostream &o, const vec4 &v){
    cout<<"["<<v[0]<<"  "<<v[1]<<"  "<<v[2]<<"  "<<v[3]<<"]"<<endl;
    return o;
}

