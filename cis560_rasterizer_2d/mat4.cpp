//written by Jiawei Wang
//University of Pennsylvania
#include "mat4.h"
#include "vec4.h"
#include <iostream>
#include <assert.h>
#include <math.h>
#include <iomanip>

using namespace std;

mat4::mat4(){
    vec4 v;
    data[0] = v;
    data[0][0] = 1.0;
    data[1] = v;
    data[1][1] = 1.0;
    data[2] = v;
    data[2][2] = 1.0;
    data[3] = v;
    data[3][3] = 1.0;
}

/// Initializes the diagonal values of the matrix to diag. All other values are 0.
mat4::mat4(float diag){
    data[0][0] = diag;
    data[1][1] = diag;
    data[2][2] = diag;
    data[3][3] = diag;
}

/// Initializes matrix with each vector representing a column in the matrix
mat4::mat4(const vec4 &col0, const vec4 &col1, const vec4 &col2, const vec4& col3){
    data[0] = col0;
    data[1] = col1;
    data[2] = col2;
    data[3] = col3;
}

// copy constructor
mat4::mat4(const mat4 &m2){
    data[0] = m2.data[0];
    data[1] = m2.data[1];
    data[2] = m2.data[2];
    data[3] = m2.data[3];
}

/// Returns the values of the column at the index
vec4  mat4::operator[](unsigned int index) const{
    assert(index>=0);
    assert(index<=3);
    return data[index];
}

/// Returns a reference to the column at the index
vec4 &mat4::operator[](unsigned int index){
    assert(index>=0);
    assert(index<=3);
    return data[index];
}

/// Returns a row of the input matrix
vec4 row(const mat4 &m, unsigned int index){
    assert(index>=0);
    assert(index<=3);
    vec4 v;
    v[0] = m[0][index];
    v[1] = m[1][index];
    v[2] = m[2][index];
    v[3] = m[3][index];
    return v;
}

/// Returns the transpose of the input matrix (v_ij == v_ji)
mat4 transpose(const mat4 &m){
    mat4 m1;
    m1[0] = row(m, 0);
    m1[1] = row(m, 1);
    m1[2] = row(m, 2);
    m1[3] = row(m, 3);
    return m1;
}

/// Creates a 3-D rotation matrix.
/// Takes an angle in degrees and an axis represented by its xyz components, and outputs a 4x4 rotation matrix
/// You may choose to only handle the three cardinal axes of rotation.
mat4 mat4::rotate(float angle, float x, float y, float z){
    double PI = 3.141592653589793238463;
    float radians = PI * angle / 180;
    float sin_angle = sin(radians);
    float cos_angle = cos(radians);
    // Roatate by x axis
    if(x == 1 && y == 0 && z == 0){
        vec4 v0(1.0, 0.0, 0.0, 0.0), v1(0.0, cos_angle, sin_angle, 0.0), v2(0.0, -sin_angle, cos_angle, 0.0), v3(0.0, 0.0, 0.0, 1.0);
        mat4 m_r_x(v0, v1, v2, v3);
        return m_r_x;
    }
    if(x == 0 && y == 1 && z == 0){
        vec4 v0(cos_angle, 0.0, -sin_angle, 0.0), v1(0.0, 1.0, 0.0, 0.0), v2(sin_angle, 0.0, cos_angle, 0.0), v3(0.0, 0.0, 0.0, 1.0);
        mat4 m_r_y(v0, v1, v2, v3);
        return m_r_y;
    }
    if(x == 0 && y == 0 && z == 1){
        vec4 v0(cos_angle, sin_angle, 0.0, 0.0), v1(-sin_angle, cos_angle, 0.0, 0.0), v2(0.0, 0.0, 1.0, 0.0), v3(0.0, 0.0, 0.0, 1.0);
        mat4 m_r_z(v0, v1, v2, v3);
        return m_r_z;
    }
}

/// Takes an xyz displacement and outputs a 4x4 translation matrix
mat4 mat4::translate(float x, float y, float z){
    vec4 v0(1.0, 0.0, 0.0, 0.0), v1(0.0, 1.0, 0.0, 0.0), v2(0.0, 0.0, 1.0, 0.0), v3(x, y, z, 1.0);
    mat4 m(v0, v1, v2, v3);
    return m;
}

/// Takes an xyz scale and outputs a 4x4 scale matrix
mat4 mat4::scale(float x, float y, float z){
    vec4 v0(x, 0.0, 0.0, 0.0), v1(0.0, y, 0.0, 0.0), v2(0.0, 0.0, z, 0.0), v3(0.0, 0.0, 0.0, 1.0);
    mat4 m(v0, v1, v2, v3);
    return m;
}

/// Generates a 4x4 identity matrix
mat4 mat4::identity(){
    vec4 v0(1.0, 0.0, 0.0, 0.0), v1(0.0, 1.0, 0.0, 0.0), v2(0.0, 0.0, 1.0, 0.0), v3(0.0, 0.0, 0.0, 1.0);
    mat4 m(v0, v1, v2, v3);
    return m;
}

///----------------------------------------------------------------------
/// Operator Functions
///----------------------------------------------------------------------

/// Assign m2 to this and return it
mat4 &mat4::operator=(const mat4 &m2){
    data[0] = m2[0];
    data[1] = m2[1];
    data[2] = m2[2];
    data[3] = m2[3];
    return *this;
}

/// Test for equality
bool mat4::operator==(const mat4 &m2) const{
    for(int i = 0; i < 4; i++){
        if(data[i] != m2[i])
            return false;
    }
    return true;
}

/// Test for inequality
bool mat4::operator!=(const mat4 &m2) const{
    for(int i = 0; i < 4; i++){
        if(data[i] != m2[i])
            return true;
    }
    return false;
}

/// Element-wise arithmetic
/// e.g. += adds the elements of m2 to this and returns this (like regular +=)
///      +  returns a new matrix whose elements are the sums of this and v2
mat4 &mat4::operator+=(const mat4 &m2){
    data[0] += m2[0];
    data[1] += m2[1];
    data[2] += m2[2];
    data[3] += m2[3];
    return *this;
}

mat4 &mat4::operator-=(const mat4 &m2){
    data[0] -= m2[0];
    data[1] -= m2[1];
    data[2] -= m2[2];
    data[3] -= m2[3];
    return *this;
}

mat4 &mat4::operator*=(float c){                 // multiplication by a scalar
    data[0] *= c;
    data[1] *= c;
    data[2] *= c;
    data[3] *= c;
    return *this;
}
mat4 &mat4::operator/=(float c){                 // division by a scalar
    data[0] /= c;
    data[1] /= c;
    data[2] /= c;
    data[3] /= c;
    return *this;
}
mat4  mat4::operator+(const mat4 &m2) const{
    mat4 m;
    m[0] = data[0] + m2[0];
    m[1] = data[1] + m2[1];
    m[2] = data[2] + m2[2];
    m[3] = data[3] + m2[3];
    return m;
}
mat4  mat4::operator-(const mat4 &m2) const{
    mat4 m;
    m[0] = data[0] - m2[0];
    m[1] = data[1] - m2[1];
    m[2] = data[2] - m2[2];
    m[3] = data[3] - m2[3];
    return m;
}
mat4  mat4::operator*(float c) const{             // multiplication by a scalar
    mat4 m;
    m[0] = data[0] * c;
    m[1] = data[1] * c;
    m[2] = data[2] * c;
    m[3] = data[3] * c;
    return m;
}
mat4  mat4::operator/(float c) const{             // division by a scalar
    mat4 m;
    m[0] = data[0] / c;
    m[1] = data[1] / c;
    m[2] = data[2] / c;
    m[3] = data[3] / c;
    return m;
}

/// Scalar multiplication (c * m)
mat4 operator*(float c, const mat4 &m){
    mat4 m1;
    m1[0] = m[0] * c;
    m1[1] = m[1] * c;
    m1[2] = m[2] * c;
    m1[3] = m[3] * c;
    return m1;
}

/// Matrix multiplication (m1 * m2)
mat4 mat4::operator*(const mat4 &m2) const{
    mat4 m;
    mat4 thisT;
    thisT = transpose(*this);
    for(int i = 0;i < 4;i++){
        for(int j = 0;j < 4;j++){
            m[j][i] = dot(thisT[i], m2[j]);
        }
    }
    return m;
}

/// Matrix/vector multiplication (m * v)
/// Assume v is a column vector (ie. a 4x1 matrix)
vec4 mat4::operator*(const vec4 &v) const{
    vec4 vec;
    mat4 thisT;
    thisT = transpose(*this);
    vec[0] = dot(thisT[0], v);
    vec[1] = dot(thisT[1], v);
    vec[2] = dot(thisT[2], v);
    vec[3] = dot(thisT[3], v);
    return vec;
}

/// Vector/matrix multiplication (v * m)
/// Assume v is a row vector (ie. a 1x4 matrix)
vec4 operator*(const vec4 &v, const mat4 &m){
    vec4 vec;
    vec[0] = dot(v, m[0]);
    vec[1] = dot(v, m[1]);
    vec[2] = dot(v, m[2]);
    vec[3] = dot(v, m[3]);
    return vec;
}

/// Prints the matrix to a stream in a nice format
std::ostream &operator<<(std::ostream &o, const mat4 &m){
    cout<<"["<<m[0][0]<<"  "<<m[1][0]<<"  "<<m[2][0]<<"  "<<m[3][0]<<endl;
    cout<<" "<<m[0][1]<<"  "<<m[1][1]<<"  "<<m[2][1]<<"  "<<m[3][1]<<endl;
    cout<<" "<<m[0][2]<<"  "<<m[1][2]<<"  "<<m[2][2]<<"  "<<m[3][2]<<endl;
    cout<<" "<<m[0][3]<<"  "<<m[1][3]<<"  "<<m[2][3]<<"  "<<m[3][3]<<"]"<<endl;
    return o;
}
