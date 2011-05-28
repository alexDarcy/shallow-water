#include "Vector3.h"
#include <iostream>
#include <cstdlib>
#include <math.h>

using namespace std;

Vector3::Vector3(): x(0),y(0),z(0) {}

Vector3::Vector3(double x0,double y0,double z0): x(x0),y(y0),z(z0) {}

Vector3::Vector3(const Vector3 & v): x(v.x),y(v.y),z(v.z) {}

Vector3& Vector3::operator= (const Vector3 & v){
  x = v.x;
  y = v.y;
  z = v.z;
  return *this;
}

Vector3& Vector3::operator+= (const Vector3 & v){
  x += v.x;
  y += v.y;
  z += v.z;
  return *this;
}

Vector3& Vector3::operator-= (const Vector3 & v){
  x -= v.x;
  y -= v.y;
  z -= v.z;
  return *this;
}
Vector3& Vector3::operator*= (const double a){
  x *= a;
  y *= a;
  z *= a;
  return *this;
}
Vector3& Vector3::operator/= (const double a){
  if (a == 0)
  {
    cerr << "Division by 0" << endl;
    exit(1);
  }
  x /= a;
  y /= a;
  z /= a;
  return *this;
}

Vector3 Vector3::operator+ (const Vector3 & v){
  Vector3 tmp = *this;
  tmp += v;
  return tmp;
}

Vector3 Vector3::operator- (const Vector3 & v){
  Vector3 tmp = *this;
  tmp -= v;
  return tmp;
}

Vector3 Vector3::operator* (const double a){
  Vector3 tmp = *this;
  tmp *= a;
  return tmp;
}
Vector3 Vector3::operator/ (const double a){
  Vector3 tmp = *this;
  tmp /= a;
  return tmp;
}

Vector3 Vector3::crossProduct(const Vector3 & v){
  Vector3 tmp = *this;
  tmp.x = y*v.z - z*v.y;
  tmp.y = -x*v.z + z*v.x;
  tmp.z = x*v.y - y*v.x;

  return tmp;
}

double Vector3::norm(){
  return sqrt(x*x + y*y + z*z);
}

void Vector3::normalize(){
  *this /= (*this).norm();
}
