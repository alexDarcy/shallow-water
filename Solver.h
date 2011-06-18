#ifndef __VECTOR_3_H
#define __VECTOR_3_H

#include <iostream>
#include <math.h>
#include "Vector3.h"
#define abs(a) (((a) < 0) ? -(a) : (a))
#define Pi 3.141592654f

using namespace std;

class Solver
{
  private:
    float g;
    int nbX;
    int nbZ;
    int nbPoints;
    int nbQuads;
    float Lx;
    float Lz;

  public:
    Solver();
    Solver(float Lx,float dx, float Lz,float dz);
    void run(float t);
    Vector3 riemannX(float* qL,float* qR);
    Vector3 riemannY(float* qL,float* qR);
    float phi(float lambda);
    ~Solver();

    float* h;
    float* u;
    float* v;

};

#endif
