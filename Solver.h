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
    float dx;
    float dz;
    float dt;
    float* F1;
    float* F2;
    float* F3;
    float* G1;
    float* G2;
    float* G3;
    float limit;
    float h1;


  public:
    Solver();
    Solver(float Lx,float dx, float Lz,float dz);
    void run(float t);
    void boundary();
    Vector3 riemannX(Vector3& qL,Vector3& qR);
    Vector3 riemannY(Vector3& qL,Vector3& qR);
    float phi(float lambda);
    ~Solver();

    float* h;
    float* u;
    float* v;
    float* q1;
    float* q2;
    float* q3;
    

};

#endif
