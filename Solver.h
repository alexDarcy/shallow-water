#ifndef __VECTOR_3_H
#define __VECTOR_3_H


#include <math.h>
#include "Vector3.h"
#define abs(a) (((a) < 0) ? -(a) : (a))

using namespace std;

class Solver
{
  private:
    float g;

  public:
    Solver();
    Vector3 riemannX(float* qL,float* qR);
    Vector3 riemannY(float* qL,float* qR);
    float phi(float lambda);
    ~Solver();
};

#endif
