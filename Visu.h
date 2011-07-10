#ifndef __VISU
#define __VISU

#include <iostream>
#include "Vector3.h"

using namespace std;

class Visu
{
  private:
    int nbX;
    int nbZ;
    int nbPoints;
    float Lx;
    float Lz;
    float dx;
    float dz;


  public:

    Visu(float* &H, float LX,float nX, float LZ,float nZ);

    ~Visu();
    void computeColors();
    void updateNormals();
    void computeIndices();
    void updateHeight(float* &H);

    float* vertices;
    float* normals;
    int* indices;
    float* colors;
    int nbQuads;
};

#endif


