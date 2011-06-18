#include "Visu.h"

Visu::Visu(float* &H, float LX,float nX, float LZ,float nZ)
{
  Lx = LX;
  Lz = LZ;
  nbX = nX;
  nbZ = nZ;
  nbPoints = nX*nZ;
  nbQuads = (nX-1)*(nZ-1);

  float dx = Lx / (nX-1);
  float dz = Lz / (nZ-1);

  vertices = new float[3*nbPoints];
  normals = new float[3*nbPoints];
  indices = new int[nbQuads*4];
  colors = new float[nbPoints*3];

  int i = 0;
  for (float z = 0; z <= Lz ; z += dz)
  {
    for (float x = 0; x <= Lx ; x += dx)
    {
      vertices[3*i] = x;
      vertices[3*i+1] = H[i];
      vertices[3*i+2] = z;
      i++;
    }
  }

  computeColors();
  computeIndices();
  updateNormals();
}

Visu::~Visu()
{
  delete[] vertices;
  delete[] normals;
  delete[] indices;
  delete[] colors;
}

void Visu::updateNormals()
{
  Vector3 v1;
  Vector3 v2;
  Vector3 v3;
  for (int i = 0; i < nbQuads ; i++)
  {
    v1.x = vertices[3*i+3] - vertices[3*i];
    v1.y = vertices[3*i+4] - vertices[3*i+1];
    v1.z = vertices[3*i+5] - vertices[3*i+2];

    v2.x = vertices[3*i+nbX*3] - vertices[3*i];
    v2.y = vertices[3*i+nbX*3+1] - vertices[3*i+1];
    v2.z = vertices[3*i+nbX*3+2] - vertices[3*i+2];
    v3 = v2.crossProduct(v1);
    v3.normalize();

    /* compute only the left of the quad */
    normals[3*i] = v3.x;
    normals[3*i+1] = v3.y;
    normals[3*i+2] = v3.z;
    normals[3*i+nbX*3] = v3.x;
    normals[3*i+nbX*3+1] = v3.y;
    normals[3*i+nbX*3+2] = v3.z;
  }

  int i = nbQuads;
  normals[3*i+3] = normals[3*i];
  normals[3*i+4] = normals[3*i+1];
  normals[3*i+5] = normals[3*i+2];

  normals[3*i+nbX*3+3] = normals[3*i+nbX*3];
  normals[3*i+nbX*3+4] = normals[3*i+nbX*3+1];
  normals[3*i+nbX*3+5] = normals[3*i+nbX*3+2];


}

void Visu::computeColors()
{
  for (int i = 0; i < nbPoints*3 ; i+= 3)
  {
    colors[i] = 0;
    colors[i+1] = 0;
    colors[i+2] = 1.0f;
  }
}

void Visu::computeIndices()
{
  int k = 0;
  for (int i = 0; i < nbQuads; i++)
  {
    if (i > 0 && i % (nbX-1) == 0)
      k++;
    indices[4*i] = k;
    indices[4*i+1] = k+1;
    indices[4*i+2] = nbX+1+k;
    indices[4*i+3] = nbX+k;

    k++;
  }
}

void Visu::updateHeight(float* &H)
{
  for (int i = 0; i < nbPoints; i++)
  {
    vertices[3*i+1] = H[i];
  }

  updateNormals();
}
