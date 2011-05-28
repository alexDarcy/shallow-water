#include <stdlib.h>
#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include "Vector3.h"

#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)
#define Pi 3.141592654f


using namespace std;

/* Mesh data */
int nbX = 31;
int nbZ = 31;
double sizeX = 5; /* start from 0 */
double sizeZ = 5;
double stepX = sizeX / (nbX-1);
double stepZ = sizeZ / (nbZ-1);

double r = 10.0f; /* radius for the camera */
double y_0 = 0.0f;
double theta = 0; /* angle to 0y */
double phi = 0; /* angle to 0z */

/* camera position */
Vector3 cameraPos(sizeX/2 + r*cos(phi)*sin(theta),r*sin(phi),
    sizeZ/2 + r*cos(phi)*cos(theta));
Vector3 verticale(0,1,0); 
double xOrigin = -1;
double yOrigin = -1;
Vector3 cameraDir(sizeX/2,0,sizeZ/2);
double speed = 0.0001;

/* Data display */
int nbQuads = (nbX-1)*(nbZ-1);
GLfloat* vertices;
GLfloat* normals;
GLfloat* colors;
bool showPoints = false;
double t;
double tmp;

/* Height of the waves */
double height(double x, double z, double t)
{
  //x++;
  //t++;
  //z++;
 // return 0.5;//0.5*sin(x*8*Pi);//*sin(z*8*Pi);
  return 0.02*(cos(x*10*Pi - t)+cos(z*10*Pi - t));//*sin(z*8*Pi);
}

/* Display the scene */
void display(void)
{
  t = glutGet (GLUT_ELAPSED_TIME) / 1000.;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  gluLookAt( cameraPos.x, cameraPos.y,cameraPos.z, /* only the camera move */
      cameraDir.x, cameraDir.y, cameraDir.z,
      verticale.x,verticale.y,verticale.z);

  // Draw ground
  glColor3f (1, 0.9, 0.7);
  glBegin(GL_QUADS);
  glVertex3f(-5.0f, -1.0f, -5.0f);
  glVertex3f(-5.0f, -1.0f,  5.0f);
  glVertex3f( 5.0f, -1.0f,  5.0f);
  glVertex3f( 5.0f, -1.0f, -5.0f);
  glEnd();

  int i = 0;
  for (double z = 0; z < sizeZ ; z += stepZ)
  {
    for (double x = 0; x <sizeX ; x += stepX)
    {
      vertices[i] = x;
      vertices[i+1] = height(x, z, t);
      vertices[i+2] = z;

      vertices[i+3] = x + stepX;
      vertices[i+4] = height(x + stepX, z, t);
      vertices[i+5] = z;

      vertices[i+6] = x + stepX;
      vertices[i+7] = height(x + stepX, z + stepZ, t);
      vertices[i+8] = z + stepZ;

      vertices[i+9]  = x;
      vertices[i+10] = height(x, z + stepZ, t);
      vertices[i+11] = z + stepZ;

      i+= 12;
    }
  }
  Vector3 v1;
  Vector3 v2;
  Vector3 v3;
  i = 0;
  for (int i = 0; i < nbQuads*3*4 ; i+= 12)
  {
    v1.x = vertices[i+3] - vertices[i];
    v1.y = vertices[i+4] - vertices[i+1];
    v1.z = vertices[i+5] - vertices[i+2];

    v2.x = vertices[i+6] - vertices[i];
    v2.y = vertices[i+7] - vertices[i+1];
    v2.z = vertices[i+8] - vertices[i+2];
    v3 = v2.crossProduct(v1);
    v3.normalize();

    for (int j = i; j < i+12 ; j+=3)
    {
      normals[j] = v3.x;
      normals[j+1] = v3.y;
      normals[j+2] = v3.z;
    }
  }
  for (int i = 0; i < nbQuads*3*4 ; i+= 3)
  {
    colors[i] = 0;
    colors[i+1] = 0;
    colors[i+2] = 1.0f;
  }

  if (showPoints)
    glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
  else
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

  /* Draw cube */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, vertices);
  glNormalPointer(GL_FLOAT, 0, normals);
  glColorPointer(3, GL_FLOAT, 0, colors);
  glDrawArrays(GL_QUADS, 0, nbQuads*4);

  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  //glColor3f (1, 0, 0);
  //glBegin (GL_LINES);
  //for (int i = 0; i < nbQuads*4*3; i+= 3)
  //{
  //  glVertex3dv (&vertices[i]);
  //  glVertex3d (vertices[i] + normals[i] ,
  //              vertices[i + 1]  + normals[i + 1] ,
  //              vertices[i + 2] + normals[i + 2] );
  //}

  glEnd ();
  
  /* End */
  glFlush();
  glutSwapBuffers();

  /* Update again and again */
  glutPostRedisplay();
}

// When the window is created or resized
void ReshapeFunc(int width, int height)
{
  glMatrixMode(GL_PROJECTION);

  glLoadIdentity();
  gluPerspective(20, width / (double) height, 5, 15);
  glViewport(0, 0, width, height);

  glMatrixMode(GL_MODELVIEW);
  glutPostRedisplay();
}

// When a key is hit
void keyHit(unsigned char key, int x, int y)
{
  /* useless */
  tmp = x;
  tmp = y;

  switch (key)
  {
    case 'p' :
      showPoints = !showPoints;
      break;
    case 'q' :
      exit(0);
      break;
  }
}

/* When mouse is clicked */
void mouseClicked(int button, int state, int x, int y)
{
  if (GLUT_LEFT_BUTTON == button)
  {
    if (state == GLUT_UP)
    {
      xOrigin = -1;
      yOrigin = -1;
    }
    else
    {
      xOrigin = x;
      yOrigin = y;
    }
  }

  tmp = x;
  tmp = y;
}

/* When the mouse is moved */
void mouseMoved(int x, int y)
{
  Vector3 v(-1,0,0);
  if (xOrigin >= 0) 
  {
    theta += (x - xOrigin)*speed;
    theta = max(theta, -90);
    theta = min(theta, 90);

    phi += (y - yOrigin)*speed;
    phi = max(phi, -90);
    phi = min(phi, 90);

    cameraPos.y = r*sin(phi);
    cameraPos.z = sizeZ/2 + r*cos(phi)*cos(theta);
    cameraPos.x = sizeX/2 + r*cos(phi)*sin(theta);
    verticale = v.crossProduct(cameraPos);
    verticale.normalize();
  }
  tmp = x;
  tmp = y;

}

int main(int argc, char* argv[])
{
  vertices = new GLfloat[nbQuads*4*3];
  normals = new GLfloat[nbQuads*4*3];
  colors = new GLfloat[nbQuads*4*3];

  /* Creation of the window */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(500, 500);
  glutInitWindowPosition (100, 100); 
  glutCreateWindow("Water");

  glClearColor(0, 0, 0, 0);
  /* Light */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
  GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
  GLfloat position[] = { -1.5f, 1.0f, -4.0f, 1.0f };
  
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glEnable(GL_COLOR_MATERIAL);
   glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  /* Declaration of the callbacks */
  glutDisplayFunc(&display);
  glutReshapeFunc(&ReshapeFunc);
  glutKeyboardFunc(&keyHit);
  glutMouseFunc(&mouseClicked);
  glutMotionFunc(&mouseMoved);

  glutMainLoop (); // the main loop
  return 0;
}
