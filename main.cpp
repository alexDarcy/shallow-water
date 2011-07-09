#include <stdlib.h>
#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <sys/resource.h>
#include <stdio.h>
#include "Vector3.h"
#include "Solver.h"
#include "Visu.h"

#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)


using namespace std;

/* Mesh data */
int nbX = 64;
int nbZ = 64;
int nbPoints = nbX*nbZ;
float LX = 4; /* start from 0 */
float LZ = 4;
int nbQuads = (nbX-1)*(nbZ-1);
Solver* s;
Visu* v;

GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
GLfloat positionLight[] = { 2.f, 5.0f, 2.f, 1.0f };

float r = 8.0f; /* radius for the camera */
float theta = Pi/15; /* angle around 0y */
float phi = Pi/7; /* angle around 0z */
float yInit = 1.f;

/* camera position */
Vector3 cameraPos(r*cos(phi)*sin(theta),
  r*sin(phi), r*cos(phi)*cos(theta));

Vector3 verticale(0,1,0); 
float xOrigin = -1;
float yOrigin = -1;
Vector3 cameraDir(LX/2,0,LZ/2);
float speed = 0.0001;

/* Data display */
GLfloat* colors;
GLuint* indices;
bool showLines = true;
float t;

/* Display the scene */
void display(void)
{
  t = glutGet (GLUT_ELAPSED_TIME) / 1000.;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  gluLookAt( cameraPos.x, cameraPos.y,cameraPos.z, /* only the camera move */
      cameraDir.x, cameraDir.y, cameraDir.z,
      verticale.x,verticale.y,verticale.z);


  if (t > 2.0f) // not too fast
  {
    s->run(t);
    v->updateHeight(s->q1);
  }

  if (showLines)
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 0, v->vertices);
  glNormalPointer(GL_FLOAT, 0, v->normals);
  glColorPointer(3, GL_FLOAT, 0, v->colors);
  glDrawElements(GL_QUADS, nbQuads*4, GL_UNSIGNED_INT, v->indices);


  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  glEnd ();
  glFlush();
  glutSwapBuffers();

  /* Update again and again */
  glutPostRedisplay();
}

// When the window is created 
void ReshapeFunc(int width, int height)
{
  glMatrixMode(GL_PROJECTION);

  glLoadIdentity();
  gluPerspective(40, width / (float) height, 1, 15);
  glViewport(0, 0, width, height);

  glMatrixMode(GL_MODELVIEW);
  glutPostRedisplay();
}

void updateCamera()
{
  cameraPos.x = r*cos(phi)*sin(theta);
  cameraPos.y = r*sin(phi);
  cameraPos.z = r*cos(phi)*cos(theta);
}

// When a key is hit
void keyHit(unsigned char key, int x, int y)
{
  /* useless */
  float tmp = x;
  tmp = y;

  switch (key)
  {
    case 'l' :
      showLines = !showLines;
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

    updateCamera();
    verticale = v.crossProduct(cameraPos);
    verticale.normalize();
  }
}

void changeStackSize()
{
  const rlim_t sizeCur = 64L*1024L*1024L;   // min stack size = 64 Mb
  const rlim_t sizeMax = 100L*1024L*1024L;   // min stack size = 64 Mb
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
  {
    if (rl.rlim_cur < sizeCur)
    {
      rl.rlim_max = sizeCur;
      rl.rlim_cur = sizeMax;
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)
      {
        fprintf(stderr, "Failed to change stack size : result = %d\n", result);
      }
    }
  }
}

int main(int argc, char* argv[])
{
  s = new Solver(LX,nbX,LZ,nbZ);
  v = new Visu(s->h,LX,nbX,LZ,nbZ);

  /* Creation of the window */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(500, 500);
  glutInitWindowPosition (100, 100); 
  glutCreateWindow("Water");
   
  glClearColor(0, 0, 0, 0);
  /* Light */
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, positionLight);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  //glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  GLfloat mat_specular[] = { 1.0F,1.0F,1.0F,1.0F };
  glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
  glMateriali(GL_FRONT,GL_SHININESS,30);


    /* Declaration of the callbacks */
  glutDisplayFunc(&display);
  glutReshapeFunc(&ReshapeFunc);
  glutKeyboardFunc(&keyHit);
  glutMouseFunc(&mouseClicked);
  glutMotionFunc(&mouseMoved);

  glutMainLoop (); // the main loop
  return 0;
}
