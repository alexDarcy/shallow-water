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
int nbX = 50;
int nbZ = 50;
int nbPoints = nbX*nbZ;
double LX = 5; /* start from 0 */
double LZ = 5;
int nbQuads = (nbX-1)*(nbZ-1);
Solver* s;
Visu* v;

GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
GLfloat positionLight[] = { -1.5f, 1.0f, -4.0f, 1.0f };

double r = 10.0f; /* radius for the camera */
double theta = 0; /* angle around 0y */
double phi = Pi/8; /* angle around 0z */

/* camera position */
Vector3 cameraPos(LX/2 + r*cos(phi)*sin(theta),r*sin(phi),
    LZ/2 + r*cos(phi)*cos(theta));
Vector3 verticale(0,1,0); 
double xOrigin = -1;
double yOrigin = -1;
Vector3 cameraDir(LX/2,0,LZ/2);
double speed = 0.0001;

/* Data display */
GLfloat* colors;
GLuint* indices;
bool showPoints = false;
double t;
double tmp;

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

  if (t > 2.0f) // not too fast
  {
    s->run(t);
    v->updateHeight(s->q1);
  }

  if (showPoints)
    glPolygonMode (GL_FRONT_AND_BACK, GL_POINT);
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
    cameraPos.z = LZ/2 + r*cos(phi)*cos(theta);
    cameraPos.x = LX/2 + r*cos(phi)*sin(theta);
    verticale = v.crossProduct(cameraPos);
    verticale.normalize();
  }
  tmp = x;
  tmp = y;

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
