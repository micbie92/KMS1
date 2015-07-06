#include<stdlib.h>
#include <GL/glx.h>    
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include<iostream>
#include<fstream>

using std::fstream;
using namespace std;

int N=125;
int powtorzen=700;
double x,y,z;


ifstream polozenia;

void viev(int w, int h) 
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(45.0f, 1, 0.1f, 100.0f);
	glMatrixMode(GL_MODELVIEW);
}
void renderScene(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	gluLookAt(5.0f, 5.0f, 5.0f,
		  0.0f, 0.0f, 0.0f,
                  0.0f, 1.0f, 0.0f);
	glColor4f(0,0,1,0.2);
	GLfloat lightpos[] = {0, 0,3, 0};
	glLightfv(GL_LIGHT0,GL_AMBIENT, lightpos);
	glTranslated(0, 0, 0);
        glutSolidSphere(2.3, 30, 30);
	glColor4f(1,1,1,1);
	
	
	GLfloat lightpos2[] = {3, 0, 0, 0};
	glColor4f(1,1,1,1);
	glLightfv(GL_LIGHT0,GL_AMBIENT, lightpos2);
	for(int i=0; i<N; i++)
	{
            polozenia >> x >> y >> z;
	    glTranslated(x, y, z);
            glutSolidSphere(0.04, 30, 30);
            glTranslated(-x, -y, -z);
     }
      
	glEnd();
	glutSwapBuffers();

	if(polozenia.eof()) 
     {
    	 polozenia.seekg(0);
         exit(1);
}
}
int main(int argc, char **argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Animacja");
	glutFullScreen();
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
        glClearColor(0.0f, 0.0f, 0.0f, 0.1f); 
	polozenia.open("XYZ.dat");
	glutDisplayFunc(renderScene);
	glutReshapeFunc(viev);
	glutIdleFunc(renderScene);
	glutMainLoop();

	return 0;
}
