/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include "freeglut.h"
#include <stdbool.h>
#include "solver.h"

/* macros */

#define IX(i,j) ((i)+(N+2)*(j))

/* global variables */

static int N;
static float dt, diff, visc;
static float force, source;
static int dvel = false;

static float* u, * v, * u_prev, * v_prev;
static float* dens, * dens_prev;
static bool* boundary;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;

Solver solver;
/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/


static void clear_data(void)
{
	solver.reset();
}

static int allocate_data(void)
{
	

	return (1);
}


/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display(void)
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display(void)
{
	glutSwapBuffers();
}

static void draw_velocity(void)
{
	int i, j;
	float x, y, h;

	h = 1.0f / N;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 1; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 1; j <= N; j++) {
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + solver.getVelocityU(i, j), y + solver.getVelocityV(i, j));
		}
	}

	glEnd();
}

static void draw_boundary(void)
{

}

static void draw_density(void)
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f / N;

	glBegin(GL_QUADS);

	for (i = 0; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 0; j <= N; j++) {
			y = (j - 0.5f) * h;

			d00 = solver.getDensity(i, j);
			d01 = solver.getDensity(i, j + 1);
			d10 = solver.getDensity(i + 1, j);
			d11 = solver.getDensity(i + 1, j + 1);

			glColor3f(d00, d00, d00); glVertex2f(x, y);
			glColor3f(d10, d10, d10); glVertex2f(x + h, y);
			glColor3f(d11, d11, d11); glVertex2f(x + h, y + h);
			glColor3f(d01, d01, d01); glVertex2f(x, y + h);
		}
	}

	glEnd();
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

static void get_from_UI(float* d, float* u, float* v)
{
	int i, j, size = (N + 2) * (N + 2);

	for (i = 0; i < size; i++) {
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2]) return;

	i = (int)((mx / (float)win_x) * N + 1);
	j = (int)(((win_y - my) / (float)win_y) * N + 1);

	if (i<1 || i>N || j<1 || j>N) return;

	if (mouse_down[0]) {
		u[IX(i, j)] = force * (mx - omx);
		v[IX(i, j)] = force * (omy - my);

		d[IX(i, j)] = source;
	}

	if (mouse_down[2]) {
		boundary[IX(i, j)] = true;
	}

	omx = mx;
	omy = my;

	return;
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'c':
	case 'C':
		solver.reset();
		break;

	case 'q':
	case 'Q':
		
		exit(0);
		break;

	case 'v':
	case 'V':
		dvel = !dvel;
		break;

	case 'i':
	case 'I':
		solver.addSource();
		break;
	}
}

static void mouse_func(int button, int state, int x, int y)
{
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y)
{
	mx = x;
	my = y;
}

static void reshape_func(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void idle_func(void)
{
	solver.update();

	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void)
{
	pre_display();

	if (dvel)
		draw_velocity();
	else
	{
		draw_density();
		draw_boundary();
	}

	post_display();
}


/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window(void)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Media Lab | Korea University");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main(int argc, char** argv)
{

	glutInit(&argc, argv);

	if (argc != 1 && argc != 6) {
		fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
		fprintf(stderr, "where:\n"); \
			fprintf(stderr, "\t N      : grid resolution\n");
		fprintf(stderr, "\t dt     : time step\n");
		fprintf(stderr, "\t diff   : diffusion rate of the density\n");
		fprintf(stderr, "\t visc   : viscosity of the fluid\n");
		fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
		fprintf(stderr, "\t source : amount of density that will be deposited\n");
		exit(1);
	}

	if (argc == 1) {
		N = 64;
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		force = 5.0f;
		source = 100.0f;
		fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
			N, dt, diff, visc, force, source);
	}
	else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
	}

	solver.init(N);

	printf("\n\nHow to use this demo:\n\n");
	printf("\t Add densities with the right mouse button\n");
	printf("\t Add velocities with the left mouse button and dragging the mouse\n");
	printf("\t Toggle density/velocity display with the 'v' key\n");
	printf("\t Clear the simulation by pressing the 'c' key\n");
	printf("\t Quit by pressing the 'q' key\n");

	dvel = 0;

	if (!allocate_data()) exit(1);
	clear_data();

	win_x = 512;
	win_y = 512;
	open_glut_window();

	glutMainLoop();

	exit(0);
}