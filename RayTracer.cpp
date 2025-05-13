/*==================================================================================
* COSC 363  Computer Graphics
* Department of Computer Science and Software Engineering, University of Canterbury.
*
* A basic ray tracer
* See Lab07.pdf   for details.
*===================================================================================
*/
#include <iostream>
#include <cmath>
#include "Plane.h"
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/freeglut.h>
#include "TextureBMP.h"
using namespace std;

TextureBMP texture;
const float EDIST = 40.0;
const int NUMDIV = 500;
const int MAX_STEPS = 5;
const float XMIN = -10.0;
const float XMAX = 10.0;
const float YMIN = -10.0;
const float YMAX = 10.0;

vector<SceneObject*> sceneObjects;

//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step) {
	glm::vec3 backgroundCol(0);					//Background colour = (0,0,0)
	glm::vec3 lightPos(0, 25, -80);				//Light's position inside the Cornell box
	glm::vec3 color(0);
	SceneObject* obj;

	ray.closestPt(sceneObjects);				//Compare the ray with all objects in the scene
	if(ray.index == -1) return backgroundCol;		//no intersection
	obj = sceneObjects[ray.index];
	
	// Check if the hit object is the floor (index 0) to create checkered pattern
	if (ray.index == 0)
	{
		// Create a checkered pattern for the floor
		int squareSize = 10; // Size of each square in the checkered pattern
		
		// Calculate checkerboard pattern based on x and z coordinates
		int ix = floor(ray.hit.x / squareSize);
		int iz = floor(ray.hit.z / squareSize);
		
		// If (ix+iz) is even, color is white, otherwise black
		if ((ix + iz) % 2 == 0) {
			color = glm::vec3(1, 1, 1); // White square
		} else {
			color = glm::vec3(0, 0, 0); // Black square
		}
		obj->setColor(color);
	}
	//object on which the closest point of intersection is found


	color = obj->lighting(lightPos, -ray.dir, ray.hit);						//Object's colour
	glm::vec3 LightVec = lightPos - ray.hit;
	Ray shadowRay(ray.hit, LightVec);
	shadowRay.closestPt(sceneObjects);
	float lightDist = glm::length(LightVec);
	if (shadowRay.index > -1 && shadowRay.dist < lightDist) {
		// In shadow
		color = 0.2f * obj->getColor();  // Only ambient lighting
	}
	if (obj->isReflective() && step < MAX_STEPS)
	{
		float rho = obj->getReflectionCoeff();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
		Ray reflectedRay(ray.hit, reflectedDir);
		glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
		color = color + (rho * reflectedColor);
	}
	return color;
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display() {
	float xp, yp;  //grid point
	float cellX = (XMAX - XMIN) / NUMDIV;  //cell width
	float cellY = (YMAX - YMIN) / NUMDIV;  //cell height
	glm::vec3 eye(0., 0., 0.);

	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a tiny quad.

	for (int i = 0; i < NUMDIV; i++) {	//Scan every cell of the image plane
		xp = XMIN + i * cellX;
		for (int j = 0; j < NUMDIV; j++) {
			yp = YMIN + j * cellY;

			glm::vec3 dir(xp + 0.5 * cellX, yp + 0.5 * cellY, -EDIST);	//direction of the primary ray

			Ray ray = Ray(eye, dir);

			glm::vec3 col = trace(ray, 1); //Trace the primary ray and get the colour value
			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp + cellX, yp);
			glVertex2f(xp + cellX, yp + cellY);
			glVertex2f(xp, yp + cellY);
		}
	}

	glEnd();
	glFlush();
}

//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL 2D orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize() {
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
	texture = TextureBMP("../Butterfly.bmp");
	glClearColor(0, 0, 0, 1);

	// Create a Cornell-style box
	// Floor - checkered pattern will be implemented in the trace function
	Plane *floor = new Plane(glm::vec3(-30, -15, -40),      // Point A
		                    glm::vec3(30, -15, -40),       // Point B
		                    glm::vec3(30, -15, -120),      // Point C
		                    glm::vec3(-30, -15, -120));    // Point D
	floor->setColor(glm::vec3(1, 1, 1));                // White base color for checkered pattern
	floor->setSpecularity(false);
	sceneObjects.push_back(floor);

	// Left wall - green
	Plane *leftWall = new Plane(glm::vec3(-30, -15, -40),    // Point A
		                       glm::vec3(-30, -15, -120),   // Point B
		                       glm::vec3(-30, 30, -120),    // Point C
		                       glm::vec3(-30, 30, -40));    // Point D
	leftWall->setColor(glm::vec3(0, 1, 0));              // Green
	leftWall->setSpecularity(false);
	sceneObjects.push_back(leftWall);

	// Right wall - red
	Plane *rightWall = new Plane(glm::vec3(30, -15, -40),    // Point A
		                        glm::vec3(30, 30, -40),     // Point B
		                        glm::vec3(30, 30, -120),    // Point C
		                        glm::vec3(30, -15, -120));  // Point D
	rightWall->setColor(glm::vec3(1, 0, 0));             // Red
	rightWall->setSpecularity(false);
	sceneObjects.push_back(rightWall);

	// Back wall - white
	Plane *backWall = new Plane(glm::vec3(-30, -15, -120),   // Point A
		                       glm::vec3(30, -15, -120),    // Point B
		                       glm::vec3(30, 30, -120),     // Point C
		                       glm::vec3(-30, 30, -120));   // Point D
	backWall->setColor(glm::vec3(1, 1, 1));              // White
	backWall->setSpecularity(false);
	sceneObjects.push_back(backWall);

	// Ceiling - white
	Plane *ceiling = new Plane(glm::vec3(-30, 30, -40),      // Point A
		                      glm::vec3(30, 30, -40),       // Point B
		                      glm::vec3(30, 30, -120),      // Point C
		                      glm::vec3(-30, 30, -120));    // Point D
	ceiling->setColor(glm::vec3(1, 1, 1));               // White
	ceiling->setSpecularity(false);
	sceneObjects.push_back(ceiling);

	// Front wall (optional, since we're looking through it)
	// Uncomment if you want a complete box
	/*
	Plane *frontWall = new Plane(glm::vec3(-30, -15, -40),   // Point A
		                        glm::vec3(-30, 30, -40),    // Point B
		                        glm::vec3(30, 30, -40),     // Point C
		                        glm::vec3(30, -15, -40));   // Point D
	frontWall->setColor(glm::vec3(1, 1, 1));             // White
	frontWall->setSpecularity(false);
	sceneObjects.push_back(frontWall);
	*/

	// Add a sphere inside the Cornell box for testing
	Sphere *sphere1 = new Sphere(glm::vec3(0, 0, -80), 10.0);
	sphere1->setColor(glm::vec3(0, 0, 1));   // Blue
	sphere1->setReflectivity(true, 0.8);
	sceneObjects.push_back(sphere1);
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(20, 20);
	glutCreateWindow("Raytracing");

	glutDisplayFunc(display);
	initialize();

	glutMainLoop();
	return 0;
}
