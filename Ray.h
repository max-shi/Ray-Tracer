//==================================================
// COSC363 Ray Tracer
// The ray class (header)
//==================================================

#ifndef H_RAY
#define H_RAY
#include <glm/glm.hpp>
#include <vector>
#include "SceneObject.h"

class Ray
{

public:
	glm::vec3 p0 = glm::vec3(0);
	glm::vec3 dir = glm::vec3(0,0,-1);
	glm::vec3 hit = glm::vec3(0);
	int index = -1;
	float dist = 0;

	Ray() {}


	Ray(glm::vec3 source, glm::vec3 direction)
	{
		const float RSTEP = 0.005f;
		p0 = source;
		dir = glm::normalize(direction);
		p0 = p0 + RSTEP * dir;
	}

	void closestPt(std::vector<SceneObject*>& sceneObjects);

};
#endif
