#ifndef H_SPHERE
#define H_SPHERE

#include "SceneObject.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>

class Sphere : public SceneObject {
protected:
	glm::vec3 center_;
	float     radius_;

	glm::mat4 transform_       = glm::mat4(1.0f);  // object → world
	glm::mat4 invTransform_    = glm::mat4(1.0f);  // world → object
	glm::mat4 normalTransform_ = glm::mat4(1.0f);  // transpose(invTransform_)

public:
	Sphere(glm::vec3 center, float radius)
	  : center_(center), radius_(radius)
	{
	}

	void translate(const glm::vec3& t);
	void rotate(float angleDeg, const glm::vec3& axis);
	void scale(const glm::vec3& s);

	// the two virtuals
	float intersect(glm::vec3 p0, glm::vec3 dir) override;
	glm::vec3 normal(glm::vec3 pos) override;
	virtual ~Sphere() {}
};

#endif
