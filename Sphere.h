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

	// exactly the same three matrices as in Torus
	glm::mat4 transform_       = glm::mat4(1.0f);  // object → world
	glm::mat4 invTransform_    = glm::mat4(1.0f);  // world → object
	glm::mat4 normalTransform_ = glm::mat4(1.0f);  // transpose(invTransform_)

public:
	Sphere(glm::vec3 center, float radius)
	  : center_(center), radius_(radius)
	{
		// already identity
	}

	// build up your object-space transform
	void translate(const glm::vec3& t) {
		transform_    = glm::translate(glm::mat4(1.0f), t) * transform_;
		invTransform_ = glm::inverse(transform_);
		normalTransform_ = glm::transpose(invTransform_);
	}
	void rotate(float angleDeg, const glm::vec3& axis) {
		// rotate about the sphere’s center
		glm::mat4 toOrigin   = glm::translate(glm::mat4(1.0f), -center_);
		glm::mat4 rot        = glm::rotate(glm::mat4(1.0f), glm::radians(angleDeg), axis);
		glm::mat4 fromOrigin = glm::translate(glm::mat4(1.0f), center_);
		transform_ = fromOrigin * rot * toOrigin * transform_;
		invTransform_ = glm::inverse(transform_);
		normalTransform_ = glm::transpose(invTransform_);
	}
	void scale(const glm::vec3& s) {
		glm::mat4 S = glm::scale(glm::mat4(1.0f), s);
		transform_ = S * transform_;
		invTransform_ = glm::inverse(transform_);
		normalTransform_ = glm::transpose(invTransform_);
	}

	// the two virtuals
	float intersect(glm::vec3 p0, glm::vec3 dir) override;
	glm::vec3 normal(glm::vec3 pos) override;
	virtual ~Sphere() {}
};

#endif
