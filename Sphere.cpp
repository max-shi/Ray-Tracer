#include "Sphere.h"
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>

float Sphere::intersect(glm::vec3 p0_world, glm::vec3 dir_world) {
	// Ray into object space
	glm::vec4 p0o4 = invTransform_ * glm::vec4(p0_world, 1.0f);
	glm::vec3  p0o  = glm::vec3(p0o4) / p0o4.w;

	glm::vec4 di4  = invTransform_ * glm::vec4(dir_world, 0.0f);
	float   tScale = glm::length(glm::vec3(di4));
	glm::vec3 dirO = glm::normalize(glm::vec3(di4));

	// Solve |p0o - center_|^2 = radius_^2
	glm::vec3 L = p0o - center_;
	float b = glm::dot(dirO, L);
	float c = glm::dot(L, L) - radius_*radius_;
	float disc = b*b - c;
	if (disc < 0.0f) return -1.0f;

	float sqrtD = std::sqrt(disc);
	float t1 = -b - sqrtD, t2 = -b + sqrtD;
	if (t2 < 0.0f) return -1.0f;
	float t_obj = (t1 > 0.0f ? t1 : t2) / tScale;

	// 3) Map hit back to world
	glm::vec3 hitO   = p0o + t_obj * dirO;
	glm::vec4 hitW4 = transform_ * glm::vec4(hitO, 1.0f);
	glm::vec3 hitW  = glm::vec3(hitW4) / hitW4.w;

	return t_obj;
}

glm::vec3 Sphere::normal(glm::vec3 hit_world) {
	// Hit into object space
	glm::vec4 hitO4 = invTransform_ * glm::vec4(hit_world, 1.0f);
	glm::vec3 hitO  = glm::vec3(hitO4) / hitO4.w;

	// normal in object space
	glm::vec3 nO = glm::normalize(hitO - center_);

	// back to world via transpose(inverse)
	glm::vec4 nW4 = normalTransform_ * glm::vec4(nO, 0.0f);
	return glm::normalize(glm::vec3(nW4));
}

void Sphere::translate(const glm::vec3& t) {
	transform_    = glm::translate(glm::mat4(1.0f), t) * transform_;
	invTransform_ = glm::inverse(transform_);
	normalTransform_ = glm::transpose(invTransform_);
}

void Sphere::rotate(float angleDeg, const glm::vec3& axis) {
	glm::mat4 toOrigin   = glm::translate(glm::mat4(1.0f), -center_);
	glm::mat4 rot        = glm::rotate(glm::mat4(1.0f), glm::radians(angleDeg), axis);
	glm::mat4 fromOrigin = glm::translate(glm::mat4(1.0f), center_);
	transform_ = fromOrigin * rot * toOrigin * transform_;
	invTransform_ = glm::inverse(transform_);
	normalTransform_ = glm::transpose(invTransform_);
}

void Sphere::scale(const glm::vec3& s) {
	glm::mat4 S = glm::scale(glm::mat4(1.0f), s);
	transform_ = S * transform_;
	invTransform_ = glm::inverse(transform_);
	normalTransform_ = glm::transpose(invTransform_);
}
