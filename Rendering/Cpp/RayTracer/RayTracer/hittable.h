#ifndef HITTABLEH
#define HITTABLEH

#include "ray.h"

class material;

// Struct used to record the hit position
struct hit_record {
	float t;
	vec3 p;
	vec3 normal;
	material *mat_ptr; // Point to class material
};

// Abstract function isHit : hit_record contains the hit position
class hittable {
public:
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) 
		const = 0;
};

class flip_normals : public hittable {
public:
	flip_normals(hittable *p) : ptr(p) {}

	virtual bool hit(
		const ray& r, float t_min, float t_max, hit_record& rec) const {

		if (ptr->hit(r, t_min, t_max, rec)) {
			rec.normal = -rec.normal;
			return true;
		}
		else
			return false;
	}

	hittable *ptr;
};

#endif