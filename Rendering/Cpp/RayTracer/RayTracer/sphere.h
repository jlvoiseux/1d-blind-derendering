#ifndef SPHEREH
#define SPHEREH

#include "hittable.h"

class sphere : public hittable {
public:
	sphere() {}
	sphere(vec3 cen, float r, material *m) : center(cen), radius(r), mat_ptr(m) {};
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	vec3 center;
	float radius;
	material *mat_ptr;
};

// Check if ray hits sphere
// Sphere equation : (x−Cx)2+(y−Cy)2+(z−Cz)2=R2
// Dot Product : dot((p−C),(p−C))=(x−Cx)2+(y−Cy)2+(z−Cz)2
// Sphere equation (vector form) : dot((p−C),(p−C))=R2
// Ray hits sphere (1) : dot((p(t)−C),(p(t)−C))=R2
// Ray hits sphere (2) : dot((A+t*B−C),(A+t*B−C))=R2
// Explicit form : t2dot(B,B)+2tdot(B,A−C)+dot(A−C,A−C)−R2=0
bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	vec3 oc = r.origin() - center;
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - radius * radius;
	float discriminant = b * b - a * c;
	if (discriminant > 0) {
		float temp = (-b - sqrt(discriminant)) / a;
		// If the closest hit point is within the bounds...
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mat_ptr;
			return true;
		}
		//Else, try with the farthest hit point
		temp = (-b + sqrt(discriminant)) / a;
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			// Sphere : normal = hitpoint - sphere center (normalized)
			rec.normal = (rec.p - center) / radius;
			rec.mat_ptr = mat_ptr;
			return true;
		}
	}
	return false;
}


#endif