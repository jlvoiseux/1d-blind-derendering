// RayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include<opencv2/core/core.hpp>
#include<opencv2/photo/photo.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <math.h>  
#include <omp.h>
using namespace std;
using namespace cv;


#include "vec3.h"
#include "ray.h"
#include "float.h"
#include "hittable_list.h"
#include "sphere.h"
#include "random.h"
#include "camera.h"
#include "material.h"
#include "rect.h"
#include "util.h"

// 0 : Not conv, 1 : conv
int check_type(const ray& r, hittable *world) {
	hit_record rec;
	if (world->hit(r, 0.001, FLT_MAX, rec)) {
		if (rec.mat_ptr->get_type() == 3)
		{
			return 1;
		}
		else {
			return 0;
		}
	}
}

// Define the color for the hittable passed in argument (the hittable can be a list of hittable)
vec3 color(const ray& r, hittable *world, int depth) {
	hit_record rec;
	// if hit, return shaded sphere
	if (world->hit(r, 0.001, FLT_MAX, rec)) {
		ray scattered;
		vec3 attenuation;
		vec3 emitted = rec.mat_ptr->emitted(rec.p);
		if (depth == 0 && rec.mat_ptr->get_type() == 6)
		{
			return(color(ray(rec.p, r.B), world, 0));
		}
		else if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {			
			return emitted + attenuation * color(scattered, world, depth + 1);
		}
		else {
			return emitted;
		}
	}
	// else, return background
	else {
		// Blue Sky
		vec3 unit_direction = unit_vector(r.direction());
		float t = 0.5*(unit_direction.y() + 1.0);
		//return (1.0 - t)*vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
		return vec3(1.0, 1.0, 1.0);
	}
}

//hittable *cornell_box() {
//	int n = 500;
//	hittable **list = new hittable*[n + 1];
//	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
//	int i = 1;
//	for (int a = -11; a < 11; a++) {
//		for (int b = -11; b < 11; b++) {
//			float choose_mat = random_double();
//			vec3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());
//			if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
//				if (choose_mat < 0.8) {  // diffuse
//					list[i++] = new sphere(center, 0.2,
//						new lambertian(vec3(random_double()*random_double(),
//							random_double()*random_double(),
//							random_double()*random_double())
//						)
//					);
//				}
//				else if (choose_mat < 0.95) { // metal
//					list[i++] = new sphere(center, 0.2,
//						new metal(vec3(0.5*(1 + random_double()),
//							0.5*(1 + random_double()),
//							0.5*(1 + random_double())),
//							0.5*random_double()));
//				}
//				else {  // glass
//					list[i++] = new sphere(center, 0.2, new dielectric(1.5));
//				}
//			}
//		}
//	}
//
//	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
//	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
//	list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
//
//	return new hittable_list(list, i);
//}

hittable *cornell_box() {
	hittable **list = new hittable*[6];
	int i = 0;
	material *red = new lambertian(vec3(0.65, 0.05, 0.05));
	material *white = new lambertian(vec3(0.73, 0.73, 0.73));
	material *green = new lambertian(vec3(0.12, 0.45, 0.15));
	material *light = new diffuse_light(vec3(15, 15, 15));
	material *blue_light = new diffuse_light(vec3(0, 0, 15));
	material *blue_light_hidden = new diffuse_light_hidden(vec3(0, 0, 30));
	material *red_light = new diffuse_light(vec3(30, 0, 0));
	material *red_light_hidden = new diffuse_light_hidden(vec3(30, 0, 0));
	material *mirror = new metal(vec3(0.8, 0.8, 0.8), 0.0f);
	material *conv_white = new conv(vec3(0.73, 0.73, 0.73));

	list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 1110, green));
	list[i++] = new yz_rect(0, 555, 0, 555, -555, red);
	//list[i++] = new xz_rect(213, 343, 227, 332, 554, light);
	list[i++] = new flip_normals(new xz_rect(-555, 1110, 0, 555, 555, white));
	list[i++] = new xz_rect(-555, 1110, 0, 555, 0, white);
	list[i++] = new flip_normals(new xy_rect(-555, 1110, -555, 1110, 555, mirror));
	list[i++] = new xy_rect(0, 20, 200, 240, 550, blue_light);
	list[i++] = new xy_rect(-20, -0, 200, 240, 550, red_light);

	return new hittable_list(list, i);
}

// Draw a PPM image of dimensions nx and ny
int main() {

	int ns = 20; // Number of samples for each pixel : (aliasing) Min value : 1
	bool do_rt_update = true;
	bool display_progress = true;

	bool compute_incoming_luminance = true;
	bool convolve_wall = false;
	bool redraw_wall = false;

	const int nx = 1280;
	const int ny = 720;

	const int thread_number = omp_get_num_procs();

	//Array allocation
	Mat image_array = Mat::zeros(ny, nx, CV_8UC3);
	Mat seen_by_wall = Mat::zeros(ny, nx, CV_8UC3);

	// Define geometry
	hittable *world = cornell_box();

	// Camera coordinates
	vec3 lower_left_corner(-8.0f, -4.5f, -1.0f);
	vec3 origin(0.0f, 0.0f, 0.0f);
	// Scale
	vec3 horizontal(16.0f, 0.0f, 0.0f);
	vec3 vertical(0.0f, 9.0f, 0.0f);

	vec3 lookfrom(278, 278, -800);
	vec3 lookat(278, 278, 0);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	float vfov = 40.0;

	camera cam(lookfrom, lookat, vec3(0, 1, 0), vfov, float(nx) / float(ny), aperture, dist_to_focus);

	int displayed_progress = 0;
	if (display_progress) {
	std:cout << "  " << displayed_progress << " %" << "\r";
	}


	// 1. Create image and extract wall incoming luminance
	if (compute_incoming_luminance)
	{
		omp_set_num_threads(thread_number);
		#pragma omp parallel for schedule(dynamic) 
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				vec3 pixel_color(0, 0, 0);
				// Aliasing appears due to geometry borders. We send ns random rays per pixel, some rays will hit the geometry, others not, and we average.
				bool wall_flag = false;
				for (int s = 0; s < ns; s++) {
					float u = float(i + random_double()) / float(nx);
					float v = float(j + random_double()) / float(ny);
					ray r = cam.get_ray(u, v);
					if (check_type(r, world) == 1)
					{
						pixel_color += color(r, world, 0);
						wall_flag = true;
						break;
					}
					pixel_color += color(r, world, 0);
				}
				if (wall_flag)
				{
					pixel_color = vec3(sqrt(pixel_color[0] <= 1.0f ? pixel_color[0] : 1.0f), sqrt(pixel_color[1] <= 1.0f ? pixel_color[1] : 1.0f), sqrt(pixel_color[2] <= 1.0f ? pixel_color[2] : 1.0f)); // gamma_correction
					int ir = int(255.99f*pixel_color[0]);
					int ig = int(255.99f*pixel_color[1]);
					int ib = int(255.99f*pixel_color[2]);
					image_array.at<Vec3b>(ny - j - 1, i) = *(new Vec3b(ib, ig, ir));
					seen_by_wall.at<Vec3b>(ny - j - 1, i) = *(new Vec3b(ib, ig, ir));
				}
				else
				{
					pixel_color /= float(ns);
					pixel_color = vec3(sqrt(pixel_color[0] <= 1.0f ? pixel_color[0] : 1.0f), sqrt(pixel_color[1] <= 1.0f ? pixel_color[1] : 1.0f), sqrt(pixel_color[2] <= 1.0f ? pixel_color[2] : 1.0f)); // gamma_correction
					int ir = int(255.99f*pixel_color[0]);
					int ig = int(255.99f*pixel_color[1]);
					int ib = int(255.99f*pixel_color[2]);
					image_array.at<Vec3b>(ny - j - 1, i) = *(new Vec3b(ib, ig, ir));
				}

				float current_progress = float(j) / float(ny - 1) * 100;
				if (omp_get_thread_num() == 0) {
					if (floor(current_progress) > displayed_progress) {
						displayed_progress = floor(current_progress);
						if (display_progress) {
							cout << "  " << displayed_progress << " %" << "\r";
						}
					}
					if (do_rt_update) {
						cv::imshow("Image", image_array);
						cv::waitKey(1);
					}
				}
			}
		}

		cv::imshow("Image", image_array);
		cv::imwrite("rendered.png", image_array);
		cv::imwrite("wall.png", seen_by_wall);
		cv::waitKey(0);
	}
	

	// 2. Read incoming luminance, fill holes and convolve
	if (convolve_wall)
	{
		// Read image
		Mat wall_luminance = cv::imread("wall.png");
		// Convert to grayscale to identify black pixels
		Mat wall_luminance_gray;
		cv::cvtColor(wall_luminance, wall_luminance_gray, COLOR_BGR2GRAY);
		// Build mask
		Mat mask = cv::Mat::zeros(wall_luminance_gray.size(), CV_8UC1);
		mask = (wall_luminance_gray == 0);
		// Fill holes
		Mat wall_luminance_filled;
		cv::inpaint(wall_luminance, mask, wall_luminance_filled, 5, cv::INPAINT_TELEA);

		// Create kernel
		Mat X, Y, cHg1, cHg2, cHg3, cHg4, Hg, H;
		int range1 = 100;
		int range2 = 100;
		meshgridTest(cv::Range(-(range1-1)/2, (range1-1)/2), cv::Range(-(range2-1)/2, (range2-1)/2), X, Y);
		X.convertTo(X, CV_32F);
		Y.convertTo(Y, CV_32F);
		double minX, maxX, minY, maxY;
		cv::minMaxLoc(X, &minX, &maxX);
		cv::minMaxLoc(Y, &minY, &maxY);
		X = M_PI * X / maxX - X / maxX;
		Y = M_PI * Y / maxY - Y / maxY;
		cHg1 = mCos(X);
		cHg2 = mCos(Y);
		cv::pow(cHg1, 2, cHg3);
		cv::pow(cHg2, 2, cHg4);
		Hg = cHg3 + cHg4;
		H = Hg / cv::sum(Hg)[0];	
		// Convolve each channel
		Mat final_wall;
		cv::filter2D(wall_luminance_filled, final_wall, wall_luminance_filled.depth(), H);
		cv::imshow("Image", final_wall);
		cv::imwrite("wall_convolved.png", final_wall);
		cv::waitKey(0);
	}

	if (redraw_wall)
	{
		Mat final_values = Mat::zeros(ny, nx, CV_8UC3);
		Mat world_values = cv::imread("rendered.png");
		Mat wall_values = cv::imread("dec.png");
		omp_set_num_threads(thread_number);
		#pragma omp parallel for schedule(dynamic) 
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				bool wall_flag = false;
				for (int s = 0; s < ns; s++) {
					float u = float(i + random_double()) / float(nx);
					float v = float(j + random_double()) / float(ny);
					ray r = cam.get_ray(u, v);
					if (check_type(r, world) == 1)
					{
						wall_flag = true;
					}
					break;

				}
				if (wall_flag)
				{
					final_values.at<Vec3b>(ny - j - 1, i) = wall_values.at<Vec3b>(ny - j - 1, i);
				}
				else
				{
					final_values.at<Vec3b>(ny - j - 1, i) = world_values.at<Vec3b>(ny - j - 1, i);
				}

				float current_progress = float(j) / float(ny - 1) * 100;
				if (omp_get_thread_num() == 0) {
					if (floor(current_progress) > displayed_progress) {
						displayed_progress = floor(current_progress);
						if (display_progress) {
							cout << "  " << displayed_progress << " %" << "\r";
						}
					}
					if (do_rt_update) {
						cv::imshow("Image", final_values);
						cv::waitKey(1);
					}
				}
			}
		}

		cv::imshow("Image", final_values);
		cv::imwrite("rendered_final.png", final_values);
		cv::waitKey(0);
	}
	

	return 0;
}


