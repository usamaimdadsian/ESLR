#pragma once
#include<vector>
#include<math.h>

float norm_no_sqrt(float x, float y);

float ang_of_vec(float vx1, float vy1, float vx2, float vy2);

float norm_v3(float* v3);

float normM(float x, float y);

float det(float a, float b, float c, float d);

float dot(float a, float b, float c, float d);

float cos_vec(float* veci, float* vecj);

float cos_vec3(float* veci, float* vecj);

float dot_v2(float x1, float y1, float x2, float y2);

float point_2_line_dis_3D(float* x_0, float* x_1, float* x_2);

float dot_v3(float* v1, float* v2);

float point_2_line_dis(float x, float y, float* l);

float point_2_line_dis(float* pt, float* linef);

bool twoLines_intersec(float* pt1, float* pt2, float* tl1, float* tl2, float intersecratio);

void mult_3_3_1(float* o, float* e, float* res_3_1);

void mult_3_4_4(float* CM, float* x, float* res_3_1);

float norm_v2(float* v2);

void cross_v3(float* v1, float* v2, float* v3);

void cross_v2(float* line, float* linef);

void mult_3_3_3(float* o, float* e, float* res_3_1);

void M_divide_b(float* M, float* b, float* v);

void norm_by_v3(float* v);

void norm_by_v4(float* v4);

void Bresenham(int x1,
	int y1,
	int const x2,
	int const y2,
	std::vector<int>& xx, std::vector<int>& yy);

void Bresenham(int x1,
	int y1,
	int const x2,
	int const y2,
	int* xx, int* yy, int& xy_size);

float max_2(float a, float b);

float min_2(float a, float b);

float max_4(float a, float b, float c, float d);

bool ID_in_array(int* ids, int t_size, int num);

bool ID_in_array(std::vector<short> ids, int num);

