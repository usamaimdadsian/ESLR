#pragma once
#include<vector>
#include<math.h>
#include "BasicMath.h"

float ang_of_vec(float vx1, float vy1, float vx2, float vy2)
{
	float dot = vx1 * vx2 + vy1 * vy2; // dot product
	float det = vx1 * vy2 - vy1 * vx2; // determinant
	return atan2(det, dot);
}

float norm_v3(float* v3)
{
	return sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
}

float normM(float x, float y)
{
	return sqrt(x * x + y * y);
}

float norm_no_sqrt(float x, float y)
{
	return (x * x + y * y);
}

float det(float a, float b, float c, float d)
{
	return a * d - b * c;
}

float dot(float a, float b, float c, float d)
{
	return a * c + b * d;
}

float cos_vec3(float* veci, float* vecj)
{
	return (veci[0] * vecj[0] + veci[1] * vecj[1] + veci[2] * vecj[2]) /
		norm_v3(veci) /
		norm_v3(vecj);
}

float cos_vec(float* veci, float* vecj)
{

	return (veci[0] * vecj[0] + veci[1] * vecj[1]) /
		(sqrt(veci[0] * veci[0] + veci[1] * veci[1]) *
			sqrt(vecj[0] * vecj[0] + vecj[1] * vecj[1]));
}

float dot_v2(float x1, float y1, float x2, float y2)
{

	return x1 * x2 + y1 * y2;
}

float point_2_line_dis(float x, float y, float* l)
{
	return abs(x * l[0] + y * l[1] + l[2]) / sqrt(l[0] * l[0] + l[1] * l[1]);
}

// 3D point to line distance
//d   =   (|(x_2-x_1)x(x_1-x_0)|)/(|x_2-x_1|)
//https://stackoverflow.com/questions/19878441/point-line-distance-calculation
float point_2_line_dis_3D(float* x_0, float* x_1, float* x_2)
{
	float x2_x1[3], x1_x0[3], cv3[3];

	x2_x1[0] = x_2[0] - x_1[0];
	x2_x1[1] = x_2[1] - x_1[1];
	x2_x1[2] = x_2[2] - x_1[2];

	x1_x0[0] = x_1[0] - x_0[0];
	x1_x0[1] = x_1[1] - x_0[1];
	x1_x0[2] = x_1[2] - x_0[2];

	cross_v3(x2_x1, x1_x0, cv3);

	return norm_v3(cv3) / norm_v3(x2_x1);
}

//2D point to line distance
float point_2_line_dis(float* pt, float* linef)
{
	return abs(pt[0] * linef[0] + pt[1] * linef[1] + linef[2]) /
		sqrt(linef[0] * linef[0] + linef[1] * linef[1]);
}

bool twoLines_intersec(float* pt1, float* pt2, float* tl1, float* tl2, float intersecratio)
{
	bool dottl_p1, dottl_p2, dotp_tl1, dotp_tl2;

	dottl_p1 = dot_v2(tl1[0] - pt1[0], tl1[1] - pt1[1], tl2[0] - pt1[0], tl2[1] - pt1[1]) <= 0;
	dottl_p2 = dot_v2(tl1[0] - pt2[0], tl1[1] - pt2[1], tl2[0] - pt2[0], tl2[1] - pt2[1]) <= 0;
	dotp_tl1 = dot_v2(pt1[0] - tl1[0], pt1[1] - tl1[1], pt2[0] - tl1[0], pt2[1] - tl1[1]) <= 0;
	dotp_tl2 = dot_v2(pt1[0] - tl2[0], pt1[1] - tl2[1], pt2[0] - tl2[0], pt2[1] - tl2[1]) <= 0;

	if (!(dottl_p1 || dottl_p2 || dotp_tl1 || dotp_tl2))
		return false;

	if ((dottl_p1 && dottl_p2) || (dotp_tl1 && dotp_tl2))
		return true;

	// pt1 is included
	if (dottl_p1)
	{
		float t_len;
		float l1_len = sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) + (pt1[1] - pt2[1]) * (pt1[1] - pt2[1]));
		float l2_len = sqrt((tl1[0] - tl2[0]) * (tl1[0] - tl2[0]) + (tl1[1] - tl2[1]) * (tl1[1] - tl2[1]));
		t_len = sqrt((pt1[0] - tl2[0]) * (pt1[0] - tl2[0]) + (pt1[1] - tl2[1]) * (pt1[1] - tl2[1]));

		if (t_len / l1_len < intersecratio && t_len / l2_len < intersecratio)
			return 0;
	}

	if (dottl_p2)
	{
		float t_len;
		float l1_len = sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) + (pt1[1] - pt2[1]) * (pt1[1] - pt2[1]));
		float l2_len = sqrt((tl1[0] - tl2[0]) * (tl1[0] - tl2[0]) + (tl1[1] - tl2[1]) * (tl1[1] - tl2[1]));

		t_len = sqrt((pt2[0] - tl1[0]) * (pt2[0] - tl1[0]) + (pt2[1] - tl1[1]) * (pt2[1] - tl1[1]));

		if (t_len / l1_len < intersecratio && t_len / l2_len < intersecratio)
			return 0;
	}

	return true;
}

void mult_3_3_1(float* o, float* e, float* res_3_1)
{
	res_3_1[0] = e[0] * o[0] + e[1] * o[1] + e[2] * o[2];
	res_3_1[1] = e[0] * o[3] + e[1] * o[4] + e[2] * o[5];
	res_3_1[2] = e[0] * o[6] + e[1] * o[7] + e[2] * o[8];
}

void mult_3_4_4(float* CM, float* x, float* res_3_1)
{
	res_3_1[0] = CM[3] + CM[0] * x[0] + CM[1] * x[1] + CM[2] * x[2];
	res_3_1[1] = CM[7] + CM[4] * x[0] + CM[5] * x[1] + CM[6] * x[2];
	res_3_1[2] = CM[11] + CM[8] * x[0] + CM[9] * x[1] + CM[10] * x[2];
}

float norm_v2(float* v2)
{
	return sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
}

void cross_v3(float* v1, float* v2, float* v3)
{
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

}

void cross_v2(float* line, float* linef)
{
	linef[0] = line[1] - line[3];
	linef[1] = line[2] - line[0];
	linef[2] = line[0] * line[3] -
		line[1] * line[2];
}

bool ID_in_array(int* ids, int t_size, int num)
{
	for (int i = 0; i < t_size; i++)
	{
		if (num == ids[i])
			return true;
	}
	return false;
}

bool ID_in_array(std::vector<short> ids, int num)
{
	for (int i = 0; i < ids.size(); i++)
	{
		if (num == ids[i])
			return true;
	}
	return false;
}

void mult_3_3_3(float* o, float* e, float* res_3_1)
{
	res_3_1[0] = e[0] * o[0] + e[1] * o[1] + e[2] * o[2];
	res_3_1[1] = e[0] * o[3] + e[1] * o[4] + e[2] * o[5];
	res_3_1[2] = e[0] * o[6] + e[1] * o[7] + e[2] * o[8];
}

void M_divide_b(float* M, float* b, float* v)
{
	v[0] = -(M[1] * M[5] * b[2] - M[2] * M[4] * b[2] -
		M[1] * M[8] * b[1] + M[2] * M[7] * b[1] +
		M[4] * M[8] * b[0] - M[5] * M[7] * b[0]) /
		(M[0] * M[4] * M[8] - M[0] * M[5] * M[7] -
			M[1] * M[3] * M[8] + M[1] * M[5] * M[6] +
			M[2] * M[3] * M[7] - M[2] * M[4] * M[6]);

	v[1] = (M[0] * M[5] * b[2] - M[2] * M[3] * b[2] -
		M[0] * M[8] * b[1] + M[2] * M[6] * b[1] +
		M[3] * M[8] * b[0] - M[5] * M[6] * b[0]) /
		(M[0] * M[4] * M[8] - M[0] * M[5] * M[7] -
			M[1] * M[3] * M[8] + M[1] * M[5] * M[6] +
			M[2] * M[3] * M[7] - M[2] * M[4] * M[6]);

	v[2] = -(M[0] * M[4] * b[2] - M[1] * M[3] * b[2] -
		M[0] * M[7] * b[1] + M[1] * M[6] * b[1] +
		M[3] * M[7] * b[0] - M[4] * M[6] * b[0]) /
		(M[0] * M[4] * M[8] - M[0] * M[5] * M[7] -
			M[1] * M[3] * M[8] + M[1] * M[5] * M[6] +
			M[2] * M[3] * M[7] - M[2] * M[4] * M[6]);

}

float dot_v3(float* v1, float* v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void norm_by_v4(float* v4)
{
	v4[0] = v4[0] / v4[3];
	v4[1] = v4[1] / v4[3];
	v4[2] = v4[2] / v4[3];
	v4[3] = 1;
}

void norm_by_v3(float* v)
{
	v[0] = v[0] / v[2];
	v[1] = v[1] / v[2];
	v[2] = 1;
}

void Bresenham(int x1,
	int y1,
	int const x2,
	int const y2,
	std::vector<int>& xx, std::vector<int>& yy)
{
	xx.clear();
	yy.clear();
	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;
	xx.push_back(x1);
	yy.push_back(y1);
	//plot(x1, y1);

	if (delta_x >= delta_y)
	{
		// error may go below zero
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2)
		{
			// reduce error, while taking into account the corner case of error == 0
			if ((error > 0) || (!error && (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}
			// else do nothing

			error += delta_y;
			x1 += ix;

			//plot(x1, y1);
			xx.push_back(x1);
			yy.push_back(y1);
		}
	}
	else
	{
		// error may go below zero
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2)
		{
			// reduce error, while taking into account the corner case of error == 0
			if ((error > 0) || (!error && (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			// else do nothing

			error += delta_x;
			y1 += iy;

			xx.push_back(x1);
			yy.push_back(y1);
		}
	}
}

void Bresenham(int x1,
	int y1,
	int const x2,
	int const y2,
	int* xx, int* yy, int& xy_size)
{
	xy_size = 0;

	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;
	xx[xy_size] = x1;
	yy[xy_size] = y1;
	xy_size++;

	if (delta_x >= delta_y)
	{
		// error may go below zero
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2)
		{
			// reduce error, while taking into account the corner case of error == 0
			if ((error > 0) || (!error && (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}
			// else do nothing

			error += delta_y;
			x1 += ix;

			xx[xy_size] = x1;
			yy[xy_size] = y1;
			xy_size++;
		}
	}
	else
	{
		// error may go below zero
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2)
		{
			// reduce error, while taking into account the corner case of error == 0
			if ((error > 0) || (!error && (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			// else do nothing
			error += delta_x;
			y1 += iy;

			xx[xy_size] = x1;
			yy[xy_size] = y1;
			xy_size++;
		}
	}
}

float max_2(float a, float b)
{
	if (a < b)
		return b;
	else
		return a;
}

float min_2(float a, float b)
{
	if (a > b)
		return b;
	else
		return a;
}

float max_4(float a, float b, float c, float d)
{
	return max_2(max_2(a, b), max_2(c, d));
}

