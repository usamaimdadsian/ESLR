#include<vector>
#include "math.h"
#include <fstream>
#include <iostream>

# include"BasicMath.h"

#include"SingleImage.h"
struct LP
{
	int i, j;
	float cx, cy, x1, y1, x2, y2;
};

void pt2FarestEnd(float x, float y, float* lines, int i, int size, float& nx, float& ny)
{
	float x1, y1, x2, y2;

	int ind = i * 7;
	x1 = lines[ind];
	y1 = lines[ind + 1];

	x2 = lines[ind + 2];
	y2 = lines[ind + 3];

	if (norm_no_sqrt(x1 - x, y1 - y) > norm_no_sqrt(x2 - x, y2 - y))
	{
		nx = x1;
		ny = y1;
	}
	else
	{
		nx = x2;
		ny = y2;
	}

}

void pt2NearestEnd(float x, float y, float* lines, int i, int size, float& nx, float& ny)
{
	float x1, y1, x2, y2;

	int ind;
	ind = 7 * i;
	x1 = lines[ind];
	y1 = lines[ind + 1];

	x2 = lines[ind + 2];
	y2 = lines[ind + 3];

	if (norm_no_sqrt(x1 - x, y1 - y) < norm_no_sqrt(x2 - x, y2 - y))
	{
		nx = x1;
		ny = y1;
	}
	else
	{
		nx = x2;
		ny = y2;
	}
}

bool LineLineIntersect(float x1, float y1, //Line 1 start
	float x2, float y2, //Line 1 end
	float x3, float y3, //Line 2 start
	float x4, float y4, //Line 2 end
	float& ixOut, float& iyOut) //Output 
{
	//http://mathworld.wolfram.com/Line-LineIntersection.html

	float x1mx2 = x1 - x2;
	float x3mx4 = x3 - x4;
	float y1my2 = y1 - y2;
	float y3my4 = y3 - y4;

	float denom = det(x1mx2, y1my2, x3mx4, y3my4);
	if (denom == 0.0)//Lines don't seem to cross
	{
		ixOut = 0;
		iyOut = 0;
		return false;
	}

	float detL1 = det(x1, y1, x2, y2);
	float detL2 = det(x3, y3, x4, y4);
	float xnom = det(detL1, x1mx2, detL2, x3mx4);
	float ynom = det(detL1, y1my2, detL2, y3my4);

	ixOut = xnom / denom;
	iyOut = ynom / denom;
	return true; //All OK
}

bool pt2LineDis(float x, float y, float* lines, int index, float dist_2)
{

	int ind = index * 7;

	if (norm_no_sqrt(lines[ind] - x, lines[ind + 1] - y) <= dist_2)
		return true;

	if (norm_no_sqrt(lines[ind + 2] - x, lines[ind + 3] - y) <= dist_2)
		return true;

	return false;
}

int boolPtinLine(float x, float y, float* lines, int i, int j, int size)
{
	// the i-th line
	float ax, ay, bx, by;

	int ind = 7 * i;
	ax = lines[ind] - x;
	ay = lines[ind + 1] - y;
	bx = lines[ind + 2] - x;
	by = lines[ind + 3] - y;

	if (ax * bx + ay * by < 0)
		return i;

	// the j-th line
	ind = 7 * j;
	ax = lines[ind] - x;
	ay = lines[ind + 1] - y;
	bx = lines[ind + 2] - x;
	by = lines[ind + 3] - y;

	if (ax * bx + ay * by < 0)
		return j;

	return 0;
}

void crossPT(float* lines, int size, float costhre, float dist, std::vector<LP>& lps)
{
	int i, j, k;
	float veca[2];
	float vecb[2];
	LP lp;

	float norma, normb, cost, cx, cy;
	float dist_2 = dist * dist;
	float x3, y3, vx3, vy3;

	int ptri, ptrj;

	for (i = 0; i < size; i++)
	{
		ptri = i * 7;

		veca[0] = lines[ptri] - lines[ptri + 2];
		veca[1] = lines[ptri + 1] - lines[ptri + 3];
		norma = normM(veca[0], veca[1]);

		for (j = i + 1; j < size; j++)
		{

			ptrj = j * 7;
			vecb[0] = lines[ptrj] - lines[ptrj + 2];
			vecb[1] = lines[ptrj + 1] - lines[ptrj + 3];

			normb = normM(vecb[0], vecb[1]);
			cost = dot(veca[0], veca[1], vecb[0], vecb[1]) / (norma * normb);

			if (abs(cost) >= costhre)
				continue;

			// check if there exist line intersection 
			if (!
				LineLineIntersect(lines[ptri], lines[ptri + 1], //Line 1 start
					lines[ptri + 2], lines[ptri + 3], //Line 1 end
					lines[ptrj], lines[ptrj + 1], //Line 1 start
					lines[ptrj + 2], lines[ptrj + 3], //Line 2 end
					cx, cy)
				)
				continue;

			// check if the intersection is in one line segment
			k = boolPtinLine(cx, cy, lines, i, j, size);
			if (k)
			{
				// take the endpoint of k-th line
				lp.i = i;
				lp.j = j;
				lp.cx = cx;
				lp.cy = cy;

				// calculate v3;
				if (k == i)
					k = j;
				else
					k = i;

				pt2NearestEnd(cx, cy, lines, k, size, x3, y3);
				vx3 = x3 - cx;
				vy3 = y3 - cy;

				if (normM(vx3, vy3) >= dist)
					continue;

				pt2FarestEnd(cx, cy, lines, k, size, x3, y3);
				lp.x1 = x3;
				lp.y1 = y3;

				if (k == i)
					k = j;
				else
					k = i;

				pt2FarestEnd(cx, cy, lines, k, size, x3, y3);
				lp.x2 = x3;
				lp.y2 = y3;

				lps.push_back(lp);
				continue;
			}

			//check the distance to endpoints
			if (!(pt2LineDis(cx, cy, lines, i, dist_2) && pt2LineDis(cx, cy, lines, j, dist_2)))
				continue;
			//then add it 

			lp.i = i;
			lp.j = j;
			lp.cx = cx;
			lp.cy = cy;

			pt2FarestEnd(cx, cy, lines, i, size, x3, y3);
			lp.x1 = x3;
			lp.y1 = y3;
			pt2FarestEnd(cx, cy, lines, j, size, x3, y3);
			lp.x2 = x3;
			lp.y2 = y3;
			lps.push_back(lp);
		}

	}

}

void callCrossPt(cv::Mat& lp_Mf, float* lines, int size, float costhre, float dist)
{

	std::vector<LP>lps;
	crossPT(lines, size, costhre, dist, lps);

	int lpsize = lps.size();

	cv::Mat lp_Mf_ = cv::Mat(lpsize, 8, CV_32FC1);
	float* outlp = (float*)lp_Mf_.data;

	float a[2], b[2], vec1[2], vec2[2];

	int id1, id2;
	float* l1, * l2;

	int ind = 0;

	for (int i = 0; i < lpsize; i++)
	{
		ind = i * 8;
		outlp[ind + 0] = lps.at(i).i;
		outlp[ind + 1] = lps.at(i).j;
		outlp[ind + 2] = lps.at(i).cx;
		outlp[ind + 3] = lps.at(i).cy;

		a[0] = lps.at(i).x1 - lps.at(i).cx;
		a[1] = lps.at(i).y1 - lps.at(i).cy;

		b[0] = lps.at(i).x2 - lps.at(i).cx;
		b[1] = lps.at(i).y2 - lps.at(i).cy;

		outlp[ind + 4] = acos(cos_vec(a, b));

		id1 = lps.at(i).i;
		id2 = lps.at(i).j;

		l1 = lines + id1 * 7;
		l2 = lines + id2 * 7;

		vec1[0] = l1[2] - l1[0];
		vec1[1] = l1[3] - l1[1];

		vec2[0] = l2[2] - l2[0];
		vec2[1] = l2[3] - l2[1];

		outlp[ind + 5] = ang_of_vec(vec1[0], vec1[1], vec2[0], vec2[1]);


		if (vec1[0] * a[0] + vec1[1] * a[1] > 0)
			outlp[ind + 6] = 1;
		else
			outlp[ind + 6] = 0;

		if (vec2[0] * b[0] + vec2[1] * b[1] > 0)
			outlp[ind + 7] = 1;
		else
			outlp[ind + 7] = 0;
	}

	lp_Mf_.copyTo(lp_Mf);
}