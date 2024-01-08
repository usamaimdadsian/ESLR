#include "PairMatch.h"
#include "BasicMath.h"
#include "Parameters.h"
# define PI 3.14159265
# define PI_6 0.523598
# define PI_2 1.570796

float angDiff(float a, float b)
{
	float d = (b - a);

	// constrain d in[0, 2 * pi]
	float two_pi = 2 * PI;
	d = d - floor(d / two_pi) * two_pi;

	// constrain d in[-pi, pi]
	if (d > PI)
		d = d - two_pi;

	return d;
}

int checkAng(float* a1, float* a2, bool& bvalid, bool& bchange)
{
	bvalid = 0;
	bchange = 0;

	if (abs(a1[0] - a2[0]) > PI_6)
	{
		bvalid = 0;
		return 0;

	}
	else if (angDiff(a1[1], a2[1]) > PI_2)
	{
		if (a1[2] == a2[3] && a1[3] == a2[2])
		{
			bvalid = 1;
			bchange = 1;
			return 0;
		}
	}
	else
	{
		if (a1[2] == a2[2] && a1[3] == a2[3])
		{
			bvalid = 1;
			return 0;
		}
	}

	return 0;

}

void depthRange(float* C, float* CM, float* CN,
	float* itpm, float* principal_r,
	float mindepth, float maxdepth,
	float& x1, float& y1, float& x2, float& y2)
{
	//pv = M\itpm;
	float pv[3];
	float pt2d1[3];
	float pt3d1[4];
	float pt2d2[3];
	float pt3d2[4];
	float vec12[2];

	float x, y, z;
	x = itpm[0];
	y = itpm[1];
	z = itpm[2];

	pv[0] = (CM[5] * CM[10] * x - CM[6] * CM[9] * x - CM[1] * CM[10] * y + CM[2] * CM[9] * y + CM[1] * CM[6] * z - CM[2] * CM[5] * z) /
		(CM[0] * CM[5] * CM[10] - CM[0] * CM[6] * CM[9] - CM[1] * CM[4] * CM[10] + CM[1] * CM[6] * CM[8] + CM[2] * CM[4] * CM[9] - CM[2] * CM[5] * CM[8]);

	pv[1] = -(CM[4] * CM[10] * x - CM[6] * CM[8] * x - CM[0] * CM[10] * y + CM[2] * CM[8] * y + CM[0] * CM[6] * z - CM[2] * CM[4] * z) /
		(CM[0] * CM[5] * CM[10] - CM[0] * CM[6] * CM[9] - CM[1] * CM[4] * CM[10] + CM[1] * CM[6] * CM[8] + CM[2] * CM[4] * CM[9] - CM[2] * CM[5] * CM[8]);

	pv[2] = (CM[4] * CM[9] * x - CM[5] * CM[8] * x - CM[0] * CM[9] * y + CM[1] * CM[8] * y + CM[0] * CM[5] * z - CM[1] * CM[4] * z) /
		(CM[0] * CM[5] * CM[10] - CM[0] * CM[6] * CM[9] - CM[1] * CM[4] * CM[10] + CM[1] * CM[6] * CM[8] + CM[2] * CM[4] * CM[9] - CM[2] * CM[5] * CM[8]);

	float norm_pv = norm_v3(pv);
	pv[0] = pv[0] / norm_pv;
	pv[1] = pv[1] / norm_pv;
	pv[2] = pv[2] / norm_pv;

	// attention principal_r should be normalized 
	float cos_v = (principal_r[0] * pv[0] + principal_r[1] * pv[1] + principal_r[2] * pv[2]);

	float dis1 = mindepth / cos_v;

	pt3d1[0] = C[0] + dis1 * pv[0];
	pt3d1[1] = C[1] + dis1 * pv[1];
	pt3d1[2] = C[2] + dis1 * pv[2];
	pt3d1[3] = 1;

	mult_3_4_4(CN, pt3d1, pt2d1);
	norm_by_v3(pt2d1);

	float dis2 = maxdepth / cos_v;

	pt3d2[0] = C[0] + dis2 * pv[0];
	pt3d2[1] = C[1] + dis2 * pv[1];
	pt3d2[2] = C[2] + dis2 * pv[2];
	pt3d2[3] = 1;

	mult_3_4_4(CN, pt3d2, pt2d2);
	norm_by_v3(pt2d2);

	// 对两个点向外扩展10个像元
	vec12[0] = pt2d2[0] - pt2d1[0];
	vec12[1] = pt2d2[1] - pt2d1[1];
	float norm_vec12 = norm_v2(vec12);
	vec12[0] = vec12[0] / norm_vec12;
	vec12[1] = vec12[1] / norm_vec12;

	// 赋值返回 此处检查出来有问题 pt2d1 应该是减去 10*vec
	x1 = pt2d1[0] - 10 * vec12[0];
	y1 = pt2d1[1] - 10 * vec12[1];

	x2 = pt2d2[0] + 10 * vec12[0];
	y2 = pt2d2[1] + 10 * vec12[1];
}

bool checkEpipolar(float* lpm, float* lpn, float* epl, float max_error)
{
	// epipolar line dist check

	float dis = point_2_line_dis(lpn[2], lpn[3], epl);
	return dis <= max_error;

}

void lineCross(float* line, float* linef)
{
	linef[0] = line[1] - line[3];
	linef[1] = line[2] - line[0];
	linef[2] = line[0] * line[3] -
		line[1] * line[2];
}

void makeCoffe(float* A, float* e, float* P, float* L, float* M, float* b, int cr)
{
	int ind = cr * 3;

	M[ind] = -P[0] * (e[0] * L[0] + e[1] * L[1] + e[2] * L[2]);
	M[ind + 1] = -P[1] * (e[0] * L[0] + e[1] * L[1] + e[2] * L[2]);
	M[ind + 2] = -e[0] * L[0] - e[1] * L[1] - e[2] * L[2];

	b[cr] = A[2] * L[0] + A[5] * L[1] + A[8] * L[2] +
		P[0] * (A[0] * L[0] + A[3] * L[1] + A[6] * L[2]) +
		P[1] * (A[1] * L[0] + A[4] * L[1] + A[7] * L[2]);

}

void obtainH(float* A, float* e,
	float* l11, float* l12,
	float* l21f, float* l22f,
	float* M, float* b, float* H)

{
	makeCoffe(A, e, l11, l21f, M, b, 0);
	makeCoffe(A, e, l11 + 2, l21f, M, b, 1);
	makeCoffe(A, e, l12, l22f, M, b, 2);

	// calculate v
	float v0, v1, v2;

	v0 = -(M[1] * M[5] * b[2] - M[2] * M[4] * b[2] -
		M[1] * M[8] * b[1] + M[2] * M[7] * b[1] +
		M[4] * M[8] * b[0] - M[5] * M[7] * b[0]) /
		(M[0] * M[4] * M[8] - M[0] * M[5] * M[7] -
			M[1] * M[3] * M[8] + M[1] * M[5] * M[6] +
			M[2] * M[3] * M[7] - M[2] * M[4] * M[6]);

	v1 = (M[0] * M[5] * b[2] - M[2] * M[3] * b[2] -
		M[0] * M[8] * b[1] + M[2] * M[6] * b[1] +
		M[3] * M[8] * b[0] - M[5] * M[6] * b[0]) /
		(M[0] * M[4] * M[8] - M[0] * M[5] * M[7] -
			M[1] * M[3] * M[8] + M[1] * M[5] * M[6] +
			M[2] * M[3] * M[7] - M[2] * M[4] * M[6]);

	v2 = -(M[0] * M[4] * b[2] - M[1] * M[3] * b[2] -
		M[0] * M[7] * b[1] + M[1] * M[6] * b[1] +
		M[3] * M[7] * b[0] - M[4] * M[6] * b[0]) /
		(M[0] * M[4] * M[8] - M[0] * M[5] * M[7] -
			M[1] * M[3] * M[8] + M[1] * M[5] * M[6] +
			M[2] * M[3] * M[7] - M[2] * M[4] * M[6]);

	// H=A-e*v';
	H[0] = A[0] - e[0] * v0;
	H[1] = A[1] - e[0] * v1;
	H[2] = A[2] - e[0] * v2;

	H[3] = A[3] - e[1] * v0;
	H[4] = A[4] - e[1] * v1;
	H[5] = A[5] - e[1] * v2;;

	H[6] = A[6] - e[2] * v0;
	H[7] = A[7] - e[2] * v1;
	H[8] = A[8] - e[2] * v2;
}

void findHomographyGo(float* inter_range,
	float* lpsm, int lpsm_size, float* lpsn, int lpsn_size, short* plMap, float* lm, float* ln,
	float* CM, float* CN, float* F, float* A, float* e, float* C,
	int imr, int imc, cv::Mat& match_score, float error_max, float max_ang)
{
	float buffer = round(error_max + 0.5);
	float best_ang = 999;

	int idl11 = 0, idl12 = 0, idl21 = 0, idl22 = 0, ind_matches = 0, bestID;

	float* l11, * l12, * l21, * l22, * lpm, * lpn, * angCellM, * angCellIND;

	float itpm[3];
	float el[3];
	float principal_r[3];
	float addvec[2];

	float M_3_3[9];
	float b_3[3];
	float H_3_3[9];
	float l21f[3], l22f[3];
	float p1[3], p2[3], p1_[3], p2_[3];
	float vec1[2], vec2[2];
	float best_H[9];

	float l_e1[3], l_e2[3], l_int_e1[3], l_int_e2[3];

	float* matchesL = new float[lpsn_size]();

	for (int i = 0; i < lpsn_size; i++)
		matchesL[i] = 0;

	cv::Mat match_score_;
	cv::Mat match_score_row(1, 11, CV_32FC1);

	bool bvalid, bchange;

	std::vector<int>xx, yy;

	float x1, y1, x2, y2, ang1;
	float mindepth, maxdepth, pr_normal;

	principal_r[0] = CM[8];
	principal_r[1] = CM[9];
	principal_r[2] = CM[10];
	pr_normal = norm_v3(principal_r);

	principal_r[0] = principal_r[0] / pr_normal;
	principal_r[1] = principal_r[1] / pr_normal;
	principal_r[2] = principal_r[2] / pr_normal;

	short lpsn_id;

	for (int i = 0; i < lpsm_size; i++)
	{
		//if (i == 82)
		//std::getchar();
		lpm = lpsm + i * 8;
		angCellM = lpm + 4;

		idl11 = lpm[0];
		idl12 = lpm[1];

		l11 = lm + idl11 * 7;
		l12 = lm + idl12 * 7;

		itpm[0] = lpm[2];
		itpm[1] = lpm[3];
		itpm[2] = 1;

		// compute the epipolar line
		mult_3_3_3(F, itpm, el);

		// compute the depth
		mindepth = inter_range[i * 2];
		maxdepth = inter_range[i * 2 + 1];

		//std::cout << " " << mindepth << " " << maxdepth;
		//std::cout << std::endl;

		// obtain search line ends

		depthRange(C, CM, CN,
			itpm, principal_r,
			mindepth, maxdepth,
			x1, y1, x2, y2);

		///////breshman all pts//////////////////
		if (abs(x1 - x2) > abs(y1 - y2))
		{
			addvec[0] = 0;
			addvec[1] = 1;
		}
		else
		{
			addvec[0] = 1;
			addvec[1] = 0;
		}

		//obtain breshman points
		Bresenham(x1, y1, x2, y2, xx, yy);

		bestID = -1;
		best_ang = 999;

		for (int mm = 0; mm < xx.size(); mm++)
			for (int k = -buffer; k <= buffer; k++)
			{

				int sx = xx.at(mm) + k * addvec[0];
				int sy = yy.at(mm) + k * addvec[1];

				if (sx <= 0 || sx >= imc || sy <= 0 || sy >= imr)
					continue;

				lpsn_id = plMap[sy * imc + sx];

				if (lpsn_id >= 0)
					continue;

				lpsn_id = -lpsn_id - 1;

				if (lpsn_id >= lpsn_size)
					continue;

				lpn = lpsn + lpsn_id * 8;
				angCellIND = lpn + 4;

				//check epipolar for junction
				if (!checkEpipolar(lpm, lpn, el, error_max))
					continue;

				// angle alignment check
				checkAng(angCellM, angCellIND, bvalid, bchange);

				if (!bvalid)
					continue;

				if (bchange)
				{
					idl21 = lpn[1];
					idl22 = lpn[0];
				}

				else
				{
					idl21 = lpn[0];
					idl22 = lpn[1];
				}

				l21 = ln + idl21 * 7;
				l22 = ln + idl22 * 7;

				// cross line
				lineCross(l21, l21f);
				lineCross(l22, l22f);

				// further check the epipolar geometry

				// first line
				p1[0] = l11[0];
				p1[1] = l11[1];
				p1[2] = 1;

				p2[0] = l11[2];
				p2[1] = l11[3];
				p2[2] = 1;

				mult_3_3_3(F, p1, l_e1);
				mult_3_3_3(F, p2, l_e2);

				cross_v3(l_e1, l21f, l_int_e1);
				cross_v3(l_e2, l21f, l_int_e2);

				norm_by_v3(l_int_e1);
				norm_by_v3(l_int_e2);

				if (!twoLines_intersec(l_int_e1, l_int_e2, l21, l21 + 2, line_2_line_intersec))
					continue;

				// second line
				p1[0] = l12[0];
				p1[1] = l12[1];
				p1[2] = 1;

				p2[0] = l12[2];
				p2[1] = l12[3];
				p2[2] = 1;

				mult_3_3_3(F, p1, l_e1);
				mult_3_3_3(F, p2, l_e2);

				cross_v3(l_e1, l22f, l_int_e1);
				cross_v3(l_e2, l22f, l_int_e2);

				norm_by_v3(l_int_e1);
				norm_by_v3(l_int_e2);

				if (!twoLines_intersec(l22, l22 + 2, l_int_e1, l_int_e2, line_2_line_intersec))
					continue;

				// obtain homography 
				obtainH(A, e,
					l11, l12,
					l21f, l22f,
					M_3_3, b_3, H_3_3);

				// check homography 
				// first line
				p1[0] = l11[0];
				p1[1] = l11[1];
				p1[2] = 1;

				p2[0] = l11[2];
				p2[1] = l11[3];
				p2[2] = 1;

				mult_3_3_3(H_3_3, p1, p1_);
				norm_by_v3(p1_);
				mult_3_3_3(H_3_3, p2, p2_);
				norm_by_v3(p2_);

				vec1[0] = l21[2] - l21[0];
				vec1[1] = l21[3] - l21[1];

				vec2[0] = p2_[0] - p1_[0];
				vec2[1] = p2_[1] - p1_[1];
				//ang_vec = std::abs(angofvec(vec1, vec2));
				//std::cout <<"ang_vec"<<" "<< ang_vec << std::endl;

				if (vec1[0] * vec2[0] + vec1[1] * vec2[1] < 0)
					continue;

				//second line
				p1[0] = l12[0];
				p1[1] = l12[1];
				p1[2] = 1;

				p2[0] = l12[2];
				p2[1] = l12[3];
				p2[2] = 1;

				mult_3_3_3(H_3_3, p1, p1_);
				norm_by_v3(p1_);
				mult_3_3_3(H_3_3, p2, p2_);
				norm_by_v3(p2_);

				vec1[0] = l22[2] - l22[0];
				vec1[1] = l22[3] - l22[1];

				vec2[0] = p2_[0] - p1_[0];
				vec2[1] = p2_[1] - p1_[1];

				ang1 = abs(ang_of_vec(vec1[0], vec1[1], vec2[0], vec2[1]));

				if (ang1 > max_ang || best_ang < ang1)
					continue;

				best_ang = ang1;

				for (int kk = 0; kk < 9; kk++)
					best_H[kk] = H_3_3[kk];

				bestID = lpsn_id;

			}

		if (bestID == -1)
			continue;

		if (matchesL[bestID] != 0 && matchesL[bestID] < best_ang)
			continue;

		// store the best
		matchesL[bestID] = best_ang;
		match_score_row.at<float>(0, 0) = itpm[0];
		match_score_row.at<float>(0, 1) = itpm[1];

		for (int kk = 0; kk < 9; kk++)
			match_score_row.at<float>(0, 2 + kk) = best_H[kk];

		match_score_.push_back(match_score_row);

	}

	match_score_.copyTo(match_score);
	delete[] matchesL;
}

void findHomography(cv::Mat inter_range_i,
	cv::Mat lpsm_Mf, cv::Mat lpsn_Mf, cv::Mat plMap_Mf, cv::Mat lm_Mf, cv::Mat ln_Mf,
	cv::Mat CM_Mf, cv::Mat CN_Mf, cv::Mat F_Mf, cv::Mat Ae_Mf, cv::Mat& matches_M_float,
	int imr, int imc,
	float error_max,
	float ang_max)
{
	cv::Mat M_ = CM_Mf.reshape(0, 3);
	cv::Mat M = M_.rowRange(0, 3).clone();
	M = M.colRange(0, 3).clone();
	cv::Mat p4 = M_.colRange(3, 4).clone();

	cv::Mat inv_M;
	cv::invert(M, inv_M);

	cv::Mat C = -inv_M * p4;

	cv::Mat A = Ae_Mf.colRange(0, 9).clone();
	cv::Mat e = Ae_Mf.colRange(9, 12).clone();

	findHomographyGo((float*)inter_range_i.data,
		(float*)lpsm_Mf.data, lpsm_Mf.rows, (float*)lpsn_Mf.data, lpsn_Mf.rows,
		(short*)plMap_Mf.data, (float*)lm_Mf.data, (float*)ln_Mf.data,
		(float*)CM_Mf.data, (float*)CN_Mf.data, (float*)F_Mf.data,
		(float*)A.data, (float*)e.data, (float*)C.data,
		imr, imc, matches_M_float, error_max, ang_max);

}


