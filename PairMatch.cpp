#include "PairMatch.h"
#include "BasicMath.h"
#include "Parameters.h"
#include "IO.h"
#include"Epipolar.h"

void get3DLine(cv::Mat CM, cv::Mat CN, cv::Mat F, float* l1, float* l2, float* p_3d_1, float* p_3d_2)
{
	cv::Mat pt1 = cv::Mat(3, 1, CV_32FC1);
	cv::Mat pt2 = cv::Mat(3, 1, CV_32FC1);
	cv::Mat pt3 = cv::Mat(3, 1, CV_32FC1);
	cv::Mat pt4 = cv::Mat(3, 1, CV_32FC1);

	pt1.at<float>(0, 0) = l1[0];
	pt1.at<float>(1, 0) = l1[1];
	pt1.at<float>(2, 0) = 1;

	cv::Mat ep1 = F * pt1;

	pt2.at<float>(0, 0) = l1[2];
	pt2.at<float>(1, 0) = l1[3];
	pt2.at<float>(2, 0) = 1;

	cv::Mat ep2 = F * pt2;

	pt3.at<float>(0, 0) = l2[0];
	pt3.at<float>(1, 0) = l2[1];
	pt3.at<float>(2, 0) = 1;

	pt4.at<float>(0, 0) = l2[2];
	pt4.at<float>(1, 0) = l2[3];
	pt4.at<float>(2, 0) = 1;

	cv::Mat l2f = pt3.cross(pt4);

	cv::Mat pt5 = l2f.cross(ep1);
	cv::Mat pt6 = l2f.cross(ep2);

	pt5 = pt5 / pt5.at<float>(2, 0);
	pt6 = pt6 / pt6.at<float>(2, 0);

	std::vector<cv::Point2f> pt1Vec;
	std::vector<cv::Point2f> pt2Vec;

	pt1Vec.push_back(cv::Point2f(pt1.at<float>(0, 0), pt1.at<float>(1, 0)));
	pt1Vec.push_back(cv::Point2f(pt2.at<float>(0, 0), pt2.at<float>(1, 0)));

	pt2Vec.push_back(cv::Point2f(pt5.at<float>(0, 0), pt5.at<float>(1, 0)));
	pt2Vec.push_back(cv::Point2f(pt6.at<float>(0, 0), pt6.at<float>(1, 0)));

	cv::Mat pnt3D;

	cv::triangulatePoints(CM, CN, pt1Vec, pt2Vec, pnt3D);

	p_3d_1[0] = pnt3D.at<float>(0, 0);
	p_3d_1[1] = pnt3D.at<float>(1, 0);
	p_3d_1[2] = pnt3D.at<float>(2, 0);
	p_3d_1[3] = pnt3D.at<float>(3, 0);

	p_3d_2[0] = pnt3D.at<float>(0, 1);
	p_3d_2[1] = pnt3D.at<float>(1, 1);
	p_3d_2[2] = pnt3D.at<float>(2, 1);
	p_3d_2[3] = pnt3D.at<float>(3, 1);

}

cv::Mat get3DLine(cv::Mat CM, cv::Mat CN, cv::Mat F, cv::Mat l1, cv::Mat l2, float& ep_ang)
{
	cv::Mat pt1 = cv::Mat(3, 1, CV_32FC1);
	cv::Mat pt2 = cv::Mat(3, 1, CV_32FC1);
	cv::Mat pt3 = cv::Mat(3, 1, CV_32FC1);
	cv::Mat pt4 = cv::Mat(3, 1, CV_32FC1);

	pt1.at<float>(0, 0) = l1.at<float>(0, 0);
	pt1.at<float>(1, 0) = l1.at<float>(0, 1);
	pt1.at<float>(2, 0) = 1;

	cv::Mat ep1 = F * pt1;

	pt2.at<float>(0, 0) = l1.at<float>(0, 2);
	pt2.at<float>(1, 0) = l1.at<float>(0, 3);
	pt2.at<float>(2, 0) = 1;

	cv::Mat ep2 = F * pt2;

	pt3.at<float>(0, 0) = l2.at<float>(0, 0);
	pt3.at<float>(1, 0) = l2.at<float>(0, 1);
	pt3.at<float>(2, 0) = 1;

	pt4.at<float>(0, 0) = l2.at<float>(0, 2);
	pt4.at<float>(1, 0) = l2.at<float>(0, 3);
	pt4.at<float>(2, 0) = 1;

	cv::Mat l2f = pt3.cross(pt4);

	cv::Mat pt5 = l2f.cross(ep1);
	cv::Mat pt6 = l2f.cross(ep2);

	pt5 = pt5 / pt5.at<float>(2, 0);
	pt6 = pt6 / pt6.at<float>(2, 0);

	std::vector<cv::Point2f> pt1Vec;
	std::vector<cv::Point2f> pt2Vec;

	pt1Vec.push_back(cv::Point2f(pt1.at<float>(0, 0), pt1.at<float>(1, 0)));
	pt1Vec.push_back(cv::Point2f(pt2.at<float>(0, 0), pt2.at<float>(1, 0)));

	pt2Vec.push_back(cv::Point2f(pt5.at<float>(0, 0), pt5.at<float>(1, 0)));
	pt2Vec.push_back(cv::Point2f(pt6.at<float>(0, 0), pt6.at<float>(1, 0)));

	cv::Mat pnt3D;

	cv::triangulatePoints(CM, CN, pt1Vec, pt2Vec, pnt3D);

	// calculate epipolar angle
	float k1 = -l2f.at<float>(0, 0) / l2f.at<float>(1, 0);
	float k2 = -ep1.at<float>(0, 0) / ep1.at<float>(1, 0);

	float ang1 = atan(abs((k2 - k1) / (1 + k1 * k2)));
	k2 = -ep2.at<float>(0, 0) / ep2.at<float>(1, 0);
	float ang2 = atan(abs((k2 - k1) / (1 + k1 * k2)));

	ep_ang = min_2(ang1, ang2);

	return  pnt3D;

}

void guidedMatching(float* lines1, float* lines_range, int* lines1_knn, float* lines2, cv::Mat line2_map,
	int lsize1, int lsize2, float* homos, cv::Mat  CM, cv::Mat  CN, cv::Mat  F, int knn_num, float dist,
	int imr, int imc, short* match4Second)
{

	float* l1, * l2, * H, * F_ptr, * CM_ptr;
	float p1[3], p2[3], p1_[3], p2_[3], vec1[2], vec2[2], p3[3], p4[3];
	float epl1[3], epl2[3], l_f1[3], l_f2[3];
	float epp1[3], epp2[3];
	int addvec[2];
	float p_3d_1[4], p_3d_2[4];

	std::vector<int>xx, yy;

	F_ptr = (float*)F.data;
	CM_ptr = (float*)CM.data;

	int* id_candidate = new int[knn_num];
	float* score_candidate = new float[knn_num];

	int* searched_ID = new int[1000];
	int searched_size;
	bool buffer_control = 0;

	int sx, sy, lid, buffer, l1_id, homos_id, lastID, bestID, candidate_size = 0;
	float bestDis, maxdis, dis1, dis2, max_depth1, min_depth1, max_depth2, min_depth2, pjl_length, vec_cos;

	buffer = round(dist + 0.5);

	for (int i = 0; i < lsize1; i++)
	{
		l1_id = i;
		l1 = lines1 + l1_id * 7;

		for (int mm = 0; mm < knn_num; mm++)
		{
			id_candidate[mm] = 0;
			score_candidate[mm] = 0;
		}

		candidate_size = 0;

		min_depth1 = lines_range[i * 4];
		max_depth1 = lines_range[i * 4 + 1];

		min_depth2 = lines_range[i * 4 + 2];
		max_depth2 = lines_range[i * 4 + 3];

		for (int j = 0; j < knn_num; j++)
		{
			lastID = 0;
			bestID = 0;
			bestDis = 0;

			homos_id = lines1_knn[i * knn_num + j];
			H = homos + homos_id * 11 + 2;

			p1[0] = l1[0];
			p1[1] = l1[1];
			p1[2] = 1;

			p2[0] = l1[2];
			p2[1] = l1[3];
			p2[2] = 1;

			// mapping with homography
			mult_3_3_3(H, p1, p1_);
			mult_3_3_3(H, p2, p2_);
			norm_by_v3(p1_);
			norm_by_v3(p2_);

			cross_v3(p1_, p2_, l_f1);
			//
			vec1[0] = p2_[0] - p1_[0];
			vec1[1] = p2_[1] - p1_[1];
			pjl_length = norm_v2(vec1);

			if (pjl_length > 3 * l1[6] || pjl_length < l1[6] / 3)
				continue;

			if (vec1[0] > vec1[1])
			{
				addvec[0] = 0;
				addvec[1] = 1;
			}
			else
			{
				addvec[0] = 1;
				addvec[1] = 0;
			}

			Bresenham(round(p1_[0]), round(p1_[1]), round(p2_[0]), round(p2_[1]), xx, yy);

			searched_size = 0;
			for (int mm = 0; mm < xx.size(); mm += 6)
				for (int k = -buffer; k <= buffer; k++)
				{
					buffer_control = !buffer_control;
					if (buffer_control)
						continue;

					sx = xx[mm] + k * addvec[0];
					sy = yy[mm] + k * addvec[1];

					if (sx <= 0 || sx >= imc || sy <= 0 || sy >= imr)
						continue;

					lid = line2_map.at<short>(sy, sx) - 1;

					if (lid < 0 || lid >= lsize2 || lastID == lid || lid == bestID)
						continue;
					lastID = lid;

					if (ID_in_array(searched_ID, searched_size, lid))
						continue;

					searched_ID[searched_size] = lid;
					searched_size++;

					l2 = lines2 + lid * 7;

					p3[0] = l2[0];
					p3[1] = l2[1];
					p3[2] = 1;

					p4[0] = l2[2];
					p4[1] = l2[3];
					p4[2] = 1;

					vec2[0] = p4[0] - p3[0];
					vec2[1] = p4[1] - p3[1];

					cross_v3(p3, p4, l_f2);

					// check direction
					// make sure they are in the same directiom
					if (vec1[0] * vec2[0] + vec1[1] * vec2[1] < 0) 
						continue;

					
					// check epipolar line 
					mult_3_3_3(F_ptr, p1, epl1);
					mult_3_3_3(F_ptr, p2, epl2);

					cross_v3(l_f2, epl1, epp1);
					cross_v3(l_f2, epl2, epp2);
					norm_by_v3(epp1);
					norm_by_v3(epp2);

					// two lines intersection;
					if (!twoLines_intersec(epp1, epp2, p3, p4, line_2_line_intersec))
						continue;
					// distance check 
					dis1 = point_2_line_dis(p3, l_f1);
					if (dis1 > dist)
						continue;

					dis2 = point_2_line_dis(p4, l_f1);
					if (dis2 > dist)
						continue;

					// intersection check 
					if (!twoLines_intersec(p3, p4, p1_, p2_, line_2_line_intersec))
						continue;

					//get the 3D line and depth check 
					get3DLine(CM, CN, F, l1, l2, p_3d_1, p_3d_2);

					norm_by_v4(p_3d_1);
					float w1 = CM_ptr[11]
						+ CM_ptr[8] * p_3d_1[0]
						+ CM_ptr[9] * p_3d_1[1]
						+ CM_ptr[10] * p_3d_1[2];

					if (w1 < min_depth1 || w1 > max_depth1)
						continue;

					norm_by_v4(p_3d_2);
					float w2 = CM_ptr[11]
						+ CM_ptr[8] * p_3d_2[0]
						+ CM_ptr[9] * p_3d_2[1]
						+ CM_ptr[10] * p_3d_2[2];

					if (w2 < min_depth1 || w2 > max_depth1)
						continue;

					maxdis = exp(-max_2(dis2, dis1) / (2 * dist));

					if (bestDis >= maxdis)
						continue;

					bestID = lid;
					bestDis = maxdis;

				}

			if (bestID == 0)
				continue;

			bool founded = 0;
			for (int mm = 0; mm < candidate_size; mm++)
				if (id_candidate[mm] == bestID)
				{
					score_candidate[mm] += bestDis;
					founded = 1;
					break;
				}

			if (!founded)
			{
				score_candidate[candidate_size] = bestDis;
				id_candidate[candidate_size] = bestID;
				candidate_size++;
			}

		}

		if (candidate_size == 0)
			continue;

		//find max score
		maxdis = 0;
		int min_id = -1;

		for (int mm = 0; mm < candidate_size; mm++)
		{
			if (score_candidate[mm] > maxdis)
			{
				maxdis = score_candidate[mm];
				min_id = id_candidate[mm];
			}
		}

		// !!!!!!!!!!!!!!!!!!!attention this plus 1 !!!!!!!!!!!!!!!!!!!!!!!!
		match4Second[i] = min_id + 1;
	}

	delete[] id_candidate;
	delete[] score_candidate;
	delete[] searched_ID;
}

void line2KDtree(cv::Mat lines_Mf, cv::Mat homo_Mf, cv::Mat& inter_knn_Mi, int support_H_num)
{
	// construct kdtree
	cv::flann::Index flannIndex(homo_Mf.colRange(0, 2), cv::flann::KDTreeIndexParams());

	// store amd query
	cv::Mat inter_2_pt3index_dist;

	cv::Mat inter_knn_Mi_;

	flannIndex.knnSearch(lines_Mf.colRange(4, 6).clone(), inter_knn_Mi_,
		inter_2_pt3index_dist, support_H_num, cv::flann::SearchParams());

	inter_knn_Mi_.copyTo(inter_knn_Mi);
}

void matSkew(cv::Mat M1f, cv::Mat& skew_Mf)
{
	cv::Mat skew_Mf_ = cv::Mat(3, 3, CV_32FC1);

	skew_Mf_.at<float>(0, 0) = 0;
	skew_Mf_.at<float>(0, 1) = -M1f.at<float>(2, 0);
	skew_Mf_.at<float>(0, 2) = M1f.at<float>(1, 0);

	skew_Mf_.at<float>(1, 0) = M1f.at<float>(2, 0);
	skew_Mf_.at<float>(1, 1) = 0;
	skew_Mf_.at<float>(1, 2) = -M1f.at<float>(0, 0);

	skew_Mf_.at<float>(2, 0) = -M1f.at<float>(1, 0);
	skew_Mf_.at<float>(2, 1) = M1f.at<float>(0, 0);
	skew_Mf_.at<float>(2, 2) = 0;

	skew_Mf_.copyTo(skew_Mf);
}

void resconstruct3DL(cv::Mat CM, cv::Mat CN, cv::Mat F,
	cv::Mat match4Second, cv::Mat lines1_Mf, cv::Mat lines2_Mf,
	cv::Mat& pt3Ds_,
	ScoreRecorder* merge_tool, int task_id,float epipolar_ang)
{
	float ep_angle;
	cv::Mat CM_ = CM.reshape(0, 3);
	cv::Mat CN_ = CN.reshape(0, 3);
	cv::Mat F_ = F.reshape(0, 3);

	cv::Mat pt3Ds = cv::Mat::zeros(lines1_Mf.rows, 10, CV_32FC1);

	cv::Mat pt3D;

	int l3d_ind = 0;

	for (int i = 0; i < match4Second.rows; i++)
	{
		int ind = match4Second.at<short>(i, 0);

		if (ind == 0)
			continue;

		l3d_ind++;
	}

	merge_tool->initial_tasks(task_id, l3d_ind);

	l3d_ind = 0;
	for (int i = 0; i < match4Second.rows; i++)
	{
		int ind = match4Second.at<short>(i, 0);

		if (ind == 0)
			continue;

		ind = ind - 1;// remember ind must minors 1

		if (lines1_Mf.at<float>(i, 6) < lines2_Mf.at<float>(ind, 6))
		{
			cv::Mat F_trans;
			cv::transpose(F_, F_trans);
			pt3D = get3DLine(CN_, CM_, F_trans,
				lines2_Mf.rowRange(ind, ind + 1).clone(),
				lines1_Mf.rowRange(i, i + 1).clone(),
				ep_angle);
		}
		else
			pt3D = get3DLine(CM_, CN_, F_,
				lines1_Mf.rowRange(i, i + 1),
				lines2_Mf.rowRange(ind, ind + 1),
				ep_angle);

		pt3Ds.at<float>(i, 0) = pt3D.at<float>(0, 0) / pt3D.at<float>(3, 0);
		pt3Ds.at<float>(i, 1) = pt3D.at<float>(1, 0) / pt3D.at<float>(3, 0);
		pt3Ds.at<float>(i, 2) = pt3D.at<float>(2, 0) / pt3D.at<float>(3, 0);
		pt3Ds.at<float>(i, 3) = 1;

		pt3Ds.at<float>(i, 4) = pt3D.at<float>(0, 1) / pt3D.at<float>(3, 1);
		pt3Ds.at<float>(i, 5) = pt3D.at<float>(1, 1) / pt3D.at<float>(3, 1);
		pt3Ds.at<float>(i, 6) = pt3D.at<float>(2, 1) / pt3D.at<float>(3, 1);
		pt3Ds.at<float>(i, 7) = 1;

		pt3Ds.at<float>(i, 8) = ep_angle;
		pt3Ds.at<float>(i, 9) = l3d_ind;

		merge_tool->add_ori_id(task_id, l3d_ind, i,
			log(ep_angle / (2 * epipolar_ang)));

		l3d_ind++;
	}

	pt3Ds.copyTo(pt3Ds_);
}

void createMap(int imr, int imc, cv::Mat inter_lines_Mf, cv::Mat lines_Mf,
	cv::Mat& inter_map_)
{
	cv::Mat inter_Mat_ushort = cv::Mat::zeros(imr, imc, CV_16SC1);
	std::vector<int> xx, yy;
	int x1, y1, x2, y2;

	float* lines_p = (float*)lines_Mf.data;
	short* inter_Mat_p = (short*)inter_Mat_ushort.data;
	float* inter_lines_p = (float*)inter_lines_Mf.data;

	int ind = 0;

	for (int i = 0; i < lines_Mf.rows; i++)
	{
		ind = i * lines_Mf.cols;
		x1 = round(lines_p[ind]);
		y1 = round(lines_p[ind + 1]);
		x2 = round(lines_p[ind + 2]);
		y2 = round(lines_p[ind + 3]);

		Bresenham(x1, y1, x2, y2, xx, yy);

		for (int j = 0; j < xx.size(); j++)
		{
			if (xx.at(j) < 0 || yy.at(j) < 0 || xx.at(j) >= imc || yy.at(j) >= imr)
				continue;

			inter_Mat_p[yy.at(j) * inter_Mat_ushort.cols + xx.at(j)] = i + 1;

			//std::cout << sparse_inter_map.ref<short>(yy.at(j), xx.at(j)) << " " << inter_Mat_p[yy.at(j) * inter_Mat_ushort.cols + xx.at(j)] << std::endl;
			//std::getchar();
		}
	}

	int px, py;
	for (int i = 0; i < inter_lines_Mf.rows; i++)
	{
		ind = i * inter_lines_Mf.cols;

		px = inter_lines_p[ind + 2];
		py = inter_lines_p[ind + 3];

		if (px < 0 || py < 0 || px >= imc || py >= imr)
			continue;

		inter_Mat_p[py * inter_Mat_ushort.cols + px] = -i - 1;
	}

	inter_Mat_ushort.copyTo(inter_map_);
}

void matchPair(std::vector<std::string>image_names, std::string input_folder,
	cv::Mat imidx_Mf, cv::Mat imsizes_Mf,
	cv::Mat cameras_Mf,
	int match_ind,
	ScoreRecorder* merge_tool,
	float epipolar_ang,
	float error_max,
	float ang_max)
{
	int mid1 = imidx_Mf.at<float>(match_ind, 0);
	int mid2 = imidx_Mf.at<float>(match_ind, 1);

	cv::Mat CM = cameras_Mf.rowRange(mid1, mid1 + 1).clone().reshape(0, 3);
	cv::Mat CN = cameras_Mf.rowRange(mid2, mid2 + 1).clone().reshape(0, 3);

	cv::Mat F, Ae;
	computeFAE(CM, CN, F, Ae);

	Ae = Ae.mul(-1);
	
	int rm2 = imsizes_Mf.at<float>(mid2, 0);
	int cm2 = imsizes_Mf.at<float>(mid2, 1);

	// load lines and intersections
	cv::Mat lines_range1, lines_range2, l2l_range1, l2l_range2;

	cv::Mat lines1_Mf, lines2_Mf, l2l_1_Mf, l2l_2_Mf;
	cv::Mat inter_range_i, line_range_i;

	std::string name1 = input_folder + "lines\\" + image_names.at(mid1) + ".m";
	std::string name2 = input_folder + "lines\\" + image_names.at(mid2) + ".m";

	readMat(name1, lines_range1);
	readMat(name2, lines_range2);

	lines1_Mf = lines_range1.colRange(0, 7).clone();
	lines2_Mf = lines_range2.colRange(0, 7).clone();
	line_range_i = lines_range1.colRange(7, 11).clone();

	name1 = input_folder + "l2l\\" + image_names.at(mid1) + ".m";
	name2 = input_folder + "l2l\\" + image_names.at(mid2) + ".m";

	readMat(name1, l2l_range1);
	readMat(name2, l2l_range2);

	l2l_1_Mf = l2l_range1.colRange(0, 8).clone();
	l2l_2_Mf = l2l_range2.colRange(0, 8).clone();
	inter_range_i = l2l_range1.colRange(8, 10).clone();

	cv::Mat l2l_Mf;

	createMap(rm2, cm2, l2l_2_Mf, lines2_Mf,
		l2l_Mf);
	
	cv::Mat homos;
	findHomography(inter_range_i,
		l2l_1_Mf, l2l_2_Mf, l2l_Mf, lines1_Mf, lines2_Mf,
		CM, CN, F, Ae, homos, rm2, cm2,
		error_max, ang_max);

	std::cout << imidx_Mf.rows - match_ind << std::endl;
	cv::Mat match4Second = cv::Mat::zeros(lines1_Mf.rows, 1, CV_16SC1);
	cv::Mat pt3Ds, matches_r2l;

	// save null mat
	if (homos.rows <= support_homo_num * 3)
	{
		saveMat(input_folder + "matches\\" + std::to_string(match_ind) + ".m",
			match4Second);

		saveMat(input_folder + "l3ds\\" + std::to_string(match_ind) + ".m",
			pt3Ds);
		return;
	}

	cv::Mat line_knn_Mi;
	line2KDtree(lines1_Mf, homos, line_knn_Mi, support_homo_num);

	// match for single lines
	float dist = error_max;

	guidedMatching((float*)lines1_Mf.data, (float*)line_range_i.data, (int*)line_knn_Mi.data, (float*)lines2_Mf.data,
		l2l_Mf, lines1_Mf.rows, lines2_Mf.rows, (float*)homos.data, CM, CN,
		F, line_knn_Mi.cols, dist,
		rm2, cm2, (short*)match4Second.data);

	resconstruct3DL(CM, CN, F,
		match4Second, lines1_Mf, lines2_Mf, pt3Ds, merge_tool, match_ind, epipolar_ang);

	//save mat
	saveMat(input_folder + "matches\\" + std::to_string(match_ind) + ".m",
		match4Second);

	saveMat(input_folder + "l3ds\\" + std::to_string(match_ind) + ".m",
		pt3Ds);

}
