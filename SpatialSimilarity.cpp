# include "SpatialSimilarity.h"
#include "BasicMath.h"
#include <iostream>
#include <fstream>
#include <mutex>
#include "ThreadPool.h"
#include "Parameters.h"
#include "IO.h"



void validPair(std::string input_folder, cv::Mat imidx_Mf, std::vector<std::vector<bool>>& valid_pair)
{

	std::string valid_file_adr = input_folder + "valid_pair.m";
	cv::Mat pair_M;
	readMat(valid_file_adr, pair_M);

	for (int i = 0; i < imidx_Mf.rows; i++)
		for (int j = 0; j < imidx_Mf.rows; j++)
			valid_pair[i][j] = pair_M.at<ushort>(i, j);

}

void lineMap(int imr, int imc, cv::Mat lines_Mf,
	cv::Mat* inter_map_)
{
	cv::Mat inter_Mat_ushort = cv::Mat::zeros(imr, imc, CV_16SC1);
	std::vector<int> xx, yy;
	int x1, y1, x2, y2;

	float* lines_p = (float*)lines_Mf.data;
	short* inter_Mat_p = (short*)inter_Mat_ushort.data;

	int ind = 0;

	int ind_l = 0;

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

	inter_Mat_ushort.copyTo(*inter_map_);

}

void matchesRight2Left(cv::Mat matches, cv::Mat lines2, cv::Mat* matches_r2l_)
{

	int p2_size = lines2.rows;

	int lid_2;

	cv::Mat matches_r2l = cv::Mat::zeros(p2_size, 1, CV_16SC1);

	for (int i = 0; i < matches.rows; i++)
	{
		lid_2 = matches.at<short>(i, 0) - 1;

		if (lid_2 < 0)
			continue;

		matches_r2l.at<short>(lid_2, 0) = i + 1;

	}

	matches_r2l.copyTo(*matches_r2l_);
}

void addAndRelease(std::vector< pair_m>* l3d_pair_vec,
	std::vector<cv::Mat>* matches_vec,
	std::vector<cv::Mat>* pt3d_vec,
	std::vector<cv::Mat>* matches_r2l_vec,
	int im11, int im12, int im21, int im22,
	int stereo_i,
	int stereo_j,
	int lid11,
	int lid12
)
{
	// add pure relations -> im21
	if ((im11 == im21)
		&& (*matches_vec)[stereo_j].at<short>(lid11, 0) > 0)
	{
		(*l3d_pair_vec).push_back({
			(int)(*pt3d_vec)[stereo_i].at<float>(lid11, 9),
			(int)(*pt3d_vec)[stereo_j].at<float>(lid11, 9),
				0 });
	}

	else if (im12 == im21 && (*matches_vec)[stereo_j].at<short>(lid12, 0) > 0)
	{
		(*l3d_pair_vec).push_back({ (int)(*pt3d_vec)[stereo_i].at<float>(lid11, 9),
			(int)(*pt3d_vec)[stereo_j].at<float>(lid12, 9),
			0 });
	}

	else if (im11 == im22)
	{
		int lid_21_ = (*matches_r2l_vec)[stereo_j].at<short>(lid11, 0) - 1;

		if (lid_21_ >= 0 && (*matches_vec)[stereo_j].at<short>(lid_21_, 0) > 0)
		{
			(*l3d_pair_vec).push_back({
			(int)(*pt3d_vec)[stereo_i].at<float>(lid11, 9),
			(int)(*pt3d_vec)[stereo_j].at<float>(lid_21_, 9),
				0 });
		}
	}

	if (im12 == im22)
	{
		int lid_21_ = (*matches_r2l_vec)[stereo_j].at<short>(lid12, 0) - 1;
		if (lid_21_ >= 0 && (*matches_vec)[stereo_j].at<short>(lid_21_, 0) > 0)
		{
			(*l3d_pair_vec).push_back({
			(int)(*pt3d_vec)[stereo_i].at<float>(lid11, 9),
			(int)(*pt3d_vec)[stereo_j].at<float>(lid_21_, 9),
				0 });
		}
	}
}

std::vector< pair_m> checkStereos1(int stereo_i, int stereo_j,
	std::vector<cv::Mat>* matches_vec,
	std::vector<cv::Mat>* matches_r2l_vec,
	std::vector<cv::Mat>* pt3d_vec,
	std::vector<cv::Mat>* map_vec,
	std::vector<cv::Mat>* lines_vec,
	float* imidx_Mf_ptr,
	cv::Mat cameras_Mf,
	cv::Mat imsizes_Mf,
	float dist,
	int buffer
)
{
	cv::Mat x3i, y3i, x3j, y3j, v3i, v3j, x1, y1, lf, x2, y2, lf_2;
	int lid11, lid12, lid21, lid22, lid_21_;

	float min_line3D_ang = maximum_ang;
	float min_line3D_cos = cos(min_line3D_ang);
	float cos_v3,
		dis1, dis2,
		max_dis, score;

	std::vector< pair_m> l3d_pair_vec;

	if ((*matches_vec)[stereo_j].rows == 0)
		return l3d_pair_vec;

	int im11 = imidx_Mf_ptr[stereo_i * 2];
	int im12 = imidx_Mf_ptr[stereo_i * 2 + 1];
	int im21 = imidx_Mf_ptr[stereo_j * 2];
	int im22 = imidx_Mf_ptr[stereo_j * 2 + 1];

	int rm1 = imsizes_Mf.at<float>(im21, 0);
	int cm1 = imsizes_Mf.at<float>(im21, 1);

	cv::Mat Pm = cameras_Mf.rowRange(im21, im21 + 1).reshape(0, 3);
	cv::Mat Pn = cameras_Mf.rowRange(im22, im22 + 1).reshape(0, 3);

	float* x1_p, * y1_p, * lf_p, * lineb1;
	float angle_3dl;
	int bestID, lastID, sx, sy;

	bool check_1st = true;

	for (int j = 0; j < (*matches_vec)[stereo_i].rows; j++)
	{
		lid11 = j;
		lid12 = (*matches_vec)[stereo_i].at<short>(j, 0) - 1;

		if (lid12 < 0)
			continue;

		// read 3D line
		x3i = (*pt3d_vec)[stereo_i].rowRange(j, j + 1).colRange(0, 4).reshape(0, 4);
		y3i = (*pt3d_vec)[stereo_i].rowRange(j, j + 1).colRange(4, 8).reshape(0, 4);
		v3i = x3i - y3i;

		if (im11 == im21)
		{
			lid21 = lid11;
			check_1st = false;
		}
		else if (im12 == im21)
		{
			lid21 = lid12;
			check_1st = false;
		}
		else if (im11 == im22)
		{
			lid21 = (*matches_r2l_vec)[stereo_j].at<short>(lid11, 0) - 1;
			check_1st = true;
		}
		else if (im12 == im22)
		{
			lid21 = (*matches_r2l_vec)[stereo_j].at<short>(lid12, 0) - 1;
			check_1st = true;
		}
		else
			continue;

		if (lid21 < 0 || lid21 >= (*matches_vec)[stereo_j].rows)continue;

		lid22 = (*matches_vec)[stereo_j].at<short>(lid21, 0) - 1;

		if (lid22 < 0 || lid22 >= (*lines_vec).at(im22).rows)
			continue;

		addAndRelease
		(
			&l3d_pair_vec,
			matches_vec,
			pt3d_vec,
			matches_r2l_vec,
			im11, im12, im21, im22,
			stereo_i,
			stereo_j,
			lid11,
			lid12
		);

		// check the line 3d angle
		x3j = (*pt3d_vec)[stereo_j].rowRange(lid21, lid21 + 1).colRange(0, 4);
		y3j = (*pt3d_vec)[stereo_j].rowRange(lid21, lid21 + 1).colRange(4, 8);
		v3j = x3j - y3j;
		cos_v3 = cos_vec3((float*)v3i.data, (float*)v3j.data);

		if (cos_v3 < min_line3D_cos)
			continue;//cos(pi/180*5)

		if (check_1st == true)
		{
			x1 = Pm * x3i;
			x1_p = (float*)x1.data;
			norm_by_v3(x1_p);
			y1 = Pm * y3i;
			y1_p = (float*)y1.data;
			norm_by_v3(y1_p);
			lf = x1.cross(y1);
			lf_p = (float*)lf.data;

			lineb1 = (float*)((*lines_vec)[im21].data) + lid21 * 11;
		}
		else
		{
			x1 = Pn * x3i;
			x1_p = (float*)x1.data;
			norm_by_v3(x1_p);
			y1 = Pn * y3i;
			y1_p = (float*)y1.data;
			norm_by_v3(y1_p);
			lf = x1.cross(y1);
			lf_p = (float*)lf.data;

			lineb1 = (float*)((*lines_vec)[im22].data) + lid22 * 11;
		}

		// distance check
		dis1 = point_2_line_dis(lineb1, lf_p);
		if (dis1 > dist)
			continue;

		dis2 = point_2_line_dis(lineb1 + 2, lf_p);
		if (dis2 > dist)
			continue;

		// intersection check 
		if (!twoLines_intersec(x1_p, y1_p, lineb1, lineb1 + 2, line_2_line_intersec))
			continue;

		max_dis = max_2(dis1, dis2);
		angle_3dl = acos(cos_v3);
		score = max_2
		(
			exp(-angle_3dl / (2 * min_line3D_ang)),
			exp(-max_dis / (2 * dist))
		);

		l3d_pair_vec.push_back({
		(int)(*pt3d_vec)[stereo_i].at<float>(lid11, 9),
		(int)(*pt3d_vec)[stereo_j].at<float>(lid21, 9),
		score });
	}

	return l3d_pair_vec;
}

std::vector< pair_m> checkStereos2(int stereo_i, int stereo_j,
	std::vector<cv::Mat>* matches_vec,
	std::vector<cv::Mat>* matches_r2l_vec,
	std::vector<cv::Mat>* pt3d_vec,
	std::vector<cv::Mat>* map_vec,
	std::vector<cv::Mat>* lines_vec,
	float* imidx_Mf_ptr,
	cv::Mat cameras_Mf,
	cv::Mat imsizes_Mf,
	float dist,
	int buffer
)
{
	cv::Mat x3i, y3i, x3j, y3j, v3i, v3j, x1, y1, lf, x2, y2, lf_2;
	int lid11, lid12, lid21, lid22;

	float best_score, cos_v3,
		dis1, dis2, dis3, dis4,
		max_dis, score, angle_3dl;

	std::vector< pair_m> l3d_pair_vec;

	std::vector<int> xx, yy;

	if ((*matches_vec)[stereo_j].rows == 0)
		return l3d_pair_vec;

	int im11 = imidx_Mf_ptr[stereo_i * 2];
	int im12 = imidx_Mf_ptr[stereo_i * 2 + 1];
	int im21 = imidx_Mf_ptr[stereo_j * 2];
	int im22 = imidx_Mf_ptr[stereo_j * 2 + 1];

	int rm1 = imsizes_Mf.at<float>(im21, 0);
	int cm1 = imsizes_Mf.at<float>(im21, 1);

	cv::Mat Pm = cameras_Mf.rowRange(im21, im21 + 1).reshape(0, 3);
	cv::Mat Pn = cameras_Mf.rowRange(im22, im22 + 1).reshape(0, 3);

	float* x1_p, * y1_p, * lf_p, * lineb1, * lineb2, * lf_p_another, * x2_p, * y2_p;
	float vec1[2];
	int addvec[2];

	int bestID, sx, sy;

	std::vector<short>searched_ID;
	bool buffer_control = 0;

	float min_line3D_ang = maximum_ang;
	float min_line3D_cos = cos(maximum_ang);

	for (int j = 0; j < (*matches_vec)[stereo_i].rows; j++)
	{
		lid11 = j;
		lid12 = (*matches_vec)[stereo_i].at<short>(j, 0) - 1;

		if (lid12 < 0)
			continue;

		// read 3D line
		x3i = (*pt3d_vec)[stereo_i].rowRange(j, j + 1).colRange(0, 4).reshape(0, 4);
		y3i = (*pt3d_vec)[stereo_i].rowRange(j, j + 1).colRange(4, 8).reshape(0, 4);
		v3i = x3i - y3i;

		// check the range
		x1 = Pm * x3i;
		x1_p = (float*)x1.data;
		norm_by_v3(x1_p);

		if (x1_p[0] <= 0 || x1_p[0] >= cm1 || x1_p[1] <= 0 || x1_p[1] >= rm1)
			continue;

		y1 = Pm * y3i;
		y1_p = (float*)y1.data;
		norm_by_v3(y1_p);
		if (y1_p[0] <= 0 || y1_p[0] >= cm1 || y1_p[1] <= 0 || y1_p[1] >= rm1)
			continue;

		lf = x1.cross(y1);
		lf_p = (float*)lf.data;

		vec1[0] = x1_p[0] - y1_p[0];
		vec1[1] = x1_p[1] - y1_p[1];

		// search in the map
		if (abs(vec1[0]) > abs(vec1[1]))
		{
			addvec[0] = 0;
			addvec[1] = 1;
		}
		else
		{
			addvec[0] = 1;
			addvec[1] = 0;
		}

		Bresenham(x1_p[0], x1_p[1], y1_p[0], y1_p[1], xx, yy);

		bestID = -1;
		best_score = 0;
		max_dis = dist;
		searched_ID.clear();

		for (int mm = 0; mm < xx.size(); mm += 5)
		{
			for (int nn = -buffer; nn <= buffer; nn++)
			{
				buffer_control = !buffer_control;
				if (buffer_control)
					continue;

				sx = xx[mm] + nn * addvec[0];
				sy = yy[mm] + nn * addvec[1];

				if (sx <= 0 || sx >= cm1 || sy <= 0 || sy >= rm1)
					continue;

				lid21 = (*map_vec)[im21].at<short>(sy, sx) - 1;
				if (lid21 < 0 || lid21 >= (*matches_vec)[stereo_j].rows)
					continue;

				lid22 = (*matches_vec)[stereo_j].at<short>(lid21, 0) - 1;
				if (lid22 < 0 || lid22 >= (*lines_vec).at(im22).rows)
					continue;// there is no match

				if (ID_in_array(searched_ID, lid21))
					continue;
				searched_ID.push_back(lid21);

				// check the line 3d angle
				x3j = (*pt3d_vec)[stereo_j].rowRange(lid21, lid21 + 1).colRange(0, 4);
				y3j = (*pt3d_vec)[stereo_j].rowRange(lid21, lid21 + 1).colRange(4, 8);
				v3j = x3j - y3j;
				cos_v3 = cos_vec3((float*)v3i.data, (float*)v3j.data);
				if (cos_v3 < min_line3D_cos)
					continue;//cos(pi/180*5)

				lineb1 = (float*)((*lines_vec)[im21].data) + lid21 * 11;

				// distance check
				dis1 = point_2_line_dis(lineb1, lf_p);
				if (dis1 > max_dis)
					continue;

				dis2 = point_2_line_dis(lineb1 + 2, lf_p);
				if (dis2 > max_dis)
					continue;

				// intersection check 
				if (!twoLines_intersec(x1_p, y1_p, lineb1, lineb1 + 2, line_2_line_intersec))
					continue;

				// then check the second line segment
				x2 = Pn * x3i;
				y2 = Pn * y3i;
				lf_2 = x2.cross(y2);

				lf_p_another = (float*)lf_2.data;
				x2_p = (float*)x2.data;
				y2_p = (float*)y2.data;
				norm_by_v3(x2_p);
				norm_by_v3(y2_p);

				lineb2 = (float*)(*lines_vec).at(im22).data + lid22 * 11;

				// distance check
				dis3 = point_2_line_dis(lineb2, lf_p_another);
				if (dis3 > max_dis)
					continue;

				dis4 = point_2_line_dis(lineb2 + 2, lf_p_another);
				if (dis4 > max_dis)
					continue;

				// intersection check 
				if (!twoLines_intersec(x2_p, y2_p, lineb2, lineb2 + 2, 0.5))
					continue;

				max_dis = max_2(max_2(dis1, dis2), max_2(dis3, dis4));
				angle_3dl = acos(cos_v3);
				score = max_2
				(
					exp(-angle_3dl / (2 * min_line3D_ang)),
					exp(-max_dis / (2 * dist))
				);

				if (best_score >= score)
					continue;

				bestID = lid21;
				best_score = score;
			}
		}

		if (bestID == -1)
		{
			continue;
		}

		l3d_pair_vec.push_back({
			(int)(*pt3d_vec)[stereo_i].at<float>(lid11, 9),
			(int)(*pt3d_vec)[stereo_j].at<float>(bestID, 9),
			best_score });
	}

	return l3d_pair_vec;
}

void merge_lines(std::vector<std::string>image_names, std::string input_folder,
	cv::Mat imidx_Mf, cv::Mat imsizes_Mf,
	cv::Mat cameras_Mf,
	float dist,
	ScoreRecorder* merge_tool)
{

	ThreadPool executor{ 20 };
	std::vector<std::future<std::vector<pair_m>>> match_task;
	std::vector<std::future<void>> read_task;

	int task_size = imidx_Mf.rows;
	float buffer = round(dist + 0.5);

	std::vector<std::vector<bool>>valid_pair(task_size);
	for (int i = 0; i < task_size; i++)
		valid_pair[i].resize(task_size);

	validPair(input_folder, imidx_Mf, valid_pair);

	int current_map_size = 0, current_line3D_size = 0,
		current_line2D_size = 0, current_match_size = 0;

	std::vector<cv::Mat>matches_vec(task_size);

	std::vector<cv::Mat>matches_r2l_vec(task_size);

	std::vector<cv::Mat>map_vec(image_names.size());

	std::vector<cv::Mat>pt3d_vec(task_size);

	std::vector<cv::Mat>lines_vec(image_names.size());

	int im11, im12, im21, im22, rm1, cm1;

	std::string name1, name2;

	float* imidx_Mf_ptr = (float*)imidx_Mf.data;

	for (int i = 0; i < imidx_Mf.rows; i++)
	{
		std::cout << imidx_Mf.rows - i << std::endl;

		im11 = imidx_Mf_ptr[i * 2];
		im12 = imidx_Mf_ptr[i * 2 + 1];

		for (int j = i; j < imidx_Mf.rows; j++)
		{
			if (valid_pair[i][j] == 0 && j != i)
				continue;

			im21 = imidx_Mf_ptr[j * 2];
			im22 = imidx_Mf_ptr[j * 2 + 1];

			// load lines
			if (lines_vec[im21].rows == 0)
			{
				name1 = image_names.at(im21);
				readMat(input_folder + "lines\\" + name1 + ".m", lines_vec[im21]);
				current_line2D_size++;
			}

			if (lines_vec[im22].rows == 0)
			{
				name1 = image_names.at(im22);
				readMat(input_folder + "lines\\" + name1 + ".m", lines_vec[im22]);
				current_line2D_size++;
			}

			// load 3d lines and matches
			if (pt3d_vec[j].rows == 0)
			{
				readMat(input_folder + "l3ds\\" + std::to_string(j) + ".m", pt3d_vec[j]);
				readMat(input_folder + "matches\\" + std::to_string(j) + ".m", matches_vec[j]);
				read_task.emplace_back(executor.commit(matchesRight2Left,
					matches_vec[j],
					lines_vec[im22], &matches_r2l_vec[j]));

				current_line3D_size++;
				current_match_size++;
			}

			// load maps
			if (map_vec[im21].rows == 0)
			{
				//create map
				rm1 = imsizes_Mf.at<float>(im21, 0);
				cm1 = imsizes_Mf.at<float>(im21, 1);
				read_task.emplace_back(executor.commit(lineMap,
					rm1, cm1, lines_vec[im21],
					&map_vec[im21]));
				current_map_size++;
			}
		}
		//wait
		for (auto& result : read_task)
			result.get();
		read_task.clear();

		// match for stereos

		for (int j = i + 1; j < imidx_Mf.rows; j++)
		{
			if (valid_pair[i][j] == 0)
				continue;

			im21 = imidx_Mf_ptr[j * 2];
			im22 = imidx_Mf_ptr[j * 2 + 1];

			if (im11 == im21 || im12 == im21 || im11 == im22 || im12 == im22)
			{
				match_task.emplace_back(executor.commit(checkStereos1,
					i, j,
					&matches_vec,
					&matches_r2l_vec,
					&pt3d_vec,
					&map_vec,
					&lines_vec,
					imidx_Mf_ptr,
					cameras_Mf,
					imsizes_Mf,
					dist,
					buffer
				));
			}
			else
			{
				match_task.emplace_back(executor.commit(checkStereos2,
					i, j,
					&matches_vec,
					&matches_r2l_vec,
					&pt3d_vec,
					&map_vec,
					&lines_vec,
					imidx_Mf_ptr,
					cameras_Mf,
					imsizes_Mf,
					dist,
					buffer
				));
			}
		}

		int cc = 0;
		for (int j = i + 1; j < imidx_Mf.rows; j++)
		{
			if (valid_pair[i][j] == 0)
				continue;

			merge_tool->add_connections(i, j, (&((std::vector< pair_m>)match_task[cc].get())));
			cc++;
		}
		match_task.clear();

		// release something for heavey task
		for (int mm = 0; mm <= i && current_map_size >= max_map_size; mm++)
		{
			if (map_vec[mm].rows != 0)
			{
				map_vec[mm].release();
				current_map_size--;
			}

		}

		for (int mm = 0; mm < i && current_line3D_size >= max_line3D_size; mm++)
		{
			if (pt3d_vec[mm].rows != 0)
			{
				pt3d_vec[mm].release();
				current_line3D_size--;
				matches_vec[mm].release();
				current_match_size--;
			}
		}

		for (int mm = 0; mm < i && current_line2D_size >= max_line2D_size; mm++)
		{
			im21 = imidx_Mf_ptr[mm * 2];
			im22 = imidx_Mf_ptr[mm * 2 + 1];

			if (lines_vec[im21].rows != 0)
			{
				lines_vec[im21].release();
				current_line2D_size--;
			}
			if (lines_vec[im22].rows != 0)
			{
				lines_vec[im22].release();
				current_line2D_size--;
			}
		}

	}

}