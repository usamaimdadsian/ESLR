#include "SingleImage.h"
#include <thread>
#include <math.h>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/types_c.h>
#include "lsd.hpp"
#include "BasicMath.h"
#include <fstream>
#include "Parameters.h"

#include "opencv2/imgproc.hpp"
#include "opencv2/ximgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"

#include"IO.h"

using namespace cv::ximgproc;

void points2KDtree(cv::Mat inter_lines_Mf, cv::Mat& inter_knn_Mi, cv::Mat lines_Mf,
	cv::Mat& lines_knn_Mi_1_, cv::Mat& lines_knn_Mi_2_,
	int support_pt_num1, int support_pt_num2, std::string pt2_file_name, float imr, float imc)
{
	std::string pt_file_adr;
	std::ifstream pt_file;
	std::string pt_line;

	cv::Mat pt2_Mf;
	cv::Mat pt2_idpt3__Mf;
	cv::Mat per_row_pt(1, 2, CV_32FC1);
	cv::Mat per_row_id(1, 1, CV_32FC1);

	cv::Mat pt2;
	readMat(pt2_file_name+".pt2", pt2);

	pt2_Mf = pt2.colRange(0, 2);
	pt2_idpt3__Mf = pt2.colRange(2,3);

	

	float imr_2 = imr / 2.0;
	float imc_2 = imc / 2.0;

	for (int i = 0; i < pt2_Mf.rows; i++)
	{
		pt2_Mf.at<float>(i, 0) = pt2_Mf.at<float>(i, 0) + imc_2;
		pt2_Mf.at<float>(i, 1) = pt2_Mf.at<float>(i, 1) + imr_2;
	}

	
	// construct kdtree
	cv::flann::Index flannIndex(pt2_Mf, cv::flann::KDTreeIndexParams());

	// intersection KD tree
	cv::Mat inter_2_pt3index_M;
	cv::Mat inter_2_pt3index_dist;
	cv::Mat inter_pts = inter_lines_Mf.colRange(2, 4).clone();

	flannIndex.knnSearch(inter_pts, inter_2_pt3index_M,
		inter_2_pt3index_dist, support_pt_num1, cv::flann::SearchParams
		());

	cv::Mat  inter_knn_Mi_ = cv::Mat(inter_pts.rows, support_pt_num1, CV_32SC1);

	int* inter_knn_Mi_ptr = (int*)inter_knn_Mi_.data;
	int* inter_2_pt3index_M_ptr = (int*)inter_2_pt3index_M.data;

	int ind_inter = 0;
	for (int i = 0; i < inter_pts.rows; i++)
	{
		for (int j = 0; j < support_pt_num1; j++)
		{

			inter_knn_Mi_ptr[ind_inter] =
				pt2_idpt3__Mf.at<float>(inter_2_pt3index_M_ptr[ind_inter], 0) - 1;
			//-1 is aligned to main_vsfm.m file
			ind_inter++;
		}

	}

	

	inter_knn_Mi_.copyTo(inter_knn_Mi);

	// line segment KD tree
	cv::Mat line_2_pt3index_1;
	cv::Mat line_2_pt3index_dist1;
	cv::Mat line_pts1 = lines_Mf.colRange(0, 2).clone();
	flannIndex.knnSearch(line_pts1, line_2_pt3index_1,
		line_2_pt3index_dist1, support_pt_num2, cv::flann::SearchParams());

	cv::Mat line_2_pt3index_2;
	cv::Mat line_2_pt3index_dist2;
	cv::Mat line_pts2 = lines_Mf.colRange(2, 4).clone();
	flannIndex.knnSearch(line_pts2, line_2_pt3index_2,
		line_2_pt3index_dist2, support_pt_num2, cv::flann::SearchParams());

	//store
	cv::Mat lines_knn_Mi_1 = cv::Mat(lines_Mf.rows, support_pt_num2, CV_32SC1);
	cv::Mat lines_knn_Mi_2 = cv::Mat(lines_Mf.rows, support_pt_num2, CV_32SC1);

	int* lines_knn_Mi_1_ptr = (int*)lines_knn_Mi_1.data;
	int* lines_knn_Mi_2_ptr = (int*)lines_knn_Mi_2.data;

	int* line_2_pt3index_1_ptr = (int*)line_2_pt3index_1.data;
	int* line_2_pt3index_2_ptr = (int*)line_2_pt3index_2.data;

	int  ind_ptr = 0;
	for (int i = 0; i < lines_Mf.rows; i++)
		for (int j = 0; j < support_pt_num2; j++)
		{

			lines_knn_Mi_1_ptr[ind_ptr] =
				pt2_idpt3__Mf.at<float>(line_2_pt3index_1_ptr[ind_ptr], 0) - 1;

			lines_knn_Mi_2_ptr[ind_ptr] =
				pt2_idpt3__Mf.at<float>(line_2_pt3index_2_ptr[ind_ptr], 0) - 1;
			//-1 is aligned to vsfm.m file
			ind_ptr = ind_ptr + 1;
		}

	lines_knn_Mi_1.copyTo(lines_knn_Mi_1_);
	lines_knn_Mi_2.copyTo(lines_knn_Mi_2_);
}

void readLsdLines(cv::Mat& lines_Mf, std::string line_file_adr)
{

	std::ifstream line_file;
	line_file.open(line_file_adr.c_str());
	std::string line_line;

	float cx, cy, x1, y1, x2, y2, length;

	cv::Mat line_Mf_;
	cv::Mat per_row(1, 7, CV_32FC1);
	float* per_row_p = (float*)per_row.data;
	while (std::getline(line_file, line_line))
	{
		std::istringstream iss_imidx(line_line);

		iss_imidx >> x1;
		iss_imidx >> y1;

		iss_imidx >> x2;
		iss_imidx >> y2;

		length = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

		cx = (x1 + x2) / 2;
		cy = (y1 + y2) / 2;

		per_row_p[0] = x1;
		per_row_p[1] = y1;
		per_row_p[2] = x2;
		per_row_p[3] = y2;
		per_row_p[4] = cx;
		per_row_p[5] = cy;
		per_row_p[6] = length;

		line_Mf_.push_back(per_row);
	}

	line_file.close();

	line_Mf_.copyTo(lines_Mf);
}

void maxMinDepth(int support_num, int lp_index, float* spacepoints,
	int* inter_knn, float* CM,
	float& mindepth, float& maxdepth)
{
	mindepth = 99999;
	maxdepth = 0;
	float w;
	int ind, indpt3;
	ind = lp_index * support_num;

	//std::cout << lp_index << " ";
	for (int i = 0; i < support_num; i++)
	{

		//std::cout << inter_knn[ind + i] << " ";
		indpt3 = inter_knn[ind + i] * 3;

		w = CM[11]
			+ CM[8] * spacepoints[indpt3]
			+ CM[9] * spacepoints[indpt3 + 1]
			+ CM[10] * spacepoints[indpt3 + 2];

		if (maxdepth < w)
			maxdepth = w;

		if (mindepth > w)
			mindepth = w;

	}

}

void maxMinPt(int feature_size, int support_num, float* spacepoints,
	int* inter_knn, float* CM, float* out_range)
{
	float mindepth, maxdepth;
	float w;
	int  ind, indpt3;
	for (int i = 0; i < feature_size; i++)
	{
		mindepth = 99999;
		maxdepth = 0;

		ind = i * support_num;

		//std::cout << lp_index << " ";
		for (int j = 0; j < support_num; j++)
		{

			//std::cout << inter_knn[ind + i] << " ";
			indpt3 = inter_knn[ind + j] * 3;

			w = CM[11]
				+ CM[8] * spacepoints[indpt3]
				+ CM[9] * spacepoints[indpt3 + 1]
				+ CM[10] * spacepoints[indpt3 + 2];

			if (maxdepth < w)
				maxdepth = w;

			if (mindepth > w)
				mindepth = w;

		}

		out_range[i * 2] = mindepth;
		out_range[i * 2 + 1] = maxdepth;

	}
}

float siftPixelAng(cv::Mat CM, float imr1, float imc1, float shift_p)
{
	float M[] = { CM.at<float>(0,0),CM.at<float>(0,1),CM.at<float>(0,2),
				  CM.at<float>(1,0),CM.at<float>(1,1),CM.at<float>(1,2),
				  CM.at<float>(2,0),CM.at<float>(2,1),CM.at<float>(2,2) };

	imr1 = imr1 / 2;;
	imc1 = imc1 / 2;;

	float v1[3];
	float v2[3];

	float size_v1[] = { imc1,imr1,1 };
	float size_v2[] = { size_v1[0] + shift_p,size_v1[1],1 };

	M_divide_b(M, size_v1, v1);
	M_divide_b(M, size_v2, v2);

	float cos_v = cos_vec3(v1, v2);
	return  1 - cos_v * cos_v;
}

void dect_lsd_lines(cv::Mat img, cv::Mat& lines_Mf_, float scale)
{
	cv::Mat img_double;
	img.convertTo(img_double, CV_32FC1);

	int line_num;
	float* img_ptr = (float*)img_double.data;
	float* line_lsd_ptr = lsd_scale(&line_num, img_ptr, img_double.cols, img_double.rows, 1.0) - 1;
	// note -1 is minored for iteration

	// add the lines to opencv Mat
	cv::Mat lines_Mf = cv::Mat(line_num, 7, CV_32FC1);
	float* lines_Mf_ptr = (float*)lines_Mf.data - 1;
	// note -1 is minored for iteration
	float x1, y1, x2, y2, cx, cy, length;
	for (int j = 0; j < line_num; j++)
	{
		x1 = *(++line_lsd_ptr) / scale;
		y1 = *(++line_lsd_ptr) / scale;
		x2 = *(++line_lsd_ptr) / scale;
		y2 = *(++line_lsd_ptr) / scale;
		cx = *(++line_lsd_ptr) / scale;
		cy = *(++line_lsd_ptr) / scale;
		++line_lsd_ptr;

		length = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
		length = std::sqrt(length);

		*(++lines_Mf_ptr) = x1;
		*(++lines_Mf_ptr) = y1;
		*(++lines_Mf_ptr) = x2;
		*(++lines_Mf_ptr) = y2;
		*(++lines_Mf_ptr) = cx;
		*(++lines_Mf_ptr) = cy;
		*(++lines_Mf_ptr) = length;
	}
	lines_Mf.copyTo(lines_Mf_);
}

void processImage(cv::Mat CM_Mf, std::vector<cv::Mat> camera_Rts, std::vector<float> cams_focals,
	cv::Mat space_points_Mf, cv::Mat imsizes_Mf_, int i, 
	std::vector<std::string>image_names, std::string input_folder,
	float costhre, float dist, int inter_support_num, int maxwidth)
{

	//1 read images
	std::cout << image_names.size() - i << " images to process\n";
	std::string image_name = image_names.at(i);
	cv::Mat img = cv::imread(input_folder + "/"+image_name, 0);

	// scale images
	// check image size
	int max_dim = std::max(img.rows, img.cols);
	unsigned int new_width = img.cols;
	unsigned int new_height = img.rows;

	imsizes_Mf_.at<float>(i, 0) = img.rows;
	imsizes_Mf_.at<float>(i, 1) = img.cols;

	// check if exist
	/*
	std::fstream _file;
	_file.open(input_folder + "/ELSR/lines/" + image_name + ".m", std::ios::in);
	if (_file)
	{
		_file.close();
		return;
	}
	_file.close();
	*/

	cv::Mat imgResized;
	float scale;
	if (maxwidth > 0 && max_dim > maxwidth)
	{
		// rescale
		scale = float(maxwidth) / float(max_dim);

		scale = scale * 0.8; // for LSD
		cv::resize(img, imgResized, cv::Size(), scale, scale);
	}
	else
	{
		scale = 0.8; // for LSD
		cv::resize(img, imgResized, cv::Size(), scale, scale);
	}

	// minum line length for intersection line
	//2 detect lines
	cv::Mat lines_Mf;
	dect_lsd_lines(imgResized, lines_Mf, scale);

	//3 detect intersection
	clock_t start = clock();
	cv::Mat inter_lines_Mf;

	callCrossPt(inter_lines_Mf, (float*)lines_Mf.data, lines_Mf.rows, costhre, dist);

	//4 load points and query knn points for intersections
	cv::Mat inter_knn_Mi, line_knn_Mi_1, line_knn_Mi_2;
	points2KDtree(inter_lines_Mf, inter_knn_Mi, lines_Mf, line_knn_Mi_1, line_knn_Mi_2,
		inter_support_num, inter_support_num * 2, input_folder +"/ELSR/"+ image_name, img.rows, img.cols);

	//5 compute the point depth for line and inter

	//5.1 calculate cameras
	cv::Mat K=cv::Mat::zeros(3, 3, CV_32FC1);
	K.at<float>(0, 0) = cams_focals[i];
	K.at<float>(1, 1) = cams_focals[i];
	K.at<float>(0, 2) = img.cols / 2.0;
	K.at<float>(1, 2) = img.rows / 2.0;
	K.at<float>(2, 2) = 1;

	cv::Mat CM_i = K * camera_Rts[i];
	

	for (int k = 0; k < 12; k++)
		CM_Mf.at<float>(i,k)=((float*)CM_i.data)[k];

	cv::Mat inter_max_min = cv::Mat(inter_lines_Mf.rows, 2, CV_32FC1);
	maxMinPt(inter_lines_Mf.rows, inter_knn_Mi.cols, (float*)space_points_Mf.data,
		(int*)inter_knn_Mi.data, (float*)CM_i.data, (float*)inter_max_min.data);

	cv::Mat line_max_min_1 = cv::Mat(lines_Mf.rows, 2, CV_32FC1);
	cv::Mat line_max_min_2 = cv::Mat(lines_Mf.rows, 2, CV_32FC1);

	maxMinPt(lines_Mf.rows, line_knn_Mi_1.cols, (float*)space_points_Mf.data,
		(int*)line_knn_Mi_1.data, (float*)CM_i.data, (float*)line_max_min_1.data);

	maxMinPt(lines_Mf.rows, line_knn_Mi_2.cols, (float*)space_points_Mf.data,
		(int*)line_knn_Mi_2.data, (float*)CM_i.data, (float*)line_max_min_2.data);

	//sift pixels 
	float sin_beta_2 = siftPixelAng(CM_i, img.rows, img.cols, depth_shift_pixel);
	float m_depth;

	for (int j = 0; j < line_max_min_1.rows; j++)
	{
		m_depth = line_max_min_1.at<float>(j, 0) * sin_beta_2;
		line_max_min_1.at<float>(j, 0) -= m_depth;
		line_max_min_1.at<float>(j, 1) += m_depth;

		m_depth = line_max_min_2.at<float>(j, 0) * sin_beta_2;
		line_max_min_2.at<float>(j, 0) -= m_depth;
		line_max_min_2.at<float>(j, 1) += m_depth;
	}

	cv::Mat line_max_min;
	hconcat(line_max_min_1, line_max_min_2, line_max_min);
	cv::Mat lines_ranges;
	hconcat(lines_Mf, line_max_min, lines_ranges);

	cv::Mat inter_ranges;
	hconcat(inter_lines_Mf, inter_max_min, inter_ranges);
	
	saveMat(input_folder + "/ELSR/lines/" + image_name + ".m", lines_ranges);
	saveMat(input_folder + "/ELSR/l2l/" + image_name + ".m", inter_ranges);
	//saveMat(input_folder + "maps\\" + image_name + ".m", search_map_Ms);
}
