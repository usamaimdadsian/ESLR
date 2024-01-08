#pragma once
#include <boost/filesystem.hpp>
#include <opencv2/core.hpp>
#include <iostream>
#include "IO.h"
int read_VisualSfM(std::string inputFolder, std::string nvmFile,
	std::vector<std::string>& image_names,
	std::vector<float>& cams_focals,
	std::vector<cv::Mat>& cams_RT,
	cv::Mat& points_space3D,
	cv::Mat& imidx_Mf_,
	int knn_image);
