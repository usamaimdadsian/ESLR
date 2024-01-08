#pragma once
#pragma once
#include <opencv2/core.hpp>
#include <opencv2/opencv.hpp>
#include "ScoreRecorder.h"

void merge_lines(std::vector<std::string>image_names, std::string input_folder,
	cv::Mat imidx_Mf, cv::Mat imsizes_Mf,
	cv::Mat cameras_Mf,
	float dist,
	ScoreRecorder* merge_tool);

