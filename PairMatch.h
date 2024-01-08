#pragma once
#include <opencv2/core.hpp>
#include <opencv2/opencv.hpp>
#include "ScoreRecorder.h"

void findHomography(cv::Mat inter_range_i,
	cv::Mat lpsm_Mf, cv::Mat lpsn_Mf, cv::Mat plMap_Mf, cv::Mat lm_Mf, cv::Mat ln_Mf,
	cv::Mat CM_Mf, cv::Mat CN_Mf, cv::Mat F_Mf, cv::Mat Ae_Mf, cv::Mat& matches_M_float,
	int imr, int imc,
	float error_max,
	float ang_max);

void matchPair(std::vector<std::string>image_names, std::string input_folder,
	cv::Mat imidx_Mf, cv::Mat imsizes_Mf,
	cv::Mat cameras_Mf,
	int match_ind,
	ScoreRecorder* merge_tool,
	float epipolar_ang,
	float error_max,
	float ang_max);
