#pragma once
#include <opencv2/core.hpp>
void computeFAE(cv::Mat cameras_all_Mf, cv::Mat imidx_Mf, cv::Mat F_Mf, cv::Mat AE_Mf, int impi);

void computeFAE(cv::Mat P1, cv::Mat P2, cv::Mat& F_Mf, cv::Mat& AE_Mf);