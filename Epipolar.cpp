#include"Epipolar.h"

#include<iostream>
void skwe_mat3(cv::Mat M1f, cv::Mat& skew_Mf)
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



void computeFAE(cv::Mat cameras_all_Mf, cv::Mat imidx_Mf, cv::Mat F_Mf, cv::Mat AE_Mf, int impi)
{
	int id1, id2;
	id1 = imidx_Mf.at<float>(impi, 0);
	id2 = imidx_Mf.at<float>(impi, 1);

	cv::Mat P1 = cameras_all_Mf.rowRange(id1, id1 + 1).clone().reshape(0, 3);
	cv::Mat P2 = cameras_all_Mf.rowRange(id2, id2 + 1).clone().reshape(0, 3);

	cv::Mat P1_inv;
	cv::invert(P1, P1_inv, cv::DECOMP_SVD);

	cv::Mat w, u, vt, e;
	//svd  Decomp for null space
	cv::SVD::compute(P1, w, u, vt,cv::SVD::FULL_UV);
	//cv::SVDecomp(P1, w, u, vt);

	cv::Mat C1 = vt.rowRange(3, 4).clone().t();
	
	cv::Mat x = P2 * C1;

	skwe_mat3(x, e);

	cv::Mat F = e * (P2 * P1_inv);

	float* f = (float*)F.data;


	for (int mm = 0; mm < 9; mm++)
		F_Mf.at<float>(impi, mm) = f[mm];

	//AE
	cv::SVD::compute(F, w, u, vt, cv::SVD::FULL_UV);
	cv::Mat eprime = u.colRange(2, 3).clone();
	cv::Mat e_prime_cross;
	skwe_mat3(eprime, e_prime_cross);

	cv::Mat A= e_prime_cross * F;

	float* a = (float*)A.data;
	float* eprime_ = (float*)eprime.data;

	for (int mm=0;mm<9;mm++)
		AE_Mf.at<float>(impi, mm) = a[mm];
	for (int mm = 0; mm < 3; mm++)
		AE_Mf.at<float>(impi, mm+9) = eprime_[mm];
}

void computeFAE(cv::Mat P1, cv::Mat P2,cv::Mat& F_Mf, cv::Mat& AE_Mf)
{
	

	cv::Mat P1_inv;
	cv::invert(P1, P1_inv, cv::DECOMP_SVD);

	cv::Mat w, u, vt, e;
	//svd  Decomp for null space
	cv::SVD::compute(P1, w, u, vt, cv::SVD::FULL_UV);
	//cv::SVDecomp(P1, w, u, vt);

	cv::Mat C1 = vt.rowRange(3, 4).clone().t();

	cv::Mat x = P2 * C1;

	skwe_mat3(x, e);

	cv::Mat F = e * (P2 * P1_inv);
	F.copyTo(F_Mf);

	//AE
	cv::SVD::compute(F, w, u, vt, cv::SVD::FULL_UV);
	cv::Mat eprime = u.colRange(2, 3).clone();
	cv::Mat e_prime_cross;
	skwe_mat3(eprime, e_prime_cross);

	cv::Mat A = e_prime_cross * F;
	cv::Mat AE = cv::Mat::zeros(1, 12,CV_32FC1);

	float* a = (float*)A.data;
	float* eprime_ = (float*)eprime.data;

	for (int mm = 0; mm < 9; mm++)
		AE.at<float>(0, mm) = a[mm];
	for (int mm = 0; mm < 3; mm++)
		AE.at<float>(0, mm + 9) = eprime_[mm];

	AE.copyTo(AE_Mf);
}