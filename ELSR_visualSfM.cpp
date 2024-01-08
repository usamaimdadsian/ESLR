//std
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <ctime>
#include <sys/stat.h>
//mine
#include "SingleImage.h"
#include "PairMatch.h"
#include "SpatialSimilarity.h"
#include "ScoreRecorder.h"
#include "Parameters.h"
#include "ThreadPool.h"
#include"read_VisualSfM.h"
#include"IO.h"



int main()
{
	clock_t start, end;
	start = clock();

	std::string inputFolder = "F:\\line3D_test_images\\P25";
	std::string nvmFile = "\\res.nvm";
	int knn_image=3;

	std::vector<std::string> image_names;
	std::vector<float> cams_focals;
	std::vector<cv::Mat> cams_RT;
	cv::Mat points_space3D;
	cv::Mat imidx_Mf;

	// read .nvm of VisualSfM
	read_VisualSfM(inputFolder, nvmFile,
		image_names,
		cams_focals,
		cams_RT,
		points_space3D,
		imidx_Mf,
		knn_image);
	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "process visual sfm "
		<< endtime << "s" << std::endl;

	//main process
	cv::Mat  imsizes_Mf = cv::Mat(image_names.size(), 2, CV_32FC1);
	cv::Mat cameras_Mf = cv::Mat(image_names.size(), 12, CV_32FC1);
	

	//1 Process single images,run in multiple threads 
	ThreadPool executor{16};
	std::vector<std::future<void>> results;
	for (int i = 0; i < image_names.size(); i++)
		results.emplace_back(executor.commit(processImage,
			cameras_Mf, cams_RT, cams_focals,
			points_space3D, imsizes_Mf, i,
			image_names, inputFolder,
			intersect_cos, intersect_dist, support_pt_num, max_image_width));
	for (auto& result : results)
		result.get();

	results.clear();

	

	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "done in " << endtime << " s" << std::endl;
	points_space3D.release();


	//2  match with imidx_Mf
	ScoreRecorder merge_tool(imidx_Mf.rows, min_epipolar_ang, min_pairs_num);
	
	std::cout << "2->double image process->\n";

	for (int i = 0; i < imidx_Mf.rows; i++)//imidx_Mf.rows
	{
		results.emplace_back(executor.commit(matchPair,
			image_names, inputFolder+"/ELSR/",
			imidx_Mf, imsizes_Mf,
			cameras_Mf, i,
			&merge_tool,
			min_epipolar_ang,
			pixel_2_line_dis,
			maximum_ang));
	}
	for (auto& result : results)
		result.get();
	results.clear();

	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "done in " << endtime << " s" << std::endl;

	std::cout << "3->merge lines for stereos->\n";
	merge_lines(image_names, inputFolder + "/ELSR/",
		imidx_Mf, imsizes_Mf,
		cameras_Mf,
		pixel_2_line_dis,
		&merge_tool);
	
	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "done in " << endtime << " s" << std::endl;

	merge_tool.vote_lines();

	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "total run time: " << endtime << " s" << std::endl;
	std::cout << "resconstruct 3D lines: " << merge_tool.ks_store.size() << std::endl;

	write2obj(inputFolder + "/ELSR/",
		"ELSR.obj",
		merge_tool.ks_store,
		merge_tool.match_ids_store,
		imidx_Mf.rows);

	std::getchar();
}




