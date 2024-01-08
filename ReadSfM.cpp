#include"read_VisualSfM.h"

void analysis_match(cv::Mat mscores, cv::Mat& imidx_Mf, int knn, int minconnect)
{

	int ind;


	cv::Mat per_pair(1, 2, CV_32FC1);

	for (int i = 0; i < mscores.cols; i++)
	{
		cv::Mat sortID;
		cv::sortIdx(mscores.row(i), sortID, cv::SORT_DESCENDING);
		for (int j = 0; j < knn; j++)
		{
			ind = sortID.at<int>(0, j);
			if (mscores.at<float>(i, ind) < minconnect)
				break;
			per_pair.at<float>(0, 0) = i;
			per_pair.at<float>(0, 1) = ind;
			imidx_Mf.push_back(per_pair);

			mscores.at<float>(ind, i) = 0;
		}
	}
}

void analysis_pair(cv::Mat imidx_Mf, cv::Mat mscores, cv::Mat& pairs_double_, int mincommon)
{
	int ind1, ind2, ind3, ind4;

	for (int i = 0; i < mscores.cols; i++)
		mscores.at<float>(i, i) = 9999;

	cv::Mat pairs_double = cv::Mat::zeros(imidx_Mf.rows, imidx_Mf.rows, CV_16UC1);

	for (int i = 0; i < imidx_Mf.rows; i++)
	{
		ind1 = imidx_Mf.at<float>(i, 0);
		ind2 = imidx_Mf.at<float>(i, 1);

		for (int j = i + 1; j < imidx_Mf.rows; j++)
		{

			ind3 = imidx_Mf.at<float>(j, 0);
			ind4 = imidx_Mf.at<float>(j, 1);

			if (mscores.at<float>(ind1, ind3) < mincommon ||
				mscores.at<float>(ind1, ind4) < mincommon ||
				mscores.at<float>(ind2, ind3) < mincommon ||
				mscores.at<float>(ind2, ind4) < mincommon)
				continue;

			pairs_double.at<ushort>(i, j) = 1;
			pairs_double.at<ushort>(j, i) = 1;

		}

	}

	pairs_double.copyTo(pairs_double_);

}

int read_VisualSfM(std::string inputFolder, std::string nvmFile,
	std::vector<std::string>& image_names,
	std::vector<float>& cams_focals,
	std::vector<cv::Mat>& cams_RT,
	cv::Mat& points_space3D,
	cv::Mat& imidx_Mf_,
	int knn_image)
{
	nvmFile = inputFolder + nvmFile;
	// check if NVM file exists
	boost::filesystem::path nvm(nvmFile);
	if (!boost::filesystem::exists(nvm))
	{
		std::cerr << "NVM file " << nvmFile << " does not exist!" << std::endl;
		return -1;
	}

	// create output directory
	boost::filesystem::create_directory
	(boost::filesystem::path(inputFolder + "/ELSR/"));

	boost::filesystem::create_directory
	(boost::filesystem::path(inputFolder + "/ELSR/l2l"));

	boost::filesystem::create_directory
	(boost::filesystem::path(inputFolder + "/ELSR/l3ds"));

	boost::filesystem::create_directory
	(boost::filesystem::path(inputFolder + "/ELSR/matches"));

	boost::filesystem::create_directory
	(boost::filesystem::path(inputFolder + "/ELSR/lines"));


	// read NVM file
	std::ifstream nvm_file;
	nvm_file.open(nvmFile.c_str());

	std::string nvm_line;
	std::getline(nvm_file, nvm_line); // ignore first line...
	std::getline(nvm_file, nvm_line); // ignore second line...

	// read number of images
	std::getline(nvm_file, nvm_line);
	std::stringstream nvm_stream(nvm_line);
	unsigned int num_cams;
	nvm_stream >> num_cams;

	if (num_cams == 0)
	{
		std::cerr << "No aligned cameras in NVM file!" << std::endl;
		return -1;
	}

	// read camera data (sequentially)
	image_names.resize(num_cams);
	cams_focals.resize(num_cams);
	cams_RT.resize(num_cams);
	cams_focals.resize(num_cams);
	cv::Mat imidx_Mf;
	


	cv::Mat R(3, 3, CV_32FC1);
	cv::Mat C(3, 1, CV_32FC1);

	for (unsigned int i = 0; i < num_cams; ++i)
	{
		std::getline(nvm_file, nvm_line);

		// image filename
		std::string filename;

		// focal_length,quaternion,center,distortion
		double focal_length, qx, qy, qz, qw;
		double Cx, Cy, Cz, dist;

		nvm_stream.str("");
		nvm_stream.clear();
		nvm_stream.str(nvm_line);
		nvm_stream >> filename >> focal_length >> qw >> qx >> qy >> qz;
		nvm_stream >> Cx >> Cy >> Cz >> dist;

		image_names[i] = filename;
		cams_focals[i] = focal_length;

		// rotation amd translation
		R.at<float>(0, 0) = 1.0 - 2.0 * qy * qy - 2.0 * qz * qz;
		R.at<float>(0, 1) = 2.0 * qx * qy - 2.0 * qz * qw;
		R.at<float>(0, 2) = 2.0 * qx * qz + 2.0 * qy * qw;

		R.at<float>(1, 0) = 2.0 * qx * qy + 2.0 * qz * qw;
		R.at<float>(1, 1) = 1.0 - 2.0 * qx * qx - 2.0 * qz * qz;
		R.at<float>(1, 2) = 2.0 * qy * qz - 2.0 * qx * qw;

		R.at<float>(2, 0) = 2.0 * qx * qz - 2.0 * qy * qw;
		R.at<float>(2, 1) = 2.0 * qy * qz + 2.0 * qx * qw;
		R.at<float>(2, 2) = 1.0 - 2.0 * qx * qx - 2.0 * qy * qy;

		C.at<float>(0, 0) = Cx;
		C.at<float>(1, 0) = Cy;
		C.at<float>(2, 0) = Cz;

		cv::Mat t = -R * C;

		cv::Mat Rt;
		cv::hconcat(R, t, Rt);
		cams_RT[i] = Rt;
	}

	// read number of images
	std::getline(nvm_file, nvm_line); // ignore line...
	std::getline(nvm_file, nvm_line);
	nvm_stream.str("");
	nvm_stream.clear();
	nvm_stream.str(nvm_line);
	unsigned int num_points;
	nvm_stream >> num_points;

	// read features (for image similarity calculation)
	cv::Mat pos3D(1, 3, CV_32FC1);

	std::vector<cv::Mat>points2D_N(num_cams);
	cv::Mat mscores = cv::Mat::zeros(num_cams, num_cams, CV_32FC1);

	std::vector<uint>cam_IDs;

	for (unsigned int i = 0; i < num_points; ++i)
	{
		// 3D position
		std::getline(nvm_file, nvm_line);
		std::istringstream iss_point3D(nvm_line);
		double px, py, pz, colR, colG, colB;
		iss_point3D >> px >> py >> pz;
		iss_point3D >> colR >> colG >> colB;

		pos3D.at<float>(0, 0) = px;
		pos3D.at<float>(0, 1) = py;
		pos3D.at<float>(0, 2) = pz;

		points_space3D.push_back(pos3D);

		//  num of views for each 3D points
		unsigned int num_views;
		iss_point3D >> num_views;

		unsigned int camID, siftID;
		float posX, posY;

		for (unsigned int j = 0; j < num_views; ++j)
		{
			iss_point3D >> camID >> siftID;
			iss_point3D >> posX >> posY;

			cam_IDs.push_back(camID);

			pos3D.at<float>(0, 0) = posX;
			pos3D.at<float>(0, 1) = posY;
			pos3D.at<float>(0, 2) = i+1;

			points2D_N[camID].push_back(pos3D);
		}

		for (int ii = 0; ii < cam_IDs.size(); ii++)
			for (int jj = ii + 1; jj < cam_IDs.size(); jj++)
			{
				mscores.at<float>(cam_IDs[ii], cam_IDs[jj])++;
				mscores.at<float>(cam_IDs[jj], cam_IDs[ii])++;
			}

		cam_IDs.clear();
	}

	// store 2D points
	for (int i = 0; i < num_cams; i++)
		saveMat(inputFolder + "/ELSR/"
			+ image_names[i] + ".pt2",
			points2D_N[i]);

	//std::cout << points2D_N[0].rowRange(0, 10);
	//analysis match pair
	cv::Mat  pairs_double;
	int mincommon = 100;
	analysis_match(mscores.clone(), imidx_Mf, knn_image, mincommon);

	analysis_pair(imidx_Mf, mscores, pairs_double, 200);

	
	saveMat(inputFolder + "/ELSR/"
		+ "valid_pair.m",
		pairs_double);
	nvm_file.close();

	//std::cout << pairs_double << std::endl;

	imidx_Mf.copyTo(imidx_Mf_);

	return 1;
}

