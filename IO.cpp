#include "IO.h"

void saveMat(std::string save_name, cv::Mat m)
{
	std::ofstream ofs(save_name, std::ios::out | std::ios::binary);
	boost::archive::binary_oarchive oa(ofs);
	oa << m;
	ofs.close();
}

void readMat(std::string name, cv::Mat& mat)
{

	std::ifstream inf1(name, std::ios::in | std::ios::binary);
	{
		boost::archive::binary_iarchive ia(inf1);
		ia >> mat;
	}
}

void write2obj(std::string input_folder, std::string outname, std::vector<unsigned int>& ks,
	std::vector<unsigned int>& ids, int pairsize)
{
	std::ofstream outfile;
	outfile.open(input_folder+outname);
	// print the result
	cv::Mat pt3ds;
	std::string name;

	int i, j;

	for (int m = 0; m < pairsize; m++)
	{
		name = input_folder + "l3ds\\" + std::to_string(m) + ".m";

		readMat(name, pt3ds);

		for (int n = 0; n < ks.size(); n++)
		{
			if (ks.at(n) != m)
				continue;

			i = m;
			j = ids.at(n);

			outfile << "v " << pt3ds.at<float>(j, 0) << " " << pt3ds.at<float>(j, 1) << " "
				<< pt3ds.at<float>(j, 2) << "\n";

			outfile << "v " << pt3ds.at<float>(j, 4) << " " << pt3ds.at<float>(j, 5) << " "
				<< pt3ds.at<float>(j, 6) << "\n";

		}

	}

	for (int i0; i < ks.size(); i++)
		outfile << "l " << i * 2 - 1 << " " << i * 2 << "\n";

	outfile.close();
}