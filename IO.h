#pragma once
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <opencv2/core.hpp>
#include <fstream>

BOOST_SERIALIZATION_SPLIT_FREE(cv::Mat)
namespace boost {
	namespace serialization {

		/** Serialization support for cv::Mat */
		template<class Archive>
		void save(Archive& ar, const cv::Mat& m, const unsigned int version)
		{
			size_t elem_size = m.elemSize();
			size_t elem_type = m.type();

			ar& m.cols;
			ar& m.rows;
			ar& elem_size;
			ar& elem_type;

			const size_t data_size = m.cols * m.rows * elem_size;
			ar& boost::serialization::make_array(m.ptr(), data_size);
		}

		/** Serialization support for cv::Mat */
		template <class Archive>
		void load(Archive& ar, cv::Mat& m, const unsigned int version)
		{
			int cols, rows;
			size_t elem_size, elem_type;

			ar& cols;
			ar& rows;
			ar& elem_size;
			ar& elem_type;

			m.create(rows, cols, elem_type);

			size_t data_size = m.cols * m.rows * elem_size;
			ar& boost::serialization::make_array(m.ptr(), data_size);
		}

	}
}

void saveMat(std::string save_name, cv::Mat m);
void readMat(std::string name, cv::Mat& mat);
void write2obj(std::string input_folder, std::string outname, std::vector<unsigned int>& ks,
	std::vector<unsigned int>& ids, int pairsize);