#pragma once
#include<vector>

struct pair_m
{
	int l3d_i, l3d_j;
	float score;
};

struct TASK
{
	unsigned int l3d_size;
	std::vector<unsigned int> ori_ID;
	std::vector<float>epipolar_ang;
	std::vector<int>consistent_size;
	std::vector<float> scores;
	std::vector<size_t> sort_ids;

	unsigned int sort_id_add = 0;

	std::vector< std::vector<unsigned int>> dst_k;
	std::vector< std::vector<unsigned int>> dst_i;

	bool valid_epipolar(float minmum_ang)
	{
		return epipolar_ang[sort_ids[sort_id_add]] > minmum_ang;
	}

	bool valid_aligned(int minmum_align)
	{
		return consistent_size[sort_ids[sort_id_add]] >= minmum_align;
	}

	float max_score()
	{
		float score_;
		while (sort_id_add < l3d_size)
		{
			score_ = scores[sort_ids[sort_id_add]];
			if (score_ == 0)
				sort_id_add++;
			else
				return score_;
		}

		return 0;
	}

	unsigned int max_id()
	{
		return sort_ids[sort_id_add];
	}

	unsigned int max_ori_id()
	{
		return ori_ID[sort_ids[sort_id_add]];
	}

	void next_score()
	{
		scores[sort_ids[sort_id_add]] = 0;
		sort_id_add++;
	}

	unsigned int ori_id(unsigned int l3d_id)
	{
		return ori_ID[sort_ids[sort_id_add]];
	}

};

class ScoreRecorder
{
	unsigned int sum_l3d_size = 0;

	std::vector<TASK> task_vec;

	void release_conflicts(unsigned int stereo_id, unsigned int l3d_id);

	void max_across_view();

	float max_score = 0;

	float current_score = 0;

	unsigned int max_stereo = 0;

	float minmum_epipolar;

	int minmum_align;

public:

	std::vector<unsigned int> ks_store;

	std::vector<unsigned int> match_ids_store;

	ScoreRecorder(int task_size, float minmum_epipolar_, int minmum_align_);

	void initial_tasks(int task_id, int l3d_size);

	void add_ori_id(int task_id, int l3d_id, int ori_id, float ep_ang);

	void add_connections(int stereo_i, int stereo_j, std::vector<pair_m>* l3d_pair_vec);

	void vote_lines();
};