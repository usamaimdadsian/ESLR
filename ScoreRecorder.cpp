#include "ScoreRecorder.h"

#include <iostream>
#include <vector>
#include <numeric>      
#include <algorithm>   

using namespace std;

template <typename T>
vector<size_t> sort_indexes(const vector<T>& v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

ScoreRecorder::ScoreRecorder(int task_size, float minmum_epipolar_, int minmum_align_)
{
	task_vec.resize(task_size);
	minmum_epipolar = minmum_epipolar_;
	minmum_align = minmum_align_;
}

void ScoreRecorder::initial_tasks(int task_id, int l3d_size)
{
	task_vec[task_id].l3d_size = l3d_size;
	task_vec[task_id].ori_ID.resize(l3d_size);
	task_vec[task_id].scores.resize(l3d_size);
	task_vec[task_id].epipolar_ang.resize(l3d_size);

	task_vec[task_id].consistent_size.resize(l3d_size);

	for (int i = 0; i < l3d_size; i++)
	{
		task_vec[task_id].scores[i] = 0;
		task_vec[task_id].consistent_size[i] = 0;
	}

	task_vec[task_id].dst_k.resize(l3d_size);
	task_vec[task_id].dst_i.resize(l3d_size);
}

void ScoreRecorder::add_ori_id(int task_id, int l3d_id, int ori_id, float ep_ang)
{

	task_vec[task_id].ori_ID[l3d_id] = ori_id;
	task_vec[task_id].ori_ID[l3d_id] = ori_id;
	task_vec[task_id].epipolar_ang[l3d_id] = ep_ang;

}

void ScoreRecorder::add_connections(int stereo_i, int stereo_j, std::vector<pair_m>* l3d_pair_vec)
{
	for (int i = 0; i < (*l3d_pair_vec).size(); i++)
	{
		task_vec[stereo_i].dst_k[(*l3d_pair_vec)[i].l3d_i].push_back(stereo_j);
		task_vec[stereo_i].dst_i[(*l3d_pair_vec)[i].l3d_i].push_back((*l3d_pair_vec)[i].l3d_j);

		task_vec[stereo_j].dst_k[(*l3d_pair_vec)[i].l3d_j].push_back(stereo_i);
		task_vec[stereo_j].dst_i[(*l3d_pair_vec)[i].l3d_j].push_back((*l3d_pair_vec)[i].l3d_i);

		if ((*l3d_pair_vec)[i].score > 0)
		{
			task_vec[stereo_i].scores[(*l3d_pair_vec)[i].l3d_i] +=
				(*l3d_pair_vec)[i].score *
				task_vec[stereo_i].epipolar_ang[(*l3d_pair_vec)[i].l3d_i];

			task_vec[stereo_j].scores[(*l3d_pair_vec)[i].l3d_j] +=
				(*l3d_pair_vec)[i].score *
				task_vec[stereo_j].epipolar_ang[(*l3d_pair_vec)[i].l3d_j];

			task_vec[stereo_j].consistent_size[(*l3d_pair_vec)[i].l3d_j]++;
			task_vec[stereo_i].consistent_size[(*l3d_pair_vec)[i].l3d_i]++;

		}

	}
}

void ScoreRecorder::max_across_view()
{
	max_score = 0;
	for (int i = 0; i < task_vec.size(); i++)
	{
		current_score = task_vec[i].max_score();

		if (current_score > max_score)
		{
			max_score = current_score;
			max_stereo = i;
		}
	}

}

void ScoreRecorder::release_conflicts
(unsigned int stereo_id, unsigned int l3d_id)
{
	if (task_vec[stereo_id].scores[l3d_id] == 0)
		return;// has been released

	task_vec[stereo_id].scores[l3d_id] = 0;

	for (int i = 0; i < task_vec[stereo_id].dst_i[l3d_id].size(); i++)
	{
		release_conflicts(
			task_vec[stereo_id].dst_k[l3d_id][i],
			task_vec[stereo_id].dst_i[l3d_id][i]
		);
	}
}

void ScoreRecorder::vote_lines()
{
	for (int i = 0; i < task_vec.size(); i++)
		task_vec[i].sort_ids = sort_indexes(task_vec[i].scores);

	while (1)
	{
		max_across_view();
		if (max_score == 0)
			break;
		//1 add to 
		if (task_vec[max_stereo].valid_epipolar(minmum_epipolar) &&
			task_vec[max_stereo].valid_aligned(minmum_align))
		{
			ks_store.push_back(max_stereo);
			match_ids_store.push_back(task_vec[max_stereo].max_ori_id());
		}

		//2 release other conflicts
		release_conflicts(max_stereo, task_vec[max_stereo].max_id());

	}
}