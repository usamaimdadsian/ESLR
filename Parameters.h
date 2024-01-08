#pragma once

#ifndef pixel_2_line_dis
#define pixel_2_line_dis 2
#endif 

#ifndef maximum_ang
#define maximum_ang 0.0872664625997 
#endif 

#ifndef intersect_dist
#define intersect_dist 10
#endif 

// for acceleration in finding the intersection
#ifndef intersect_cos
#define intersect_cos 0.8660
#endif 

#ifndef support_pt_num
#define support_pt_num 10
#endif 

#ifndef support_homo_num
#define support_homo_num 10
#endif 

#ifndef line_2_line_intersec
#define line_2_line_intersec 0.5
#endif 

#ifndef min_epipolar_ang
#define min_epipolar_ang 0.087266
#endif 

#ifndef depth_shift_pixel
#define depth_shift_pixel 10.0
#endif 

#ifndef min_pairs_num
#define min_pairs_num 2
#endif 

/*
  One can scale the image by change this threshold.
  it is set as a large number to avoid scale.
 */
#ifndef max_image_width
#define max_image_width 999999 
#endif 

/***********************************************************/
/*The following parameters are used for managing the memory*/
/***********************************************************/
#ifndef max_map_size
#define max_map_size 100
#endif 

#ifndef max_line3D_size
#define max_line3D_size 200
#endif 

#ifndef max_line2D_size
#define max_line2D_size 200
#endif 

