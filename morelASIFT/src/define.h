#ifndef _DEFINE_H_
#define _DEFINE_H_

//#include <cv.h>

#define ABS(x) ( ( (x) < 0 )? -(x) : (x) )

//#define MAX(a, b) ( a < b ? b : a)

typedef struct _t_patch_{
	int dx;
	int dy;
	int mag;
	int ori;
	int ori_bin;
	int desc_bin_fraction;
	int local_x;
	int local_y;
} t_patch;

typedef struct _t_keypoint_param {
	float tilt;
	float angle;
} KeyPoint_param;

#endif
