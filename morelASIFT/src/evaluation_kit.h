/*
*	Data: 2015.09.25
*	Author: Joohyuk Yum
*	Revision:
*
*/

#ifndef __EVALUATION_KIT_H__
#define __EVALUATION_KIT_H__

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "define.h"

class evaluation_kit {
public:
	evaluation_kit(char *fname, int mode, int io);
	~evaluation_kit();

	/*
	*	DOT_SIFT: Dot SIFT file contains information for a evaluation about repeatability and matching score.
	*	DOT_DESC: Dot DESC file contains results of SIFT algorithm in the form of SIFT hardware output.
	 *	DOT_SIFT: Dot ASIFT file contains information for a evaluation about repeatability and matching score.
	*/
	enum { DOT_SIFT, DOT_DESC, DOT_ASIFT };
	/*
	*	READ: the file pointer member is set for read
	*	WRITE: the file pointer member is set for write
	*/
	enum { READ, WRITE };

	int initial(char *fname, int mode, int io);
	int operator() (std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors);
	int operator() (std::vector<cv::KeyPoint> &keypoints, std::vector<KeyPoint_param> &paramKeypoints, cv::Mat &descriptors);
	void operator() (const cv::Mat &img, const std::vector<cv::KeyPoint> &keypoints, const std::vector<KeyPoint_param> &paramKeypoints, const cv::Mat &descriptors);

protected:
	int printData(const std::vector<cv::KeyPoint> &keypoints, const cv::Mat &descriptors);
	int printData(const std::vector<cv::KeyPoint> &keypoints, const std::vector<KeyPoint_param> &paramKeypoints, const cv::Mat &descriptors);
	int getData(std::vector<cv::KeyPoint>& keypoints, cv::Mat &descriptors);
	void getellipse(float radius, float angle, float tilt, float &ell_a, float &ell_b, float &ell_c);
	int solveQuadEquation(float a, float b, float c, float &sol_a, float &sol_b);

private:
	FILE *fp;
	int MODE;
	int IO;
};


#endif
