/*
 *	Data: 2015.09.25
 *	Author: Joohyuk Yum
 *	Revision:
 *
 */

#include "evaluation_kit.h"

/* PUBLIC members */
evaluation_kit::evaluation_kit(char *fname, int mode, int io) : MODE(mode), IO(io) {
	if (mode == DOT_SIFT || mode == DOT_ASIFT) {
		if (io == READ)
			fp = std::fopen(fname, "r");
		else
			fp = std::fopen(fname, "w");
	}
	else { //DOT_DESC
		if (io == READ)
			fp = std::fopen(fname, "rb");
		else
			fp = std::fopen(fname, "wb");
	}

	if (fp == NULL)
		std::cout << "[ERROR] The file '" << fname << "' doesn't exist" << std::endl;
}

evaluation_kit::~evaluation_kit() {
	if (fp != NULL)
		std::fclose(fp);
}

int evaluation_kit::initial(char *fname, int mode, int io) {
	if (fp != NULL)
		std::fclose(fp);

	if (mode == DOT_SIFT || mode == DOT_ASIFT) {
		if (io == READ)
			fp = std::fopen(fname, "rt");
		else
			fp = std::fopen(fname, "wt");
	}
	else { //DOT_DESC
		if (io == READ)
			fp = std::fopen(fname, "rb");
		else
			fp = std::fopen(fname, "wb");
	}

	if (fp == NULL) {
		std::cout << "[ERROR] The file '" << fname << "' doesn't exist" << std::endl;
		return -1;
	}
	else
		return 1;
}

int evaluation_kit::operator() (std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors) {
	if (fp == NULL) {
		std::cout << "[ERROR] The file doesn't exist" << std::endl;
		return 0;
	}
	if (MODE == DOT_ASIFT) {
		std::cout << "[ERROR] evaluation_kit::operator() (std::vector<cv::KeyPoint> &, cv::Mat &) doesn't support a instance initialized for 'asift'" << std::endl;
		return 0;
	}

	int ret(0);
	if (IO == WRITE)
		ret = printData(keypoints, descriptors);
	else if(IO == READ)
		ret = getData(keypoints, descriptors);
	return ret;
}

int evaluation_kit::operator() (std::vector<cv::KeyPoint> &keypoints, std::vector<KeyPoint_param> &paramKeypoints, cv::Mat &descriptors) {
	if (fp == NULL) {
		std::cout << "[ERROR] The file doesn't exist" << std::endl;
		return 0;
	}
	int ret(0);
	if (IO == WRITE)
		ret = printData(keypoints, paramKeypoints, descriptors);
	else if (IO == READ)
		std::cout << "[ERROR] This operation is not supported yet" << std::endl;

	return ret;
}

void evaluation_kit::operator() (const cv::Mat &img, const std::vector<cv::KeyPoint> &keypoints, const std::vector<KeyPoint_param> &paramKeypoints, const cv::Mat &descriptors) {
	const cv::Scalar BLACK(0, 0, 0), WHITE(255, 255, 255);
	const float PI_PER_180 = CV_PI / 180.0;
	//cv::Mat display(img.clone());
	int step = (int)keypoints.size() / 200;
	for (int it_kps(0); it_kps < (int)keypoints.size(); it_kps = it_kps + step) {
		cv::Mat display(img.clone());
		cv::circle(display, keypoints[it_kps].pt, 2, WHITE, 2);
		float rad = keypoints[it_kps].size / 2.0;
		float tilt = paramKeypoints[it_kps].tilt;
		float angle = paramKeypoints[it_kps].angle * PI_PER_180;
		float ell_a, ell_b, ell_c;
		getellipse(rad, angle, tilt, ell_a, ell_b, ell_c);
		printf("rad(%f), tilt(%f), angle(%f)\n", rad, tilt, paramKeypoints[it_kps].angle);
		//display ellipse
		//bool condition = paramKeypoints[it_kps].angle == 108 && paramKeypoints[it_kps].tilt == 2;
		bool condition = 1;
		if (condition) {
			for (int it_deg(-180); it_deg < 180; it_deg++) {
				float radian = it_deg * PI_PER_180;
				float px = rad * cos(radian);
				float py = rad * sin(radian) * paramKeypoints[it_kps].tilt;
				float qx = cos(angle)*px - sin(angle)*py;
				float qy = sin(angle)*px + cos(angle)*py;
				cv::Point pt(qx + keypoints[it_kps].pt.x, qy + keypoints[it_kps].pt.y);
				cv::circle(display, pt, 1.0, BLACK, 1.0);
			}
			for (int ell_x(-rad); ell_x < rad; ell_x++) {
				float a = ell_a * ell_x * ell_x - 1;
				float b = ell_b * ell_x;
				float c = ell_c;
				float ell_y1, ell_y2;
				if (solveQuadEquation(c, b, a, ell_y1, ell_y2)) {
					cv::Point pt1(ell_x + keypoints[it_kps].pt.x, ell_y1 + keypoints[it_kps].pt.y);
					cv::Point pt2(ell_x + keypoints[it_kps].pt.x, ell_y2 + keypoints[it_kps].pt.y);
					cv::circle(display, pt1, 1.0, WHITE, 1.0);
					cv::circle(display, pt2, 1.0, WHITE, 1.0);
				}
			}
			for (int ell_y(-rad); ell_y < rad; ell_y++) {
				float a = ell_a;
				float b = ell_b * ell_y;
				float c = ell_c * ell_y * ell_y - 1;
				float ell_x1, ell_x2;
				if (solveQuadEquation(a, b, c, ell_x1, ell_x2)) {
					cv::Point pt1(ell_x1 + keypoints[it_kps].pt.x, ell_y + keypoints[it_kps].pt.y);
					cv::Point pt2(ell_x2 + keypoints[it_kps].pt.x, ell_y + keypoints[it_kps].pt.y);
					cv::circle(display, pt1, 1.0, WHITE, 1.0);
					cv::circle(display, pt2, 1.0, WHITE, 1.0);
				}
			}
		}
		cv::imshow("evaluation_kit::ellipse", display);
		cv::waitKey();
	}
	//cv::imshow("evaluation_kit::ellipse", display);
	//cv::waitKey();
}

//void evaluation_kit::operator() (const cv::Mat &img, const std::vector<cv::KeyPoint> &keypoints, const std::vector<KeyPoint_param> &paramKeypoints, const cv::Mat &descriptors) {
//	const cv::Scalar BLACK(0, 0, 0), WHITE(255, 255, 255);
//	const float PI_PER_180 = CV_PI / 180.0;
//	cv::Mat display(img.clone());
//	int step = (int)keypoints.size() / 100;
//	for (int it_kps(0); it_kps < (int)keypoints.size(); it_kps = it_kps + step) {
//		cv::circle(display, keypoints[it_kps].pt, 2, WHITE, 2);
//		float rad = keypoints[it_kps].size / 2.0;
//		float tilt = paramKeypoints[it_kps].tilt;
//		float angle = paramKeypoints[it_kps].angle * PI_PER_180;
//		float ell_a, ell_b, ell_c;
//		getellipse(rad, angle, tilt, ell_a, ell_b, ell_c);
//		//display ellipse
//		//bool condition = paramKeypoints[it_kps].angle == 108 && paramKeypoints[it_kps].tilt == 2;
//		bool condition = 1;
//		if (condition) {
//			for (int it_deg(-180); it_deg < 180; it_deg++) {
//				float radian = it_deg * PI_PER_180;
//				float px = rad * cos(radian);
//				float py = rad * sin(radian) * paramKeypoints[it_kps].tilt;
//				float qx = cos(angle)*px - sin(angle)*py;
//				float qy = sin(angle)*px + cos(angle)*py;
//				cv::Point pt(qx + keypoints[it_kps].pt.x, qy + keypoints[it_kps].pt.y);
//				cv::circle(display, pt, 1.0, BLACK, 1.0);
//			}
//			//for (int ell_x(-rad); ell_x < rad; ell_x++) {
//			//	float a = ell_a * ell_x * ell_x - 1;
//			//	float b = ell_b * ell_x;
//			//	float c = ell_c;
//			//	float ell_y1, ell_y2;
//			//	if (solveQuadEquation(c, b, a, ell_y1, ell_y2)) {
//			//		cv::Point pt1(ell_x + keypoints[it_kps].pt.x, ell_y1 + keypoints[it_kps].pt.y);
//			//		cv::Point pt2(ell_x + keypoints[it_kps].pt.x, ell_y2 + keypoints[it_kps].pt.y);
//			//		cv::circle(display, pt1, 1.0, BLACK, 1.0);
//			//		cv::circle(display, pt2, 1.0, BLACK, 1.0);
//			//	}
//			//}
//			//for (int ell_y(-rad); ell_y < rad; ell_y++) {
//			//	float a = ell_a;
//			//	float b = ell_b * ell_y;
//			//	float c = ell_c * ell_y * ell_y - 1;
//			//	float ell_x1, ell_x2;
//			//	if (solveQuadEquation(a, b, c, ell_x1, ell_x2)) {
//			//		cv::Point pt1(ell_x1 + keypoints[it_kps].pt.x, ell_y + keypoints[it_kps].pt.y);
//			//		cv::Point pt2(ell_x2 + keypoints[it_kps].pt.x, ell_y + keypoints[it_kps].pt.y);
//			//		cv::circle(display, pt1, 1.0, BLACK, 1.0);
//			//		cv::circle(display, pt2, 1.0, BLACK, 1.0);
//			//	}
//			//}
//		}
//	}
//	cv::imshow("evaluation_kit::ellipse", display);
//	cv::waitKey();
//}

/* PROTECTED members */
int evaluation_kit::printData(const std::vector<cv::KeyPoint> &keypoints, const cv::Mat &descriptors) {
	int keypoint_size = (int)keypoints.size();
	if (!keypoint_size)
		return 0;

	if (MODE == DOT_SIFT) {

		fprintf(fp, "128\n");
		fprintf(fp, "%d\n", keypoint_size);

		for (int l = 0; l < (int)keypoints.size(); l++) {
			fprintf(fp, "%.2f ", keypoints[l].pt.x);
			fprintf(fp, "%.2f ", keypoints[l].pt.y);

			double radius = (keypoints[l].size / 2.0) * (keypoints[l].size / 2.0);
			double one_over_radius = 1.0 / radius;
			fprintf(fp, "%.7f 0 %.7f ", one_over_radius, one_over_radius);

			for (int k = 0; k<128; k++)
				fprintf(fp, "%d ", (int)descriptors.at<float>(l, k));
			fprintf(fp, "\n");
		}
	}
	else if (MODE == DOT_DESC) {

		int tmpInt = keypoint_size;
		fwrite(&tmpInt, 4, 1, fp);

		for (int l = 0; l < (int)keypoints.size(); l++) {
			int keypoint_x = keypoints[l].pt.x;
			tmpInt = (double)keypoint_x / pow(2.0, keypoints[l].octave);
			fwrite(&tmpInt, 2, 1, fp);

			int keypoint_y = keypoints[l].pt.y;
			tmpInt = (double)keypoint_y / pow(2.0, keypoints[l].octave);
			fwrite(&tmpInt, 2, 1, fp);

			int keypoint_scale = keypoints[l].class_id - 1;
			tmpInt = keypoint_scale;
			fwrite(&tmpInt, 2, 1, fp);

			int keypoint_octave = keypoints[l].octave;
			tmpInt = keypoint_octave;
			fwrite(&tmpInt, 2, 1, fp);

			for (int k = 0; k < 128; k++) {
				int keypoint_descriptor = descriptors.at<float>(l, k);
				tmpInt = keypoint_descriptor;
				fwrite(&tmpInt, 1, 1, fp);
			}
		}
	}
	else
		return 0;
	return 1;
}

int evaluation_kit::printData(const std::vector<cv::KeyPoint> &keypoints, const std::vector<KeyPoint_param> &paramKeypoints, const cv::Mat &descriptors) {
	int keypoint_size = (int)keypoints.size();
	if (!keypoint_size)
		return 0;

	if (MODE == DOT_ASIFT) {
		/* FILE format
		*	[The size of a descriptor vector]
		*	[The number of keypoints]
		*	[kp.x] [kp.y] [ell_a] [ell_b] [ell_c] [tilt] [angle] [descriptor vector] ... \n
		*	[kp.x] [kp.y] [ell_a] [ell_b] [ell_c] [tilt] [angle] [descriptor vector] ... \n
		*	...
		*/
		const float PI_PER_180 = CV_PI / 180.0;
		fprintf(fp, "128\n");
		fprintf(fp, "%d\n", keypoint_size);

		for (int l = 0; l < (int)keypoints.size(); l++) {
			fprintf(fp, "%.2f ", keypoints[l].pt.x);
			fprintf(fp, "%.2f ", keypoints[l].pt.y);

			float rad = keypoints[l].size / 2.0;
			float tilt = paramKeypoints[l].tilt;
			float angle = paramKeypoints[l].angle * PI_PER_180;
			float ell_a, ell_b, ell_c;
			getellipse(rad, angle, tilt, ell_a, ell_b, ell_c);

			fprintf(fp, "%.7f %.7f %.7f ", ell_a, ell_b, ell_c);

			fprintf(fp, "%.7f %.7f ", paramKeypoints[l].tilt, paramKeypoints[l].angle);

			for (int k = 0; k < 128; k++) {
				fprintf(fp, "%d ", (int)descriptors.at<float>(l, k));
			}
			fprintf(fp, "\n");
		}
	}
	else if (MODE == DOT_SIFT) {
		const float PI_PER_180 = CV_PI / 180.0;
		fprintf(fp, "128\n");
		fprintf(fp, "%d\n", keypoint_size);

		for (int l = 0; l < (int)keypoints.size(); l++) {
			fprintf(fp, "%.2f ", keypoints[l].pt.x);
			fprintf(fp, "%.2f ", keypoints[l].pt.y);

			float rad = keypoints[l].size / 2.0;
			float tilt = paramKeypoints[l].tilt;
			float angle = paramKeypoints[l].angle * PI_PER_180;
			float ell_a, ell_b, ell_c;
			getellipse(rad, angle, tilt, ell_a, ell_b, ell_c);

			fprintf(fp, "%.7f %.7f %.7f ", ell_a, ell_b, ell_c);

			for (int k = 0; k < 128; k++) {
				fprintf(fp, "%d ", (int)descriptors.at<float>(l, k));
			}
			fprintf(fp, "\n");
		}
	}
	else if (MODE == DOT_DESC) {
		int tmpInt = keypoint_size;
		fwrite(&tmpInt, 4, 1, fp);

		for (int l = 0; l < (int)keypoints.size(); l++) {
			int keypoint_x = keypoints[l].pt.x;
			int tmpInt1 = (double)keypoint_x / pow(2.0, keypoints[l].octave);
			fwrite(&tmpInt1, 2, 1, fp);

			int keypoint_y = keypoints[l].pt.y;
			int tmpInt2 = (double)keypoint_y / pow(2.0, keypoints[l].octave);
			fwrite(&tmpInt2, 2, 1, fp);

			int keypoint_scale = keypoints[l].class_id - 1;
			int tmpInt3 = keypoint_scale;
			fwrite(&tmpInt3, 2, 1, fp);

			int keypoint_octave = keypoints[l].octave;
			int tmpInt4 = keypoint_octave;
			fwrite(&tmpInt4, 2, 1, fp);

			for (int k = 0; k<128; k++) {
				int keypoint_descriptor = descriptors.at<float>(l, k);
				int tmpInt5 = keypoint_descriptor;
				fwrite(&tmpInt5, 1, 1, fp);
			}
		}
	}
	else
		return 0;

	return 1;
}

int evaluation_kit::getData(std::vector<cv::KeyPoint>& keypoints, cv::Mat &descriptors) {
	if (MODE == DOT_SIFT) {
		std::cout << "[ERROR] 'getData from a DOT_SIFT file is not supported yet" << std::endl;
		return 0;
	}
	else if (MODE == DOT_DESC) {
		int nKps;
		fread(&nKps, sizeof(int), 1, fp);

		keypoints.resize(nKps);
		descriptors.create((int)keypoints.size(), 128, cv::DataType<float>::type);

		for (int iterKps = 0; iterKps < nKps; iterKps++) {
			int tmpInt;
			fread(&tmpInt, sizeof(int), 1, fp);
			keypoints[iterKps].pt.y = (int)(tmpInt & 0xffff0000) >> 16;
			keypoints[iterKps].pt.x = (int)(tmpInt & 0x0000ffff);

			fread(&tmpInt, sizeof(int), 1, fp);
			keypoints[iterKps].octave = (tmpInt & 0xffff0000) >> 16;
			keypoints[iterKps].size = tmpInt & 0x0000ffff;          // scale index
			keypoints[iterKps].angle = -1;

			keypoints[iterKps].pt.y = (float)keypoints[iterKps].pt.y * pow(2.0, keypoints[iterKps].octave);
			keypoints[iterKps].pt.x = (float)keypoints[iterKps].pt.x * pow(2.0, keypoints[iterKps].octave);

			for (int iterDesc = 0; iterDesc < 128; iterDesc = iterDesc + 4) {
				fread(&tmpInt, sizeof(int), 1, fp);
				descriptors.at<float>(iterKps, iterDesc + 3) = (float)((tmpInt & 0xff000000) >> 24);
				descriptors.at<float>(iterKps, iterDesc + 2) = (float)((tmpInt & 0x00ff0000) >> 16);
				descriptors.at<float>(iterKps, iterDesc + 1) = (float)((tmpInt & 0x0000ff00) >> 8);
				descriptors.at<float>(iterKps, iterDesc + 0) = (float)((tmpInt & 0x000000ff));
			}
		}
	}
	return 1;
}

//(ell_a) * x^2 + (2*ell_b) * x*y + (ell_c) * y^2 = 1
void evaluation_kit::getellipse(float radius, float angle, float tilt, float &ell_a, float &ell_b, float &ell_c) {
	float cx = cos(angle);
	float sx = sin(angle);
	float rx = radius;
	float ry = radius * tilt;
	ell_a = (cx*cx / (rx*rx)) + (sx*sx / (ry*ry));
	ell_b = 2 * cx*sx*(1 / (rx*rx) - 1 / (ry*ry));
	ell_c = (sx*sx / (rx*rx)) + (cx*cx / (ry*ry));
}

int evaluation_kit::solveQuadEquation(float a, float b, float c, float &sol_a, float &sol_b) {
	float d = b*b - 4 * a*c;
	if (d < 0.0)
		return 0;
	else {
		sol_a = (-b + std::sqrt(d)) / (2 * a);
		sol_b = (-b - std::sqrt(d)) / (2 * a);
		return 1;
	}
}
