
#include <iostream>
#include <fstream>
#include "match_util.h"

using namespace std;
using namespace cv;

void determine_good_match_by_brute_force(const Mat &descriptors_object, const Mat &descriptors_scene, vector<DMatch> &good_matches) {
    good_matches.clear();
		
	Ptr<DescriptorMatcher> descriptorMatcher;
	descriptorMatcher = new BFMatcher(NORM_L2);

	vector<vector<DMatch> > matches;
    descriptorMatcher->knnMatch(descriptors_object, descriptors_scene, matches, 2);
	
	double ratio = 0.8;
//	 double ratio = 0.60;
	
	for( int i = 0; i < descriptors_object.rows; i++ ) {
		if(matches[i][0].distance < ratio*matches[i][1].distance)
			good_matches.push_back( matches[i][0]);
		//good_matches.push_back(matches[i][0]);
    }
}

void determine_good_match_by_brute_force(const Mat &descriptors_object, const Mat &descriptors_scene, const vector<KeyPoint> &lkps, const vector<KeyPoint> &rkps, vector<DMatch> &good_matches) {		
	double ratio = 0.6;
    int search_thr = 20;

	good_matches.clear();
	Ptr<DescriptorMatcher> descriptorMatcher;
	descriptorMatcher = new BFMatcher(NORM_L2);

	vector<vector<DMatch> > matches;
    descriptorMatcher->knnMatch(descriptors_object, descriptors_scene, matches, search_thr);
	
	for( int i = 0; i < descriptors_object.rows; i++ ) {
		for(int k = 1; k < search_thr; k++) {
			if(rkps[matches[i][0].trainIdx].response == rkps[matches[i][k].trainIdx].response)
				if(matches[i][0].distance < ratio*matches[i][k].distance)
					good_matches.push_back( matches[i][0]);
		}
    }
}

void matched_keypoint_to_mat(const vector<DMatch> &good_matches, const vector<KeyPoint> &keypoints_object, const vector<KeyPoint> &keypoints_scene, Mat &obj_mat, Mat &scene_mat) {
    obj_mat.create((int)good_matches.size(), 1, CV_32FC2);
    scene_mat.create((int)good_matches.size(), 1, CV_32FC2);

    for(int i=0;i<good_matches.size();i++) {
        obj_mat.at<Vec2f>(i, 0) = keypoints_object[good_matches[i].queryIdx].pt;
        scene_mat.at<Vec2f>(i, 0) = keypoints_scene[good_matches[i].trainIdx].pt;
    }
}

void draw_transformed_box(const Mat &corners, const Point2f &offset, Mat &img, const Scalar &color, int thick) {
    for(int i = 0; i < 4; i++) {
        Point2f st = Point2f(corners.at<Vec2f>(i, 0)[0], corners.at<Vec2f>(i, 0)[1]) + offset;
        Point2f ed = Point2f(corners.at<Vec2f>((i+1)%4, 0)[0], corners.at<Vec2f>((i+1)%4, 0)[1]) + offset;

		 cv::line(img, st, ed, color, thick);
    }
}

void match_keypoint(const Mat &lhs_image, const vector<KeyPoint> &lhs_kps, const Mat &rhs_image, const vector<KeyPoint> &rhs_kps) {
	vector<DMatch> good_matches;

	clock_t start;
	clock_t end;
	
	start = clock();

	good_matches.resize(lhs_kps.size());
	
	for( int i = 0; i < lhs_kps.size(); i++ ) {
		good_matches[i].queryIdx = i;
		good_matches[i].trainIdx = i;
    }
	
	cout<<"N_of_matches : "<<good_matches.size()<<endl;

	Mat org_kp_mat, comp_kp_mat;
	Mat H;
	Mat matchImg = rhs_image.clone();

	if(good_matches.size() > 3) {
		matched_keypoint_to_mat(good_matches, lhs_kps, rhs_kps,	org_kp_mat, comp_kp_mat);
		H = findHomography( org_kp_mat, comp_kp_mat, CV_RANSAC );
		//-- Get the corners from the image_1 ( the object to be "detected" )
		std::vector<Point2f> obj_corners(4);
		obj_corners[0] = cvPoint(0,0);
		obj_corners[1] = cvPoint( lhs_image.cols, 0 );
		obj_corners[2] = cvPoint( lhs_image.cols, lhs_image.rows );
		obj_corners[3] = cvPoint( 0, lhs_image.rows-20 );

		Mat obj_corners_pt_mat(obj_corners);
		Mat scene_corners_pt_mat;
		cv::perspectiveTransform( obj_corners_pt_mat, scene_corners_pt_mat, H);
		
		//-- Draw lines between the corners (the mapped object in the scene - image_2 )
		//draw_transformed_box(scene_corners_pt_mat, Point2f( lhs_image.cols, 0), matchImg, Scalar(0, 255, 0), 2);
		draw_transformed_box(scene_corners_pt_mat, Point2f( 0, 0), matchImg, Scalar(0, 255, 0), 3);	
	}
	
	end = clock();
	printf("keypoint_match : %lf ms elapsed\n", double(end - start));

	//imshow("match", matchImg);
}

void find_transform(const cv::Mat &lhs_image, const std::vector<cv::KeyPoint> &lhs_kps, const cv::Mat &lhs_desc, const cv::Mat &rhs_image, const std::vector<cv::KeyPoint> &rhs_kps, const cv::Mat &rhs_desc) {
 	vector<DMatch> good_matches;
	Mat matchImg;
	Mat matchImg_t1;
	Mat matchImg_t2;

	clock_t start;
	clock_t end;

	start = clock();
	determine_good_match_by_brute_force(rhs_desc, lhs_desc, good_matches);
	end = clock();
	printf("determine_good_match_by_brute_force : %lf ms elapsed\n", double(end - start));

	//ofstream match_print;
	//match_print.open("match_result.csv");
	//
	//for(int i = 0; i < good_matches.size(); i++) {
	//	match_print << lhs_kps[good_matches[i].queryIdx].pt.x << "," << lhs_kps[good_matches[i].queryIdx].pt.y << "," << lhs_kps[good_matches[i].queryIdx].response << ",";
	//	match_print << rhs_kps[0][good_matches[i].trainIdx].pt.x << "," << rhs_kps[0][good_matches[i].trainIdx].pt.y << "," << rhs_kps[0][good_matches[i].trainIdx].response << endl;
	//}
	//match_print.close();
	
	start = clock();
	
//	int max_row = lhs_image.rows > rhs_image.rows ? lhs_image.rows : rhs_image.rows;
//	matchImg.create(max_row, lhs_image.cols + rhs_image.cols, lhs_image.type());
//	lhs_image.copyTo(matchImg(cvRect(0, 0, lhs_image.cols, lhs_image.rows)));
//	rhs_image.copyTo(matchImg(cvRect(lhs_image.cols, 0, rhs_image.cols, rhs_image.rows)));
//
//	for(int i = 0; i < (int)lhs_kps.size(); i++) {
//		circle(matchImg, Point(lhs_kps[i].pt.x, lhs_kps[i].pt.y), 1, Scalar(0,255,255), 1);
//	}
//	for(int i = 0; i < (int)rhs_kps.size(); i++) {
//		circle(matchImg, Point(rhs_kps[i].pt.x + lhs_image.cols, rhs_kps[i].pt.y), 1, Scalar(0,255,255), 1);
//	}
//
//	for(int i = 0; i < good_matches.size(); i++) {
//		cv::line(matchImg, Point(lhs_kps[good_matches[i].queryIdx].pt.x, lhs_kps[good_matches[i].queryIdx].pt.y)
//			, Point(rhs_kps[good_matches[i].trainIdx].pt.x + lhs_image.cols, rhs_kps[good_matches[i].trainIdx].pt.y), Scalar(0, 255, 0), 1);
//	}
	drawMatches(rhs_image, rhs_kps, lhs_image, lhs_kps, good_matches, matchImg);
	lhs_image.copyTo(matchImg_t1);
	matchImg.copyTo(matchImg_t2);

	end = clock();
	printf("drawMatches : %lf ms elapsed\n", double(end - start));
	
	cout<<"N_of_matches : "<<good_matches.size()<<endl;
	
	Mat org_kp_mat, comp_kp_mat;
	Mat H;
	
	if(good_matches.size() > 3) {
		start = clock();
		matched_keypoint_to_mat(good_matches, rhs_kps, lhs_kps, org_kp_mat, comp_kp_mat);
		//matched_keypoint_to_mat(good_matches, lhs_kps, rhs_kps[0], org_kp_mat, comp_kp_mat);
		end = clock();
		printf("matched_keypoint_to_mat : %lf ms elapsed\n", double(end - start));
		
		start = clock();
		H = findHomography( org_kp_mat, comp_kp_mat, CV_RANSAC );
		end = clock();
		printf("findHomography : %lf ms elapsed\n", double(end - start));
		
		//-- Get the corners from the image_1 ( the object to be "detected" )
		std::vector<Point2f> obj_corners(4);
		obj_corners[0] = cvPoint(0,0);
		obj_corners[1] = cvPoint( rhs_image.cols, 0 );
		obj_corners[2] = cvPoint( rhs_image.cols, rhs_image.rows );
		obj_corners[3] = cvPoint( 0, rhs_image.rows );

		Mat obj_corners_pt_mat(obj_corners);
		Mat scene_corners_pt_mat;
		start = clock();
		perspectiveTransform( obj_corners_pt_mat, scene_corners_pt_mat, H);
		end = clock();
		printf("perspectiveTransform : %lf ms elapsed\n", double(end - start));
	
		//-- Draw lines between the corners (the mapped object in the scene - image_2 )
		//draw_transformed_box(scene_corners_pt_mat, Point2f(0, 0), matchImg_t1, Scalar(0, 255, 0), 4);
		//draw_transformed_box(scene_corners_pt_mat, Point2f( rhs_image[0].cols, 0), matchImg, Scalar(0, 255, 0), 4);
		//cv::putText(matchImg, full_string, Point(10, 30), FONT_HERSHEY_PLAIN , 2, Scalar(0,255,0), 3);
	}

	imwrite("match.bmp", matchImg);
	cv::Mat scaledMatchImg;
	float scl(1.0);
	cv::resize(matchImg, scaledMatchImg, cv::Size(), scl, scl, 1);
	imshow("match", scaledMatchImg);
	waitKey(10);
}
