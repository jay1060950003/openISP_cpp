#pragma once
#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <vector>

#include <opencv2/core.hpp>

int median(int, int, int);
cv::Mat padding(cv::Mat&, uchar*);
std::vector<cv::Mat> split_bayer(cv::Mat&, std::string&);
cv::Mat reconstruct_bayer(std::vector<cv::Mat>&, std::string&);
float cnf_fade(float&, char);
float sum_cvMat(cv::Mat);

#endif
