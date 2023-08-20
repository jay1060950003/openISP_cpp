#pragma once
#ifndef _MODULE_H_
#define _MODULE_H_
#include <iostream>
#include <cmath>
#include <iomanip>

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core/utils/logger.hpp>

#include "common.h"

const float sdr_max_value = 255.0;

void dpc(cv::Mat&, uchar);
void blc(cv::Mat&, std::string, float, float);
void aaf(cv::Mat&, cv::Mat&);
void wbgc(cv::Mat&, std::string, std::vector<float>, ushort);
void bnr(cv::Mat&, std::string, uchar);
void cfa(cv::Mat&, std::string);
void ccm(cv::Mat&, cv::Mat&);
void gc(cv::Mat&, float, ushort);
void csc(cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&, float[3][4]);
void nlm(cv::Mat&, uchar, uchar, uchar, uchar);
void bnf(cv::Mat&, uchar, uchar, ushort clip);
void ee(cv::Mat&, cv::Mat&, char[3][5], ushort);
void bcc(cv::Mat&, uchar, uchar, ushort);
void yuv2rgb(cv::Mat&, cv::Mat&, cv::Mat&, cv::Mat&);
#endif