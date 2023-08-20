#include "common.h"

int median(int a, int b, int c) {
	return a >= b ? (b >= c ? b : (a >= c ? c : a)) : (a >= c ? a : (b >= c ? c : b));
}

cv::Mat padding(cv::Mat& img, uchar* pads)
{
	cv::Mat img_pad = cv::Mat(cv::Size(img.cols + pads[2] + pads[3], img.rows + pads[0] + pads[1]), img.type());
	cv::copyMakeBorder(img, img_pad, pads[0], pads[1], pads[2], pads[3], cv::BORDER_REFLECT);
	return img_pad;
}

std::vector<cv::Mat> split_bayer(cv::Mat& bayer_array, std::string& bayer_pattern) {
	cv::Mat R = cv::Mat(bayer_array.size() / 2, bayer_array.type());
	cv::Mat Gr = cv::Mat(bayer_array.size() / 2, bayer_array.type());
	cv::Mat Gb = cv::Mat(bayer_array.size() / 2, bayer_array.type());
	cv::Mat B = cv::Mat(bayer_array.size() / 2, bayer_array.type());
	ushort* p_R, * p_Gr, * p_Gb, * p_B, * p_bayer1, * p_bayer2;
	for (int i = 0; i < bayer_array.rows / 2; ++i) {
		p_R = R.ptr<ushort>(i);
		p_Gr = Gr.ptr<ushort>(i);
		p_Gb = Gb.ptr<ushort>(i);
		p_B = B.ptr<ushort>(i);
		p_bayer1 = bayer_array.ptr<ushort>(2 * i);
		p_bayer2 = bayer_array.ptr<ushort>(2 * i + 1);
		for (int j = 0; j < bayer_array.cols / 2; ++j) {
			if (bayer_pattern == "gbrg") {
				p_R[j] = p_bayer1[2 * j + 1];
				p_Gr[j] = p_bayer2[2 * j + 1];
				p_Gb[j] = p_bayer1[2 * j];
				p_B[j] = p_bayer2[2 * j];
			}
			else if (bayer_pattern == "rggb") {
				p_R[j] = p_bayer1[2 * j];
				p_Gr[j] = p_bayer2[2 * j];
				p_Gb[j] = p_bayer1[2 * j + 1];
				p_B[j] = p_bayer2[2 * j + 1];
			}
			else if (bayer_pattern == "bggr") {
				p_R[j] = p_bayer2[2 * j + 1];
				p_Gr[j] = p_bayer1[2 * j + 1];
				p_Gb[j] = p_bayer2[2 * j];
				p_B[j] = p_bayer1[2 * j];
			}
			else if (bayer_pattern == "grbg") {
				p_R[j] = p_bayer2[2 * j];
				p_Gr[j] = p_bayer1[2 * j];
				p_Gb[j] = p_bayer2[2 * j + 1];
				p_B[j] = p_bayer1[2 * j + 1];
			}
			else {
				throw "bayer pattern is not declared!\n";
			}
		}
	}
	return std::vector<cv::Mat> { R,Gr,Gb,B };
}

cv::Mat reconstruct_bayer(std::vector<cv::Mat>& sub_arrays, std::string& bayer_pattern) {
	cv::Mat bayer_array = cv::Mat(sub_arrays[0].size() * 2, sub_arrays[0].type());
	ushort* p_R, * p_Gr, * p_Gb, * p_B, * p_bayer1, * p_bayer2;
	for (int i = 0; i < sub_arrays[0].rows; ++i) {
		p_R = sub_arrays[0].ptr<ushort>(i);
		p_Gr = sub_arrays[1].ptr<ushort>(i);
		p_Gb = sub_arrays[2].ptr<ushort>(i);
		p_B = sub_arrays[3].ptr<ushort>(i);
		p_bayer1 = bayer_array.ptr<ushort>(2 * i);
		p_bayer2 = bayer_array.ptr<ushort>(2 * i + 1);
		for (int j = 0; j < sub_arrays[0].cols; ++j) {
			if (bayer_pattern == "gbrg") {
				p_bayer1[2 * j + 1] = p_R[j];
				p_bayer2[2 * j + 1] = p_Gr[j];
				p_bayer1[2 * j] = p_Gb[j];
				p_bayer2[2 * j] = p_B[j];
			}
			else if (bayer_pattern == "rggb") {
				p_bayer1[2 * j] = p_R[j];
				p_bayer2[2 * j] = p_Gr[j];
				p_bayer1[2 * j + 1] = p_Gb[j];
				p_bayer2[2 * j + 1] = p_B[j];
			}
			else if (bayer_pattern == "bggr") {
				p_bayer2[2 * j + 1] = p_R[j];
				p_bayer1[2 * j + 1] = p_Gr[j];
				p_bayer2[2 * j] = p_Gb[j];
				p_bayer1[2 * j] = p_B[j];
			}
			else if (bayer_pattern == "grbg") {
				p_bayer2[2 * j] = p_R[j];
				p_bayer1[2 * j] = p_Gr[j];
				p_bayer2[2 * j + 1] = p_Gb[j];
				p_bayer1[2 * j + 1] = p_B[j];
			}
			else {
				throw "bayer pattern is not declared!\n";
			}
		}
	}
	return bayer_array;
}
