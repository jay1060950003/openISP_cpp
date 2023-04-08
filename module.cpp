#include "module.h"

void dpc(cv::Mat& src, uchar thres)
{
	clock_t begin = clock();
	uchar pads[4] = { 2, 2, 2, 2 };
	cv::Mat src_p = padding(src, pads);
	ushort* p0, * p1, * p2, * p3, * p4, * p5, * p6, * p7, * p8, * p;
	int grad_h1 = 0, grad_h2 = 0, grad_h3 = 0, grad_v1 = 0, grad_v2 = 0, grad_v3 = 0;
	int grad_45_1 = 0, grad_45_2 = 0, grad_45_3 = 0, grad_135_1 = 0, grad_135_2 = 0, grad_135_3 = 0;
	int grad_h = 0, grad_v = 0, grad_45 = 0, grad_135 = 0;
	std::vector<int> gradient(4, 0);
	int grad_sum = 0;
	for (int i = 0; i < src_p.rows - 4; ++i) {
		std::cout << "\r" << "DPC: ";
		std::cout << std::setw(6) << std::fixed << std::setprecision(2) << (float)i / (src_p.rows - 2) * 100 << "%";
		p0 = src_p.ptr<ushort>(i + 2);
		p1 = src_p.ptr<ushort>(i);
		p2 = src_p.ptr<ushort>(i);
		p3 = src_p.ptr<ushort>(i);
		p4 = src_p.ptr<ushort>(i + 2);
		p5 = src_p.ptr<ushort>(i + 2);
		p6 = src_p.ptr<ushort>(i + 4);
		p7 = src_p.ptr<ushort>(i + 4);
		p8 = src_p.ptr<ushort>(i + 4);
		p = src.ptr<ushort>(i);
		for (int j = 0; j < src_p.cols - 4; j++) {
			grad_h1 = abs(p1[j] + p3[j + 4] - 2 * p2[j + 2]);
			grad_h2 = abs(p4[j] + p5[j + 4] - 2 * p0[j + 2]);
			grad_h3 = abs(p6[j] + p8[j + 4] - 2 * p7[j + 2]);
			grad_v1 = abs(p1[j] + p4[j] - 2 * p6[j]);
			grad_v2 = abs(p2[j + 2] + p7[j + 2] - 2 * p0[j + 2]);
			grad_v3 = abs(p3[j + 4] + p8[j + 4] - 2 * p5[j + 4]);
			grad_45_1 = 2 * abs(p2[j + 2] - p4[j]);
			grad_45_2 = abs(p3[j + 4] + p6[j] - 2 * p0[j + 2]);
			grad_45_3 = 2 * abs(p5[j + 4] - p7[j + 2]);
			grad_135_1 = 2 * abs(p2[j + 2] -p5[j + 4]);
			grad_135_2 = abs(p1[j] + p8[j + 4] - 2 * p0[j + 2]);
			grad_135_3 = 2 * abs(p4[j] - p7[j + 2]);

			grad_h = median(grad_h1, grad_h2, grad_h3);
			grad_v = median(grad_v1, grad_v2, grad_v3);
			grad_45 = median(grad_45_1, grad_45_2, grad_45_3);
			grad_135 = median(grad_135_1, grad_135_2, grad_135_3);
			gradient = { grad_h, grad_v, grad_45, grad_135 };
			auto minPosition = std::min_element(gradient.begin(), gradient.end());
			if (minPosition == gradient.begin() && grad_h2 > 4*(grad_h1+grad_h3))
			{
				if (abs(p4[j] - p0[j + 2]) < abs(p0[j + 2] - p5[j + 4])) {
					p[j] = p4[j] + (p2[j + 2] + p7[j + 2] - p1[j] - p6[j]) / 2;
				}
				else {
					p[j] = p5[j + 4] + (p2[j + 2] + p7[j + 2] - p3[j + 4] - p8[j + 4]) / 2;
				}
			}
			if (minPosition == gradient.begin()+1 && grad_v2 > 4 * (grad_v1 + grad_v3))
			{
				if (abs(p2[j + 2] - p0[j + 2]) < abs(p0[j + 2] - p7[j + 2])) {
					p[j] = p2[j + 2] + (p4[j] + p5[j + 4] - p1[j] - p3[j + 4]) / 2;
				}
				else {
					p[j] = p7[j + 2] + (p4[j] + p5[j + 4] - p6[j] - p8[j + 4]) / 2;
				}
			}
			if (minPosition == gradient.begin()+2)
			{
				grad_sum = abs(grad_135_1 - grad_135_2) + abs(grad_135_1 - grad_135_3) + abs(grad_135_2 - grad_135_3);
				if (grad_sum > 100) {
					if (grad_45_2 > 3 * (grad_45_1 + grad_45_3) && grad_135_2 > 3 * (grad_135_1 + grad_135_3)) {
						if (abs(p3[j + 4] - p0[j + 2]) < abs(p0[j + 2] - p6[j])) {
							p[j] = p3[j + 4] + (p4[j] + p7[j + 2] - p2[j + 2] - p5[j + 4]) / 2;
						}
						else {
							p[j] = p6[j] - (p4[j] + p7[j + 2] - p2[j + 2] - p5[j + 4]) / 2;
						}
					}
				}
				else {
					if (grad_45_2 > 3 * (grad_45_1 + grad_45_3)) {
						if (abs(p3[j + 4] - p0[j + 2]) < abs(p0[j + 2] - p6[j])) {
							p[j] = p3[j + 4] + (p4[j] + p7[j + 2] - p2[j + 2] - p5[j + 4]) / 2;
						}
						else {
							p[j] = p6[j] - (p4[j] + p7[j + 2] - p2[j + 2] - p5[j + 4]) / 2;
						}
					}
				}
			}
			if (minPosition == gradient.begin()+3)
			{
				grad_sum = abs(grad_45_1 - grad_45_2) + abs(grad_45_1 - grad_45_3) + abs(grad_45_2 - grad_45_3);
				if (grad_sum > 100) {
					if (grad_135_2 > 3 * (grad_135_1 + grad_135_3) && grad_45_2 > 3 * (grad_45_1 + grad_45_3)) {
						if (abs(p1[j] - p0[j + 2]) < abs(p0[j + 2] - p8[j + 4])) {
							p[j] = p1[j] + (p5[j + 4] + p6[j] - p2[j + 2] - p4[j]) / 2;
						}
						else {
							p[j] = p8[j + 4] - (p5[j + 4] + p6[j] - p2[j + 2] - p4[j]) / 2;
						}
					}
				}
				else {
					if (grad_135_2 > 3 * (grad_135_1 + grad_135_3)) {
						if (abs(p1[j] - p0[j + 2]) < abs(p0[j + 2] - p8[j + 4])) {
							p[j] = p1[j] + (p5[j + 4] + p6[j] - p2[j + 2] - p4[j]) / 2;
						}
						else {
							p[j] = p8[j + 4] - (p5[j + 4] + p6[j] - p2[j + 2] - p4[j]) / 2;
						}
					}
				}
			}

		}
	}
	clock_t end = clock();
	std::cout << "\r" << "DPC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void blc(cv::Mat& src, std::string bayer_pattern, float black_level, float white_level)
{
	clock_t begin = clock();
	std::vector<cv::Mat> bayer4 = split_bayer(src, bayer_pattern);
	int pixel_R = 0, pixel_Gr = 0, pixel_Gb = 0, pixel_B = 0;
	ushort* p_R, * p_Gr, * p_Gb, * p_B;
	for (int i = 0; i < bayer4[0].rows; ++i) {
		p_R = bayer4[0].ptr<ushort>(i);
		p_Gr = bayer4[1].ptr<ushort>(i);
		p_Gb = bayer4[2].ptr<ushort>(i);
		p_B = bayer4[3].ptr<ushort>(i);
		for (int j = 0; j < bayer4[0].cols; ++j) {
			pixel_R = p_R[j] - black_level;
			pixel_B = p_B[j] - black_level;
			pixel_Gr = p_Gr[j] - black_level;
			pixel_Gb = p_Gb[j] - black_level;

			p_R[j] = ((pixel_R < 0) ? 0 : (pixel_R > white_level ? white_level : pixel_R));
			p_Gr[j] = ((pixel_Gr < 0) ? 0 : (pixel_Gr > white_level ? white_level : pixel_Gr));
			p_Gb[j] = ((pixel_Gb < 0) ? 0 : (pixel_Gb > white_level ? white_level : pixel_Gb));
			p_B[j] = ((pixel_B < 0) ? 0 : (pixel_B > white_level ? white_level : pixel_B));
		}
	}
	src = reconstruct_bayer(bayer4, bayer_pattern);
	clock_t end = clock();
	std::cout << "BLC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void aaf(cv::Mat& src, cv::Mat& kernel)
{
	clock_t begin = clock();
	cv::Mat origin_img = src.clone();
	cv::filter2D(origin_img, src, CV_16UC1, kernel);
	clock_t end = clock();
	std::cout << "AAF Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms."<< std::endl;
}

void wbgc(cv::Mat& src, std::string bayer_pattern, std::vector<float> gain, ushort clip)
{
	gain = { 1 / gain[0] * 1024, 1 / gain[1] * 1024, 1 / gain[2] * 1024 };
	clock_t begin = clock();
	std::vector<cv::Mat> rggb = split_bayer(src, bayer_pattern);
	rggb[0] *= (gain[0] / 1024);
	rggb[1] *= (gain[1] / 1024);
	rggb[2] *= (gain[1] / 1024);
	rggb[3] *= (gain[2] / 1024);
	src = reconstruct_bayer(rggb, bayer_pattern);
	ushort* p_img;
	for (int i = 0; i < src.rows; ++i) {
		p_img = src.ptr<ushort>(i);
		for (int j = 0; j < src.cols; ++j) {
			p_img[j] = (p_img[j] < 0) ? 0 : (p_img[j] > clip ? clip : p_img[j]);
		}
	}
	clock_t end = clock();
	std::cout << "AWB Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void bnr(cv::Mat& src, std::string bayer_pattern, uchar ksize) {
	clock_t begin = clock();
	std::vector<cv::Mat> rggb = split_bayer(src, bayer_pattern);
	medianBlur(rggb[0], rggb[0], ksize);
	medianBlur(rggb[1], rggb[1], ksize);
	medianBlur(rggb[2], rggb[2], ksize);
	medianBlur(rggb[3], rggb[3], ksize);
	src = reconstruct_bayer(rggb, bayer_pattern);
	clock_t end = clock();
	std::cout  << "BNF Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void cfa(cv::Mat& src, std::string bayer_pattern)
{
	clock_t begin = clock();
	if (bayer_pattern == "gbrg") {
		cv::cvtColor(src, src, cv::COLOR_BayerGB2RGB);
	}
	else if (bayer_pattern == "rggb") {
		cv::cvtColor(src, src, cv::COLOR_BayerRG2RGB);
	}
	else if (bayer_pattern == "bggr") {
		cv::cvtColor(src, src, cv::COLOR_BayerBG2RGB);
	}
	else if (bayer_pattern == "grbg") {
		cv::cvtColor(src, src, cv::COLOR_BayerGR2RGB);
	}
	else {
		throw "bayer pattern is not declared!\n";
	}
	clock_t end = clock();
	std::cout << "CFA Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void ccm(cv::Mat& src, cv::Mat& ccm)
{
	clock_t begin = clock();
	cv::Mat img0 = src.clone();
	cv::Vec3s* p, * p0;
	for (int i = 0; i < src.rows; ++i) {
		p0 = img0.ptr<cv::Vec3s>(i);
		p = src.ptr<cv::Vec3s>(i);
		for (int j = 0; j < src.cols; ++j) {
			for (int k = 0; k < 3; ++k) {
				p[j][k] = ccm.ptr<float>(k)[0] * p0[j][0] + ccm.ptr<float>(k)[1] * p0[j][1] + ccm.ptr<float>(k)[2] * p0[j][2];
			}
		}
	}
	clock_t end = clock();
	std::cout << "CCM Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void gc(cv::Mat& src, float gamma, ushort clip)
{
	clock_t begin = clock();
	std::vector<int> LUT{};
	for (int i = 0; i < clip + 1; ++i) {
		float lx = std::pow(((float)i / clip), (float)gamma) * (float)sdr_max_value;
		LUT.push_back((int)lx);
	}

	cv::Vec3s* p;
	for (int i = 0; i < src.rows; ++i) {
		p = src.ptr<cv::Vec3s>(i);
		for (int j = 0; j < src.cols; ++j) {
			for (int k = 0; k < 3; ++k) {
				p[j][k] = ((p[j][k] < 0) ? ushort(LUT[0]) : (p[j][k] > clip ? ushort(LUT[clip]) : ushort(LUT[p[j][k]])));
			}
		}
	}
	clock_t end = clock();
	std::cout << "GaC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void csc(cv::Mat& src, cv::Mat& y, cv::Mat& u, cv::Mat& v, float YUV420[3][4])
{
	clock_t begin = clock();
	cv::Vec3s* p;
	float* p_YUV420;
	ushort* p_y, * p_u, * p_v;
	for (int i = 0; i < src.rows; ++i) {
		p = src.ptr<cv::Vec3s>(i);
		p_y = y.ptr<ushort>(i);
		p_u = u.ptr<ushort>(i);
		p_v = v.ptr<ushort>(i);
		for (int j = 0; j < src.cols; ++j) {
			p_y[j] = (p[j][0] * YUV420[0][0] + p[j][1] * YUV420[0][1] + p[j][2] * YUV420[0][2] + YUV420[0][3]);
			p_u[j] = (p[j][0] * YUV420[1][0] + p[j][1] * YUV420[1][1] + p[j][2] * YUV420[1][2] + YUV420[1][3]);
			p_v[j] = (p[j][0] * YUV420[2][0] + p[j][1] * YUV420[2][1] + p[j][2] * YUV420[2][2] + YUV420[2][3]);
		}
	}

	clock_t end = clock();
	std::cout << "CSC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void nlm(cv::Mat& src, uchar ds, uchar Ds, uchar h, uchar clip)
{
	clock_t begin = clock();
	ushort center_y = 0, center_x = 0, start_y = 0, start_x = 0, pixel_pad = 0;
	double sweight = .0, average = .0, wmax = .0, w = .0, dist = .0, pixel = 0;
	ushort* p_src, *p_p, *p_p0, * p_neighbor, * p_center_w;

	uchar pads[4] = { Ds, Ds, Ds, Ds };
	cv::Mat src_p = padding(src, pads);
	for (int y = 0; y < src_p.rows - 2 * Ds; ++y) {
		std::cout << "\r" << "NLM: ";
		std::cout << std::setw(6) << std::fixed << std::setprecision(2) << (float)y / (src_p.rows - 2) * 100 << "%";
		center_y = y + Ds;
		p_src = src.ptr<ushort>(y);
		p_p0 = src_p.ptr<ushort>(center_y);
		for (int x = 0; x < src_p.cols - 2 * Ds; ++x) {
			center_x = x + Ds;
			//calWeights
			sweight = .0, average = .0, wmax = .0, dist = .0;
			for (int m = 0; m < 2 * Ds + 1 - 2 * ds - 1; ++m) {
				start_y = center_y - Ds + ds + m;
				p_p = src_p.ptr<ushort>(start_y);
				for (int n = 0; n < 2 * Ds + 1 - 2 * ds - 1; ++n) {
					start_x = center_x - Ds + ds + n;
					pixel_pad = p_p[start_x];
					if (m != center_y || n != center_x) {
						for (int i = 0; i < 2 * ds + 1; ++i) {
							p_neighbor = src_p.ptr<ushort>(start_y - ds + i);
							p_center_w = src_p.ptr<ushort>(center_y - ds + i);
							for (int j = 0; j < 2 * ds + 1; ++j) {
								dist += (((p_neighbor[start_x - ds + j] - p_center_w[center_x - ds + j]) ^ 2) / (2 * ds + 1) ^ 2);
							}
						}
						w = exp(-dist / (h ^ 2));
						if (w > wmax) wmax = w;
						sweight += w;
						average += (w * pixel_pad);
					}
				}
			}
			average += (wmax * p_p0[center_x]);
			sweight += wmax;
			pixel = average / sweight;
			pixel = (pixel <= 0) ? 0 : (pixel >= clip ? clip : pixel);
			p_src[x] = (ushort)pixel;
		}
	}
	clock_t end = clock();
	std::cout << "\r"  << "NLM Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void bnf(cv::Mat& src, uchar sigmoid_s, uchar sigmoid_p, ushort clip)
{
	clock_t begin = clock();

	uchar pads[4] = { 2, 2, 2, 2 };
	cv::Mat src_p = padding(src, pads);
	int sum_weight = 0, sum_imgA = 0;
	ushort pixel_center = 0;
	double dive = .0, w = .0;
	ushort* p, * p_p;

	for (int y = 0; y < src_p.rows - 4; ++y) {
		std::cout << "\r" << "BNF: ";
		std::cout << std::setw(6) << std::fixed << std::setprecision(2) << (float)y / (src_p.rows - 4) * 100 << "%";
		p = src.ptr<ushort>(y);
		for (int x = 0; x < src_p.cols - 4; ++x) {
			pixel_center = src_p.ptr<ushort>(y + 2)[x + 2];
			sum_imgA = 0;
			sum_weight = 0;
			dive = .0;
			for (int i = 0; i < 5; ++i) {
				p_p = src_p.ptr<ushort>(y + i);
				for (int j = 0; j < 5; ++j) {
					w = exp(-((abs(p_p[x + j] - pixel_center)) ^ 2) / 2 / sigmoid_p ^ 2) * exp(-((2 - i) ^ 2 + (2 - j) ^ 2) / 2 / sigmoid_s);
					sum_imgA += p_p[x + j] * w;
					sum_weight += w;
				}
			}
			dive = sum_imgA / sum_weight;
			dive = (dive < 0) ? 0 : (dive > clip ? clip : dive);
			p[x] = (ushort)dive;
		}
	}
	clock_t end = clock();
	std::cout << "\r" << "BNF Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void ee(cv::Mat& src, cv::Mat& edgemap, char edge_filter[3][5], ushort clip)
{
	 clock_t begin = clock();

	uchar pads[4] = { 1,1,2,2 };
	cv::Mat src_p = padding(src, pads);
	double em_img, ee_img, tmp_em_img;
	ushort* p_src, * p_edgemap;
	for (int y = 0; y < src_p.rows - 2; ++y) {
		std::cout << "\r" << "EEH: ";
		std::cout << std::setw(6) << std::fixed << std::setprecision(2) << (float)y / (src_p.rows - 4) * 100 << "%";
		p_edgemap = edgemap.ptr<ushort>(y);
		p_src = src.ptr<ushort>(y);
		for (int x = 0; x < src_p.cols - 4; ++x) {
			em_img = 0.0;
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 5; ++j) {
					em_img += src_p.ptr<ushort>(y + i)[x + j] * edge_filter[i][j];
				}
			}
			em_img = em_img / 8;
			ee_img = src_p.ptr<ushort>(y + 1)[x + 2] + em_img;
			ee_img = (ee_img < 0) ? 0 : (ee_img > clip ? clip : ee_img);
			p_edgemap[x] = (ushort)em_img;
			p_src[x] = (ushort)ee_img;
		}
	}
	clock_t end = clock();
	std::cout << "\r" << "EEH Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void bcc(cv::Mat& src, uchar brightness, uchar contrast, ushort bcc_clip)
{
	clock_t begin = clock();
	contrast = contrast / (2 ^ 5);
	ushort* p;
	int pixel;
	for (int y = 0; y < src.rows; ++y) {
		p = src.ptr<ushort>(y);
		for (int x = 0; x < src.cols; ++x) {
			pixel = p[x];
			pixel = pixel + brightness + (pixel - 127) * contrast;
			pixel = (pixel < 0) ? 0 : (pixel > bcc_clip ? bcc_clip : pixel);
			p[x] = (ushort)pixel;
		}
	}
	clock_t end = clock();
	std::cout << "BCC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void yuv2rgb(cv::Mat& src, cv::Mat& y, cv::Mat& u, cv::Mat& v)
{
	clock_t begin = clock();
	ushort* p_y, * p_u, * p_v;
	cv::Vec3s* p;
	for (int i = 0; i < src.rows; ++i){
		p = src.ptr<cv::Vec3s>(i);
		p_y = y.ptr<ushort>(i);
		p_u = u.ptr<ushort>(i);
		p_v = v.ptr<ushort>(i);
		for (int j = 0; j < src.cols; ++j){
			p[j][0] = cv::saturate_cast<ushort>(p_y[j] + 1.402 * (p_v[j] - 128));									// R
			p[j][1] = cv::saturate_cast<ushort>(p_y[j] - 0.34413 * (p_u[j] - 128) - 0.71414 * (p_v[j] - 128));		// G
			p[j][2] = cv::saturate_cast<ushort>(p_y[j] + 1.772 * (p_u[j] - 128));									// B
		}
	}
	clock_t end = clock();
	std::cout << "Y2R Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;

}
