#include "module.h"

void dpc(cv::Mat& src, uchar thres)
{
	clock_t begin = clock();
	uchar pads[4] = { 2, 2, 2, 2 };
	cv::Mat src_p = padding(src, pads);
	ushort* p0, * p1, * p2, * p3, * p4, * p5, * p6, * p7, * p8, * p;
	for (int i = 0; i < src_p.rows - 4; ++i) {
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
			if ((p1[j] - p0[j + 2]) > thres && (p2[j + 2] - p0[j + 2]) > thres && (p3[j + 4] - p0[j + 2]) > thres
				&& (p4[j] - p0[j + 2]) > thres && (p5[j + 4] - p0[j + 2]) > thres && (p6[j] - p0[j + 2]) > thres
				&& (p7[j + 2] - p0[j + 2]) > thres && (p8[j + 4] - p0[j + 2]) > thres) {
				std::vector<int> gradient = { abs(2 * p0[j + 2] - p2[j + 2] - p7[j + 2]),abs(2 * p0[j + 2] - p4[j] - p5[j + 4]),
											 abs(2 * p0[j + 2] - p1[j] - p8[j + 4]), abs(2 * p0[j + 2] - p3[j + 4] - p6[j]) };
				auto minPosition = std::min_element(gradient.begin(), gradient.end());
				if (minPosition == gradient.begin()) {
					p[j] = (p2[j + 2] + p7[j + 2] + 1) / 2;
				}
				else if (minPosition == gradient.begin() + 1)
				{
					p[j] = (p4[j] + p5[j + 4] + 1) / 2;
				}
				else if (minPosition == gradient.begin() + 2)
				{
					p[j] = (p1[j] + p8[j + 4] + 1) / 2;
				}
				else {
					p[j] = (p3[j + 4] + p6[j] + 1) / 2;
				}
			}
			else {
				p[j] = p0[j + 2];
			}
		}
	}
	clock_t end = clock();
	std::cout << "DPC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
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

void cnf(cv::Mat& src, std::string bayer_pattern, std::vector<float> gain, uchar thres)
{
	gain = { 1 / gain[0] * 1024, 1 / gain[1] * 1024, 1 / gain[2] * 1024 };
	clock_t begin = clock();
	std::vector<cv::Mat> rggb = split_bayer(src, bayer_pattern);
	uchar pads[4] = { 4, 4, 4, 4 };
	cv::Mat src_p = padding(src, pads);
	std::vector<cv::Mat> rggb_pad = split_bayer(src_p, bayer_pattern);
	ushort* R, * B, chroma_corrected = 0;
	float avg_r = 0, avg_g = 0, avg_b = 0, y = 0, fade1 = 0, fade = 0, damp_factor_r = 0, damp_factor_b = 0, max_avg = 0;

	if (gain[0] <= 1024)	damp_factor_r = 256;
	else if (gain[0] <= 1229)	damp_factor_r = 128;
	else damp_factor_r = 77;

	if (gain[2] <= 1024)	damp_factor_b = 256;
	else if (gain[2] <= 1229)	damp_factor_b = 128;
	else damp_factor_b = 77;
	for (int i = 0; i < rggb[0].rows; ++i) {
		std::cout << "\r" << "CNF: ";
		std::cout << std::setw(8) << std::fixed << std::setprecision(2) << (float)i / rggb[0].rows * 100 << "%";
		R = rggb[0].ptr<ushort>(i);
		B = rggb[3].ptr<ushort>(i);
		for (int j = 0; j < rggb[0].cols; ++j) {
			avg_r = sum_cvMat(rggb_pad[0](cv::Rect(j, i, 5, 5))) / 25.0;
			avg_g = sum_cvMat(rggb_pad[1](cv::Rect(j, i, 5, 5))) / 50.0 + sum_cvMat(rggb_pad[2](cv::Rect(j, i, 5, 5))) / 50.0;
			avg_b = sum_cvMat(rggb_pad[3](cv::Rect(j, i, 5, 5))) / 25.0;

			y = (306.0 * avg_r + 601.0 * avg_g + 117.0 * avg_b) / 1024.0;
			fade1 = cnf_fade(y, 'y');

			// r_noise
			if ((R[j] - avg_g) > thres && (R[j] - avg_b) > thres && (avg_r - avg_g) > thres && (avg_r - avg_b) < thres) {
				if (avg_g > avg_b) {
					max_avg = avg_g;
				}
				else {
					max_avg = avg_b;
				}
				chroma_corrected = max_avg + (damp_factor_r * R[j] - max_avg) / 256;
				fade = fade1 * cnf_fade(avg_r, 'r');
				R[j] = ushort(fade * chroma_corrected + (1 - fade) * R[j]);
			}
			// b_noise
			if ((B[j] - avg_g) > thres && (B[j] - avg_r) > thres && (avg_b - avg_g) > thres && (avg_b - avg_r) < thres) {
				if (avg_g > avg_r) {
					max_avg = avg_g;
				}
				else {
					max_avg = avg_r;
				}
				chroma_corrected = max_avg + (damp_factor_b * B[j] - max_avg) / 256;
				fade = fade1 * cnf_fade(avg_b, 'b');
				B[j] = ushort(fade * chroma_corrected + (1 - fade) * B[j]);
			}
		}
	}
	src = reconstruct_bayer(rggb, bayer_pattern);
	clock_t end = clock();
	std::cout << "\r" << "CNF Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
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

void bnf(cv::Mat& src, float bnf_dw[5][5], uchar* bnf_rw, uchar* bnf_rthres, ushort clip)
{
	clock_t begin = clock();

	uchar pads[4] = { 2, 2, 2, 2 };
	cv::Mat src_p = padding(src, pads);
	int sum_weight = 0, sum_imgA = 0;
	ushort pixel_center = 0, rdiff = 0;
	double dive = .0;
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
					rdiff = abs(p_p[x + j] - pixel_center);
					if (rdiff >= bnf_rthres[0]) {
						rdiff = bnf_rw[0];
					}
					else{
						if (rdiff < bnf_rthres[0] && rdiff > bnf_rthres[1]) {
							rdiff = bnf_rw[1];
						}
						else{
							if (rdiff < bnf_rthres[0] && rdiff > bnf_rthres[2]){
								rdiff = bnf_rw[2];
							}
							else{
								rdiff = bnf_rw[3];
							}
						}
					}
					sum_imgA += (p_p[x+j] * rdiff * bnf_dw[i][j]);
					sum_weight += (rdiff * bnf_dw[i][j]);
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

double emlut(double val, uchar* thres, uchar* gain, char* clip)
{
	double lut = .0;
	if (val < -thres[1]) {
		lut = gain[1] * val;
	}
	else{
		if (val < -thres[0] && val > -thres[1]) {
			lut = .0;
		}
		else
		{
			if (val < thres[0] && val > -thres[1]) {
				lut = gain[0] * val;
			}
			else
			{
				if (val > thres[0] && val < thres[1]) {
					lut = .0;
				}
				else
				{
					if (val > thres[1])	lut = gain[1] * val;
				}
			}
		}
	}
	lut = MAX(clip[0], MIN(lut / 256, clip[1]));
	return lut;
}

void ee(cv::Mat& src, cv::Mat& edgemap, char edge_filter[3][5], uchar* ee_thres, uchar* ee_gain, char* ee_emclip, ushort clip)
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
			em_img = em_img / 8.0;
			ee_img = src_p.ptr<ushort>(y + 1)[x + 2] + emlut(em_img, ee_thres, ee_gain, ee_emclip);
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

void fcs(cv::Mat& src_u, cv::Mat& src_v, cv::Mat& edgemap, uchar fcs_delta_min , uchar fcs_delta_max, ushort clip)
{
	clock_t begin = clock();
	double threshold_delta = fcs_delta_max - fcs_delta_min;
	threshold_delta = (threshold_delta < 1e-6) ? 1e-6 : threshold_delta;
	int slope = -(65536 / threshold_delta);

	ushort* p_u, * p_v, * p_edgemap;
	ushort pixel_u = 0, pixel_v = 0, pixel_edgemap = 0;
	int gain_map = 0;
	for (int y = 0; y < src_u.rows; ++y) {
		p_u = src_u.ptr<ushort>(y);
		p_v = src_v.ptr<ushort>(y);
		p_edgemap = edgemap.ptr<ushort>(y);
		for (int x = 0; x < src_u.cols; ++x) {
			pixel_u = p_u[x];
			pixel_v = p_v[x];
			gain_map = slope * (abs(p_edgemap[x]) - fcs_delta_max);
			gain_map = gain_map < 0 ? 0 : (gain_map > 65536 ? 65536 : gain_map);
			pixel_u = (pixel_u - 128) * gain_map / 65536 + 128;
			pixel_v = (pixel_v - 128) * gain_map / 65536 + 128;
			pixel_u = (pixel_u < 0) ? 0 : (pixel_u > clip ? clip : pixel_u);
			pixel_v = (pixel_v < 0) ? 0 : (pixel_v > clip ? clip : pixel_v);
			p_u[x] = pixel_u;
			p_v[x] = pixel_v;
		}
	}
	clock_t end = clock();
	std::cout << "FCS Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
}

void hsc(cv::Mat& src_u, cv::Mat& src_v, uchar hue, ushort saturatin, ushort clip)
{
	clock_t begin = clock();
	std::vector<ushort> lut_sin{}, lut_cos{};
	for (int i = 0; i < 360; ++i) {
		lut_sin.push_back(round(sin(i * 3.1415926 / 180) * 256.0));
		lut_cos.push_back(round(cos(i * 3.1415926 / 180) * 256.0));
	}
	ushort* p_u, * p_v;
	ushort pixel_u, pixel_v;
	for (int y = 0; y < src_u.rows; ++y) {
		p_u = src_u.ptr<ushort>(y);
		p_v = src_v.ptr<ushort>(y);
		for (int x = 0; x < src_u.cols; ++x) {
			pixel_u = (p_u[x] - 128) * lut_cos[hue] + (p_v[x] - 128) * lut_sin[hue] + 128;
			pixel_v = (p_v[x] - 128) * lut_cos[hue] + (p_u[x] - 128) * lut_sin[hue] + 128;
			pixel_u = saturatin * (p_u[x] - 128) / 256.0 + 128;
			pixel_v = saturatin * (p_v[x] - 128) / 256.0 + 128;
			pixel_u = pixel_u < 0 ? 0 : (pixel_u > clip ? clip : pixel_u);
			pixel_v = pixel_v < 0 ? 0 : (pixel_v > clip ? clip : pixel_v);
			p_u[x] = pixel_u;
			p_v[x] = pixel_v;
		}
	}
	clock_t end = clock();
	std::cout << "HSC Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC * 1000 << " ms." << std::endl;
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
