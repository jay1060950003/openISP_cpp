#include <iostream>
#include <vector>
#include <ctime>

#include "libraw.h"
#include "module.h"

int main() {
    cv::utils::logging::setLogLevel(cv::utils::logging::LOG_LEVEL_ERROR);   // 只输出错误日志
    clock_t begin = clock();

    const char* file = "IMG_1.dng";
    // open dng file
    LibRaw* iProcessor = new LibRaw;    
    iProcessor->open_file(file);
    iProcessor->unpack();

    // raw2Mat and cut2roi
    cv::Mat raw = cv::Mat(iProcessor->imgdata.sizes.raw_height, iProcessor->imgdata.sizes.raw_width, CV_16UC1, iProcessor->imgdata.rawdata.raw_image);
    raw = raw(cv::Rect(iProcessor->imgdata.sizes.left_margin, iProcessor->imgdata.sizes.top_margin, iProcessor->imgdata.sizes.width, iProcessor->imgdata.sizes.height));

    // init parameter
    uchar dpc_thres = 30;
    std::string BAYER = "rggb";
    float black = iProcessor->imgdata.color.dng_levels.dng_black;
    float white = iProcessor->imgdata.color.dng_levels.dng_whitelevel[0];
    cv::Mat aaf_kernel = (cv::Mat_<float>(5, 5) <<
        1, 0, 1, 0, 1,
        0, 0, 0, 0, 0,
        1, 0, 8, 0, 1,
        0, 0, 0, 0, 0,
        1, 0, 1, 0, 1) / 16;
    std::vector<float> parameter(iProcessor->imgdata.color.dng_levels.asshotneutral, iProcessor->imgdata.color.dng_levels.asshotneutral + 3);
    ushort CLIP = 4095;
    uchar cnf_thres = 60;
    uchar bnr_size = 3;
    cv::Mat ccmatrix = cv::Mat(3, 4, CV_32F, iProcessor->imgdata.color.cmatrix)(cv::Rect(0, 0, 3, 3));
    float gamma = 0.42;
    float RGB2YUV420[3][4] = {
        {0.257, 0.504, 0.098, 16},
        {-.148, -.291, 0.439, 128 },
        {0.439, -.368, -.071, 128},
    };
    uchar nlm_dw = 1;
    uchar nlm_Dw = 3;
    uchar nlm_thres = 10;
    float bnf_dw[5][5] = {
        {0.5, 0.75,  2, 0.75, 0.5 },
        {0.75,   4,  8,    4, 0.75},
        {2,      4, 64,     4,    2},
        {0.75,   4,  8,    4, 0.75},
        {0.5, 0.75,  2, 0.75, 0.5 },
    };
    uchar bnf_rw[4] = { 0, 8, 16, 32 };
    uchar bnf_rthres[3] = { 128, 32, 8 };
    char edge_filter[3][5] = {
        {-1, 0, -1, 0, -1},
        {-1, 0, 8,  0, -1},
        {-1, 0, -1, 0, -1} };
    uchar ee_gain[2] = { 32, 128 };
    uchar ee_thres[2] = { 32, 64 };
    char ee_emclip[2] = { -64, 64 };
    uchar brightness = 10;
    uchar contrast = 10;
    uchar fcs_delta_min = 8;
    uchar fcs_delta_max = 32;
    uchar hue = 128;
    ushort saturation = 256;

    dpc(raw, dpc_thres);                        //DPC
    blc(raw, BAYER, black, white);              //BLC
    aaf(raw, aaf_kernel);                       //AAF
    wbgc(raw, BAYER, parameter, CLIP);          //WBGC
    //cnf(raw, BAYER, parameter, cnf_thres);      //CNF
    bnr(raw, BAYER, bnr_size);
    cfa(raw, BAYER);                            //CFA
    ccm(raw, ccmatrix);                         //CCM
    gc(raw, gamma, CLIP);                       //GC
    cv::Mat raw_ = raw.clone();
    raw_.convertTo(raw_, CV_8UC3);
    cv::imwrite("result_raw.png", raw_);

    cv::Mat y = cv::Mat(cv::Size(raw.cols, raw.rows), CV_16U);
    cv::Mat u = cv::Mat(cv::Size(raw.cols, raw.rows), CV_16U);
    cv::Mat v = cv::Mat(cv::Size(raw.cols, raw.rows), CV_16U);
    csc(raw, y, u, v, RGB2YUV420);              //CSC
    nlm(y, nlm_dw, nlm_Dw, nlm_thres, CLIP);    //NLM
    
    bnf(y, bnf_dw, bnf_rw, bnf_rthres, CLIP);   //BNF
    cv::Mat edge_map = cv::Mat(cv::Size(y.cols, y.rows), CV_16U);
    ee(y, edge_map, edge_filter, ee_thres, ee_gain, ee_emclip, CLIP);           //EE


    bcc(y, brightness, contrast, CLIP);                                         //BCC
    //fcs(u, v, edge_map, fcs_delta_min, fcs_delta_max, CLIP);                    //FCS，效果不好
    hsc(u, v, hue, saturation, CLIP);                                           //HSC
    yuv2rgb(raw, y, u, v);                                                      //2RGB
    raw.convertTo(raw, CV_8UC3);
    cv::imwrite("result.png", raw);
    y.convertTo(y, CV_8UC1);
    cv::imwrite("result_y.png", y);
    u.convertTo(u, CV_8UC1);
    cv::imwrite("result_u.png", u);
    v.convertTo(v, CV_8UC1);
    cv::imwrite("result_v.png", v);
    //edge_map.convertTo(edge_map, CV_8UC1);
    //cv::imwrite("result_edgemap.png", edge_map);

    clock_t end = clock();    
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "ISP Done! Elapsed " << double(end - begin) / CLOCKS_PER_SEC << " s." << std::endl;
    // recycle LibRaw and delete ptr
    iProcessor->recycle();
    delete iProcessor;
}