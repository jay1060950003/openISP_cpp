# openISP_cpp
Open Image Signal Processor
You can found the description in [Practice record](https://jay1060950003.github.io/2022/11/02/isp_pipeline%E7%BB%83%E4%B9%A0%E5%AE%9E%E5%BD%95/)

## Introduction

As the name implies, openISP_cpp is a C++ implementation of the [openISP](https://github.com/cruxopen/openISP) project.

Here is the running time in my i5-12500H machine with the 4032x3024 input Bayer array:

|Module             |openISP(1920x1080) |cpp(4032x3024)|
|:-----------------:|:------:|:----------:|
|DPC                |20.57s  |17.8s       |
|BLC                |11.75s  |0.38s       |
|AAF                |16.87s  |0.52s       |
|AWB                |7.54s   |0.41s       |
|CNF                |73.99s  |-----       |
|BNF                |------  |1.16s       |
|CFA                |40.71s  |0.05s       |
|CCM                |56.85s  |0.93s       |
|GAC                |25.71s  |0.61s       |
|CSC                |60.32s  |0.34s       |
|NLM                |1600.95s|12.55s      |
|BNF                |801.24s |8.18s       |
|EEH                |68.60s  |1.65s       |
|FCS                |25.07s  |-----       |
|HSC                |56.34s  |-----       |
|End-to-end pipeline|2894.41s|47.99s       |

## Usage

1. Clone this repo
2. Compile (I use vs2019 which configured `opencv4.6.0`)
3. Adjust the path and parameters of the code, than you can do the pipeline

**Note:**
- I used my `iPhone 11` to get the `raw dng file`. You can try to get the raw data in dng format with your own device. It is very simple on Android phones. On an iPhone you may need `Pro device `or apps like `Procam` or `Halide`
- To process the dng file, libraw is used to get the data

## Thanks to

[libraw](https://github.com/LibRaw/LibRaw): Library for reading and processing of RAW digicam images

[openISP](https://github.com/cruxopen/openISP): Open Image Signal Processor (openISP)

[fast-openISP](https://github.com/QiuJueqin/fast-openISP): Fast Open Image Signal Processor (fast-openISP)