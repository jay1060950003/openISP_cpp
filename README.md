# openISP_cpp
Open Image Signal Processor

## Introduction

As the name implies, openISP_cpp is a C++ implementation of the [openISP](https://github.com/cruxopen/openISP) project.

Here is the running time in my i5-12500H machine with the 4032x3024 input Bayer array:

|Module             |openISP(1920x1080) |cpp(4032x3024)|
|:-----------------:|:------:|:----------:|
|DPC                |20.57s  |0.37s       |
|BLC                |11.75s  |0.42s       |
|AAF                |16.87s  |0.59s       |
|AWB                |7.54s   |0.45s       |
|CNF                |73.99s  |2.62s       |
|CFA                |40.71s  |0.07s       |
|CCM                |56.85s  |1.08s       |
|GAC                |25.71s  |0.68s       |
|CSC                |60.32s  |0.35s       |
|NLM                |1600.95s|12.34s       |
|BNF                |801.24s |3.62s       |
|EEH                |68.60s  |1.65s       |
|FCS                |25.07s  |0.08s       |
|HSC                |56.34s  |0.22s       |
|End-to-end pipeline|2894.41s|28.28s       |

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