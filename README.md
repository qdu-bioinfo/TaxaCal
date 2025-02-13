# TaxaCal

## Contents

- [Introduction](#Introduction)

- [Package](#Package)

- [Installation](#Installation)

- [Example dataset](#Example-dataset)

- [Training and Calibration](#Training-and-Calibration)

## Introduction

In this study, we present TaxaCal (Taxonomic Calibrator), a machine learning algorithm designed to calibrate species-level taxonomy profiles in 16S amplicon data using a two-tier correction strategy. TaxaCal effectively reduces biases in amplicon sequencing, mitigating discrepancies between microbial profiles derived from 16S and WGS. Moreover, TaxaCal enables seamless cross-platform comparisons between these two sequencing approaches, significantly improving disease detection in 16S-based microbiome data.

## Package

scikit-learn==1.6.1

scipy===1.15.1

numpy==2.2.2

## Installation

```
cd TaxaCal-main
sh init.sh
```

## Example dataset

We provide a demonstration dataset stored in the /data directory. This dataset includes: (1) training_16s_genus.csv, which contains the genus-level relative abundance of 20 training 16S samples; (2) training_wgs_genus.csv, which provides the genus-level relative abundance of the corresponding 20 training WGS samples; (3) training_wgs_species.csv, which contains the species-level relative abundance of the same 20 training WGS samples; and (4) test_16s.csv, which represents the genus-level relative abundance of the 16S sample to be calibrated.

The format of the taxonomic relative abundance tables: 
```
Sample	g_01	g_02	g_03	g_04	g_05	g_06	g_07	g_08	g_09
Sample1	0.1	0	0.3	0.1	0.1	0.1	0.1	0	0.2
Sample2	0.3	0.1	0.1	0	0.1	0.2	0	0.1	0.1
Sample3	0	0.2	0.1	0.3	0	0	0.4	0	0

...
SampleN	0	0.1	0.2	0.4	0	0	0.3	0	0
```
you can set the custom training sample path in the config section under the '/bin' folder.

## Training and Calibration

For this demo dataset, you can quickly complete the training and calibration process by running run.sh in the '/TaxaCal' folder to obtain a calibrated 16S species relative abundance file.

```
cd TaxaCal-main
chmod a+x run.sh
./run.sh
```

### The framework of TaxaCal algorithm

![image](https://github.com/qdu-bioinfo/TaxaCal/blob/main/img.png)


## Result

The calibrated 16S species abundance file will be saved in the specified output directory after the script has finished running.
