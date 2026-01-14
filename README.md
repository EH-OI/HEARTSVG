# HEARTSVG
HEARTSVG: a high-efficiency and robust method for identifying spatial expression patterns in large-scale spatial transcriptomic data
We proposed HEARTSVG, which builds upon serial correlation tests for rapidly and precisely detecting SVGs without relying on the underlying data-generative model or pre-defining spatial patterns.Furthermore, we developed the software HEARTSVG to provide an SVG auto-clustering module for predictions of spatial domain and visualizations, allowing users to discover novel biological phenomena. 

To be able to use the R package HEARTSVG successfully, we suggest manually installing the following packages and their dependent packages first: spatstat, Seurat, plotly, plot3D, TTR, stringr, fpp2, seasonal, cluster, poolr, reshape2, data.table, plyr, dplyr, ggplot2, MASS, EnvStats, BreakPoints, gprofiler2

```{r,warning=F,message=F,results='hide',error=F}
library(spatstat);library(Seurat);library(plot3D);library(plotly);library(TTR);
library(stringr);library(fpp2);library(seasonal);library(cluster);library(poolr);
library(reshape2);library(data.table);library(plyr);library(dplyr);library(ggplot2);
library(MASS);library(EnvStats);library(BreakPoints);library(gprofiler2)

```
## Required R modules
```{r ,warning=F}
R (>= 3.5.0)
```

## Dependencies
```{r ,warning=F}
fpp2 (>=2.4)
seasonal (>=1.9.0)
cluster (>=2.1.3)
poolr (>=1.1.1)
Seurat (4.1.1)
spatstat (>=2.3.4)
ggplot2 (>=3.4.2)
EnvStats (>=2.7.0)
BreakPoints (>=1.2)
plotly (>=4.10.0)
plot3D (>=1.4)
gprofiler2 (>=0.2.1)


```

## Installation

```{r}
library(devtools)
devtools::install_github('https://github.com/cz0316/HEARTSVG',force = T)
library(HEARTSVG)
```

## Operating systems tested on (HEARTSVG 1.1.0) :
macOS Ventura 13.4.1 <br>
Windows 10 <br>


## 2023/12 Update
We optimized the code to achieve reductions in runtime and memory usage, and the optimized function is named **'heartsvg_fast'**. However, due to an unknown error, we are unable to compile the 'heartsvg_fast' function into the R package. As a result, we are providing it separately (file_name: [heartsvg_fast.R](https://github.com/cz0316/HEARTSVG/blob/main/heartsvg_fast.R)). We are currently investigating the cause of the error.




## 0 Load data

We use the ST data from the Wu et al. study as an example to introduce the usage of the R package HEARTSVG. The article and data download link are: https://doi.org/10.1158/2159-8290.CD-21-0316; http://www.cancerdiversity.asia/scCRLM/.


```{r ,warning=F}
stc2.dt=readRDS(file = '~/HEARTSVG/data/HEARTSVG_demo.RDS')
dim(stc2.dt)
## [1] 4174 15429
```
If your data is of type seuratObject, we provide two functions, **filter.gene** and **seurat.trans** , for converting the data into a format suitable for HEARTSVG.
We use the function **filter.gene** to filter the non-expressed genes and low expression proportion genes. Then we use the function **seurat.trans** to transform the data type, from SeuratObject to a data.frame.

```{r}
## stc2=filter.gene(stc2,thr=0.01) # thr=0.01 represents the minimum expression proportion of each gene across all cells.
## dim(stc2)
## [1] 8812 4174

## stc2.dt=seurat.trans(stc2)
## dim(stc2.dt)
## [1] 4174 8814
```


## 1 Find SVG

The input data of the function **heartsvg** is a data.frame with dimensions n*(p+2) contains gene coordinates and gene expression counts. 

The first 2 columns should represent coordinates and their column names must be 'row' and 'col', respectively. The remaining p columns correspond to the expression levels of p genes. There are n rows, representing n spots.
```{r}
demo=stc2.dt[,1:1002]
demo[1:6,1:6]
##                    row col SAMD11 NOC2L HES4 ISG15
## AAACAACGAATAGTTC-1   0  16      0     0    1     1
## AAACAAGTATCTCCCA-1  50 102      0     2    0     2
## AAACAATCTACTAGCA-1   3  43      0     1    0     0
## AAACACCAATAACTGC-1  59  19      1     0   11     6
## AAACAGAGCGACTCCT-1  14  94      0     0    1     2
## AAACAGCTTTCAGAAG-1  43   9      0     0    0     2
```


```{r}
## find SVG
result=heartsvg(demo,scale=T)
head(result)
##       gene pval p_adj rank
## 729  RPS27    0     0    1
## 604   PIGR    0     0    2
## 746 S100A6    0     0    3
## 185  CSRP1    0     0    4
## 725  RPL11    0     0    5
## 732   RPS8    0     0    6
```

Regarding the parameter **'scale'**, the default value is T. This parameter is used to mitigate the impact of extreme values caused by technical factors. 

If you are using technologies such as 10X Visium, which produce data with a moderate level of sparsity and some noise, we recommend setting **scale=T**.

Conversely, if you are using technologies such as HDST, which generate highly sparse data with very low expression levels, we recommend setting **scale=F**. 

Additionally, if you have already performed noise reduction, normalization, or other preprocessing steps on your ST data, we also recommend setting **scale=F**.


## 2 Visualization
Different visualizations of top SVGs.

### 2.1 2D plot
```{r}
svg.2Dplot(data = demo,gene=result$gene[1:3],method = 'raw')
```

![image](https://user-images.githubusercontent.com/57090974/227206644-dfdc02ab-94c6-480d-81a3-f92cb0dcd8a3.png)
![image](https://user-images.githubusercontent.com/57090974/227206706-d284478f-4c75-43d8-9855-c99b65c694be.png)
![image](https://user-images.githubusercontent.com/57090974/227206753-c0825a0e-8c85-4beb-8c73-fc0971ca9469.png)



### 2.2 Marginal expression plot of SVG

The marginal expression of genes are obtained by the semi-pooling process.

```{r}
svg.Mplot(data=demo,gene=result$gene[1:3],method = 'raw')
```

![image](https://user-images.githubusercontent.com/57090974/227206902-da9df34a-803b-4e73-a503-b3e10cc5d23f.png)
![image](https://user-images.githubusercontent.com/57090974/227206963-d036d030-827c-4791-8df3-79aa6a7a9e19.png)
![image](https://user-images.githubusercontent.com/57090974/227207038-98177ad3-46d6-48c8-8f82-9bfb5145f8e4.png)



### 2.3 3D plot

#### 2.3.1 3D scatter plot

```{r}
svg.3Dplot(data=demo,gene = c('RPS8'),type = 'scatter')

```

![newplot](https://user-images.githubusercontent.com/57090974/227207527-221c9d5c-8a7b-47da-8bd7-1f2686dd60b0.png)


#### 2.3.1 3D surface plot

```{r}
svg.3Dplot(data=demo,gene=c('RPS27'),type = 'surface')

```

![image](https://user-images.githubusercontent.com/57090974/227207699-86fa6205-3293-40fe-bbc0-3c74d8b3345e.png)


## 3 Prediction of spatial domains

We offers an auto-clustering module for SVGs that facilitates further analyses, including predicting spatial domains, conducting functional analyses, and visualization based on the set of SVG files.  
We categorize SVGs based on their similarity in spatial expression and adopt a data-driven approach to estimate the number of spatially coherent regions in both SVG expression and histology.

### 3.1 clustering for SVGs
```{r}
clust_demo=svg.clust(data=demo,svg=result$gene[1:500],method = 'h')
table(clust_demo$cluster)
## 
##  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17 
## 403  16  59   6   1   1   2   1   3   1   1   1   1   1   1   1   1 

```

```{r}
km_demo=svg.clust(data=demo,svg=result$gene[1:500],method = 'k',n.c=6)
table(km_demo$cluster)
## 
##   1   2   3   4   5   6 
## 114 139 123  27  68  29
```


### 3.2 visualization of spatial domains

```{r}
svg.2Dplot(data=demo,gene=clust_demogene[clustdemogene[clust_democluster==1],method = 'mean',title = 'domain-1')
```
![image](https://user-images.githubusercontent.com/57090974/227209204-cb369b7f-273b-4568-bc21-5477bb2ed4f0.png)


```{r}
svg.2Dplot(data=demo,gene=clust_demogene[clustdemogene[clust_democluster==2],method = 'mean',title = 'domain-2')
```

cl1.en=svg.enrichment(svg=clust_demo$gene[clust_demo$cluster==1])

head(cl1.en)


## 4 Enrichment

We perform enrichment analysis of SVG sets based on the package gprofiler2.

```{r}
cl1.en=svg.enrichment(svg=clust_demogene[clustdemogene[clust_democluster==1])

head(cl1.en)
##     query significant      p_value term_size query_size intersection_size
## 1 query_1        TRUE 5.007490e-08      5861        405               174
## 2 query_1        TRUE 2.599282e-07      4088        405               132
## 3 query_1        TRUE 5.064712e-07      6348        405               181
## 4 query_1        TRUE 5.795109e-07      5657        405               166
## 5 query_1        TRUE 1.015270e-06      6446        405               182
## 6 query_1        TRUE 2.018340e-06      2664        405                95
##   precision     recall    term_id source
## 1 0.4296296 0.02968777 GO:0048856  GO:BP
## 2 0.3259259 0.03228963 GO:0042221  GO:BP
## 3 0.4469136 0.02851292 GO:0048518  GO:BP
## 4 0.4098765 0.02934418 GO:0048522  GO:BP
## 5 0.4493827 0.02823456 GO:0032502  GO:BP
## 6 0.2345679 0.03566066 GO:0010033  GO:BP
##                                   term_name effective_domain_size source_order
## 1          anatomical structure development                 21128        13923
## 2                      response to chemical                 21128        10676
## 3 positive regulation of biological process                 21128        13622
## 4   positive regulation of cellular process                 21128        13626
## 5                     developmental process                 21128         8073
## 6             response to organic substance                 21128         3939
##                              parents
## 1                         GO:0032502
## 2                         GO:0050896
## 3             GO:0008150, GO:0050789
## 4 GO:0009987, GO:0048518, GO:0050794
## 5                         GO:0008150
## 6                         GO:0042221

```

![Rplot](https://user-images.githubusercontent.com/57090974/227210739-74e236d1-ab4b-4032-9b92-50565af266ad.png)


## Citation
Yuan, X., Ma, Y., Gao, R. et al. HEARTSVG: a fast and accurate method for identifying spatially variable genes in large-scale spatial transcriptomics. Nat Commun 15, 5700 (2024). https://doi.org/10.1038/s41467-024-49846-1


doi: https://doi.org/10.1038/s41467-024-49846-1

---
# Morphology-aware SVG Detection using Spatial Neighborhood Features

## Improving HEARTSVG: A Proxy Feature Approach for Spatial Transcriptomics
## 프록시 추가를 통한 HEARTSVG의 한계점 보완 제안

Project by: Taeyi Kim (김태이)
Base Method: HEARTSVG

📌 Project Overview

이 프로젝트는 대규모 공간 전사체(Spatial Transcriptomics) 데이터 분석 툴인 HEARTSVG를 심층 분석하고, 기존 알고리즘의 한계를 보완하는 새로운 방법론(Morphology-aware Proxy Features)을 제안 및 검증한 연구입니다.

연구는 다음 세 단계로 진행되었습니다:

Paper Review: HEARTSVG 논문 분석 및 기존 방법론 비교

Code Analysis: 알고리즘 로직 분석 및 데모 데이터 구동

Improvement Proposal (Core): 한계점 도출 및 Proxy Feature(조직밀집도, 엔트로피) 제안

1. Background & Limitations

1.1 HEARTSVG Overview

HEARTSVG는 비모수적(Distribution-free) 통계 검정 방식을 사용하여 대규모 공간 전사체 데이터에서 **공간 변이 유전자(SVG)**를 빠르고 정확하게 식별하는 알고리즘입니다.

1.2 Problem Definition (한계점)

HEARTSVG는 공간 좌표($x, y$)와 발현량($e$)만을 변수로 사용합니다. 이로 인해 다음과 같은 한계가 존재합니다.

형태학적 정보(Morphology) 부재: 조직의 밀도나 경계선 같은 병리학적 특징을 반영하지 못함.

단순 구조적 마커 편향: 단순히 세포가 많이 뭉쳐있는 곳의 유전자(예: 미토콘드리아 유전자, Housekeeping gene)가 최상위 랭크(Rank 1~10)를 차지하는 경향이 있음.

Disease-relevant Gene 누락: 정작 암의 증식이나 전이와 관련된 핵심 유전자는 발현량이 낮거나 국소적이어서 하위권으로 밀려나는 문제 발생.

2. Proposed Method: Proxy Features

H&E 염색 이미지와 같은 외부 데이터 없이, 오직 좌표와 발현 데이터만으로 조직의 형태학적 특징을 추정하는 두 가지 **'Proxy Feature'**를 고안하여 기존 알고리즘에 가중치로 적용했습니다.

🛠️ Method 1: Local Density Score (조직 밀집도)

가설: 암 조직은 세포 분열이 활발하여 정상 조직보다 세포 밀도가 높다. 따라서 밀도가 높은 영역의 유전자는 암의 증식/생존과 직결될 것이다.

구현: $k$-NN ($k=10$) 알고리즘을 이용해 이웃 간 평균 거리의 역수로 밀도 점수 산출.


$$Density \propto \frac{1}{\text{mean(distance to k-neighbors)}}$$

🛠️ Method 2: Neighborhood Entropy (이웃 다양성)

가설: 암세포와 정상 세포가 맞닿는 **'종양 경계선'**이나 **'미세환경(TME)'**은 다양한 세포가 혼재되어 있다. 복잡성이 높은 곳의 유전자는 면역 반응/침윤과 관련이 깊다.

구현: 비지도 학습(K-means)으로 가상의 라벨을 부여한 뒤, 주변 이웃($k=20$)의 Shannon Entropy 계산.

3. Experiment Results

논문에서 사용된 대장암(Colorectal Cancer, CRC) 데이터를 사용하여 Baseline(기존)과 제안 모델(Proxy 적용)을 비교 분석했습니다.

📊 1. Baseline (HEARTSVG Original)

Top Rank: MT-ATP6, MT-CO2 (미토콘드리아 유전자), B2M (6위)

해석: 단순히 세포가 존재하는 위치나 구조적 특징만을 포착함. 병리학적 의미가 부족함.

📊 2. Improvement with Density Proxy

Top Rank: LGR5 (4위), UHRF1 (2위)

성과:

기존 3,042위였던 대장암 줄기세포 마커 **LGR5**를 4위로 급상승시킴.

암세포 증식 마커 UHRF1을 2위로 식별.

의의: 암세포가 가장 빽빽하게 뭉쳐있는 **'종양 핵심부(Tumor Core)'**를 정확히 타격함.

📊 3. Improvement with Entropy Proxy

Top Rank: IGLV6-57 (6위), COL10A1, MMP7

성과:

면역 유전자(IGLV)와 기질/침윤 효소(COL, MMP)가 상위권 진입.

조직 전반에 균일한 B2M은 6위 $\to$ 5,982위로 하락 (단순 배경 신호 필터링).

의의: 암세포와 주변 환경이 상호작용하는 **'종양 미세환경(TME) 및 경계선'**을 성공적으로 포착함.

🖼️ Visualization Comparison

Baseline (Noise/Structure)

Density Proxy (Tumor Core)

Entropy Proxy (Boundary/TME)







MT-ATP6 (Rank 1)

LGR5 (Rank 4)

DSCR8 / IGLV6-57

4. Conclusion

본 연구는 고가의 이미지 데이터 처리 과정 없이, **좌표 기반의 수치적 프록시(Density, Entropy)**만으로도 HEARTSVG의 생물학적 민감도를 획기적으로 개선할 수 있음을 입증했습니다.

Density Proxy: 암의 성장과 증식(Core) 규명 특화

Entropy Proxy: 암의 전이와 면역 반응(Boundary) 규명 특화

이는 향후 대규모 공간 전사체 데이터를 저비용으로 신속하게 분석하여 질병의 단계와 특성을 입체적으로 파악하는 새로운 분석 프레임워크가 될 수 있습니다.

