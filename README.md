# NMIBCClassifier
 
<!-- badges: start -->
<!-- badges: end -->

A transcriptomic model for predicting non-muscle-invasive bladder cancer (NMIBC) consensus subtypes. This approach utilizes a single-sample classifier based on the Pearson nearest-centroid method, which evaluates the correlation between the transcriptomic profile of an individual sample and the representative profiles of the three NMIBC consensus groups (cluster 1,2,3) as described in Liu et al (to be published).

## Citation
You can cite ...

## Installation

You may install this package with devtools:

``` r
library(devtools)
devtools::install_github("redmele/NMIBCClassifier", build_vignettes = TRUE)
library(NMIBCClassifier)
```

## Usage

``` r
classify_nmibc(data, cor_cut = 0.1)
```
Where `data` is a dataframe with rows represent unique genes and columns represent samples. Gene identifiers (rownames) can be provided as HUGO gene symbols.RNA-seq data should be log-transformed, for example using log2(normalized counts + 1).

`cor_cut`  is a numeric value specifying the minimum Pearson correlation threshold for classification.Samples with a maximum correlation below this threshold will remain unclassified, and their classification results will be set to NA. Default value is 0.1.

## Example

```{r example}
library(NMIBCClassifier)

data("dejongA")

result <- classify_nmibc(dejongA,cor_cut = 0.1)
head(result)
#     consensus_cluster correlation_to_consensus_cluster_1 correlation_to_consensus_cluster_2
#R1                   1                  0.622375078555726                  0.344249168079062
#R10                  1                  0.661031953629024                  0.430581247870673
#R100                 1                  0.609151076253045                  0.363211102037651
#R105                 1                   0.65241067422918                  0.468824850632178
#R108                 1                   0.62776099581035                  0.390511502972139
#R11                  1                  0.643480455642961                  0.354965986849694
#     correlation_to_consensus_cluster_3 normalized_difference               p_value
#R1                    0.470076987240111     0.244704674983188  2.18947891281627e-90
#R10                   0.569648494614646     0.138243633326179 1.24201813563067e-105
#R100                   0.50630296106019     0.168838436312852  1.19687055940123e-85
#R105                   0.59785028133111    0.0836289089882427 4.88183486773281e-102
#R108                  0.527453480451906     0.159786154329262  2.21243686288274e-92
#R11                    0.51128412713991     0.205439539528768  1.94927108685008e-98
```

The classifier returns a datafram with 6 columns:

`consensus_cluster` provides the predicted consensus class label for each sample. If the maximum correlation for a sample is below the specified cor_cut threshold, the classification result is set to NA, indicating low confidence in the prediction.

`p_value` represents the p-value associated with the Pearson correlation between the sample and its nearest centroid. Smaller p-values indicate stronger statistical significance for the correlation.

`normalized_difference` measures how distinct a sample is from other consensus classes. It ranges from 0 to 1, where 0 indicates that the sample is too close to other consensus classes to be confidently assigned a single class, and 1 indicates that the sample is highly representative of its assigned class and well-separated from others. The normalized difference is calculated as:

$$
\frac{\text{correlation to nearest centroid} - \text{correlation to second nearest centroid}}{\text{correlation to nearest centroid}}
$$

The remaining columns `correlation_to_consensus_cluster_X`, provide the Pearson correlation values between each sample and each consensus class, where `X` represents the class number. These values help to understand the relationship between the sample and all consensus classes.
