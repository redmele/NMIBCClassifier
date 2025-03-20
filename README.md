# NMIBCClassifier
 
<!-- badges: start -->
<!-- badges: end -->

A transcriptomic model for predicting non-muscle-invasive bladder cancer (NMIBC) consensus subtypes. This approach utilizes a single-sample classifier based on the Pearson nearest-centroid method, which evaluates the correlation between the transcriptomic profile of an individual sample and the representative profiles of the four NMIBC consensus clusters (CMC1-4) as described in Liu et al (to be published).

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
classify_nmibc(data, cor_cut = 0.2)
```
Where `data` is a dataframe with rows represent unique genes and columns represent samples. Gene identifiers (rownames) can be provided as HUGO gene symbols.RNA-seq data should be log-transformed, for example using log2(normalized counts + 1).

`cor_cut`  is a numeric value specifying the minimum Pearson correlation threshold for classification.Samples with a maximum correlation below this threshold will remain unclassified, and their classification results will be set to NA. Default value is 0.2.

## Example

```{r example}
library(NMIBCClassifier)

data("dejongA")

result <- classify_nmibc(dejongA,cor_cut = 0.1)
head(result)
     consensus_cluster correlation_to_consensus_cluster_1 correlation_to_consensus_cluster_2
R1                   3                  0.615020049284025                  0.709953064725233
R10                  3                  0.618341583708948                  0.701746159370815
R100                 3                  0.595296487937614                  0.711287595982398
R105                 3                    0.6608163932423                  0.728995427667581
R108                 3                  0.615290770820966                  0.711151843594884
R11                  3                  0.613669640349605                   0.71223493452466
     correlation_to_consensus_cluster_3 correlation_to_consensus_cluster_4 normalized_difference
R1                    0.761157907307797                  0.730567653528214    0.0401891033199404
R10                   0.739127244744341                  0.718764561738891    0.0275496312038855
R100                  0.748335812411318                  0.711522549210876    0.0491935072328304
R105                  0.754232178835922                  0.740848808975043    0.0177443633889165
R108                  0.756398672762588                  0.715757671391593    0.0537296042873293
R11                   0.749577035998255                  0.719597275535167    0.0399955695323056
                   p_value
R1   2.76647108611937e-263
R10  1.79137685318766e-240
R100 1.01868170552663e-249
R105 7.47683498835475e-256
R108 3.76523257100649e-258
R11  5.38293473229588e-251
```

The classifier returns a datafram with 7 columns:

`consensus_cluster` provides the predicted consensus class label for each sample. If the maximum correlation for a sample is below the specified cor_cut threshold, the classification result is set to NA, indicating low confidence in the prediction.

`p_value` represents the p-value associated with the Pearson correlation between the sample and its nearest centroid. Smaller p-values indicate stronger statistical significance for the correlation.

`normalized_difference` measures how distinct a sample is from other consensus classes. It ranges from 0 to 1, where 0 indicates that the sample is too close to other consensus classes to be confidently assigned a single class, and 1 indicates that the sample is highly representative of its assigned class and well-separated from others. The normalized difference is calculated as:

$$
\frac{\text{correlation to nearest centroid} - \text{correlation to second nearest centroid}}{\text{correlation to nearest centroid}}
$$

The remaining columns `correlation_to_consensus_cluster_X`, provide the Pearson correlation values between each sample and each consensus class, where `X` represents the class number. These values help to understand the relationship between the sample and all consensus classes.
