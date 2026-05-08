# NMIBCClassifier
 
<!-- badges: start -->
<!-- badges: end -->

A transcriptomic model for predicting non-muscle-invasive bladder cancer (NMIBC) consensus subtypes. This approach utilizes a single-sample classifier based on the Pearson nearest-centroid method, which evaluates the correlation between the transcriptomic profile of an individual sample and the representative profiles of the four NMIBC integrated molecular clusters (IMC1-4) as described in Liu et al (to be published).

## Citation
You can cite "An integrated Molecular Classification System Identifies Distinct Prognostic Clusters of Non-muscle-invasive Bladder Cancer".

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

classify_nmibc(dejongA)

The classifier returns a datafram with 7 columns:

`IMC` provides the predicted IMC label for each sample. If the maximum correlation for a sample is below the specified cor_cut threshold, the classification result is set to NA, indicating low confidence in the prediction.

`p_value` represents the p-value associated with the Pearson correlation between the sample and its nearest centroid. Smaller p-values indicate stronger statistical significance for the correlation.

`normalized_difference` measures how distinct a sample is from other consensus classes. It ranges from 0 to 1, where 0 indicates that the sample is too close to other consensus classes to be confidently assigned a single class, and 1 indicates that the sample is highly representative of its assigned class and well-separated from others. The normalized difference is calculated as:

$$
\frac{\text{correlation to nearest centroid} - \text{correlation to second nearest centroid}}{\text{correlation to nearest centroid}}
$$

The remaining columns `correlation_to_IMC_X`, provide the Pearson correlation values between each sample and each IMC, where `X` represents the class number. These values help to understand the relationship between the sample and all IMCs.
