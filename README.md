## Machine learning models for classification

In a first instance, we are placing here a basic workflow for machine learning classification of binary problems 
based on a set of genomic features.
The first examples will be based on the R library `caret`, and will address the classification of taxonomic groups 
based on genomic features extracted with the tool GBRAP ([Vischioni et al. 2022](https://www.biorxiv.org/content/10.1101/2021.09.21.461110v2)).

In this first approach we use the **RFE (recursive feature elimination) algorithm** to select a subset of features 
to be later used in the final classification model.

A second script [classify_probiotics.R](workflow/classify_probiotics.R) is used to tell probiotics from non probiotics 
(an example of **one-class classification** problem), based on a set of genomic features.
