## Machine learning models for classification

In a first instance, we are placing here a basic workflow for machine learning classification of binary problems 
based on a set of genomic features.
The first examples will be based on the R library `caret`, and will address the classification of taxonomic groups 
based on genomic features extracted with the tool GBRAP ([Yaddehige et al. 2025](https://doi.org/10.1093/molbev/msaf114)).

In this first approach we use the **RFE (recursive feature elimination) algorithm** to select a subset of features 
to be later used in the final classification model.

A second script [classify_probiotics.R](workflow/classify_probiotics.R) is used to tell probiotics from non probiotics 
(an example of **one-class classification** problem), based on a set of genomic features.

### ML for probiotics discovery

1. training/test split: [make_test.R](https://github.com/filippob/ml_classification/blob/main/support_scripts/make_test.R)
2. model tuning:
    -  [twoclass-svm-tuning.r](https://github.com/filippob/ml_classification/blob/main/model_scripts/twoclass-svm-tuning.r)
    -  [oneclass-svm-tuning.r](https://github.com/filippob/ml_classification/blob/main/model_scripts/oneclass-svm-tuning.r)
3. i) **model evaluation** and ii) **final predictions** on unlabelled data:
   - [twoclass-svm-predict.r](https://github.com/filippob/ml_classification/blob/main/model_scripts/twoclass-svm-predict.r)
   - [oneclass-predict.r](https://github.com/filippob/ml_classification/blob/main/model_scripts/oneclass-predict.r)
4. parse results: use the script `support_scripts/parse_results.r` ([here](https://github.com/filippob/ml_classification/blob/main/support_scripts/parse_results.r)) to collect all results from the different model runs:
    - results are saved separately per problem (twoclass, oneclass classification) and method (SVM)
    - the different replicates show the variability of results
