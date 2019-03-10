# ProteoML
Machine learning based Classification and Prediction of Intrinsically Disordered Proteins/Regions in Human Proteome

With applications in Genetics, Bioinformatics and Pharmacy, our model can be used for expediting research in the cure of neurodegenerative disorders like Alzheimer’s, Parkinson’s and cancers like Leukemia, Sarcoma, along with bypassing complex lab techniques, all in the field of Health Care.
The first protein structures and their sequence, after complex isolation and purification, were discovered and studied using crystallography methods initially and in more recent years, using spectroscopy. Until quite recently, proteins were thought to have fixed 3D structures. But recent advances have shown that some proteins (IDPs) and some regions of proteins (IDRs) also have a disordered structure that have both important functional and harmful implications. Their ability to transiently shift between phases, along with flexible multi-domain assemblies, have come to be known for important physiological functions as well as characteristic of certain diseased states.
Using DisProt as the database for IDPs and IDRs, which contains only a fraction of all proteins and regions of the human proteome whose data was obtained from UniProt - IDPs are ~800 in # and IDRs are ~2000 in number compared to the entire human proteome which is ~70000 in number - we have tried to predict which among the ~70000 might have the features for disorderly-ness, based on 11 characteristics. Prediction was done on the human proteome data set with an accuracy of ~82%. This model helps to obtain these characteristics of proteins in diseases much faster, thereby helping to build the IDP/IDR database.

# Architecture

Random Forest Classifier
4000 estimators
Feature selection metric - Gini index
Pruning on depth - 11

# Results

Accuracy = ~82%
sensitivity: 0.8767334360554699
specificity: 0.7688098495212038
Matthew's correlation coefficient: 0.6459675624258796
