# SDNE-MDA
Previous studies indicated that miRNA plays an important role in human biological processes especially in the field of diseases. However, constrained by biotechnology, only a small part of the miRNA-disease associations has been verified by biological experiment. This impel that more and more researchers pay attention to develop efficient and high-precision computational methods for predicting the potential miRNA-disease associations. Based on the assumption that molecules are related to each other in human physiological processes, we developed a novel structural deep network embedding model (SDNE-MDA) for predicting miRNA-disease association using molecular associations network. Specifically, the SDNE-MDA model first integrating miRNA attribute information by Chao Game Representation (CGR) algorithm and disease attribute information by disease semantic similarity. Secondly, we extract feature by structural deep network embedding from the heterogeneous molecular associations network. Then, a comprehensive feature descriptor is constructed by combining attribute information and behavior information. Finally, Convolutional Neural Network (CNN) is adopted to train and classify these feature descriptors. In the five-fold cross validation experiment, SDNE-MDA achieved AUC of 0.9447 with the prediction accuracy of 87.38% on the HMDD v3.0 dataset. To further verify the performance of SDNE-MDA, we contrasted it with different feature extraction models and classifier models. Moreover, the case studies with three important human diseases, including Breast Neoplasms, Kidney Neoplasms, Lymphoma were implemented by the proposed model. As a result, 47, 46 and 46 out of top-50 predicted disease-related miRNAs have been confirmed by independent databases. These results anticipate that SDNE-MDA would be a reliable computational tool for predicting potential miRNA-disease associations.