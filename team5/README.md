# Deep learning models of ribosome expressions



List of Python scripts (.py)

0_onehot_position: one-hot encoding for the original sequential data and positional encoding for transformer<br>
input file:SorghumInput.10000line.csv.zip -- targeted DNA sequence need to be mannual set in 0_onehot_position. <br>
output file:  numpy format file 

1_Onehot_cds_CNN: using one-hot encoded cds sequences to predict Ribo-seq TPM with a CNN model <br>
input file: numpy format file from 0_onehot_position for one-hot encoded cds sequences<br>
output file: trained models and figure of training process and R-squared of training sets and validation sets

2_Onehot_CNN_mRNA_Dense: using one-hot encoded cds sequences and RNA-seq TPM to predict Ribo-seq TPM with a CNN-dense mixed model <br>
input file:  numpy format file from 0_onehot_position for one-hot encoded cds sequences and RNA-seq TPM<br>
output file:  trained models and figure of training process and R-squared of training sets and validation sets

3_Embbed_cds_CNN_mRNA_Dense:using embedded cds sequences and RNA-seq TPM to predict Ribo-seq TPM with a CNN-dense mixed model <br>
input file:  numpy format file from 0_onehot_position for embedded cds sequences and RNA-seq TPM<br>
output file:  trained models and figure of training process and R-squared of training sets and validation sets

4_Onehot_trans: using one-hot encoded cds sequences to predict Ribo-seq TPM with a transformer model <br>
input file: numpy format file from 0_onehot_position for one-hot encoded cds sequences <br>
output file:  trained models and figure of training process and R-squared of training sets and validation sets


Ref 1. [Building deep learning networks for predicting protein expression from DNA sequences using small sample datasets](https://www.nature.com/articles/s41467-022-34902-5)

Ref 2. [Riboformer](https://github.com/lingxusb/Riboformer.git)
