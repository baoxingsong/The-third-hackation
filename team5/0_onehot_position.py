import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import LabelEncoder, OneHotEncoder


def datadeal_generator(file_path):
    data = pd.read_csv(file_path, header=None)
    data = data.dropna()
    output_labels = data.iloc[:, 2].values.astype(float)
    np.save("./Ribo_median.npy", output_labels)
    def string_to_array(my_string):
        my_array = np.array(list(my_string))
        return my_array

    label_encoder = LabelEncoder()
    label_encoder.fit(np.array(['A', 'C', 'G', 'T', 'N']))
    fit_array = np.array([0, 1, 2, 4])
    fit_array = fit_array.reshape(4, 1)

    def one_hot_encoder(my_array):
        integer_encoded = label_encoder.transform(my_array)
        onehot_encoder = OneHotEncoder(sparse=False, dtype=int, handle_unknown="ignore")
        integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
        onehot_encoder.fit(fit_array)
        onehot_encoded = onehot_encoder.transform(integer_encoded)
        return onehot_encoded

    for index, row in data.iterrows():
        #编码指定列
        col_encode = row[1]
        arr = string_to_array(col_encode)
        encoded_arr = one_hot_encoder(arr)
        yield encoded_arr
    del data


def process_data(file_path,out_path):
    array_generator = datadeal_generator(file_path=file_path)
    print(type(array_generator))
    array_list = list(array_generator)  # 调用生成器以获取所有编码数组
    del array_generator
    stacked_array = np.stack(array_list, axis=2)
    stacked_array = np.transpose(stacked_array, (2, 0, 1))
    print(stacked_array.shape)

    def positional_encoding(matrix_shape):
        L, d_model = matrix_shape
        pos_enc = np.zeros((L, d_model))
        for pos in range(L):
            for i in range(0, d_model, 2):
                pos_enc[pos, i] = np.sin(pos / (10000 ** ((2 * i) / d_model)))
                pos_enc[pos, i + 1] = np.cos(pos / (10000 ** ((2 * (i + 1)) / d_model)))
        return pos_enc

    encoded_matrices = []
    for matrix in stacked_array:
        L, d_model = matrix.shape
        pos_enc_matrix = positional_encoding((L, d_model))
        encoded_matrix = matrix + pos_enc_matrix
        encoded_matrices.append(encoded_matrix)
    del stacked_array
    # 创建一个新的列表来存储转换后的数组
    encoded_matrices_float32 = []
    for matrix in encoded_matrices:
        encoded_matrix_float32 = matrix.astype(np.float32)
        encoded_matrices_float32.append(encoded_matrix_float32)
    del encoded_matrices
    # 将列表转换为NumPy数组
    encoded_tensor_float32 = np.stack(encoded_matrices_float32)
    np.save(out_path, encoded_tensor_float32)
    del encoded_tensor_float32
# 调用处理函数
process_data(file_path='./SorghumInput.csv',out_path='./up_position.npy')



