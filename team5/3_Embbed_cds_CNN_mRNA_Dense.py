import os
import random
import pickle
import scipy
import scipy.stats
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Flatten, Conv2D, MaxPooling2D, MultiHeadAttention, Add, \
    BatchNormalization, Dropout,Multiply,Embedding
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

data1 = np.load("./CDS_embbed.npy")
data1 = np.expand_dims(data1, axis=-1)
data2 = np.load("./TPM.npy")
label = np.load("./Ribo_median.npy")


# 按照8:2的比例划分训练集和验证集
X1_train, X1_val,X2_train, X2_val, y_train, y_val = train_test_split(data1,data2, label, test_size=0.2, random_state=42)


def r_square(y_true, y_pred):
    from keras import backend as K
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res / (SS_tot + K.epsilon()))


def build_model(input_shape1, input_shape2):
    inputs1 = Input(shape=input_shape1)
    inputs2 = Input(shape=input_shape2)

    # 添加嵌入层，词向量维度设为8
    x = Embedding(input_dim=12646, output_dim=8)(inputs1)

    x = tf.expand_dims(x, axis=-1)

    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2, 1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Flatten()(x)
    x = Dense(1024, activation='relu')(x)

    x_mRNA = Dense(1024, activation='relu')(inputs2)
    multiplied = Multiply()([x, x_mRNA])
    x = Dense(512, activation='relu')(multiplied)
    x = Dense(256, activation='relu')(x)
    outputs = Dense(1, activation='softplus')(x)

    model = Model([inputs1, inputs2], outputs)
    return model



input_shape1 = (12646)
input_shape2 = (1)

model = build_model(input_shape1, input_shape2)
model.summary()


optimizer = Adam(lr=0.001,
                 beta_1=0.9,
                 beta_2=0.999,
                 epsilon=1e-08,
                 decay=0.0)
model.compile(loss='mean_squared_error',
              optimizer=optimizer,
              metrics=[r_square])

checkpoint_path = "./training_entire"
earlystop = EarlyStopping(monitor='val_loss',
                          mode='min',
                          verbose=1,
                          patience=20)

modelcheckpoint = ModelCheckpoint(filepath=checkpoint_path,
                                  monitor='val_loss',
                                  verbose=1,
                                  save_best_only=True,
                                  mode='auto')

epoch = 80
batch_size = 4

history = model.fit([X1_train, X2_train], y_train, epochs=epoch, batch_size=batch_size,
                    validation_data=([X1_val, X2_val], y_val),
                    callbacks=[modelcheckpoint, earlystop],
                    verbose=1,
                    shuffle=True)


# 获取训练过程中的损失值历史
train_loss = history.history['loss']
val_loss = history.history['val_loss']

# 生成 epoch 数组，用于横坐标
epochs = range(1, len(train_loss) + 1)

# 绘制训练集和验证集的损失值变化曲线，使用不同颜色但相同形状的线条
plt.figure()
plt.plot(epochs, train_loss, 'b-', label='Training loss')  # 'b-' 表示蓝色实线
plt.plot(epochs, val_loss, 'g-', label='Validation loss')  # 'g-' 表示绿色实线
plt.title('Training and Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss (MSE)')
plt.legend()
# 保存图像
plt.savefig('./figure/CNN_dense_2.png')