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
    BatchNormalization, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

physical_devices = tf.config.list_physical_devices('GPU')
tf.config.set_visible_devices(physical_devices[1], 'GPU')
tf.config.experimental.set_memory_growth(physical_devices[1], True)

data = np.load("./up_position.npy")
data = np.expand_dims(data, axis=-1)

data.shape

label = np.load("./Ribo_median.npy")
label

# 按照8:2的比例划分训练集和验证集
X_train, X_val, y_train, y_val = train_test_split(data, label, test_size=0.2, random_state=42)

# 检查数据形状
print(f'X_train shape: {X_train.shape}')
print(f'y_train shape: {y_train.shape}')
print(f'X_val shape: {X_val.shape}')
print(f'y_val shape: {y_val.shape}')



def r_square(y_true, y_pred):
    from keras import backend as K
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res / (SS_tot + K.epsilon()))



def build_model(input_shape):
    inputs = Input(shape=input_shape)

    x1 = MultiHeadAttention(num_heads=4, key_dim=96)(inputs, inputs)
    x = Add()([inputs, x1])
    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2,1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x2 = MultiHeadAttention(num_heads=4, key_dim=96)(x, x)
    x = Add()([x, x2])
    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2,1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x3 = MultiHeadAttention(num_heads=4, key_dim=96)(x, x)
    x = Add()([x, x3])
    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2,1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x4 = MultiHeadAttention(num_heads=4, key_dim=96)(x, x)
    x = Add()([x, x4])
    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2,1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x5 = MultiHeadAttention(num_heads=4, key_dim=96)(x, x)
    x = Add()([x, x5])
    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2,1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x6 = MultiHeadAttention(num_heads=4, key_dim=96)(x, x)
    x = Add()([x, x6])
    x = Conv2D(1024, (4, 4), activation='relu', padding="same")(x)
    x = MaxPooling2D((2,1))(x)
    x = BatchNormalization()(x)
    x = Dropout(0.25)(x)

    x = Flatten()(x)
    x = Dense(1024, activation='relu')(x)
    x = Dense(512, activation='relu')(x)
    x = Dense(256, activation='relu')(x)
    outputs = Dense(1, activation='softplus')(x)

    model = Model(inputs, outputs)
    return model


input_shape = (2000, 4, 1)
model = build_model(input_shape)
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


epoch=80
batch_size=4
history = model.fit(X_train, y_train, epochs=epoch, batch_size=batch_size, validation_data=(X_val, y_val),
                    callbacks=[modelcheckpoint, earlystop],
                    verbose=1,
                    shuffle=True)


# 预测训练集和验证集
train_predictions = model.predict(X_train)
val_predictions = model.predict(X_val)

# 计算 R² 值
train_r2 = r2_score(y_train, train_predictions)
val_r2 = r2_score(y_val, val_predictions)

print(f"Training R²: {train_r2}")
print(f"Validation R²: {val_r2}")

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
plt.savefig('./figure/trans_1.png')
