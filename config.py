# -*- coding: utf8 -*-
from easydict import EasyDict as edict
import json
import numpy as np
import os

configPara = edict()

configPara.if_continue_train=True
configPara.test=True

configPara.nn1=128
configPara.nn2=128
configPara.rate=4

configPara.lr_init = 2e-4
configPara.beta1 = 0.5
configPara.n_epoch =100

configPara.L1_lambda=1
configPara.tvDiff_lambda=200
configPara.loop_lambda=1
configPara.gan_lambda=0.01

configPara.test_freq=20
configPara.save_model_freq=10

configPara.if_aug = True

if configPara.if_aug:
    configPara.scale = 5
else:
    configPara.scale = 0
    
configPara.if_scale = True

nn1=configPara.nn1
nn2=configPara.nn2
scale=configPara.scale
tvDiff_lambda=configPara.tvDiff_lambda
configPara.type = 4
TYPE = configPara.type          # 0 : 训练 ，  1: 测试plut   ， 2：测试sigsbee  3: 测试plut 10  ， 4：测试sigsbee20

configPara.samples_save_dir = "train_result/samples/"
if TYPE == 0:
    configPara.test_save_dir="train_result/test_result"           # 训练用
elif TYPE == 1:
    configPara.test_save_dir="train_result/test_result_snr"       # 测试pluto用
elif TYPE == 2:
    configPara.test_save_dir="train_result/test_result_sigsbee"       # 测试sigsbee
else:
    configPara.test_save_dir="train_result/test_result_snr_10"
configPara.checkpoint_dir = "train_result/checkpoint/"
configPara.buffer_dir="train_result/buffer/"

configPara.train_image_path="../01-2D-cnn-denoise-liu/train_Data/"

if TYPE == 0:
    configPara.test_image_path="../01-2D-cnn-denoise-liu/test_data_txt/input/"           # 训练用
    configPara.test_image_path_target = "../01-2D-cnn-denoise-liu/test_data_txt/target/"
elif TYPE == 1:
    configPara.test_image_path="../01-2D-cnn-denoise-liu/test_data_txt_snr/input/"       # 测试pluto用
    configPara.test_image_path_target = "../01-2D-cnn-denoise-liu/test_data_txt_snr/target/"
elif TYPE == 2:
    configPara.test_image_path="../01-2D-cnn-denoise-liu/test_data_txt_sigsbee/input/"       # 测试sigsbee
    configPara.test_image_path_target = "../01-2D-cnn-denoise-liu/test_data_txt_sigsbee/target/"
elif TYPE == 3:
    configPara.test_image_path = "../01-2D-cnn-denoise-liu/test_data_txt_snr_10/input/"  # 测试pluto16用
    configPara.test_image_path_target = "../01-2D-cnn-denoise-liu/test_data_txt_snr_10/target/"



configPara.test_image_path_20="../01-2D-cnn-denoise-liu/test_data_txt_sigsbee_20/"   #测试20张sigsbee

if not os.path.exists(configPara.checkpoint_dir):
    os.makedirs(configPara.checkpoint_dir)
if not os.path.exists(configPara.test_save_dir):
    os.makedirs(configPara.test_save_dir)
if not os.path.exists(configPara.samples_save_dir):
    os.makedirs(configPara.samples_save_dir)
if not os.path.exists(configPara.buffer_dir):
    os.makedirs(configPara.buffer_dir)
