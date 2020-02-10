#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:11:27 2020

@author: peter
"""

from tslearn.utils import *
from tslearn.datasets import UCR_UEA_datasets
from tslearn.clustering import TimeSeriesKMeans

my_first_time_series = [1, 3, 4, 2]
formatted_time_series = to_time_series(my_first_time_series)
print(formatted_time_series.shape)

my_first_time_series = [1, 3, 4, 2]
my_second_time_series = [1, 2, 4, 2]
formatted_dataset = to_time_series_dataset([my_first_time_series, my_second_time_series])
print(formatted_dataset.shape)

my_third_time_series = [1, 2, 4, 2, 2]
formatted_dataset = to_time_series_dataset([my_first_time_series,
                                            my_second_time_series,
                                            my_third_time_series])
print(formatted_dataset.shape)

X_train, y_train, X_test, y_test = UCR_UEA_datasets(
        ).load_dataset("TwoPatterns")
print(X_train.shape)
print(y_train.shape)

#time_series_dataset = load_time_series_txt("path/to/your/file.txt")
#save_time_series_txt("path/to/another/file.txt", dataset_to_be_saved)

km = TimeSeriesKMeans(n_clusters=3, metric="dtw")
km.fit(X_train)