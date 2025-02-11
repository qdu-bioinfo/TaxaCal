import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.model_selection import KFold

from scipy.spatial.distance import braycurtis

'''
先对属进行KNN 挑选 选出WGS种K个属对对应的种，将该值作为标准，之后通过该值 同比例校正16S
'''


def bray_curtis_distance(x1, x2):
    """
    计算 Bray-Curtis 距离
    """
    distance = np.sum(np.abs(x1 - x2)) / np.sum(np.abs(x1 + x2))
    return distance


def find_k_nearest_neighbors(train_data, test_data, k):
    """
    使用KNN算法找到测试数据中每个样本在训练数据中的k个最近邻居
    """
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='brute', metric=bray_curtis_distance)
    nbrs.fit(train_data)
    distances, indices = nbrs.kneighbors(test_data)
    return distances, indices


def get_knn_genus(train_data, test_data, k):
    '''
    获取WGS中临近的K个属
    '''
    distances, indices = find_k_nearest_neighbors(train_data, test_data, k)
    return distances, indices


def get_species_by_genus(genus_neighbors_data, species_data):
    '''
    获取K个属所对应的种
    '''
    # 在种层级中筛选出行名和K临近的行名一样的行
    selected_rows = species_data[species_data.index.isin(genus_neighbors_data.index)]
    # 筛选出 species 表列名中包含 genus 表列名的列
    selected_columns = [col for col in species_data.columns if
                        any(col.startswith(prefix) for prefix in genus_neighbors_data.columns)]
    return selected_rows[selected_columns]


def same_rate_cal(standard, uncal_genus):
    '''
    同比例校正
    '''
    # 通过species倒推到genus。
    df_standard_genus = standard.T.groupby(standard.columns.str.rsplit(';', n=1).str[0]).sum().T
    result_df = pd.DataFrame(index=uncal_genus.index)  # 存放结果
    # 遍历属 对每个属进行同比例的校正
    for genus in uncal_genus.columns:
        sums_genus = uncal_genus[genus]
        # 首先，检查 df_standard_genus 是否包含 genus 列
        if genus in df_standard_genus.columns:
            sums_species = df_standard_genus[genus]
        else:
            sums_species = pd.Series(0, index=df_standard_genus.index)

        list_species = [item for item in standard.columns if item.startswith(genus)]
        for species_name in list_species:
            if sums_species.sum() != 0:
                result = (standard[species_name].mean() / sums_species.sum() * sums_genus.sum())
            else:
                result = 0

            row = str(uncal_genus.index[0])
            temp_df = pd.DataFrame(data=[result], index=[row], columns=[species_name])
            result_df = pd.concat([result_df, temp_df], axis=1)
    return result_df


def distance_same_samples(df_16s, df_wgs):
    all_columns = list(set(df_16s.columns) | set(df_wgs.columns))
    # 创建新的表1和表2，确保包含相同的列和相同的顺序
    new_16s = df_16s.reindex(columns=all_columns, fill_value=0)
    new_wgs = df_wgs.reindex(columns=all_columns, fill_value=0)
    # 使用rename方法将行名中的"_16S"替换为"_out",inplace=True在原始DataFrame上进行修改，而不是创建一个新的DataFrame。
    #new_16s.rename(index=lambda x: x.replace("_16S", "_output"), inplace=True)  # CRC
    # new_16s.rename(index=lambda x: x.replace("_16S", ""), inplace=True)# AGP

    # 计算每个样本的 Bray-Curtis 距离
    bray_curtis_distances = pd.DataFrame()
    index_list = new_16s.index.tolist()
    for i in range(len(index_list)):
        # 如果数据类型不一致，尝试将其转换为数值类型
        if new_16s.iloc[i].dtype == object:
            new_16s.iloc[i] = pd.to_numeric(new_16s.iloc[i], errors='coerce')
        if new_wgs.iloc[i].dtype == object:
            new_wgs.iloc[i] = pd.to_numeric(new_wgs.iloc[i], errors='coerce')
        distance = braycurtis(new_16s.iloc[i], new_wgs.iloc[i])
        df_distance = pd.DataFrame([[distance]], index=[new_wgs.index[i]])
        bray_curtis_distances = pd.concat([bray_curtis_distances, df_distance])
    return bray_curtis_distances


def evaluate_k(train_data, train_labels, k, cv=5):
    """
    使用K折交叉验证评估给定k值的kNN模型性能。

    参数：
    - train_data: 训练数据的属层级丰度矩阵
    - train_labels: 训练数据的种层级丰度矩阵
    - k: 当前评估的k值
    - cv: 交叉验证的折数

    返回：
    - 平均均方误差（MSE）
    """
    kf = KFold(n_splits=cv, shuffle=True, random_state=42)
    mse_scores = []
    df_distance_score = pd.DataFrame();
    for train_index, val_index in kf.split(train_data):
        X_train, X_val = train_data.iloc[train_index], train_data.iloc[val_index]
        y_train, y_val = train_labels.iloc[train_index], train_labels.iloc[val_index]

        #训练KNN模型
        result = KNN(X_train, y_train, X_val, k)

        # 计算距离；
        distance_16s_wgs = distance_same_samples(result,y_val)
        df_distance_score = pd.concat([df_distance_score, distance_16s_wgs])

    # 新的列名列表
    new_columns = ["k=" + str(k)]
    df_distance_score.columns = new_columns
    return df_distance_score

def KNN(df_wgs_genus,df_wgs_species,df_16s,k):
    # 获取WGS中临近的K个属
    distances, indices = get_knn_genus(df_wgs_genus, df_16s, k=k)
    cal_result = pd.DataFrame()
    for i in range(len(df_16s)):
        #print(f"样本 {i + 1} 的最近邻居索引：{indices[i]}")
        genus_neighbors_data = df_wgs_genus.iloc[indices[i]]  # k个临近属
        # 通过K个属找到对应的种
        species_neighbors_data = get_species_by_genus(genus_neighbors_data, df_wgs_species)
        # 求species 层级的k个邻居的平均值
        row_mean = species_neighbors_data.mean(axis=0)
        df_average = pd.DataFrame([row_mean], index=['average'])
        # 将上述平均值作为标准，同比例校正
        result_df = same_rate_cal(df_average, pd.DataFrame([df_16s.iloc[i]]))
        cal_result = pd.concat([cal_result, result_df])

    return cal_result


def species_cal(df_wgs_genus,df_wgs_species,df_test_16s):

    # 确保训练集和测试集具有相同的属
    df_train = df_wgs_genus.reindex(columns=df_test_16s.columns, fill_value=0)

    # 定义k值的范围
    k_values = range(1, 11)  # k从1到10
    cv = 5  # 5折交叉验证
    mse_results = pd.DataFrame();

    # 定义种层级丰度矩阵（假设为df_wgs）
    species_data = df_wgs_species

    # 评估每个k值
    for k in k_values:
        mse = evaluate_k(df_train, species_data, k, cv=cv)
        mse_results = pd.concat([mse_results, mse], axis=1)
        #print(f"k={k}, Cross-Validation MSE={mse:.4f}")

    # 计算每个k值的平均Bray-Curtis距离
    mean_avg = mse_results.mean()
    # 显示平均距离
    print(mean_avg)
    # 找到最小平均距离对应的k值
    best_k = int(mean_avg.idxmin().split('=')[-1])
    best_distance = mean_avg.min()
    print(f"最佳k值为: {best_k}")

    # 最佳k值下的矫正结果
    df_cal_result = KNN(df_train,species_data,df_test_16s,best_k)
    return df_cal_result

