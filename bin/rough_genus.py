import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import pandas as pd


def normalize(df):
    # 计算每行的和
    row_sums = df.sum(axis=1)
    # 检查是否存在和为 0 的行
    zero_sum_rows = row_sums[row_sums == 0].index
    if len(zero_sum_rows) > 0:
        #print(f"警告：以下行的和为 0，将被删除: {zero_sum_rows.tolist()}")
        df = df.drop(zero_sum_rows)  # 删除和为 0 的行
    # 归一化
    df_normalized = df.div(df.sum(axis=1), axis=0)

    return df_normalized


"""
    训练线性回归模型，用于校正 16S 数据。

    参数:- abd_16s: 16S 丰度数据 (DataFrame)
        - abd_wgs: WGS 丰度数据 (DataFrame)
        - pair_count: 样本数量
    返回:
    - models: 训练好的模型字典，键为物种名称，值为对应的线性回归模型
"""
def train_step(abd_16s, abd_wgs, pair_count):
    # 第一步：归一化 16S 和 WGS 数据
    abd_16s_norm = normalize(abd_16s)
    abd_wgs_norm = normalize(abd_wgs)

    # 初始化模型字典
    models = {}

    for species in abd_16s.columns:
        species_16s = abd_16s_norm[species].values
        species_wgs = abd_wgs_norm[species].values
        # 删除行和=0的行
        # 删除16S=0且wgs=0

        # 删除当前物种中 ref_sample == 0 且 target_sample == 0 的样本
        valid_indices = (species_16s != 0) | (species_wgs != 0)
        species_16s = species_16s[valid_indices]
        species_wgs = species_wgs[valid_indices]

        # 如果样本为空，则跳过
        if len(species_16s) == 0 or len(species_wgs) == 0:
            continue

        # 训练每个物种的模型
        X = species_16s.reshape(-1, 1)  # 16S 丰度作为自变量
        y = species_wgs  # WGS 丰度作为因变量
        model = LinearRegression()
        model.fit(X, y)
        models[species] = model  # 将训练好的模型存入字典

    return models



"""
    对单个样本的丰度数据进行校准。
    参数:- abd: 单个样本的丰度数据（一维数组或 Series）
    返回:- cal_abd: 校准后的丰度数据（一维数组）
"""
def calibrate(model_ols, abd):
    # 归一化
    abd_sum = abd.sum()
    if abd_sum <= 0:
        return abd  # 如果总和为 0，直接返回原始数据
    cal_abd = abd / abd_sum  # 归一化

    # 校准
    for species, value in cal_abd.items():
        if species not in model_ols or value == 0:
            continue  # 如果模型不存在或丰度为 0，跳过
        model = model_ols[species]  # 使用物种名称查找模型
        cal_abd[species] = max(model.predict([[value]])[0], 0)  # 预测并确保非负
        cal_abd[species] *= abd_sum  # 恢复原始总和

    # 再次归一化
    cal_abd_sum = cal_abd.sum()
    if cal_abd_sum > 0:
        cal_abd /= cal_abd_sum

    return cal_abd


"""
处理所有样本的丰度数据。
    参数:- input_abd: 输入丰度数据（DataFrame）
        - model_ols: 线性回归模型字典
        - sample_names: 样本名称列表
    返回:- output_table: 校正后的丰度数据（DataFrame）
"""
def process_samples(model_ols, input_abd):
    output_table = pd.DataFrame(columns=input_abd.columns)
    for i, row in input_abd.iterrows():  # 使用 iterrows 遍历 DataFrame
        sample_name = row.name  # 获取样本名称
        abd = row  # 获取当前样本的丰度数据
        cal_abd = calibrate(model_ols, abd)  # 校准
        output_table.loc[sample_name] = cal_abd  # 添加到输出表
    return output_table


"""
主函数，用于测试代码。
"""
def genus_cal(train_16s,train_wgs,other_16s,pair_count):

    """
    train_16s = pd.read_csv('E:/1shen/结果整理汇总/调整代码所需数据/train16S_15.csv', encoding='utf-8',
                            sep=',', header=0, index_col=0)
    train_wgs = pd.read_csv('E:/1shen/结果整理汇总/调整代码所需数据/trainWGS_15.csv', encoding='utf-8',
                            sep=',', header=0, index_col=0)
    other_16s = pd.read_csv('E:/1shen/结果整理汇总/调整代码所需数据/待校正16s.csv', encoding='utf-8',
                            sep=',', header=0, index_col=0)

    pair_count = 15
    """
    # 调用训练函数
    trained_models = train_step(train_16s, train_wgs, pair_count)

    # 调用校正函数
    df_calibrated = process_samples(trained_models, other_16s)
    #print("校正后的丰度数据：")
    #print(df_calibrated)
    # 计算每列的和，# 删除列和为 0 的列
    col_sum = df_calibrated.sum(axis=0)
    df_cleaned = df_calibrated.loc[:, col_sum != 0]
    return df_cleaned
    #df_cleaned.to_csv('E:/1shen/结果整理汇总/调整代码所需数据/校正后.csv', index=True)


# 运行主函数
