from config import parse_args
import os
import pandas as pd
from rough_genus import genus_cal
from refined_species import species_cal

def pre_data(args):
    train_WGS_genus = pd.read_csv(args.train_WGS_g, header=0, index_col=0)
    train_WGS_species = pd.read_csv(args.train_WGS_s, header=0, index_col=0)
    train_16S_genus = pd.read_csv(args.train_16S_g, header=0, index_col=0)
    test_genus = pd.read_csv(args.test, header=0, index_col=0)
    #output = pd.read_csv(args.o)
    output_path = args.o

    def somcol(df_16s,df_wgs,df_test):
        all_columns = list(set(df_16s.columns) | set(df_wgs.columns) | set(df_test.columns))
        new_16s = df_16s.reindex(columns=all_columns, fill_value=0)
        new_wgs = df_wgs.reindex(columns=all_columns, fill_value=0)
        new_test = df_test.reindex(columns=all_columns, fill_value=0)
        return new_16s,new_wgs,new_test

    new_train_16S,new_train_wgs,new_test = somcol(train_16S_genus,train_WGS_genus,test_genus)
    pair_count = len(new_train_16S)
    print("开始校正")
    df_genus_cal = genus_cal(new_train_16S,new_train_wgs,new_test,pair_count)
    df_species_cal = species_cal(new_train_wgs, train_WGS_species, df_genus_cal)
    df_species_cal.to_csv(output_path, index=True)
    print("校正完成")
    print("文件输出到", args.o)


def main():
    # 解析命令行参数
    args = parse_args()  # 这行代码会调用 parse_args 函数来解析命令行传入的参数

    # 检查输出目录是否存在
    if args.o and isinstance(args.o, str):
        if not os.path.exists(os.path.dirname(args.o)):
            os.makedirs(os.path.dirname(args.o), exist_ok=True)
    else:
        print("错误：输出路径无效")

    #if not os.path.exists(args.o):
        # 如果目录不存在，则创建该目录
        #os.makedirs(args.o)

    # 执行任务
    pre_data(args)

if __name__ == '__main__':
    main()
