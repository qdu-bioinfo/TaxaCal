import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.model_selection import KFold

from scipy.spatial.distance import braycurtis


def bray_curtis_distance(x1, x2):
    distance = np.sum(np.abs(x1 - x2)) / np.sum(np.abs(x1 + x2))
    return distance


def find_k_nearest_neighbors(train_data, test_data, k):
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='brute', metric=bray_curtis_distance)
    nbrs.fit(train_data)
    distances, indices = nbrs.kneighbors(test_data)
    return distances, indices


def get_knn_genus(train_data, test_data, k):
    distances, indices = find_k_nearest_neighbors(train_data, test_data, k)
    return distances, indices


def get_species_by_genus(genus_neighbors_data, species_data):
    '''
    Get the species corresponding to K genera
    '''

    selected_rows = species_data[species_data.index.isin(genus_neighbors_data.index)]
    # Filter out the columns whose species table names contain the genus table name
    selected_columns = [col for col in species_data.columns if
                        any(col.startswith(prefix) for prefix in genus_neighbors_data.columns)]
    return selected_rows[selected_columns]


def same_rate_cal(standard, uncal_genus):
    '''
    Same ratio correction
    '''
    df_standard_genus = standard.T.groupby(standard.columns.str.rsplit(';', n=1).str[0]).sum().T
    result_df = pd.DataFrame(index=uncal_genus.index)  # save result

    for genus in uncal_genus.columns:
        sums_genus = uncal_genus[genus]
        # check if df_standard_genus contains a genus column
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

    new_16s = df_16s.reindex(columns=all_columns, fill_value=0)
    new_wgs = df_wgs.reindex(columns=all_columns, fill_value=0)

    # Calculate the Bray-Curtis distance for each sample
    bray_curtis_distances = pd.DataFrame()
    index_list = new_16s.index.tolist()
    for i in range(len(index_list)):
        # If the data type is inconsistent, try to convert it to a numeric type
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
    Use K-fold cross validation to evaluate the performance of the kNN model for a given k value.

    Parameters:- train_data: genus-level abundance matrix of the training data
                - train_labels: species-level abundance matrix of the training data
                - k: the k value currently evaluated
                - cv: the number of folds for cross validation
    Return:
        - Mean Squared Error (MSE)

    """
    kf = KFold(n_splits=cv, shuffle=True, random_state=42)
    mse_scores = []
    df_distance_score = pd.DataFrame();
    for train_index, val_index in kf.split(train_data):
        X_train, X_val = train_data.iloc[train_index], train_data.iloc[val_index]
        y_train, y_val = train_labels.iloc[train_index], train_labels.iloc[val_index]

        result = KNN(X_train, y_train, X_val, k)

        distance_16s_wgs = distance_same_samples(result,y_val)
        df_distance_score = pd.concat([df_distance_score, distance_16s_wgs])

    # 新的列名列表
    new_columns = ["k=" + str(k)]
    df_distance_score.columns = new_columns
    return df_distance_score

def KNN(df_wgs_genus,df_wgs_species,df_16s,k):
    distances, indices = get_knn_genus(df_wgs_genus, df_16s, k=k)
    cal_result = pd.DataFrame()
    for i in range(len(df_16s)):
        genus_neighbors_data = df_wgs_genus.iloc[indices[i]]

        species_neighbors_data = get_species_by_genus(genus_neighbors_data, df_wgs_species)

        row_mean = species_neighbors_data.mean(axis=0)
        df_average = pd.DataFrame([row_mean], index=['average'])

        result_df = same_rate_cal(df_average, pd.DataFrame([df_16s.iloc[i]]))
        cal_result = pd.concat([cal_result, result_df])

    return cal_result


def species_cal(df_wgs_genus,df_wgs_species,df_test_16s):
    df_train = df_wgs_genus.reindex(columns=df_test_16s.columns, fill_value=0)
    k_values = range(1, 11)  # k from 1 to 10
    cv = 5
    mse_results = pd.DataFrame();
    species_data = df_wgs_species # Species-level abundance
    # Evaluate each value of k
    for k in k_values:
        mse = evaluate_k(df_train, species_data, k, cv=cv)
        mse_results = pd.concat([mse_results, mse], axis=1)
        #print(f"k={k}, Cross-Validation MSE={mse:.4f}")
    mean_avg = mse_results.mean()
    #print(mean_avg)
    best_k = int(mean_avg.idxmin().split('=')[-1])
    best_distance = mean_avg.min()
    print(f"The best K is: {best_k}")
    df_cal_result = KNN(df_train,species_data,df_test_16s,best_k)
    return df_cal_result

