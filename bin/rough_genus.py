import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import pandas as pd


def normalize(df):
    row_sums = df.sum(axis=1)
    zero_sum_rows = row_sums[row_sums == 0].index
    if len(zero_sum_rows) > 0:
        df = df.drop(zero_sum_rows)  # Delete rows with sum equal to 0
    # Normalization
    df_normalized = df.div(df.sum(axis=1), axis=0)
    return df_normalized


def train_step(abd_16s, abd_wgs, pair_count):
    """
        Train a linear regression model to correct 16S data.

        Parameters: - abd_16s: 16S abundance data (DataFrame)
                    - abd_wgs: WGS abundance data (DataFrame)
                    - pair_count: number of samples
        Returns:
                - models: trained model dictionary, with the key being the species name and the value being the corresponding linear regression model
    """

    abd_16s_norm = normalize(abd_16s)
    abd_wgs_norm = normalize(abd_wgs)
    models = {}

    for species in abd_16s.columns:
        species_16s = abd_16s_norm[species].values
        species_wgs = abd_wgs_norm[species].values

        # Delete samples with ref_sample == 0 and target_sample == 0 in the current species
        valid_indices = (species_16s != 0) | (species_wgs != 0)
        species_16s = species_16s[valid_indices]
        species_wgs = species_wgs[valid_indices]

        if len(species_16s) == 0 or len(species_wgs) == 0:
            continue

        X = species_16s.reshape(-1, 1)
        y = species_wgs
        model = LinearRegression()
        model.fit(X, y)
        models[species] = model

    return models



def calibrate(model_ols, abd):
    """
        Calibrate the abundance data of a single sample.
        Parameters:- abd: Abundance data of a single sample
        Returns:- cal_abd: Calibrated abundance data
    """

    abd_sum = abd.sum()
    if abd_sum <= 0:
        return abd
    cal_abd = abd / abd_sum

    for species, value in cal_abd.items():
        if species not in model_ols or value == 0:
            continue
        model = model_ols[species]
        cal_abd[species] = max(model.predict([[value]])[0], 0)
        cal_abd[species] *= abd_sum

    cal_abd_sum = cal_abd.sum()
    if cal_abd_sum > 0:
        cal_abd /= cal_abd_sum

    return cal_abd


def process_samples(model_ols, input_abd):
    output_table = pd.DataFrame(columns=input_abd.columns)
    for i, row in input_abd.iterrows():
        sample_name = row.name
        abd = row
        cal_abd = calibrate(model_ols, abd)
        output_table.loc[sample_name] = cal_abd
    return output_table



#main
def genus_cal(train_16s,train_wgs,other_16s,pair_count):
    # train
    trained_models = train_step(train_16s, train_wgs, pair_count)
    # cal
    df_calibrated = process_samples(trained_models, other_16s)

    col_sum = df_calibrated.sum(axis=0)
    df_cleaned = df_calibrated.loc[:, col_sum != 0]
    return df_cleaned

