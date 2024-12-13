# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:59:20 2022

@author: Niko Tutland
"""

import pandas as pd
from pathlib import Path
import numpy as np
import pickle
import matplotlib.pyplot as plt
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, root_mean_squared_error
from skranger.ensemble import RangerForestRegressor


def main():
    ## User inputs ########
    fia_path = Path(__file__).parent / "FIA_raw"  # where are the raw FIA data located?
    out_path = (
        Path(__file__).parent / "RF_Models"
    )  # where do you want the model objects saved?
    plot_path = out_path / "Plots"  # where do you want the plots saved?
    plot_path.mkdir(exist_ok=True)
    fia_done = True  # has an FIA_all file already been made?
    ## End user inputs ####

    ## Use 11 states from western USFS Regions that measure tree age
    states = ["WA", "OR", "CA", "ID", "NV", "UT", "AZ", "NM", "CO", "WY", "MT"]

    ## Compile all FIA data from Rocky Mountain and Pacific Northwest research stations
    FIA_all = compile_fia(fia_done, fia_path, states, out_path)

    ## Fit random forest models for each major species group
    fit_rf(FIA_all, out_path, plot_path, group=1, prop=0.5)
    fit_rf(FIA_all, out_path, plot_path, group=2, prop=0.2)
    # fit_rf(FIA_all, out_path, plot_path, 3)
    # fit_rf(FIA_all, out_path, plot_path, 4)

    ##### End main function


def fit_rf(
    data: pd.DataFrame,
    out_path: Path,
    plot_path: Path,
    group: int,
    prop: int = 1,
) -> None:
    FIA_reg = data[data["MAJOR_SPGRPCD"] == group].sample(frac=prop)
    FIA_reg = FIA_reg[
        ["AGE", "dia", "ht", "cclcd", "physclcd", "sicond", "stdage", "balive"]
    ]  # select predictors and response
    FIA_reg["cclcd"] = FIA_reg["cclcd"].astype(
        "category"
    )  # make sure canopy class is a categorical variable
    FIA_reg["physclcd"] = FIA_reg["physclcd"].astype(
        "category"
    )  # make sure physiognomic class is a categorical variable

    # Get predictors and response
    x_all = FIA_reg.drop(["AGE"], axis="columns")  # predictors for whole dataset
    y_all = FIA_reg["AGE"]  # response for whole dataset

    print("Tuning hyperparameters...")
    tune_spec = RangerForestRegressor(  # specify the model used for tuning
        seed=47, n_estimators=500  # set seed  # use 500 trees
    )
    mtry = np.arange(
        1, 8, 1
    )  # how many variables (predictors) to try at each split? test between 1 and all (7)
    min_n = np.arange(
        2, 41, 1
    )  # what is the minimum number of samples considered at a split? test between 2 and 40
    tune_grid = {  # set up tuning grid for the above hyperparameters
        "mtry": mtry,  # mtry = max_features
        "min_node_size": min_n,  # min_n = min_samples_split
    }
    tune_res = RandomizedSearchCV(  # set up v-fold cross-validation
        estimator=tune_spec,
        param_distributions=tune_grid,  # use random combinations of values in the tuning grid
        n_iter=10,  # test 10 combinations
        cv=10,  # v = 10 (10-fold cross validation)
        random_state=4747,  # set seed
        n_jobs=8,  # use parallel processing
        scoring="neg_root_mean_squared_error",  # choose best parameters based on lowest RMSE
    )
    tune_res.fit(x_all, y_all)  # do the cross-validation
    best_rmse = tune_res.best_params_  # extract the best parameters
    print(f"Chosen parameters: {best_rmse}")

    print("Evaluating tuned model using cross-validation...")
    test_spec = RangerForestRegressor(  # specify a model to evaluate on the test set
        seed=541,  # set seeed
        mtry=best_rmse.get("mtry"),  # use the tuned mtry
        min_node_size=best_rmse.get("min_node_size"),  # use the tuned min_n
        n_estimators=500,  # use 500 trees
    )
    # Define the number of folds
    n_folds = 10
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=47)  # 10-fold CV
    rmse_list = []
    r2_list = []
    pred_dict = {"Fold": [], "y_test": [], "y_pred": []}
    fold = 0
    # Iterate through each fold
    for train_index, test_index in kf.split(x_all):
        fold += 1
        # Split the data into training and testing sets
        x_train, x_test = x_all.iloc[train_index], x_all.iloc[test_index]
        y_train, y_test = y_all.iloc[train_index], y_all.iloc[test_index]
        # Fit model
        test_spec.fit(x_train, y_train.values.ravel())
        # Predict on the test set
        y_pred = test_spec.predict(x_test)
        # Calculate RMSE and R2 for this fold
        r2_test = r2_score(y_test, y_pred)
        rmse_test = root_mean_squared_error(y_test, y_pred)
        # Store the RMSE
        rmse_list.append(rmse_test)
        r2_list.append(r2_test)
        # Store the predictions
        for test_val in y_test:
            pred_dict["y_test"].append(test_val)
        for pred_val in y_pred:
            pred_dict["y_pred"].append(pred_val)
        fold_list = [fold] * len(y_test)
        for fold_id in fold_list:
            pred_dict["Fold"].append(fold_id)

    print(
        f"Average R-Squared in {n_folds}-fold CV: {np.mean(r2_list)}\n"
        f"Average RMSE in {n_folds}-fold CV: {np.mean(rmse_list)}"
    )  # how did the tuned model do?

    # Save the CV predictions and metrics for plotting in R
    metric_df = pd.DataFrame(
        {"Fold": np.linspace(1, fold, fold), "RMSE": rmse_list, "R2": r2_list}
    )
    metric_df.to_csv(plot_path / f"cv_metrics_spgrp{group}.csv", index=False)
    pred_df = pd.DataFrame(pred_dict)
    pred_df.to_csv(plot_path / f"cv_predictions_spgrp{group}.csv", index=False)

    print("Fitting final model on all data...")
    final_spec = (
        RangerForestRegressor(  # specify a the final model to train on all the data
            seed=541,  # set seed
            mtry=best_rmse.get("mtry"),  # use the tuned mtry
            min_node_size=best_rmse.get("min_node_size"),  # use the tuned min_n
            n_estimators=500,  # use 500 trees
        )
    )
    final_spec.fit(x_all, y_all)

    # Pickle the model object
    rf_file = open(out_path / f"RF_model_spgrp{group}.obj", "wb")
    pickle.dump(final_spec, rf_file)
    rf_file.close()

    print("Saving model object...")
    # Pickle the model object
    rf_file = open(out_path / f"RF_model_spgrp{group}.obj", "wb")
    pickle.dump(test_spec, rf_file)
    rf_file.close()
    print(f"Model fitting for species group {group} complete.")


def compile_fia(
    FIA_done: bool, fia_path: Path, states: list, out_path: Path
) -> pd.DataFrame:
    if FIA_done == False:
        # Get relevant species info
        ref = pd.read_csv(fia_path / "REF_SPECIES.csv")
        ref = ref[["SPCD", "SPECIES_SYMBOL", "MAJOR_SPGRPCD"]]
        ref = ref.rename(columns={"SPCD": "spcd"})
        FIA_all = pd.DataFrame()
        # i = "NC" #for testing loop
        for i in states:
            print(f"Processing {i} FIA data...")
            tree = pd.read_csv(fia_path / f"{i}_TREE.csv")  # read in tree table
            tree = tree[
                [
                    "cn",
                    "plt_cn",
                    "condid",
                    "dia",
                    "ht",
                    "cclcd",
                    "spcd",
                    "totage",
                    "bhage",
                ]
            ]  # select relevant columns
            tree = tree[
                (tree["totage"] > 0) | (tree["bhage"] > 0)
            ]  # filter trees with age measurements
            tree = tree[tree["cn"].notna()]  # cn should not be na
            tree = tree[tree["cn"] > 0]  # or zero
            tree = tree[tree["dia"].notna()]  # diameter should not be na
            tree = tree[tree["dia"] > 0]  # or zero
            tree = tree[tree["ht"].notna()]  # height should not be na
            tree = tree[tree["ht"] > 0]  # or zero
            tree = tree[tree["cclcd"].notna()]  # canopy class should not be na
            tree = tree[tree["cclcd"] > 0]  # or zero
            tree.loc[tree["totage"] > 0, "AGE"] = tree[
                "totage"
            ]  # if totage is populated, use it
            tree.loc[tree["totage"] == 0, "AGE"] = (
                tree["bhage"] + 10
            )  # if not, use bhage+10
            cond = pd.read_csv(fia_path / f"{i}_COND.csv")  # read in condtion table
            cond = cond[["plt_cn", "condid", "physclcd", "sicond", "stdage", "balive"]]
            treecond = tree.merge(cond, how="left", on=["plt_cn", "condid"])
            treecond["sicond"][treecond["sicond"] == 0] = treecond["sicond"][
                treecond["sicond"] != 0
            ].mean()
            treecond["stdage"][treecond["stdage"] == 0] = treecond["stdage"][
                treecond["stdage"] != 0
            ].mean()
            treecond["sicond"][treecond["sicond"].isna()] = treecond["sicond"][
                treecond["sicond"].notna()
            ].mean()
            treecond["stdage"][treecond["stdage"].isna()] = treecond["stdage"][
                treecond["stdage"].notna()
            ].mean()
            treecond = treecond[treecond["balive"].notna()]  # live ba should not be na
            treecond = treecond[treecond["balive"] > 0]  # or zero
            treecond = treecond[
                treecond["physclcd"].notna()
            ]  # physiognomic class should not be na
            treecond = treecond[treecond["physclcd"] > 0]  # or zero
            ## Append each state
            FIA_all = pd.concat([FIA_all, treecond])
        FIA_all = FIA_all.merge(ref, how="left", on="spcd")
        FIA_all.to_csv(out_path / "FIA_all_agereg.csv", index=False)
    else:
        FIA_all = pd.read_csv(out_path / "FIA_all_agereg.csv")
    return FIA_all


if __name__ == "__main__":
    main()
