'''

Scenario reduction algorithm

Author: Jonghwan Kwon (kwonj@anl.gov)

Note:
- main code is based on the one in https://gitlab.com/supsi-dacd-isaac/scenred

'''


import numpy as np
from scipy.spatial.distance import pdist,squareform,cdist
import matplotlib
import sys
import matplotlib.pyplot as plt
import pandas as pd
import os



def run_scenario_reduction(**kwargs):

    print ("Start ALEAF scenario reduction algorithm!")

    # set default parameters and pars kwargs
    print ("==== current setting ========================")
    setting = {'time_resolution': "Hourly",
            'num_scenarios_list': [60],
            'fixing_extreme_days_flag': True,
            'generate_input_data_flag': True,
            'generate_duration_curve_flag': False,
            'num_data_set': 3,
            'type_of_data_set': ["load_shape", "wind_shape", "solar_shape"],
            'data_location_timeseries': 'C:/ALEAF/ALEAF/aleaf/data/LC_GEP/MISO',
            'data_input_path': 'Data_Input',
            'data_output_path': 'Data_Output',
            'plot_output_path': 'Plot_Output',
            'windCapacity': 1.0,
            'solarCapacity': 1.0,
            'peakDemand': 1.0
            }

    for key in setting.keys():
        if key in kwargs:
            setting[key] = kwargs[key]
        print (str(key + "\t:  " + str(setting[key])).expandtabs(45))
    print ("=============================================")

    # Check input output folders
    if not os.path.exists(setting["data_input_path"]): os.makedirs(setting["data_input_path"])
    if not os.path.exists(setting["data_output_path"]): os.makedirs(setting["data_output_path"])
    if not os.path.exists(setting["plot_output_path"]): os.makedirs(setting["plot_output_path"])

    # Generate input data for the scenario reduction algorithm
    if setting["generate_input_data_flag"] == True:
        print ("Processing input data files")
        if setting["time_resolution"] == "Hourly":
            generate_input_data_hourly(setting["data_input_path"], setting["data_location_timeseries"], setting["windCapacity"], setting["solarCapacity"], setting["peakDemand"])
        elif setting["time_resolution"] == "Five-min":
            generate_input_data_5min(setting["data_input_path"], setting["data_location_timeseries"], setting["windCapacity"], setting["solarCapacity"], setting["peakDemand"])

    # Read input data for the scenario reduction algorithm
    if setting["time_resolution"] == "Hourly":
        load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh = read_input_data_Hourly(setting["data_input_path"])
    elif setting["time_resolution"] == "Five-min":
        load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh = read_input_data_5min(setting["data_input_path"])

    # Identify extreme points
    extreme_scenarios, extreme_datapoint = identify_extreme_points(load_shape, wind_shape, solar_shape, load_MWh,
                                                                   wind_MWh, solar_MWh, net_load_MWh)

    # Aggregate data points
    num_timestep = np.shape(load_shape.iloc[:-1, 1:].to_numpy())[0]
    num_total_scenario = np.shape(load_shape.iloc[:-1, 1:].to_numpy())[1]
    data = np.zeros((num_timestep, num_total_scenario, setting["num_data_set"]))

    flag_load_shape = flag_wind_shape = flag_solar_shape = flag_load_MWh = flag_wind_MWh = flag_solar_MWh = flag_net_load_MWh = False
    for idx in range(0, setting["num_data_set"]):
        if ("load_shape" in setting["type_of_data_set"]) & (flag_load_shape == False):
            data[:, :, idx] = load_shape.iloc[:-1, 1:].to_numpy()
            flag_load_shape = True
        elif ("wind_shape" in setting["type_of_data_set"]) & (flag_wind_shape == False):
            data[:, :, idx] = wind_shape.iloc[:-1, 1:].to_numpy()
            flag_wind_shape = True
        elif ("solar_shape" in setting["type_of_data_set"]) & (flag_solar_shape == False):
            data[:, :, idx] = solar_shape.iloc[:-1, 1:].to_numpy()
            flag_solar_shape = True
        elif ("load_MWh" in setting["type_of_data_set"]) & (flag_load_MWh == False):
            data[:, :, idx] = load_MWh.iloc[:-1, 1:].to_numpy()
            flag_load_MWh = True
        elif ("wind_MWh" in setting["type_of_data_set"]) & (flag_wind_MWh == False):
            data[:, :, idx] = wind_MWh.iloc[:-1, 1:].to_numpy()
            flag_wind_MWh = True
        elif ("solar_MWh" in setting["type_of_data_set"]) & (flag_solar_MWh == False):
            data[:, :, idx] = solar_MWh.iloc[:-1, 1:].to_numpy()
            flag_solar_MWh = True
        elif ("net_load_MWh" in setting["type_of_data_set"]) & (flag_net_load_MWh == False):
            data[:, :, idx] = net_load_MWh.iloc[:-1, 1:].to_numpy()
            flag_net_load_MWh = True


    print("Total number of cases to run: %d" % (len(setting["num_scenarios_list"])))
    idx = 1

    for list in setting["num_scenarios_list"]:

        num_scenarios = list
        print("(case %d) running a case with %d number of representative days" % (idx, num_scenarios))
        idx += 1
        # print("setting:  ", setting["fixing_extreme_days_flag"])
        if setting["fixing_extreme_days_flag"] == True:
            if num_scenarios < len(extreme_scenarios):
                print(
                    'The given number of scenarios to select is less than the number of extreme cases. Forcing it to the number of extreme cases')
                num_scenarios = len(extreme_scenarios)

        # Prepare fixed scenarios set
        extreme_set = [load_shape.columns.get_loc(c) - 1 for c in extreme_scenarios]
        extreme_set.sort()

        # Call scenarios reduction algorithm
        [selected_scenarios, selected_scenarios_prob, scenarios_tree_structure] = scenario_reduction_core(
            setting["fixing_extreme_days_flag"], extreme_set,
            np.copy(data),
            nodes=np.linspace(
                num_scenarios,
                num_scenarios,
                num_timestep,
                dtype=int))

        # Identify selected scenarios
        s_prob = selected_scenarios_prob[0, :]
        sce_reduction_result = scenarios_tree_structure[0, :]
        selected_sce = []
        all_sce_prob = []
        selected_sce_prob = []
        rep_day_input_data = []
        # rep_day_input_data.append("index,Day,Probability")
        prob_idx = 0
        for i in np.arange(np.shape(scenarios_tree_structure)[1]):
            if sce_reduction_result[i] == True:
                selected_sce.append(i + 1)
                all_sce_prob.append(s_prob[prob_idx])
                selected_sce_prob.append(s_prob[prob_idx])
                rep_day_input_data.append([prob_idx + 1, i + 1, s_prob[prob_idx]])
                prob_idx += 1
            else:
                all_sce_prob.append(0)

        rep_day_input_df = pd.DataFrame(rep_day_input_data, columns=['index', 'Day', 'Probability'])

        # pd.DataFrame(all_sce_prob).T.to_csv('Data_Outputs' + '/output_scenarios_prob_' + str(num_scenarios) + '.csv', index=False)
        # pd.DataFrame(selected_sce_prob).T.to_csv(
        #     setting["data_output_path"] + "/" + 'output_selected_scenarios_prob_' + str(num_scenarios) + '.csv',
        #     index=False)
        # pd.DataFrame(selected_sce).T.to_csv(
        #     setting["data_output_path"] + "/" +"output_selected_scenarios_" + str(num_scenarios) + '.csv', index=False)
        rep_day_input_df.to_csv(setting["data_output_path"] + "/" + "repDays_" + str(num_scenarios) + '.csv', index=False)

        # Plot duration curves
        if setting["generate_duration_curve_flag"] == True:
            plot_duration_curve(wind_shape, selected_sce_prob, selected_sce, num_scenarios, "Wind Shape (%)", setting["plot_output_path"])
            plot_duration_curve(solar_shape, selected_sce_prob, selected_sce, num_scenarios, "Solar Shape (%)", setting["plot_output_path"])
            plot_duration_curve(wind_MWh, selected_sce_prob, selected_sce, num_scenarios, "Wind Output (MWh)", setting["plot_output_path"])
            plot_duration_curve(solar_MWh, selected_sce_prob, selected_sce, num_scenarios, "Solar Output (MWh)", setting["plot_output_path"])
            plot_duration_curve(net_load_MWh, selected_sce_prob, selected_sce, num_scenarios, "Net Load (MWh)", setting["plot_output_path"])
            plot_duration_curve(load_MWh, selected_sce_prob, selected_sce, num_scenarios, "Load (MWh)", setting["plot_output_path"])

            # Plot ramp duration curves
            plot_ramp_duration_curve(net_load_MWh, selected_sce_prob, selected_sce, num_scenarios, "Net Load Ramp (MWh)", setting["plot_output_path"])

    print ("==== DONE ! ================================")


def get_dist(X,metric):
    D = squareform(pdist(X,metric))
    return D


def scenario_reduction_core(fixing_extreme_days, extreme_set, samples, **kwargs):

    '''
    This revised code is based on the one in https://gitlab.com/supsi-dacd-isaac/scenred
    '''

    T = samples.shape[0]
    n_obs = samples.shape[1]

    defaultNodes = np.ones((T,1))
    # pars kwargs
    pars = {'nodes':defaultNodes,
            'tol': 10,
            'metric':'cityblock'}
    for key in ('nodes','tol','metric'):
        if key in kwargs:
            pars[key] = kwargs[key]

    # Obtain the observation matrix, size (n*T) x n_obs, from which to compute distances. Normalize observations
    X = []
    S = samples
    for i in np.arange(np.shape(S)[2]):
        V = S[:,:,i]
        V_norm = (V-np.mean(V,1).reshape(-1,1))/(np.std(V,1)+1e-6).reshape(-1,1)
        X.append(V_norm)

    X=np.vstack(X).T

    D = get_dist(X, pars['metric'])
    D = D + np.eye(D.shape[0]) * (1 + np.max(D.ravel()))
    infty = 1e12

    D[1, :] = infty
    D[:, 1] = infty

    # generate the tolerance vector
    if all(pars['nodes'].ravel() == defaultNodes.ravel()):
        #Tol = np.tile(pars['tol'], (1, T))[0]
        Tol = np.fliplr(pars['tol'] / (1.5 ** (T - np.arange(T).reshape(1,-1)+1))).ravel()
        Tol[0] = infty
    else:
        Tol = infty * np.ones((1, T)).ravel()

    J = np.asanyarray(np.ones((T,n_obs)),bool)
    L = np.zeros((n_obs,n_obs))
    P = np.ones((T, n_obs)) / n_obs
    branches = n_obs
    for i in np.fliplr(np.arange(T).reshape(1,-1))[0]:
        delta_rel = 0
        delta_p = 0
        D_i = D
        delta_rel_2 = 0
        delta_p_2 = 0

        basic_idx = np.asanyarray(np.hstack([np.ones((1, i+1)), np.zeros((1, T - i-1))]),bool)
        sel_idx = np.tile(basic_idx, (1, samples.shape[2])).ravel()
        X_filt = X[J[i,:], :]
        X_filt = X_filt[:, sel_idx.ravel()]
        D_j = get_dist(X_filt, pars['metric'])
        D_j[np.asanyarray(np.eye(np.shape(D_j)[0]),bool)] = 0 # do not consider self - distance
        delta_max = np.min(np.sum(D_j,0))

        while (delta_rel < Tol[i]) and (branches > pars['nodes'][i]):
            D_i[~J[i,:],:] =  infty  # set distance of discarded scenarios to infinity, ignoring them
            D_i[:, ~J[i,:]] = infty
            d_s = np.sort(D_i,0)  # sort distances with respect to the non-discarded scenarios
            idx_s = np.argsort(D_i,0)
            z = d_s[0,:]*P[i,:]  # vector of weighted probabilities
            #z[z == 0] = infty  # set prob. of removed scenario to inf in order to ignore them
            z[~J[i,:]] = infty  # set prob. of removed scenario to inf in order to ignore them
            if fixing_extreme_days == True:
                z[np.array(extreme_set)] = infty    # prevent scenarios in the extreme_set being selected
            idx_rem = np.argmin(z)  # find the scenario which cause the smallest p*d-deviation when merged, and its index
            dp_min = np.min(z)

            idx_aug = idx_s[0, idx_rem]  # retrieve who's being augmented with the probability of idx_rem
            J[i, idx_rem] = False  # mark it as a removed scenarion in the current timestep
            P[i, idx_aug] = P[i, idx_rem] + P[i,idx_aug] # add the probability of the removed scenario to the closest scenario
            P[i, idx_rem] = 0  # set probability of removed scenarios to 0
            branches = np.sum(P[i,:] > 0)  # count remaining branches
            L[idx_aug, idx_rem] = 1  # keep track of which scenario has merged
            L[idx_aug, L[idx_rem,:] > 0] = 1  # the scenario who's been augmented heredit all the scenarios previously merged with the removed scenario
            S[0: i+1, idx_rem,:] =  S[0: i+1, idx_aug,:]  # make the merged scenarios equal up to the root node
            to_merge_idx = np.argwhere(L[[idx_rem],:])  # make all the scenario previously merged with the removed one equal to the one in idx_aug, up to the root node
            for j in np.arange(np.shape(to_merge_idx)[0]):
                S[: i+1, to_merge_idx[j,1],:] =  S[: i+1, idx_aug,:]

            # update the differential accuracy
            if not Tol[i]==infty:
                delta_p = delta_p + dp_min
            delta_rel = delta_p / delta_max

            delta_p_2 = delta_p_2 + dp_min
            delta_rel_2 = delta_p_2 / delta_max

        if i > 0:
            # Update available scenarios in the previous timestep
            J[i - 1,:] = J[i,:]

            # update previous timestep probabilities
            P[i - 1,:] = P[i,:]
            D[~J[i,:], ~J[i,:]] = infty

        #print('Branches t=%i: %i' %  (i, branches))

    S = S[:, J[-1,:] > 0,:]
    P = P[:, J[-1,:] > 0]

    return S, P, J


def read_input_data_5min(data_input_path):
    load_shape = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_load_Shape.csv")
    wind_shape = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_wind_Shape.csv")
    solar_shape = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_solar_Shape.csv")

    load_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_load_MWh.csv")
    wind_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_wind_MWh.csv")
    solar_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_solar_MWh.csv")

    net_load_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_5min_netload.csv")

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh


def read_input_data_Hourly(data_input_path):
    load_shape = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_load_Shape.csv")
    wind_shape = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_wind_Shape.csv")
    solar_shape = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_solar_Shape.csv")

    load_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_load_MWh.csv")
    wind_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_wind_MWh.csv")
    solar_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_solar_MWh.csv")

    net_load_MWh = pd.read_csv(data_input_path + "/" + "Input_Raw_Scenarios_60min_netload.csv")

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh


def generate_input_data_5min(data_input_path, data_location_timeseries, windCapacity, solarCapacity, peakDemand):

    timeSeriesParams_load = pd.read_csv(
        data_location_timeseries + "/" + "timeseries_data_files/Load/timeseries_load_5mins.csv")
    timeSeriesParams_wind = pd.read_csv(
        data_location_timeseries + "/" + "timeseries_data_files/WIND/timeseries_wind_5mins.csv")
    timeSeriesParams_solar = pd.read_csv(
        data_location_timeseries + "/" + "timeseries_data_files/PV/timeseries_pv_5mins.csv")

    nrows = 24 * 12 + 1
    ncols = 365
    columnNames = ["Time"]
    for i in range(1, ncols + 1):
        columnNames.append("Scenario" + str(i))
    rawScenarios = pd.DataFrame(index=range(nrows), columns=columnNames)
    rawScenarios["Time"] = np.append(np.arange(1, len(rawScenarios)), "Prob")

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_load.iloc[j:k, 1].values * peakDemand \
                                            - timeSeriesParams_wind.iloc[j:k, 1].values * windCapacity \
                                            - timeSeriesParams_solar.iloc[j:k, 1].values * solarCapacity
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_netload.csv', index=False)

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_load.iloc[j:k, 1].values * peakDemand
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_load_MWh.csv', index=False)

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_wind.iloc[j:k, 1].values * windCapacity
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_wind_MWh.csv', index=False)

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_solar.iloc[j:k, 1].values * solarCapacity
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_solar_MWh.csv', index=False)

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_load.iloc[j:k, 1].values
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_load_Shape.csv', index=False)

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_wind.iloc[j:k, 1].values
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_wind_Shape.csv', index=False)

    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_solar.iloc[j:k, 1].values
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_5min_solar_Shape.csv', index=False)


def generate_input_data_hourly(data_input_path, data_location_timeseries, windCapacity, solarCapacity, peakDemand):

    timeSeriesParams_load = pd.read_csv(
        data_location_timeseries +  "timeseries_data_files/Load/timeseries_load_hourly.csv")
    timeSeriesParams_wind = pd.read_csv(
        data_location_timeseries +  "timeseries_data_files/WIND/timeseries_wind_hourly.csv")
    timeSeriesParams_solar = pd.read_csv(
        data_location_timeseries +  "timeseries_data_files/PV/timeseries_pv_hourly.csv")

    # timeSeriesParams_load = pd.read_csv(
    #     data_location_timeseries + "/" + "timeseries_data_files/Load/timeseries_load_hourly_scenario_reduction.csv")
    # timeSeriesParams_wind = pd.read_csv(
    #     data_location_timeseries + "/" + "timeseries_data_files/WIND/timeseries_wind_hourly_scenario_reduction.csv")
    # timeSeriesParams_solar = pd.read_csv(
    #     data_location_timeseries + "/" + "timeseries_data_files/PV/timeseries_pv_hourly_scenario_reduction.csv")

    nrows = 24 + 1
    ncols = 365
    columnNames = ["Time"]
    for i in range(1, ncols + 1):
        columnNames.append("Scenario" + str(i))
    rawScenarios = pd.DataFrame(index=range(nrows), columns=columnNames)
    rawScenarios["Time"] = np.append(np.arange(1, len(rawScenarios)), "Prob")

    # Netload
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_load.iloc[j:k, 1].values * peakDemand \
                                            - timeSeriesParams_wind.iloc[j:k, 1].values * windCapacity \
                                            - timeSeriesParams_solar.iloc[j:k, 1].values * solarCapacity
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_netload.csv', index=False)

    # Load (MWh)
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_load.iloc[j:k, 1].values * peakDemand
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_load_MWh.csv', index=False)

    # Wind (MWh)
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_wind.iloc[j:k, 1].values * windCapacity
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_wind_MWh.csv', index=False)

    # Solar (MWh)
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_solar.iloc[j:k, 1].values * solarCapacity
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_solar_MWh.csv', index=False)

    # Load shape
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_load.iloc[j:k, 1].values
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_load_Shape.csv', index=False)

    # Wind shape
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_wind.iloc[j:k, 1].values
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_wind_Shape.csv', index=False)

    # Solar shape
    j = 0
    k = 0
    for i in range(1, ncols + 1):
        k = j + nrows - 1
        rawScenarios.iloc[0:nrows - 1, i] = timeSeriesParams_solar.iloc[j:k, 1].values
        j = k

    rawScenarios.to_csv(data_input_path + "/" + r'Input_Raw_Scenarios_60min_solar_Shape.csv', index=False)


def identify_extreme_points(load, wind, solar, load_MWh, wind_MWh, solar_MWh, net_load_MWh, **kwargs):
    extreme_scenarios = []
    extreme_datapoint = {}
    idx = 0

    # peak demand
    value = load_MWh.iloc[:-1, 1:].max().max()
    sce_id = load_MWh.iloc[:-1, 1:].max().idxmax()
    extreme_scenarios.append(sce_id)
    extreme_datapoint[idx] = {}
    extreme_datapoint[idx]["type"] = "peak demand"
    extreme_datapoint[idx]["value"] = value
    extreme_datapoint[idx]["scenario"] = sce_id
    idx += 1

    # peak net demand
    value = net_load_MWh.iloc[:-1, 1:].max().max()
    sce_id = net_load_MWh.iloc[:-1, 1:].max().idxmax()
    extreme_scenarios.append(sce_id)
    extreme_datapoint[idx] = {}
    extreme_datapoint[idx]["type"] = "peak net demand"
    extreme_datapoint[idx]["value"] = value
    extreme_datapoint[idx]["scenario"] = sce_id
    idx += 1

    # # peak ramp (load, pos and neg)
    # benchmark_df = load_MWh.copy(deep=True)
    # load_benchmark = benchmark_df.iloc[:-1, 1:]
    # ramp_benchmark = pd.DataFrame().reindex_like(load_benchmark)
    # for idx in range(1, ramp_benchmark.shape[0]):
    #     ramp_benchmark.iloc[idx, :] = load_benchmark.iloc[idx, :] - load_benchmark.iloc[idx - 1, :]
    # for idx in range(1, ramp_benchmark.shape[1]):
    #     ramp_benchmark.iloc[0, idx] = load_benchmark.iloc[0, idx] - load_benchmark.iloc[-1, idx - 1]
    # ramp_benchmark.iloc[0, 0] = load_benchmark.iloc[0, 0] - load_benchmark.iloc[-1, -1]
    # # ramp_benchmark = ramp_benchmark.abs()

    # value = ramp_benchmark.iloc[:, 1:].max().max()
    # sce_id = ramp_benchmark.iloc[:, 1:].max().idxmax()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[2] = {}
    # extreme_datapoint[2]["type"] = "peak ramp (pos, load)"
    # extreme_datapoint[2]["value"] = value
    # extreme_datapoint[2]["scenario"] = sce_id

    # value = ramp_benchmark.iloc[:, 1:].min().min()
    # sce_id = ramp_benchmark.iloc[:, 1:].min().idxmin()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[3] = {}
    # extreme_datapoint[3]["type"] = "peak ramp (neg, load)"
    # extreme_datapoint[3]["value"] = value
    # extreme_datapoint[3]["scenario"] = sce_id

    # # peak net ramp (pos and neg)
    # benchmark_df = net_load_MWh.copy(deep=True)
    # load_benchmark = benchmark_df.iloc[:-1, 1:]
    # ramp_benchmark = pd.DataFrame().reindex_like(load_benchmark)
    # for idx in range(1, ramp_benchmark.shape[0]):
    #     ramp_benchmark.iloc[idx, :] = load_benchmark.iloc[idx, :] - load_benchmark.iloc[idx - 1, :]
    # for idx in range(1, ramp_benchmark.shape[1]):
    #     ramp_benchmark.iloc[0, idx] = load_benchmark.iloc[0, idx] - load_benchmark.iloc[-1, idx - 1]
    # ramp_benchmark.iloc[0, 0] = load_benchmark.iloc[0, 0] - load_benchmark.iloc[-1, -1]

    # value = ramp_benchmark.iloc[:, 1:].max().max()
    # sce_id = ramp_benchmark.iloc[:, 1:].max().idxmax()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[4] = {}
    # extreme_datapoint[4]["type"] = "peak ramp (pos, netload)"
    # extreme_datapoint[4]["value"] = value
    # extreme_datapoint[4]["scenario"] = sce_id

    # value = ramp_benchmark.iloc[:, 1:].min().min()
    # sce_id = ramp_benchmark.iloc[:, 1:].min().idxmin()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[5] = {}
    # extreme_datapoint[5]["type"] = "peak ramp (neg, netload)"
    # extreme_datapoint[5]["value"] = value
    # extreme_datapoint[5]["scenario"] = sce_id

    # peak wind shape
    value = wind.iloc[:-1, 1:].max().max()
    sce_id = wind.iloc[:-1, 1:].max().idxmax()
    extreme_scenarios.append(sce_id)
    extreme_datapoint[idx] = {}
    extreme_datapoint[idx]["type"] = "peak wind shape"
    extreme_datapoint[idx]["value"] = value
    extreme_datapoint[idx]["scenario"] = sce_id
    idx += 1

    # peak solar shape
    value = solar.iloc[:-1, 1:].max().max()
    sce_id = solar.iloc[:-1, 1:].max().idxmax()
    extreme_scenarios.append(sce_id)
    extreme_datapoint[idx] = {}
    extreme_datapoint[idx]["type"] = "peak solar shape"
    extreme_datapoint[idx]["value"] = value
    extreme_datapoint[idx]["scenario"] = sce_id
    idx += 1

    # # peak wind shape (daily sum)
    # value = wind.iloc[:-1, 1:].sum().max()
    # sce_id = wind.iloc[:-1, 1:].sum().idxmax()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[idx] = {}
    # extreme_datapoint[idx]["type"] = "peak daily wind shape"
    # extreme_datapoint[idx]["value"] = value
    # extreme_datapoint[idx]["scenario"] = sce_id
    # idx += 1

    # # peak solar shape (daily sum)
    # value = solar.iloc[:-1, 1:].sum().max()
    # sce_id = solar.iloc[:-1, 1:].sum().idxmax()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[idx] = {}
    # extreme_datapoint[idx]["type"] = "peak daily solar shape"
    # extreme_datapoint[idx]["value"] = value
    # extreme_datapoint[idx]["scenario"] = sce_id
    # idx += 1

    # lowest wind shape (daily sum)
    value = wind.iloc[:-1, 1:].sum().min()
    sce_id = wind.iloc[:-1, 1:].sum().idxmin()
    extreme_scenarios.append(sce_id)
    extreme_datapoint[idx] = {}
    extreme_datapoint[idx]["type"] = "lowest daily wind shape"
    extreme_datapoint[idx]["value"] = value
    extreme_datapoint[idx]["scenario"] = sce_id
    idx += 1

    # lowest solar shape (daily sum)
    value = solar.iloc[:-1, 1:].sum().min()
    sce_id = solar.iloc[:-1, 1:].sum().idxmin()
    extreme_scenarios.append(sce_id)
    extreme_datapoint[idx] = {}
    extreme_datapoint[idx]["type"] = "lowest daily solar shape"
    extreme_datapoint[idx]["value"] = value
    extreme_datapoint[idx]["scenario"] = sce_id
    idx += 1

    # # peak wind ramp (pos and neg)
    # benchmark_df = wind.copy(deep=True)
    # load_benchmark = benchmark_df.iloc[:-1, 1:]
    # ramp_benchmark = pd.DataFrame().reindex_like(load_benchmark)
    # for idx in range(1, ramp_benchmark.shape[0]):
    #     ramp_benchmark.iloc[idx, :] = load_benchmark.iloc[idx, :] - load_benchmark.iloc[idx - 1, :]
    # for idx in range(1, ramp_benchmark.shape[1]):
    #     ramp_benchmark.iloc[0, idx] = load_benchmark.iloc[0, idx] - load_benchmark.iloc[-1, idx - 1]
    # ramp_benchmark.iloc[0, 0] = load_benchmark.iloc[0, 0] - load_benchmark.iloc[-1, -1]


    # value = ramp_benchmark.iloc[:, 1:].max().max()
    # sce_id = ramp_benchmark.iloc[:, 1:].max().idxmax()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[12] = {}
    # extreme_datapoint[12]["type"] = "peak wind ramp (pos, shape)"
    # extreme_datapoint[12]["value"] = value
    # extreme_datapoint[12]["scenario"] = sce_id

    # value = ramp_benchmark.iloc[:, 1:].min().min()
    # sce_id = ramp_benchmark.iloc[:, 1:].min().idxmin()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[13] = {}
    # extreme_datapoint[13]["type"] = "peak wind ramp (neg, shape)"
    # extreme_datapoint[13]["value"] = value
    # extreme_datapoint[13]["scenario"] = sce_id

    # # peak solar ramp (pos and neg)
    # benchmark_df = solar.copy(deep=True)
    # load_benchmark = benchmark_df.iloc[:-1, 1:]
    # ramp_benchmark = pd.DataFrame().reindex_like(load_benchmark)
    # for idx in range(1, ramp_benchmark.shape[0]):
    #     ramp_benchmark.iloc[idx, :] = load_benchmark.iloc[idx, :] - load_benchmark.iloc[idx - 1, :]
    # for idx in range(1, ramp_benchmark.shape[1]):
    #     ramp_benchmark.iloc[0, idx] = load_benchmark.iloc[0, idx] - load_benchmark.iloc[-1, idx - 1]
    # ramp_benchmark.iloc[0, 0] = load_benchmark.iloc[0, 0] - load_benchmark.iloc[-1, -1]

    # value = ramp_benchmark.iloc[:, 1:].max().max()
    # sce_id = ramp_benchmark.iloc[:, 1:].max().idxmax()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[14] = {}
    # extreme_datapoint[14]["type"] = "peak solar ramp (pos, shape)"
    # extreme_datapoint[14]["value"] = value
    # extreme_datapoint[14]["scenario"] = sce_id

    # value = ramp_benchmark.iloc[:, 1:].min().min()
    # sce_id = ramp_benchmark.iloc[:, 1:].min().idxmin()
    # extreme_scenarios.append(sce_id)
    # extreme_datapoint[15] = {}
    # extreme_datapoint[15]["type"] = "peak solar ramp (neg, shape)"
    # extreme_datapoint[15]["value"] = value
    # extreme_datapoint[15]["scenario"] = sce_id

    return extreme_scenarios, extreme_datapoint


def plot_ramp_duration_curve(dataframe_raw, selected_sce_prob, selected_sce, num_scenarios, type, output_path):

    # benchmark data
    benchmark_df = dataframe_raw.copy(deep=True)
    netload_benchmark = benchmark_df.iloc[:-1, 1:]
    ramp_benchmark = pd.DataFrame().reindex_like(netload_benchmark)
    for idx in range(1, ramp_benchmark.shape[0]):
        ramp_benchmark.iloc[idx, :] = netload_benchmark.iloc[idx, :] - netload_benchmark.iloc[idx - 1, :]
    for idx in range(1, ramp_benchmark.shape[1]):
        ramp_benchmark.iloc[0, idx] = netload_benchmark.iloc[0, idx] - netload_benchmark.iloc[-1, idx - 1]
    ramp_benchmark.iloc[0, 0] = netload_benchmark.iloc[0, 0] - netload_benchmark.iloc[-1, -1]
    ramp_benchmark = ramp_benchmark.abs()

    rdc_benchmark_np = ramp_benchmark.to_numpy()
    rdc_benchmark_data = np.concatenate([rdc_benchmark_np[c] for c in range(0, len(rdc_benchmark_np))])
    rdc_benchmark_data[::-1].sort()
    rx_benchmark = [5 + 5 * x for x in range(0, len(rdc_benchmark_data))]

    # Scenario reduction
    rdc_sce = []
    rx_sce = []
    total_ramp_MWh = 0
    adjusted_prob = [i * 365 for i in selected_sce_prob]

    for idx in range(0, len(selected_sce)):
        rdc_sce.append(ramp_benchmark.iloc[:, selected_sce[idx] - 1].to_numpy())
        rx_sce.append([5 * adjusted_prob[idx]] * len(ramp_benchmark.iloc[:, idx]))
        total_ramp_MWh += sum(ramp_benchmark.iloc[:, selected_sce[idx] - 1]) * adjusted_prob[idx]

    rdc_sce_data = np.concatenate([rdc_sce[c] for c in range(0, len(rdc_sce))])
    rx_sce_data = np.concatenate([rx_sce[c] for c in range(0, len(rx_sce))])

    zipped_lists = zip(rdc_sce_data, rx_sce_data)
    sorted_pairs = sorted(zipped_lists, reverse=True)

    tuples = zip(*sorted_pairs)
    rdc_sce_data, rx_sce_data = [list(tuple) for tuple in tuples]

    rx_sce_final = [rx_sce_data[0]]
    for idx in range(1, len(rx_sce_data)):
        time_step = rx_sce_final[idx - 1] + rx_sce_data[idx]
        rx_sce_final.append(time_step)

    gap_percent = 100 * (sum(rdc_benchmark_data) - total_ramp_MWh) / sum(rdc_benchmark_data)
    Total_ramp_gap_value = (sum(rdc_benchmark_data) / 12 - total_ramp_MWh / 12)

    NRMSD = calculate_NRMSD(rx_benchmark, rdc_benchmark_data, rx_sce_final, rdc_sce_data)

    textstr = '\n'.join((
        r'# of scenarios=%.0f' % (num_scenarios),
        r'Annual gap (percent)=%.2f' % (gap_percent),
        r'NRMSD =%.3f' % (NRMSD)),
    )

    fig, axes = plt.subplots(1, 1)
    axes.plot(rx_benchmark, rdc_benchmark_data, color="blue", label="benchmark")
    axes.plot(rx_sce_final, rdc_sce_data, color="red", label="scenario reduction")
    axes.legend()
    axes.set_title("Duration Curve Comparison:" + type)
    axes.set_xlabel("time steps")
    axes.set_ylabel(type)

    # # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #
    # # place a text box in upper left in axes coords
    axes.text(0.55, 0.80, textstr, transform=axes.transAxes, fontsize=11,
              verticalalignment='top', bbox=props)

    # plt.show()
    plt.savefig(output_path + '/ramp_duration_curve_' + type + '_' + str(num_scenarios) + '.png', transparent=True,
                dpi=160,
                bbox_inches="tight")
    plt.close()
    # print (num_scenarios, "\t", type, "\t", gap_percent, "\t", NRMSD)


def plot_duration_curve(dataframe_raw, selected_sce_prob, selected_sce, num_scenarios, type, output_path):

    # benchmark data
    benchmark_df = dataframe_raw.copy(deep=True)
    benchmark_raw = benchmark_df.iloc[:-1, 1:].to_numpy()
    benchmark_data = np.concatenate([benchmark_raw[c] for c in range(0, len(benchmark_raw))])
    benchmark_data[::-1].sort()
    benchmark_x = [5 + 5 * x for x in range(0, len(benchmark_data))]

    # Scenario reduction
    scenario = []
    scenario_x = []
    total_value = 0
    adjusted_prob = [i * 365 for i in selected_sce_prob]

    for idx in range(0, len(selected_sce)):
        scenario.append(benchmark_df.iloc[:-1, selected_sce[idx]].to_numpy())
        scenario_x.append([5 * adjusted_prob[idx]] * len(benchmark_df.iloc[:-1, idx]))
        total_value += sum(benchmark_df.iloc[:-1, selected_sce[idx]]) * adjusted_prob[idx]

    scenario_data = np.concatenate([scenario[c] for c in range(0, len(scenario))])
    scenario_x_raw = np.concatenate([scenario_x[c] for c in range(0, len(scenario_x))])

    zipped_lists = zip(scenario_data, scenario_x_raw)
    sorted_pairs = sorted(zipped_lists, reverse=True)

    tuples = zip(*sorted_pairs)
    scenario_data, scenario_x_raw = [list(tuple) for tuple in tuples]

    scenario_x_final = [scenario_x_raw[0]]
    for idx in range(1, len(scenario_x_raw)):
        time_step = scenario_x_final[idx - 1] + scenario_x_raw[idx]
        scenario_x_final.append(time_step)

    gap_percent = 100 * (sum(benchmark_data) - total_value) / sum(benchmark_data)
    gap_value = (sum(benchmark_data) / 12 - total_value / 12)

    NRMSD = calculate_NRMSD(benchmark_x, benchmark_data, scenario_x_final, scenario_data)

    textstr = '\n'.join((
        r'# of scenarios=%.0f' % (num_scenarios),
        r'Annual gap (percent)=%.2f' % (gap_percent),
        r'NRMSD =%.3f' % (NRMSD)),
    )

    fig, axes = plt.subplots(1, 1)
    axes.plot(benchmark_x, benchmark_data, color="blue", label="benchmark")
    axes.plot(scenario_x_final, scenario_data, color="red", label="scenario reduction")
    axes.legend()
    axes.set_title("Duration Curve Comparison:" + type)
    axes.set_xlabel("time steps")
    axes.set_ylabel(type)

    # # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #
    # # place a text box in upper left in axes coords
    axes.text(0.55, 0.80, textstr, transform=axes.transAxes, fontsize=11,
              verticalalignment='top', bbox=props)

    # plt.show()
    plt.savefig(output_path + '/duration_curve_' + type + '_' + str(num_scenarios) + '.png', transparent=True, dpi=160,
                   bbox_inches="tight")
    plt.close()

    # print (num_scenarios, "\t", type, "\t", gap_percent, "\t", NRMSD)


def calculate_NRMSD(benchmark_x, benchmark_data, scenario_x, scenario_data):

    from scipy import interpolate
    import math

    widths = np.array(benchmark_x)
    heights = np.array(benchmark_data)

    heights_smooth = interpolate.splrep(widths, heights)

    # Select desired width values
    width_vals = benchmark_x

    # splev returns the value of your spline evaluated at the width values.
    benchmark_heights = interpolate.splev(width_vals, heights_smooth)

    # Scenario
    widths_sce = np.array(scenario_x)
    heights_sce = np.array(scenario_data)

    heights_smooth_scenario = interpolate.splrep(widths_sce, heights_sce)

    scenario_heights = interpolate.splev(width_vals, heights_smooth_scenario)

    # Normalized root mean square deviation
    numer = math.sqrt((1.0 / len(benchmark_heights)) * sum(
        pow(scenario_heights[i] - benchmark_heights[i], 2) for i in range(0, len(benchmark_heights))))
    denom = (1.0 / len(benchmark_heights)) * sum(benchmark_heights[i] for i in range(0, len(benchmark_heights)))
    NRMSD = numer / denom

    return NRMSD

