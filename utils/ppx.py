import numpy as np
import pandas as pd
import sqlite3
from pathlib import Path
import yaml

with open("ppx_settings.yml", "r") as setfile:
    settings = yaml.load(setfile, Loader=yaml.FullLoader)

p_runs = settings["runs"].copy()
if "baseline" in list(settings.keys()):
    p_runs.remove(settings["baseline"])


def fmt_dollarscents(x):
    return "${:.2f}".format(x)

def fmt_dollars(x):
    return "${:.0f}".format(x)

def fmt_1dec(x):
    return "{:.1f}".format(x)

def fmt_3dec(x):
    return"{:.3f}".format(x)


def summarize_annual_dispatch(cons):
    wa_lambda = pd.DataFrame()
    ENS = pd.DataFrame()

    for run in settings["runs"]:
        con = cons[run]
        df = pd.read_sql_query("SELECT * FROM annual_dispatch_summary", con)

        if len(wa_lambda) == 0:
            wa_lambda = df[["wa_lambda"]]
            wa_lambda = wa_lambda.rename(columns={"wa_lambda": run})
        else:
            wa_lambda[run] = df[["wa_lambda"]]

        wa_lambda.index.names = ["Simulation year"]

        if len(ENS) == 0:
            ENS = df[["ENS"]]
            ENS = ENS.rename(columns={"ENS": run})
        else:
            ENS[run] = df[["ENS"]]

        ENS.index.names = ["Simulation year"]

    wa_lambda_D = wa_lambda.copy()
    ENS_D = ENS.copy()

    if "baseline" in list(settings.keys()):
        for run in p_runs:
            rname = f"{run}_D"
            wa_lambda_D[rname] = wa_lambda_D[run] / wa_lambda_D[settings["baseline"]]
            ENS_D[rname] = ENS_D[run] / ENS_D[settings["baseline"]]

    L_formats = {}
    E_formats = {}
    for run in settings["runs"]:
        L_formats[run] = fmt_dollarscents
        E_formats[run] = fmt_1dec
    for run in p_runs:
        rname = f"{run}_D"
        L_formats[rname] = fmt_3dec
        E_formats[rname] = fmt_3dec

    for col in list(wa_lambda_D.columns):
        wa_lambda_D[col] = wa_lambda_D[col].apply(L_formats[col])
    for col in list(ENS_D.columns):
        ENS_D[col] = ENS_D[col].apply(E_formats[col])

    dispatch_results = {
        "wa_lambda": wa_lambda_D,
        "ENS": ENS_D,
    }

    return dispatch_results


def summarize_agent_decisions(cons):
    for run in settings["runs"]:
        con = cons[run]
        decs = pd.read_sql_query("SELECT * FROM agent_decisions", con)
        decs["run"] = run
        decs = decs.drop(columns=["allowed"])

        try:
            all_decs = pd.concat([all_decs, decs], axis=0, ignore_index=True)
        except:
            all_decs = decs

    conditions = []
    for run in settings["runs"]:
        for agent in settings["agents"]:
            for y in sorted(all_decs["base_pd"].unique()):
                data = all_decs[(all_decs["run"] == run) & (all_decs["agent_id"] == agent) & (all_decs["base_pd"] == y)]
                if len(data) > 0:
                    if data["mode"].reset_index(drop=True)[0] == "ret_only":
                        condition = "Financial distress"
                    else:
                        if (data["NPV"] > 0).any():
                            condition = "Unconstrained decision"
                        else:
                            condition = "All NPVs negative"
                    conditions.append({
                        "run": run,
                        "agent": agent,
                        "base_pd": y,
                        "condition": condition,
                    })

    conditions = pd.DataFrame(conditions)

    conditions_pivot = pd.pivot_table(
        conditions,
        values = "condition",
        index = ["base_pd"],
        columns = ["run", "agent"],
        aggfunc = "sum",
    )[settings["runs"]]

    xtr = all_decs[all_decs["project_type"] == "new_xtr"]
    ret = all_decs[all_decs["project_type"] == "retirement"]

    xtr_pivot = pd.pivot_table(
        xtr,
        values = "units_to_execute",
        index = ["base_pd", "unit_type", "lag"],
        columns = ["run", "agent_id"],
        aggfunc = "sum"
    )[settings["runs"]]

    xtr_now = xtr.copy()
    xtr_now = xtr_now[xtr_now["lag"] == 0]

    xtr_exec = xtr_now.copy()
    xtr_exec = xtr_exec[xtr_exec["units_to_execute"] != 0]
    print(xtr_exec)
    if len(xtr_exec) > 0:
        xtr_exec_p = pd.pivot_table(
            xtr_exec,
            values = "units_to_execute",
            index = ["base_pd", "unit_type"],
            columns = ["run", "agent_id"],
            aggfunc = "sum",
        )[settings["runs"]]
        xtr_exec_p = xtr_exec_p.fillna(0)
    else:
        xtr_exec_p = None

    ret_pivot = pd.pivot_table(
        ret,
        values = "units_to_execute",
        index = ["base_pd", "unit_type", "ret_pd", "lag"],
        columns = ["run", "agent_id"],
        aggfunc = "sum"
    )[settings["runs"]]

    ret_exec = ret[ret["units_to_execute"] != 0]
    if len(ret_exec) > 0:
        ret_exec_p = pd.pivot_table(
            ret_exec,
            values = "units_to_execute",
            index = ["base_pd", "unit_type", "ret_pd", "lag"],
            columns = ["run", "agent_id"],
            aggfunc = "sum"
        )[settings["runs"]]
    else:
        ret_exec_p = None

    decisions = {
        "conditions": conditions_pivot,
        "xtr_pivot": xtr_pivot,
        "xtr_exec_pivot": xtr_exec_p,
        "ret_pivot": ret_pivot,
        "ret_exec_pivot": ret_exec_p,
    }

    return decisions


def summarize_agent_fss(cons):
    # Create dictionaries to facilitate run sorting
    i = 0
    runs_fwd = {}   # {run_name: i_run_name}
    runs_rev = {}   # {i_run_name: run_name}
    for run in settings["runs"]:
        new_name = f"{i}_{run}"
        runs_fwd[run] = new_name
        runs_rev[new_name] = run
        i += 1

    for run in settings["runs"]:
        con = cons[run]
        fs = pd.read_sql_query("SELECT * FROM forecasted_agent_fss", con)
        fs["run"] = runs_fwd[run]
        fs["delta_pd"] = fs["projected_pd"] - fs["base_pd"]

        try:
            all_fss = pd.concat([all_fss, fs], axis=0, ignore_index=True)
        except:
            all_fss = fs

    # Filter out distant years--not needed for diagnostics
    all_fss = all_fss[all_fss["delta_pd"] < 14].reset_index(drop=True)
    all_fss = all_fss[["run", "agent_id", "base_pd", "projected_pd", "revenue", "FCF", "moodys_score"]]

    fs_results = {}

    for agent in settings["agents"]:
        fs_a = all_fss[all_fss["agent_id"] == agent].reset_index(drop=True).drop(columns=["agent_id"])
        fs_results[agent] = fs_a

        for datum in ["revenue", "FCF", "moodys_score"]:
            dname = f"{agent}_{datum}"
            fs_results[dname] = pd.pivot_table(
                fs_a,
                values = datum,
                index = ["run", "projected_pd"],
                columns = ["base_pd"],
                aggfunc = "sum",
            )

            fs_results[dname] = fs_results[dname].rename(index=runs_rev)

    return fs_results


def write_to_excel(everything):
    # Write the results to xlsx
    xlsxwriter = pd.ExcelWriter(
        Path("/filespace/k/kebiegel/abce/utils") / f"{settings['groupname']}.xlsx"
    )

    names = {
        "wa_lambda": "Wtd Avg Lambda",
        "ENS": "Energy Not Served",
        "conditions": "Decision conditions",
        "xtr_pivot": "New xtr pivot",
        "xtr_exec_pivot": "New xtr executed",
        "ret_pivot": "Retirement pivot",
        "ret_exec_pivot": "Retirements executed",
        "201_revenue": "201 revenue",
        "201_FCF": "201 FCF",
        "201_moodys_score": "201 score",
        "202_revenue": "202 revenue",
        "202_FCF": "202 FCF",
        "202_moodys_score": "202 score",
    }

    with xlsxwriter as writer:
        for key in everything.keys():
            everything[key].to_excel(writer, sheet_name = names[key])
#        everything["wa_lambda"].to_excel(writer, sheet_name = "Wtd Avg Lambda")
#        everything["ENS"].to_excel(writer, sheet_name = "Energy Not Served")
#        everything["conditions"].to_excel(writer, sheet_name = "Decision conditions")
#        everything["xtr_pivot"].to_excel(writer, sheet_name = "New xtr pivot")
#        everything["xtr_exec_pivot"].to_excel(writer, sheet_name = "New xtr L0 pivot")
#        everything["ret_pivot"].to_excel(writer, sheet_name = "Retirement pivot")
#        everything["ret_exec_pivot"].to_excel(writer, sheet_name = "Retirements executed")
#        everything["201_revenue"].to_excel(writer, sheet_name = "201 revenue")
#        everything["201_FCF"].to_excel(writer, sheet_name = "201 FCF")
#        everything["201_moodys_score"].to_excel(writer, sheet_name = "201 score")
#        everything["202_revenue"].to_excel(writer, sheet_name = "202 revenue")
#        everything["202_FCF"].to_excel(writer, sheet_name = "202 FCF")
#        everything["202_moodys_score"].to_excel(writer, sheet_name = "202 score")



def run():
    cons = {}
    for fname in settings["runs"]:
        full_fname = Path(settings["root_dir"]) / f"{settings['common_prefix']}{fname}" / "abce_db.db"
        print(full_fname)
        cons[fname] = sqlite3.connect(full_fname)

    # Invoke all processing functions
    dispatch_results = summarize_annual_dispatch(cons)
    decision_results = summarize_agent_decisions(cons)
    fs_results = summarize_agent_fss(cons)

    # Save all the dataframes to an xlsx
    everything = {
                  "wa_lambda": dispatch_results["wa_lambda"],
                  "ENS": dispatch_results["ENS"],
                  "conditions": decision_results["conditions"],
                  "xtr_pivot": decision_results["xtr_pivot"],
                  "ret_pivot": decision_results["ret_pivot"],
                 }

    if decision_results["xtr_exec_pivot"] is not None:
        everything["xtr_exec_pivot"] = decision_results["xtr_exec_pivot"]

    if decision_results["ret_exec_pivot"] is not None:
        everything["ret_exec_pivot"] = decision_results["ret_exec_pivot"]

    for agent in settings["agents"]:
        everything[f"{agent}_revenue"] = fs_results[f"{agent}_revenue"]
        everything[f"{agent}_FCF"] = fs_results[f"{agent}_FCF"]
        everything[f"{agent}_moodys_score"] = fs_results[f"{agent}_moodys_score"]

    write_to_excel(everything)


if __name__ == "__main__":
    run()
