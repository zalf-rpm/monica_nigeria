#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from collections import defaultdict, OrderedDict
import csv
from datetime import datetime
import gc
import json
import numpy as np
import os
from pyproj import CRS, Transformer
import sqlite3
import sys
import timeit
import types
import zmq

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    "mbm-local-remote": {
        "path-to-data-dir": "data/",
        "path-to-output-dir": "out/",
        "path-to-csv-output-dir": "csv-out/"
    },
    "remoteConsumer-remoteMonica": {
        "path-to-data-dir": "./data/",
        "path-to-output-dir": "/out/out/",
        "path-to-csv-output-dir": "/out/csv-out/"
    }
}


def create_output(msg):
    cm_count_to_vals = defaultdict(dict)
    for data in msg.get("data", []):
        results = data.get("results", [])

        is_daily_section = data.get("origSpec", "") == '"daily"'

        for vals in results:
            if "CM-count" in vals:
                cm_count_to_vals[vals["CM-count"]].update(vals)
            elif is_daily_section:
                cm_count_to_vals[vals["Date"]].update(vals)

    cmcs = list(cm_count_to_vals.keys())
    cmcs.sort()
    last_cmc = cmcs[-1]
    if "Year" not in cm_count_to_vals[last_cmc]:
        cm_count_to_vals.pop(last_cmc)

    return cm_count_to_vals


def write_row_to_grids(row_col_data, row, col_0, no_of_cols, header, path_to_output_dir, setup_id):
    "write grids row by row"

    if not hasattr(write_row_to_grids, "nodata_row_count"):
        write_row_to_grids.nodata_row_count = defaultdict(lambda: 0)
        write_row_to_grids.list_of_output_files = defaultdict(list)

    def make_dict_nparr():
        return defaultdict(lambda: np.full((no_of_cols,), -9999, dtype=np.float))

    output_grids = {
        "Yield": {"data": make_dict_nparr(), "cast-to": "float", "digits": 2},
        "AbBiom": {"data": make_dict_nparr(), "cast-to": "float", "digits": 2},
        "LAI": {"data": make_dict_nparr(), "cast-to": "float", "digits": 2},
        # "tradefavg": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        # "heatredavg": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        # "frostredavg": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
    }
    output_keys = list(output_grids.keys())

    cmc_to_crop = {}

    is_no_data_row = True
    # skip this part if we write just a nodata line
    if row in row_col_data:
        no_data_cols = no_of_cols
        for col in range(col_0, col_0 + no_of_cols + 1):
            if col in row_col_data[row]:
                rcd_val = row_col_data[row][col]
                if rcd_val == -9999:
                    continue
                else:
                    no_data_cols -= 1
                    cmc_and_year_to_vals = defaultdict(lambda: defaultdict(list))
                    for cell_data in rcd_val:
                        # if we got multiple datasets per cell, iterate over them and aggregate them in the following step
                        for cm_count, data in cell_data.items():
                            for key in output_keys:
                                # store mapping cm_count to crop name for later file name creation
                                if cm_count not in cmc_to_crop and "Crop" in data:
                                    cmc_to_crop[cm_count] = data["Crop"]

                                # only further process/store data we actually received
                                if key in data:
                                    v = data[key]
                                    if isinstance(v, list):
                                        for i, v_ in enumerate(v):
                                            cmc_and_year_to_vals[(cm_count, data["Year"])][f"{key}_{i + 1}"].append(v_)
                                    else:
                                        cmc_and_year_to_vals[(cm_count, data["Year"])][key].append(v)
                                # if a key is missing, because that monica event was never raised/reached, create the empty list
                                # so a no-data value is being produced
                                else:
                                    cmc_and_year_to_vals[(cm_count, data["Year"])][key]

                    # potentially aggregate multiple data per cell and finally store them for this row
                    for (cm_count, year), key_to_vals in cmc_and_year_to_vals.items():
                        for key, vals in key_to_vals.items():
                            output_vals = output_grids[key]["data"]
                            if len(vals) > 0:
                                output_vals[(cm_count, year)][col - col_0] = sum(vals) / len(vals)
                            else:
                                output_vals[(cm_count, year)][col - col_0] = -9999

        is_no_data_row = no_data_cols == no_of_cols

    if is_no_data_row:
        write_row_to_grids.nodata_row_count[setup_id] += 1

    def write_nodata_rows(file_):
        for _ in range(write_row_to_grids.nodata_row_count[setup_id]):
            rowstr = " ".join(["-9999" for __ in range(no_of_cols)])
            file_.write(rowstr + "\n")

    # iterate over all prepared data for a single row and write row
    for key, y2d_ in output_grids.items():
        y2d = y2d_["data"]
        cast_to = y2d_["cast-to"]
        digits = y2d_.get("digits", 0)
        if cast_to == "int":
            mold = lambda x: str(int(x))
        else:
            mold = lambda x: str(round(x, digits))

        for (cm_count, year), row_arr in y2d.items():
            crop = cmc_to_crop[cm_count] if cm_count in cmc_to_crop else "none"
            crop = crop.replace("/", "").replace(" ", "")
            path_to_file = path_to_output_dir + crop + "_" + key + "_" + str(year) + "_" + str(cm_count) + ".asc"

            if not os.path.isfile(path_to_file):
                with open(path_to_file, "w") as _:
                    _.write(header)
                    write_row_to_grids.list_of_output_files[setup_id].append(path_to_file)

            with open(path_to_file, "a") as file_:
                write_nodata_rows(file_)
                rowstr = " ".join(["-9999" if int(x) == -9999 else mold(x) for x in row_arr])
                file_.write(rowstr + "\n")

    # clear the no-data row count when no-data rows have been written before a data row
    if not is_no_data_row:
        write_row_to_grids.nodata_row_count[setup_id] = 0

    # if we're at the end of the output and just empty lines are left, then they won't be written in the
    # above manner because there won't be any rows with data where they could be written before
    # so add no-data rows simply to all files we've written to before
    if is_no_data_row \
            and write_row_to_grids.list_of_output_files[setup_id] \
            and write_row_to_grids.nodata_row_count[setup_id] > 0:
        for path_to_file in write_row_to_grids.list_of_output_files[setup_id]:
            with open(path_to_file, "a") as file_:
                write_nodata_rows(file_)
        write_row_to_grids.nodata_row_count[setup_id] = 0

    if row in row_col_data:
        del row_col_data[row]


def run_consumer(leave_after_finished_run=True, server={"server": None, "port": None}):
    """collect data from workers"""

    config = {
        "mode": "mbm-local-remote",  # remote "mbm-local-remote", local "cj-local-remote"
        "port": server["port"] if server["port"] else "7777",  # local 7778,  remote 7777
        "server": server["server"] if server["server"] else "localhost",  # "login01.cluster.zalf.de",
        "timeout": 600000  # 10 minutes
    }

    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    paths = PATHS[config["mode"]]

    if not "out" in config:
        config["out"] = paths["path-to-output-dir"]
    if not "csv-out" in config:
        config["csv-out"] = paths["path-to-csv-output-dir"]

    print("consumer config:", config)

    context = zmq.Context()
    socket = context.socket(zmq.PULL)

    socket.connect("tcp://" + config["server"] + ":" + config["port"])
    socket.RCVTIMEO = config["timeout"]
    leave = False
    write_normal_output_files = False

    setup_id_to_data = defaultdict(lambda: {
        "header": None,
        "no_of_cols": None,
        "no_of_rows": None,
        "out_dir_exists": False,
        "row_col_data": defaultdict(lambda: defaultdict(list)),
        "cols@row_received": {},
        "next_row": None
    })

    def process_message(msg):
        if len(msg["errors"]) > 0:
            print("There were errors in message:", msg, "\nSkipping message!")
            return

        if not hasattr(process_message, "wnof_count"):
            process_message.wnof_count = 0
            process_message.setup_count = 0

        leave = False

        if not write_normal_output_files:
            custom_id = msg["customId"]
            setup_id = custom_id["setup_id"]
            region = custom_id["region"]
            planting = custom_id["planting"]
            nitrogen = custom_id["nitrogen"]
            crop = custom_id["crop"]

            data = setup_id_to_data[setup_id]

            row = custom_id["s_row"]
            col = custom_id["s_col"]
            no_of_cols = custom_id["no_of_s_cols"]
            no_of_rows = custom_id["no_of_s_rows"]
            row_0 = custom_id["s_row_0"]
            col_0 = custom_id["s_col_0"]
            if row not in data["cols@row_received"]:
                data["cols@row_received"][row] = 0
            if data["next_row"] is None:
                data["next_row"] = row_0
            if data["header"] is None:
                data["header"] = f"""ncols        {no_of_cols}
nrows        {no_of_rows}
xllcorner    {custom_id["b_lon_0"]}
yllcorner    {custom_id["b_lat_0"] - (no_of_rows * custom_id["s_resolution"])}
cellsize     {custom_id["s_resolution"]}
NODATA_value -9999
"""
            is_nodata = custom_id["nodata"]

            debug_msg = "received work result " + str(process_message.received_env_count) \
                        + " customId: " + str(msg.get("customId", "")) \
                        + " next row: " + str(data["next_row"]) \
                        + " cols@row to go: " + str(no_of_cols - data["cols@row_received"][row]) + "@" \
                        + str(row) + " cols_per_row: " + str(no_of_cols)
            print(debug_msg)
            # debug_file.write(debug_msg + "\n")
            if is_nodata:
                data["row_col_data"][row][col] = -9999
            else:
                data["row_col_data"][row][col].append(create_output(msg))
            data["cols@row_received"][row] += 1

            process_message.received_env_count = process_message.received_env_count + 1

            while (data["next_row"] in data["row_col_data"] and
                   data["cols@row_received"][data["next_row"]] == no_of_cols):   #\
                    #or (len(data["cols@row_received"]) > data["next_row"] and
                    #    data["cols@row_received"][data["next_row"]] == 0):

                path_to_out_dir = f"{config['out']}{setup_id}_reg-{region}_{crop}_plant-{planting}_{nitrogen}-N/"
                print(path_to_out_dir)
                if not data["out_dir_exists"]:
                    if os.path.isdir(path_to_out_dir) and os.path.exists(path_to_out_dir):
                        data["out_dir_exists"] = True
                    else:
                        try:
                            os.makedirs(path_to_out_dir)
                            data["out_dir_exists"] = True
                        except OSError:
                            print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
                            exit(1)

                write_row_to_grids(data["row_col_data"], data["next_row"], col_0, no_of_cols, data["header"],
                                   path_to_out_dir, setup_id)

                debug_msg = "wrote row: " + str(data["next_row"]) \
                            + " next_row: " + str(data["next_row"] + 1) \
                            + " rows unwritten: " + str(list(data["row_col_data"].keys()))
                print(debug_msg)
                # debug_file.write(debug_msg + "\n")

                data["next_row"] += 1  # move to next row (to be written)

                # this setup is finished
                if leave_after_finished_run and data["next_row"] > (row_0 + no_of_rows):
                    process_message.setup_count += 1

        elif write_normal_output_files:

            if msg.get("type", "") in ["jobs-per-cell", "no-data", "setup_data"]:
                # print "ignoring", result.get("type", "")
                return

            print("received work result ", process_message.received_env_count, " customId: ",
                  str(msg.get("customId", "")))

            custom_id = msg["customId"]
            setup_id = custom_id["setup_id"]
            row = custom_id["srow"]
            col = custom_id["scol"]
            # crow = custom_id.get("crow", -1)
            # ccol = custom_id.get("ccol", -1)
            # soil_id = custom_id.get("soil_id", -1)

            process_message.wnof_count += 1

            path_to_out_dir = config["out"] + str(setup_id) + "/" + str(row) + "/"
            print(path_to_out_dir)
            if not os.path.exists(path_to_out_dir):
                try:
                    os.makedirs(path_to_out_dir)
                except OSError:
                    print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
                    exit(1)

            # with open("out/out-" + str(i) + ".csv", 'wb') as _:
            with open(path_to_out_dir + "col-" + str(col) + ".csv", "w", newline='') as _:

                writer = csv.writer(_, delimiter=",")
                for data_ in msg.get("data", []):
                    results = data_.get("results", [])
                    orig_spec = data_.get("origSpec", "")
                    output_ids = data_.get("outputIds", [])

                    if len(results) > 0:
                        writer.writerow([orig_spec.replace("\"", "")])
                        for row in monica_io3.write_output_header_rows(output_ids,
                                                                       include_header_row=True,
                                                                       include_units_row=True,
                                                                       include_time_agg=False):
                            writer.writerow(row)

                        for row in monica_io3.write_output(output_ids, results):
                            writer.writerow(row)

                    writer.writerow([])

            process_message.received_env_count = process_message.received_env_count + 1

        return leave

    process_message.received_env_count = 1

    while not leave:
        try:
            # start_time_recv = timeit.default_timer()
            msg = socket.recv_json()  # encoding="latin-1"
            # elapsed = timeit.default_timer() - start_time_recv
            # print("time to receive message" + str(elapsed))
            # start_time_proc = timeit.default_timer()
            leave = process_message(msg)
            # elapsed = timeit.default_timer() - start_time_proc
            # print("time to process message" + str(elapsed))
        except zmq.error.Again as _e:
            print('no response from the server (with "timeout"=%d ms) ' % socket.RCVTIMEO)
            return
        except Exception as e:
            print("Exception:", e)
            # continue

    print("exiting run_consumer()")
    # debug_file.close()


if __name__ == "__main__":
    run_consumer()
