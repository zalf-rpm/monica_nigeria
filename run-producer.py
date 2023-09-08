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

from collections import defaultdict
import copy
import csv
from datetime import date, timedelta
import gzip
import json
import math
from netCDF4 import Dataset
import numpy as np
import os
from pyproj import CRS, Transformer
import sqlite3
import sqlite3 as cas_sq3
import sys
import time
import zmq

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    # adjust the local path to your environment
    "mbm-local-local": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        # "path-to-soil-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        "path-to-soil-dir": "/home/berg/Desktop/soil/",
        "monica-path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "mbm-local-remote": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        # "path-to-soil-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        "path-to-soil-dir": "/home/berg/Desktop/soil/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "hpc-local-remote": {
        #"path-to-climate-dir": "/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        # mounted path to archive or hard drive with climate data
        "path-to-soil-dir": "/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "remoteProducer-remoteMonica": {
        "path-to-climate-dir": "/data/",  # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-to-soil-dir": "/project/soil/global_soil_dataset_for_earth_system_modeling/",
        "path-debug-write-folder": "/out/debug-out/",
    }
}


def run_producer(server={"server": None, "port": None}):
    context = zmq.Context()
    socket = context.socket(zmq.PUSH)  # pylint: disable=no-member

    config = {
        "mode": "mbm-local-local",  # local:"cj-local-remote" remote "mbm-local-remote"
        "server-port": server["port"] if server["port"] else "6666",  # local: 6667, remote 6666
        "server": server["server"] if server["server"] else "localhost",  # "login01.cluster.zalf.de",
        "start_lat": "83.95833588",
        "end_lat": "-55.95833206",
        "start_lon": "-179.95832825",
        "end_lon": "179.50000000",
        "region": "nigeria",
        "resolution": "5min",  # 30sec,
        "path_to_dem_grid": "",
        "sim.json": "sim.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups.csv",
        "run-setups": "[3]"
    }

    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    print("config:", config)

    s_resolution = {"5min": 5 / 60., "30sec": 30 / 3600.}[config["resolution"]]
    s_res_scale_factor = {"5min": 60., "30sec": 3600.}[config["resolution"]]

    region_to_lat_lon_bounds = {
        "nigeria": {"tl": {"lat": 14.0, "lon": 2.7}, "br": {"lat": 4.25, "lon": 14.7}},
        "africa": {"tl": {"lat": 37.4, "lon": -17.55}, "br": {"lat": -34.9, "lon": 51.5}},
        "earth": {
            "5min": {"tl": {"lat": 83.95833588, "lon": -179.95832825},
                     "br": {"lat": -55.95833206, "lon": 179.50000000}},
            "30sec": {"tl": {"lat": 83.99578094, "lon": -179.99583435},
                      "br": {"lat": -55.99583435, "lon": 179.99568176}}
        }
    }

    # select paths
    paths = PATHS[config["mode"]]
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["server-port"]))

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"])

    # transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    # soil_crs_to_x_transformers = {}
    # wgs84_crs = CRS.from_epsg(4326)
    # utm32_crs = CRS.from_epsg(25832)
    # transformers[wgs84] = Transformer.from_crs(wgs84_crs, gk5_crs, always_xy=True)

    def get_lat_0_lon_0_resolution_from_grid_metadata(metadata):
        lat_0 = float(metadata["yllcorner"]) \
                    + (float(metadata["cellsize"]) * float(metadata["nrows"])) \
                    - (float(metadata["cellsize"]) / 2.0)
        lon_0 = float(metadata["xllcorner"]) + (float(metadata["cellsize"]) / 2.0)
        resolution = float(metadata["cellsize"])
        return {"lat_0": lat_0, "lon_0": lon_0, "res": resolution}

    # eco regions
    path_to_eco_grid = paths["path-to-data-dir"] + "/eco_regions/agro_eco_regions.asc"
    eco_metadata, _ = Mrunlib.read_header(path_to_eco_grid)
    eco_grid = np.loadtxt(path_to_eco_grid, dtype=int, skiprows=6)
    aer_ll0r = get_lat_0_lon_0_resolution_from_grid_metadata(eco_metadata)

    def check_for_nill_dates(mgmt):
        for key, value in mgmt.items():
            if "date" in key and value == "Nill":
                return False
        return True
    
    def mgmt_date_to_rel_date(mgmt_date):
        day_str, month_short_name = mgmt_date.split("-")
        month_str = "00"
        if month_short_name == "Jan":
            month_str = "01"
        elif month_short_name == "Feb":
            month_str = "02"
        elif month_short_name == "Mar":
            month_str = "03"
        elif month_short_name == "Apr":
            month_str = "04"
        elif month_short_name == "May":
            month_str = "05"
        elif month_short_name == "Jun":
            month_str = "06"
        elif month_short_name == "Jul":
            month_str = "07"
        elif month_short_name == "Aug":
            month_str = "08"
        elif month_short_name == "Sep":
            month_str = "09"
        elif month_short_name == "Oct":
            month_str = "10"
        elif month_short_name == "Nov":
            month_str = "11"
        elif month_short_name == "Dec":
            month_str = "12"

        return "0000-" + month_str + "-" + day_str

    # open netcdfs
    path_to_soil_netcdfs = paths["path-to-soil-dir"] + "/" + config["resolution"] + "/"
    if config["resolution"] == "5min":
        soil_data = {
            "sand": {"var": "SAND", "file": "SAND5min.nc", "conv_factor": 0.01},  # % -> fraction
            "clay": {"var": "CLAY", "file": "CLAY5min.nc", "conv_factor": 0.01},  # % -> fraction
            "corg": {"var": "OC", "file": "OC5min.nc", "conv_factor": 0.01},  # scale factor
            "bd": {"var": "BD", "file": "BD5min.nc", "conv_factor": 0.01 * 1000.0},  # scale factor * 1 g/cm3 = 1000 kg/m3
        }
    else:
        soil_data = None  # ["Sand5min.nc", "Clay5min.nc", "OC5min.nc", "BD5min.nc"]
    soil_datasets = {}
    soil_vars = {}
    for elem, data in soil_data.items():
        ds = Dataset(path_to_soil_netcdfs + data["file"], "r", format="NETCDF4")
        soil_datasets[elem] = ds
        soil_vars[elem] = ds.variables[data["var"]]

    def create_soil_profile(row, col):
        # skip first 4.5cm layer and just use 7 layers
        layers = []

        layerDepth = 8
        # find the fill value for the soil data
        for elem2 in soil_data.keys():
            for i in range(8):
                if np.ma.is_masked(soil_vars[elem2][i, row, col]):
                    if i < layerDepth:
                        layerDepth = i
                    break
                    #return None
        layerDepth -= 1

        if layerDepth < 4:
            return None
        
        for i, real_depth_cm, monica_depth_m in [(0, 4.5, 0), (1, 9.1, 0.1), (2, 16.6, 0.1), (3, 28.9, 0.1),
                                                 (4, 49.3, 0.2), (5, 82.9, 0.3), (6, 138.3, 0.6), (7, 229.6, 0.7)][1:]:
            if i <= layerDepth:
                layers.append({
                    "Thickness": [monica_depth_m, "m"],
                    "SoilOrganicCarbon": [soil_vars["corg"][i, row, col] * soil_data["corg"]["conv_factor"], "%"],
                    "SoilBulkDensity": [soil_vars["bd"][i, row, col] * soil_data["bd"]["conv_factor"], "kg m-3"],
                    "Sand": [soil_vars["sand"][i, row, col] * soil_data["sand"]["conv_factor"], "fraction"],
                    "Clay": [soil_vars["clay"][i, row, col] * soil_data["clay"]["conv_factor"], "fraction"]
                })
        return layers

    sent_env_count = 1
    start_time = time.perf_counter()

    # run calculations for each setup
    for _, setup_id in enumerate(run_setups):

        if setup_id not in setups:
            continue
        start_setup_time = time.perf_counter()

        setup = setups[setup_id]
        gcm = setup["gcm"]
        scenario = setup["scenario"]
        ensmem = setup["ensmem"]
        crop = setup["crop"]

        region = setup["region"] if "region" in setup else config["region"]
        lat_lon_bounds = region_to_lat_lon_bounds.get(region, {
            "tl": {"lat": float(config["start_lat"]), "lon": float(config["start_lon"])},
            "br": {"lat": float(config["end_lat"]), "lon": float(config["end_lon"])}
        })

        planting = setup["planting"].lower()
        nitrogen = setup["nitrogen"].lower()
        management_file = f"agro_ecological_regions_{planting}_planting_{nitrogen}_nitrogen.csv"
        # load management data
        management = Mrunlib.read_csv(paths["path-to-data-dir"] + "/eco_regions/" + management_file, key="id")

        # height data for germany
        path_to_dem_grid = paths["path-to-data-dir"] + setup["dem_asc_grid"]
        dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid, 5)
        dem_grid = np.loadtxt(path_to_dem_grid, dtype=float, skiprows=5)
        # print("read: ", path_to_dem_grid)
        dem_ll0r = get_lat_0_lon_0_resolution_from_grid_metadata(dem_metadata)

        # slope data
        path_to_slope_grid = paths["path-to-data-dir"] + "slope_0.009_4326_wgs84_nigeria.asc"
        slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid, 6)
        slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
        print("read: ", path_to_slope_grid)
        slope_ll0r = get_lat_0_lon_0_resolution_from_grid_metadata(slope_metadata)

        # read template sim.json
        with open(setup.get("sim.json", config["sim.json"])) as _:
            sim_json = json.load(_)
        # change start and end date acording to setup
        if setup["start_date"]:
            sim_json["climate.csv-options"]["start-date"] = str(setup["start_date"])
        if setup["end_date"]:
            end_year = int(setup["end_date"].split("-")[0])
            sim_json["climate.csv-options"]["end-date"] = str(setup["end_date"])

            # read template site.json
        with open(setup.get("site.json", config["site.json"])) as _:
            site_json = json.load(_)

        if len(scenario) > 0 and scenario[:3].lower() == "rcp":
            site_json["EnvironmentParameters"]["rcp"] = scenario

        # read template crop.json
        with open(setup.get("crop.json", config["crop.json"])) as _:
            crop_json = json.load(_)
            crop_json["cropRotation"][0]["worksteps"][1]["crop"][2] = crop

        crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = setup[
            "use_vernalisation_fix"] if "use_vernalisation_fix" in setup else False

        # create environment template from json templates
        env_template = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })

        c_lon_0 = -179.75
        c_lat_0 = +89.25
        c_resolution = 0.5

        s_lat_0 = region_to_lat_lon_bounds["earth"][config["resolution"]]["tl"]["lat"]
        s_lon_0 = region_to_lat_lon_bounds["earth"][config["resolution"]]["tl"]["lon"]
        b_lat_0 = lat_lon_bounds["tl"]["lat"]
        b_lon_0 = lat_lon_bounds["tl"]["lon"]

        lats_scaled = range(int(lat_lon_bounds["tl"]["lat"] * s_res_scale_factor),
                            int(lat_lon_bounds["br"]["lat"] * s_res_scale_factor) - 1,
                            -int(s_resolution * s_res_scale_factor))
        no_of_lats = len(lats_scaled)
        s_row_0 = int((s_lat_0 - (lats_scaled[0] / s_res_scale_factor)) / s_resolution)
        for lat_scaled in lats_scaled:
            lat = lat_scaled / s_res_scale_factor

            print(lat, )

            lons_scaled = range(int(lat_lon_bounds["tl"]["lon"] * s_res_scale_factor),
                                int(lat_lon_bounds["br"]["lon"] * s_res_scale_factor) + 1,
                                int(s_resolution * s_res_scale_factor))
            no_of_lons = len(lons_scaled)
            s_col_0 = int(((lons_scaled[0] / s_res_scale_factor) - s_lon_0) / s_resolution)
            for lon_scaled in lons_scaled:
                lon = lon_scaled / s_res_scale_factor
                print(lon, )

                c_col = int((lon - c_lon_0) / c_resolution)
                c_row = int((c_lat_0 - lat) / c_resolution)

                s_col = int((lon - s_lon_0) / s_resolution)
                s_row = int((s_lat_0 - lat) / s_resolution)

                # set management
                aer_col = int((lon - aer_ll0r["lon_0"]) / aer_ll0r["res"])
                aer_row = int((aer_ll0r["lat_0"] - lat) / aer_ll0r["res"])

                aer = 0
                valid_mgmt = False
                if 0 <= aer_row < int(eco_metadata["nrows"]) \
                        and 0 <= aer_col < int(eco_metadata["ncols"]):
                    aer = eco_grid[aer_row, aer_col]
                    if aer > 0 and aer in management:
                        mgmt = management[aer]
                        if check_for_nill_dates(mgmt):
                            valid_mgmt = True
                            for ws in env_template["cropRotation"][0]["worksteps"]:
                                if ws["type"] == "Sowing":
                                    ws["date"] = mgmt_date_to_rel_date(mgmt["Sowing date"])
                                    ws["PlantDensity"] = [float(mgmt["Planting density"]), "plants/m2"]
                                elif ws["type"] == "AutomaticHarvest":
                                    ws["latest-date"] = mgmt_date_to_rel_date(mgmt["Harvest date"])
                                elif ws["type"] == "Tillage":
                                    ws["date"] = mgmt_date_to_rel_date(mgmt["Tillage date"])
                                elif ws["type"] == "MineralFertilization":
                                    app_no = int(ws["application"])
                                    app_str = str(app_no) + ["st", "nd", "rd", "th"][app_no - 1]
                                    ws["date"] = mgmt_date_to_rel_date(mgmt[f"N {app_str} date"])
                                    ws["amount"] = [float(mgmt[f"N {app_str} application (kg/ha)"]), "kg"]

                soil_profile = create_soil_profile(s_row, s_col)
                if soil_profile is None or aer == 0 or aer not in management or not valid_mgmt:
                    env_template["customId"] = {
                        "setup_id": setup_id,
                        "lat": lat, "lon": lon,
                        "b_lat_0": b_lat_0, "b_lon_0": b_lon_0,
                        "s_resolution": s_resolution,
                        "s_row": s_row, "s_col": s_col,
                        "s_row_0": s_row_0, "s_col_0": s_col_0,
                        "no_of_s_cols": no_of_lons, "no_of_s_rows": no_of_lats,
                        "c_row": int(c_row), "c_col": int(c_col),
                        "env_id": sent_env_count,
                        "planting": planting,
                        "nitrogen": nitrogen,
                        "region": region,
                        "crop": crop,
                        "nodata": True
                    }
                    socket.send_json(env_template)
                    print("sent nodata env ", sent_env_count, " customId: ", env_template["customId"])
                    sent_env_count += 1
                    continue

                dem_col = int((lon - dem_ll0r["lon_0"]) / dem_ll0r["res"])
                dem_row = int((dem_ll0r["lat_0"] - lat) / dem_ll0r["res"])
                height_nn = dem_grid[dem_row, dem_col]

                slope_col = int((lon - slope_ll0r["lon_0"]) / slope_ll0r["res"])
                slope_row = int((slope_ll0r["lat_0"] - lat) / slope_ll0r["res"])
                slope = slope_grid[slope_row, slope_col]

                env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup[
                    "LeafExtensionModifier"]

                env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                if setup["elevation"]:
                    env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

                if setup["slope"]:
                    env_template["params"]["siteParameters"]["slope"] = slope / 90.0

                if setup["latitude"]:
                    env_template["params"]["siteParameters"]["Latitude"] = lat

                if setup["FieldConditionModifier"]:
                    fcms = setup["FieldConditionModifier"].split("|")
                    fcm = float(fcms[aer-1])
                    if fcm > 0:
                        env_template["cropRotation"][0]["worksteps"][1]["crop"]["cropParams"]["species"][
                            "FieldConditionModifier"] = fcm

                env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup[
                    "fertilization"]
                env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]
                env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
                env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup[
                    "WaterDeficitResponseOn"]
                env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup[
                    "EmergenceMoistureControlOn"]
                env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup[
                    "EmergenceFloodingControlOn"]

                env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
                hist_sub_path = "isimip/3b_v1.1_CMIP6/csvs/{gcm}/historical/{ensmem}/row-{crow}/col-{ccol}.csv.gz".format(
                    gcm=gcm, ensmem=ensmem, crow=c_row, ccol=c_col)
                sub_path = "isimip/3b_v1.1_CMIP6/csvs/{gcm}/{scenario}/{ensmem}/row-{crow}/col-{ccol}.csv.gz".format(
                    gcm=gcm, scenario=scenario, ensmem=ensmem, crow=c_row, ccol=c_col
                )
                if setup["incl_historical"] and scenario != "historical":
                    climate_data_paths = [
                        paths["monica-path-to-climate-dir"] + hist_sub_path,
                        paths["monica-path-to-climate-dir"] + sub_path
                    ]
                else:
                    climate_data_paths = [paths["monica-path-to-climate-dir"] + sub_path]
                env_template["pathToClimateCSV"] = climate_data_paths
                print("pathToClimateCSV:", env_template["pathToClimateCSV"])

                env_template["customId"] = {
                    "setup_id": setup_id,
                    "lat": lat, "lon": lon,
                    "s_lat_0": s_lat_0, "s_lon_0": s_lon_0,
                    "s_resolution": s_resolution,
                    "s_row": s_row, "s_col": s_col,
                    "s_row_0": s_row_0, "s_col_0": s_col_0,
                    "no_of_s_cols": no_of_lons, "no_of_s_rows": no_of_lats,
                    "c_row": int(c_row), "c_col": int(c_col),
                    "env_id": sent_env_count,
                    "planting": planting,
                    "nitrogen": nitrogen,
                    "region": region,
                    "crop": crop,
                    "nodata": False
                }

                socket.send_json(env_template)
                print("sent env ", sent_env_count, " customId: ", env_template["customId"])

                sent_env_count += 1

        stop_setup_time = time.perf_counter()
        print("Setup ", (sent_env_count - 1), " envs took ", (stop_setup_time - start_setup_time), " seconds")

    stop_time = time.perf_counter()

    # write summary of used json files
    try:
        print("sending ", (sent_env_count - 1), " envs took ", (stop_time - start_time), " seconds")
        print("exiting run_producer()")
    except Exception:
        raise


if __name__ == "__main__":
    run_producer()
