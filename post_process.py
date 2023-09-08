#!/usr/bin/python
# -*- coding: UTF-8


# read folder with ascii files
# merge every 10 years into one file
# by using 1. average
#          2. standard deviation

import os
import numpy as np
from dataclasses import dataclass
import gzip
import math
import sys
import errno

# read folder

PATHS = {
    "local": {
        "sourcepath": "./out/3_reg-nigeria_plant-early_high-N/",
        "outputpath": "./out/3_reg-nigeria_plant-early_high-N/merged/",
        "std": "std/",  # path to std images
        "avg": "avg/",  # path to avg images
    },
    "cluster": {
        "sourcepath": "/source/",
        "outputpath": "/out/merged/",
        "std": "std/",  # path to std images
        "avg": "avg/",  # path to avg images
    }
}
USER = "local"
NONEVALUE = -9999
# start and end years for each decade
start_year = 1970  # 1601  # 1970
end_year = 2020  # 2101  # 2009
yearRanges = [(y, y+9) for y in range(start_year, end_year, 10)]
crops = ["maizegrainmaize"]  # TBD use grain maize
types = ["Yield", "LAI", "AbBiom"]
# e.g. maizesilagemaize_Yield_1987_18.asc 
fileTemplateInput = "{0}_{1}_{2}_{3}.asc"  # 0 = crop, 1 = type, 2 = year, 3 = index
fileTemplateOutAvg = "{0}_{1}_{2}_{3}_avg.asc"  # 0 = crop, 1 = type, 2 = yearstart, 3 = yearend
fileTemplateOutStd = "{0}_{1}_{2}_{3}_std.asc"  # 0 = crop, 1 = type, 2 = yearstart, 3 = yearend


def main():
    # read cmd line args
    pathId = USER
    sourceFolder = ""
    outputFolder = ""
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k == "path":
                pathId = v
            if k == "source":
                sourceFolder = v
            if k == "out":
                outputFolder = v

    if not sourceFolder:
        sourceFolder = PATHS[pathId]["sourcepath"]
    if not outputFolder:
        outputFolder = PATHS[pathId]["outputpath"]

    stdFolder = os.path.join(outputFolder, PATHS[pathId]["std"])
    avgFolder = os.path.join(outputFolder, PATHS[pathId]["avg"])

    makeDir(stdFolder)
    makeDir(avgFolder)

    # read files per decade
    indexOffset = 1
    for yearRange in yearRanges:
        for crop in crops:
            for type in types:
                roundToDecimal = 0
                if type == "LAI":
                    roundToDecimal = 2
                readFilesPerDecade(crop, type, roundToDecimal, yearRange, indexOffset, sourceFolder, stdFolder,
                                   avgFolder)
        indexOffset += 10


def readFilesPerDecade(crop, type, roundToDecimal, yearRange, indexOffset, sourcepath, stdFolder, avgFolder):
    # read files per decade
    files = []
    index = indexOffset
    for year in range(yearRange[0], yearRange[1] + 1):
        file = fileTemplateInput.format(crop, type, year, index)
        index += 1
        # join path and file
        file = os.path.join(sourcepath, file)

        grid = readFile(file)
        files.append(grid)
    # calculate average
    avg = np.nanmean(files, axis=0)
    # calculate standard deviation
    std = np.nanstd(files, axis=0)
    # round values to int
    avg = np.round(avg, decimals=roundToDecimal)
    std = np.round(std, decimals=roundToDecimal)

    fmtStr = '%1.' + str(roundToDecimal) + 'f'

    # set nan values to -9999
    avg[np.isnan(avg)] = NONEVALUE
    std[np.isnan(std)] = NONEVALUE
    # write files
    fileAvg = fileTemplateOutAvg.format(crop, type, yearRange[0], yearRange[1])
    fileAvg = os.path.join(avgFolder, fileAvg)
    # writeAsciiHeader(fileAvg, readAsciiHeader(file))
    # write file with acuracy of 0 decimals (int)
    np.savetxt(fileAvg, avg, fmt=fmtStr, delimiter=' ', newline='\n', header=asciiHeaderString(readAsciiHeader(file)),
               footer='', comments='# ', encoding=None)
    fileStd = fileTemplateOutStd.format(crop, type, yearRange[0], yearRange[1])
    fileStd = os.path.join(stdFolder, fileStd)
    # writeAsciiHeader(fileStd, readAsciiHeader(file))

    np.savetxt(fileStd, std, fmt=fmtStr, delimiter=' ', newline='\n', header=asciiHeaderString(readAsciiHeader(file)),
               footer='', comments='# ', encoding=None)


def readFile(file):
    print("File:", file)
    header = readAsciiHeader(file)
    ascii_data_array = np.loadtxt(header.ascii_path, dtype=float, skiprows=6)
    # Set the nodata values to nan
    ascii_data_array[ascii_data_array == header.ascii_nodata] = np.nan
    ascii_data_array[ascii_data_array == 0] = np.nan

    return ascii_data_array


def makeDir(dir):
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


@dataclass
class AsciiHeader:
    ascii_path: str
    ascci_cols: int
    ascii_rows: int
    ascii_xll: float
    ascii_yll: float
    ascii_cs: float
    ascii_nodata: float
    image_extent: list


def readAsciiHeader(ascii_path):
    if ascii_path.endswith(".gz"):
        # Read in ascii header data
        with gzip.open(ascii_path, 'rt') as source:
            ascii_header = source.readlines()[:6]
    else:
        # Read in ascii header data
        with open(ascii_path, 'r') as source:
            ascii_header = source.readlines()[:6]

    # Read the ASCII raster header
    ascii_header = [item.strip().split()[-1] for item in ascii_header]
    ascci_cols = int(ascii_header[0])
    ascii_rows = int(ascii_header[1])
    ascii_xll = float(ascii_header[2])
    ascii_yll = float(ascii_header[3])
    ascii_cs = float(ascii_header[4])
    ascii_nodata = float(ascii_header[5])

    image_extent = [
        ascii_xll, ascii_xll + ascci_cols * ascii_cs,
        ascii_yll, ascii_yll + ascii_rows * ascii_cs]

    return AsciiHeader(ascii_path, ascci_cols, ascii_rows, ascii_xll, ascii_yll, ascii_cs, ascii_nodata, image_extent)


def asciiHeaderString(header):
    return "ncols         {0}\n".format(header.ascci_cols) + \
        "nrows         {0}\n".format(header.ascii_rows) + \
        "xllcorner     {0}\n".format(header.ascii_xll) + \
        "yllcorner     {0}\n".format(header.ascii_yll) + \
        "cellsize      {0}\n".format(header.ascii_cs) + \
        "NODATA_value  {0}\n".format(header.ascii_nodata)


def writeAsciiHeader(file, header):
    with open(file, 'w') as f:
        f.write("ncols         {0}\n".format(header.ascci_cols))
        f.write("nrows         {0}\n".format(header.ascii_rows))
        f.write("xllcorner     {0}\n".format(header.ascii_xll))
        f.write("yllcorner     {0}\n".format(header.ascii_yll))
        f.write("cellsize      {0}\n".format(header.ascii_cs))
        f.write("NODATA_value  {0}\n".format(header.ascii_nodata))


if __name__ == "__main__":
    main()
