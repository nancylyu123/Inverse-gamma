
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 10:05:46 2017

@author: yulit
"""


from __future__ import division
# import pylab
import os
import numpy
import scipy.ndimage
import logging
import sys
import datetime


ROOT_PATH = 'C:/directory'
DO_SHIFTS = False
MEASURED_TOTAL_SKIP_HEADER = 278
MEASURED_TOTAL_SKIP_FOOTER = 183
MEASURED_X_AXIS_SKIP_HEADER = 320
MEASURED_X_AXIS_SKIP_FOOTER = 181
CALCULATED_X_AXIS_SKIP_HEADER = 5
CALCULATED_X_AXIS_SKIP_FOOTER = 220
ZOOMING_RATIO = 10


# Do spline interpolation of calculate_data so that it's resolution is 0.1 mm instead of 1 mm
# Perform gamma for every measurement point.
# It finds every non zero point in measured data and searches the coresponding dose point
# In calculated dataset with a mask of size (dis*2+1). So if dis=2mm, the mask is a 5x5 numpy array.
def gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, distance=2, percent=0.02, threshold=0.05):
    logging.info('Calculating gamma, distance=%f, percent=%f, threshold=%f', distance, percent, threshold)
    measured_max = measure_data.max()
    calculated_max = calculate_data.max()  # max dose point of calculated data
    count_total = 0
    count_pass = 0
    gamma_image = numpy.zeros(shape=(measure_data.shape[0], measure_data.shape[1]))  # Create a gamma map with zeros
    for measure_row in range(0, measure_Y_axis.shape[0]):
        for measure_column in range(0, measure_X_axis.shape[0]):
            if measure_data[measure_row, measure_column] != 0 and measure_data[measure_row, measure_column] / measured_max >= threshold:
                calculate_row = numpy.where(calculate_Y_axis == measure_Y_axis[measure_row] * 10)[0][0]
                calculate_column = numpy.where(calculate_X_axis == measure_X_axis[measure_column] * 10)[0][0]
                if calculate_row - distance >= 0 and calculate_column - distance >= 0 and calculate_row + distance <= calculate_data.shape[0] - 1 and calculate_column + distance <= calculate_data.shape[1] - 1:
                    count_total = count_total + 1
                    # offset_x < 0 or offset_y < 0 should be excluded from preceding 'if' clause
                    offset_x = int((calculate_row - distance) * ZOOMING_RATIO)
                    offset_y = int((calculate_column - distance) * ZOOMING_RATIO)

                    sample = measure_data[measure_row, measure_column]
                    mask_size = int(distance * ZOOMING_RATIO * 2 + 1)
                    mask = numpy.zeros(shape=(mask_size, mask_size))
                    for x in numpy.arange(0, mask_size):
                        for y in numpy.arange(0, mask_size):
                            reference = interpolated_data[offset_x + x, offset_y + y]
                            if distance > 0:
                                radius = numpy.sqrt((x - distance * ZOOMING_RATIO)**2+(y - distance * ZOOMING_RATIO)**2) / ZOOMING_RATIO  # distance to central point
                                mask[x, y] = numpy.sqrt(((sample - reference) / calculated_max / percent)**2 + (radius / distance)**2)
                            elif distance == 0:
                                mask[x, y] = numpy.sqrt(((sample - reference) / calculated_max / percent)**2)

                    gamma_image[measure_row, measure_column] = mask.min()
                    logging.debug('Calculated row: %d, column: %d in measured data', measure_row, measure_column)
                    if gamma_image[measure_row, measure_column] <= 1:
                        count_pass = count_pass + 1
    return count_pass, count_total, gamma_image


def interpolation(calculate_data):
    size_x = calculate_data.shape[0]
    size_y = calculate_data.shape[1]
    zoom_x = ((size_x - 1) * ZOOMING_RATIO + 1) / size_x
    zoom_y = ((size_y - 1) * ZOOMING_RATIO + 1) / size_y
    return scipy.ndimage.interpolation.zoom(calculate_data, [zoom_x, zoom_y])


def inverse_gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, threshold=1, dose_distance_ratio=1):
    logging.info('Running inverse_gamma')
    for criteria in numpy.arange(1, 10, 0.1):
        distance_criteria = criteria
        percent_dose = criteria*dose_distance_ratio / 100
        pass_count, total_count, gamma_image = gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, distance_criteria, percent_dose)
        if pass_count / total_count >= threshold:
            return criteria


def inverse_gamma_fixdta(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, threshold=1, distance_criteria=1):
    logging.info('Running inverse_gamma_fixdta')
    for criteria in numpy.arange(1, 50, 0.1):
        percent_dose = criteria/100
        pass_count, total_count, gamma_image = gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, distance_criteria, percent_dose)
        if pass_count / total_count >= threshold:
            return criteria


def inverse_gamma_fixdd(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, threshold=1, percent_dose=0.03):
    logging.info('Running inverse_gamma_fixdd')
    for criteria in numpy.arange(0, 20, 0.1):
        distance_criteria = criteria
        pass_count, total_count, gamma_image = gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data, distance_criteria, percent_dose)
        if pass_count / total_count >= threshold:
            return criteria


def main(argv):
    logging.getLogger(None).setLevel(logging.INFO)

    logging.info('Entering main() function')
    file = open("output.txt", "w")
    root_path = ROOT_PATH
    for path in os.listdir(root_path):
        sub_path = root_path + '/' + path
        if os.path.isdir(sub_path):
            num_arcs = 0
            for file_path in os.listdir(sub_path):
                if file_path.startswith('arc') and file_path.endswith('.txt'):
                    num_arcs += 1
            file.write(path + "\n")
            file.write('Number of arcs: ' + str(num_arcs) + "\n")
            logging.info('Patient name: %s', path)
            logging.info('Number of arcs: %d', num_arcs)
            logging.info('Time starting processing: {:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now()))
            for arc_number in range(1, 1 + num_arcs):
                file.write('Arc ' + str(arc_number) + "\n")
                #            print('Arc ' + str(arc_number)) #print arc number
                measured_file = sub_path + '/arc' + str(arc_number) + '.txt'
                calculated_file = sub_path + '/calc' + str(arc_number) + '.snc'
                logging.info('Processing measured file %s', measured_file)
                logging.info('Processing calculated file %s', calculated_file)
                measure_raw_data = numpy.genfromtxt(measured_file, skip_header=MEASURED_TOTAL_SKIP_HEADER, skip_footer=MEASURED_TOTAL_SKIP_FOOTER)
                # Import the X and Y coordinates of measured data
                measure_X_axis = (numpy.genfromtxt(measured_file, skip_header=MEASURED_X_AXIS_SKIP_HEADER, skip_footer=MEASURED_X_AXIS_SKIP_FOOTER))[1:]
                measure_Y_axis = measure_raw_data[:, 0]
                # This is the measured dataset only without coordinates
                measure_data = measure_raw_data[:,2:]
                # Acquire array from snc file
                calculate_raw_data = numpy.genfromtxt(calculated_file, skip_header=6)
                # Import X and Y coordinates of calculated data
                calculate_X_axis = (numpy.genfromtxt(calculated_file, skip_header=CALCULATED_X_AXIS_SKIP_HEADER, skip_footer=CALCULATED_X_AXIS_SKIP_FOOTER))[1:]
                calculate_Y_axis = calculate_raw_data[:, 0]
                # This is the calculated dataset only without coordinates
                calculate_data = calculate_raw_data[:, 1:]
                interpolated_data = interpolation(calculate_data)
                no_shift_pass, no_shift_total, no_shift_gamma_image = gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data)
                # print gamma results
                file.write("Gamma pass rate: " + str(no_shift_pass / no_shift_total) + "\n")
                file.write("Points passed: " + str(no_shift_pass) + "\n")
                file.write("Points failed: " + str(no_shift_total - no_shift_pass) + "\n")
                str_inverse_gamma = str(inverse_gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data))
                str_inverse_gamma_fixdta = str(inverse_gamma_fixdta(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data))
                str_inverse_gamma_fixdd = str(inverse_gamma_fixdd(measure_data, measure_X_axis, measure_Y_axis, calculate_data, calculate_X_axis, calculate_Y_axis, interpolated_data))
                interpolated_data = None
                if DO_SHIFTS:
                    best_passrate = no_shift_pass/no_shift_total
                    best_shift_pass = 0
                    best_shift_total = 0
                    best_distance = 0
                    best_passrate_shift = 0, 0
                    best_locations = [(0,0)]
                    # Perform shifts
                    for shift_x in numpy.arange(-2, 3, 1):
                        for shift_y in numpy.arange(-2, 3, 1):
                            logging.info('Calculating gamma for shift_x=%d and shift_y=%d', shift_x, shift_y)
                            calculate_shift = scipy.ndimage.interpolation.shift(calculate_data, (shift_y, shift_x))
                            interpolated_data = interpolation(calculate_shift)
                            shift_pass, shift_total, shift_gamma_image = gamma(measure_data, measure_X_axis, measure_Y_axis, calculate_shift, calculate_X_axis, calculate_Y_axis, interpolated_data)
                            if shift_pass/shift_total > best_passrate:
                                best_passrate = shift_pass/shift_total
                                best_distance = numpy.sqrt(shift_x**2+shift_y**2)
                                best_shift_pass = shift_pass
                                best_shift_total = shift_total
                                best_passrate_shift = shift_x, shift_y
                                best_locations = []
                                best_locations.append((shift_x, shift_y))
                            elif shift_pass/shift_total == best_passrate:
                                current_distance = numpy.sqrt(shift_x**2+shift_y**2)
                                if current_distance < best_distance:
                                    best_distance = current_distance
                                    best_passrate_shift = shift_x, shift_y
                                best_locations.append((shift_x, shift_y))
                    file.write("Best gamma pass rate: " + str(best_passrate) + "\n")
                    file.write("Points passed: " + str(best_shift_pass) + "\n")
                    file.write("Points failed: " + str(best_shift_total-best_shift_pass) + "\n")
                    file.write("Best gamma pass rate shift: "+str(best_passrate_shift)+"\n")
                    file.write("Best locations: " + str(best_locations) + "\n")
                file.write("Dose distance ratio: " + str_inverse_gamma + "\n")
                file.write(str_inverse_gamma_fixdta + "%, 1 mm\n")
                file.write("3% " + str_inverse_gamma_fixdd + "mm\n")
                file.write("============================================================" + "\n")
            logging.info('Time finishing processing: {:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now()))
    file.close()


if __name__ == "__main__":
    main(sys.argv)
