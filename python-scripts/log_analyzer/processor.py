import collections
import re
import csv
import copy
import numpy as np
import os
import sys
from strings import *
from util import *

class Processor:

    def __init__(self):
        self.HARDEN_COL = 5
        self.LINENUM_COL = 1
        self.TOTAL_COL = 6
        self.current_section = 0
        self.current_line_number = -1
        self.current_token_number = 1
        self.total_initial_tripped = 0
        self.curr_line = 0
        self.total_tripped = 0
        self.isReinforced = -1
        self.total_shedding_load = 0.0
        self.total_tripped_generators = 0
        self.plan_length = 0
        self.total_replications = 0
        self.hardening_plans = { }
        self.line_map = { }
        self.matrix_map = { }
        self.config_values = { }
        self.csv_rows = []
        self.plans = [ ]
        self.encoded_matrix = [ ]

    def process(self, token_collection, outfiles, infiles, config_values, flag1):
        """
        Main processing function for entire analyzer. Orchaestrates smaller
        processing functions such as csv_process() and matrix_process().
        """
        flag2 = 1
        token_collection_list = []
        self.config_values = config_values
        if self.config_values[NUMLINES] == "":
            print("You forgot to specify the number of lines in the config file. Please try again.")
            sys.exit(0)

        for token in token_collection:
            token_collection_list.append(token)

        # if infiles:
        #     self.build_line_map(infiles[1])
        #     #self.print_line_map()

        if flag2 == 1:
            self.csv_process(token_collection_list, os.path.join('output_files', outfiles[0]), flag1)
            self.matrix_process(token_collection_list, os.path.join('output_files', outfiles[1]))

    def csv_process(self, token_collection, filename, flag1):
        """
        Main processing function for generating CSV output. Orchaestrates smaller
        processing functions for each log section (currently five + header)
        """
        current_replication = 0

        for token in token_collection:
            if flag1:
                print(token)
            if token.type == 'HEADER':
                self.csv_process_header(token, current_replication)
            elif token.type == 'COLUMN':
                pass
            elif token.type == 'DATA':

                # # # # # # # # # # # # #
                # Process Each Section  #
                # # # # # # # # # # # # #

                # Section 1
                if self.current_section == 1:                        # Initial Tripped
                    self.csv_process_section_one(token)
                # Section 2
                elif self.current_section == 2:
                    self.csv_process_section_two(token)
                # Section 3
                elif self.current_section == 3:
                    self.csv_process_section_three(token)
                # Section 4
                elif self.current_section == 4:
                    self.csv_process_section_four()
                # Section 5
                elif self.current_section == 5:
                    current_replication = self.csv_process_section_five(token)
            else:
                print("Nothing")

        self.store_and_reset(current_replication)

        #print(self.total_replications)
        #print(len(self.line_map))
        self.write_to_csv(filename)

    def csv_process_header(self, token, current_replication):
        """
        Function used to determine how to process header tokens for csv file
        output.
        """
        self.isReinforced = -1
        self.current_section = self.get_section(token)
        if self.current_section == 5:
            if len(self.plans) > 0:
                self.store_and_reset(current_replication)

    def csv_process_section_one(self, token):
        """
        Function used to determine how to process section one data tokens for
        csv file output. Section one tells us which lines are reinforced and
        which lines were inititally tripped.
        """
        if self.current_token_number > self.TOTAL_COL:
            self.current_token_number = 1
        if self.current_token_number == 1:               # LineNum token
            self.current_line_number = int(token.value)  # Save LineNum
            self.current_token_number += 1
        elif self.current_token_number == 5:             # Check if Hardened
            self.isReinforced = int(token.value)
            if self.isReinforced == 1:                   # If Line Hardened
                self.plans.append(self.current_line_number)   # Add to plan
            else:
                self.total_initial_tripped += 1          # Otherwise failed
            self.current_token_number += 1
        else:
            self.current_token_number += 1               # No match

    def csv_process_section_two(self, token):
        """
        Function used to determine how to process section two data tokens for
        csv file output. Section two tells us which lines are tripped after
        cascading power failure is complete. For csv output, we simply need to
        know how many.
        """
        next_line = int(token.line)
        if next_line > self.curr_line:
            self.curr_line = next_line
            self.total_tripped += 1

    def csv_process_section_three(self, token):
        """
        Function used to determine how to process section three data tokens for
        csv file output. Section three tell us shedding load for a given bus. We
        use this information to calculate total load shedding.
        """
        if self.current_token_number > 4:
            self.current_token_number = 1
        if self.current_token_number == 4:
            self.total_shedding_load += float(token.value)
        self.current_token_number += 1

    def csv_process_section_four(self):
        """
        Function used to determine how to process section four data tokens for
        csv file output. Section four tells us which generators tripped. For csv
        processing, we simply want to know how many.
        """
        self.total_tripped_generators += 1

    def csv_process_section_five(self, token):
        """
        Function used to determine how to process section five data tokens for
        csv file output. Section five tells us the replication index.
        """
        current_replication = int(token.value)
        if current_replication > self.total_replications:
            self.total_replications = current_replication
        return current_replication

    def matrix_process(self,token_collection, outfile):
        """
        Main processing function for generating sparse matrix output.
        Orchaestrates smaller processing functions for each log section (currently five + header).
        Currently, code is not broken up into smaller functions. That is on todo list.
        """
        self.build_matrix_map()
        numlines = int(self.config_values[NUMLINES])
        current_replication = 0
        rep_count = 0
        current_plan = 0
        keys = list(self.hardening_plans.keys())
        row = len(self.hardening_plans)*self.total_replications
        col = 3*numlines
        self.build_matrix(row, col+1)
        #print(len(keys))

        count = 1

        for token in token_collection:
            #print("HERE")
            #print(token)
            if token.type == 'HEADER':
                self.isReinforced = -1
                self.current_section = self.get_section(token)
                if self.current_section == 5:
                    self.reset()
                    rep_count += 1
                    if rep_count == self.total_replications:
                        current_plan += 1
                        rep_count = 0
                        #print(current_plan)
            elif token.type == 'COLUMN':
                pass
            elif token.type == 'DATA':

                # # # # # # # # # # # # #
                # Process Each Section  #
                # # # # # # # # # # # # #

                # Section 1
                if self.current_section == 1:                        # Initial Tripped
                    if self.current_token_number > self.TOTAL_COL:
                        self.current_token_number = 1
                    if self.current_token_number == 1:               # LineNum token
                        self.current_line_number = int(token.value)  # Save LineNum
                        self.current_token_number += 1
                        # print(current_plan)
                        key = tuple([keys[current_plan-1], current_replication])
                        # print(key)
                        #print(self.current_line_number+numlines)
                        self.update_matrix(self.matrix_map[key], 0, current_replication)
                        self.update_matrix(self.matrix_map[key], self.current_line_number+numlines, 1)
                    elif self.current_token_number == 5:             # Check if Hardened
                        self.isReinforced = int(token.value)
                        if self.isReinforced == 1:
                            line = 0
                            key = tuple([keys[current_plan-1], current_replication])
                            #print(self.matrix_map[key])
                            self.update_matrix(self.matrix_map[key], self.current_line_number, 1)
                        self.current_token_number += 1

                    else:
                        self.current_token_number += 1               # No match
                # Section 2
                elif self.current_section == 2:
                    if count == 1:
                        #print("LAST: " + str(int(token.value)+(numlines*2)))
                        self.update_matrix(self.matrix_map[key], int(token.value)+(numlines*2), 1)
                        count += 1
                    else:
                        if count == 3:
                            count = 1
                        else:
                            count += 1
                elif self.current_section == 5:
                    current_replication = int(token.value)

        self.reset()
        if str(self.config_values[DELIMITER]) == "sp":
            self.print_matrix(outfile)
        elif str(self.config_values[DELIMITER]) == ",":
            self.print_matrix_csv("output_files/matrix.csv", numlines)
        else:
            self.print_matrix(outfile)
            self.print_matrix_csv("output_files/matrix.csv", numlines)

    def build_matrix_map(self):
        """
        Helper function used to build up the matrix map. Given a particular
        hardening plan, and a replication count, it will map this key to a row
        number in the matrix.
        """
        #((1, 3), 1) : 1
        total_count = 0
        temp_count = 1
        plan = 0
        length = len(self.hardening_plans)
        total = length * self.total_replications
        keys = list(self.hardening_plans.keys())
        # key = tuple([keys[0], 1])
        # print(key)

        while total_count < total:
            if temp_count == self.total_replications+1:
                plan += 1
                temp_count = 1
            key = tuple([keys[plan], temp_count])
            self.matrix_map[key] = total_count
            total_count += 1
            temp_count += 1

        #print(self.matrix_map)

    def update_matrix(self, row, col, value):
        self.encoded_matrix[row][col] = value

    def print_matrix_csv(self, outfile, numlines):
        """
        Writes the csv header and csv_rows list to a csv file. This is the final
        output.
        """
        with open(outfile, 'w', newline='') as csvfile:
            filewriter = csv.writer(csvfile)
            column_header = self.build_matrixcsv_header(numlines)
            filewriter.writerow(column_header)
            for row in self.encoded_matrix:
                filewriter.writerow(row)

    def build_matrixcsv_header(self, numlines):
        """
        Builds the csv_header based on output requirements.
        """
        col_list = []
        matrix = self.encoded_matrix.sort()
        col_list.append("Replication")
        for i in range(numlines):
            col_list.append("HLI_"+str(i+1))

        for i in range(numlines):
            col_list.append("ITL_"+str(i+1))

        for i in range(numlines):
            col_list.append("TTL_"+str(i+1))

        return col_list


    def print_matrix(self, outfile):
        with open(outfile, "w") as f:
            for row in self.encoded_matrix:
                listline = str(row)
                listline = listline.replace("[", "")
                listline = listline.replace("]", "\n")
                templine = listline.replace(",", " ")
                f.write(templine)

    def build_matrix(self, row, col):
        self.encoded_matrix = [[0 for x in range(col)] for y in range(row)]

    def store_and_reset(self, replication):
        """
        Helper function used to store a hardening plan into the dictionary and
        reset all values required for the new processing cycle.
        """
        plan_key = self.add_plan(self.plans, replication)
        if self.plan_length < len(plan_key):
            self.plan_length = len(plan_key)
        self.add_to_csvrows(self.hardening_plans[plan_key])
        self.reset()

    def reset(self):
        """
        Resets valus to their default states.
        """
        self.total_initial_tripped = 0
        self.total_tripped = 0
        self.total_tripped_generators = 0
        self.current_line_number = 0
        self.current_token_number = 1
        self.total_shedding_load = 0.0
        self.isReinforced = 0
        self.plans = [ ]

    def add_to_csvrows(self, replication_index):
        """
        Builds a list of csv rows so that we only need to do one set of writes
        at the end of the csv processing cycle.
        """
        row = []
        for hline in self.plans:
            row.append(hline)

        row.append(str(self.total_initial_tripped))
        row.append(str(self.total_tripped))
        row.append(str(self.total_shedding_load))
        row.append(str(self.total_tripped_generators))
        row.append(str(replication_index))

        self.csv_rows.append(row)

    def write_to_csv(self, filename):
        """
        Writes the csv header and csv_rows list to a csv file. This is the final
        output.
        """
        with open(filename, 'w', newline='') as csvfile:
            filewriter = csv.writer(csvfile)
            column_header = self.build_csv_header()
            filewriter.writerow(column_header)
            for row in self.csv_rows:
                filewriter.writerow(row)

    def build_csv_header(self):
        """
        Builds the csv_header based on output requirements.
        """
        col_list = []
        for i in range(self.plan_length):
            col_list.append('LineNum')

        col_list.append('Total Initial Tripped Lines')
        col_list.append('Total Tripped Lines')
        col_list.append('Total Shedding Load Amount')
        col_list.append('Total Tripped Generators')
        col_list.append('Replication Index')

        return col_list

    def add_plan(self, plan, replication):
        """
        Deprecated. Convert plan into a tuple and store it in the dictionary along with its
        replication index.
        """
        plan_key = tuple(plan)
        self.hardening_plans[plan_key] = replication
        return plan_key

    def get_section(self, token):
        """
        Get section header string from the section_map in strings.py. Used for
        processing header tokens.
        """
        myregex = '=+'
        line = str(token.value)
        section_text = re.sub(myregex, '', line)
        section_text = section_text.strip()
        return section_table[section_text]

    def build_line_map(self, infile):
        """
        Deprecated. Reads in a line map from a given file and builds a line mapping
        of power lines using the To and From buses as keys into the hashmap.
        """
        temp_key = [ ]

        with open(infile, "r") as csvfile:
            linereader = csv.reader(csvfile, delimiter=',')
            for line in linereader:
                temp_key.append(line[1])                # Get "From" bus
                temp_key.append(line[2])                # Get "To" bus
                key = tuple(temp_key)                   # Use "From" and "To" as hashmap key
                value = line[0]                         # Use line num as value in hashmap
                self.line_map[key] = value
                temp_key = []

    def print_line_map(self):
        for entry in self.line_map:
            print(str(entry) + " = " + str(self.line_map[entry]))
