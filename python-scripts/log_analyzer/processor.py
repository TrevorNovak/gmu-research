import collections
import re
import csv
import copy
import numpy as np
from strings import *

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
        self.csv_rows = []
        self.plans = [ ]
        self.encoded_matrix = [ ]

    def process(self, token_collection, outfile, infile, flag1):
        flag2 = 1
        token_collection_list = []

        for token in token_collection:
            token_collection_list.append(token)

        if infile:
            self.build_line_map(infile)
            #self.print_line_map()


        if flag2 == 1:
            self.first_process(token_collection_list, outfile, flag1)
            self.matrix_process(token_collection_list)

    def build_line_map(self, infile):
        temp_key = [ ]

        with open(infile, "r") as csvfile:
            linereader = csv.reader(csvfile, delimiter=',')
            for line in linereader:
                temp_key.append(line[1])                # Get "From" bus
                temp_key.append(line[2])                # Get "To" bus
                key = tuple(temp_key)                   # Use "From" and "To" as hashmap key
                value = line[0]                         # Use line num as value in hashmap
                if key in self.line_map:
                    print("Current Line: " + str(value) +  " Line: " + str(self.line_map[key])+ " Key: " + str(key))

                self.line_map[key] = value
                temp_key = []

    def print_line_map(self):
        for entry in self.line_map:
            print(str(entry) + " = " + str(self.line_map[entry]))

    def first_process(self, token_collection, filename, flag1):
        current_replication = 0

        for token in token_collection:
            if flag1:
                print(token)
            if token.type == 'HEADER':
                self.isReinforced = -1
                self.current_section = self.get_section(token)
                if self.current_section == 5:
                    if len(self.plans) > 0:
                        self.store_and_reset(current_replication)
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
                    elif self.current_token_number == 5:             # Check if Hardened
                        self.isReinforced = int(token.value)
                        if self.isReinforced == 1:                   # If Line Hardened
                            self.plans.append(self.current_line_number)   # Add to plan
                        else:
                            self.total_initial_tripped += 1          # Otherwise failed
                        self.current_token_number += 1
                    else:
                        self.current_token_number += 1               # No match
                # Section 2
                elif self.current_section == 2:
                    next_line = int(token.line)
                    if next_line > self.curr_line:
                        self.curr_line = next_line
                        self.total_tripped += 1
                # Section 3
                elif self.current_section == 3:
                    if self.current_token_number > 4:
                        self.current_token_number = 1
                    if self.current_token_number == 4:
                        self.total_shedding_load += float(token.value)
                    self.current_token_number += 1
                # Section 4
                elif self.current_section == 4:
                    self.total_tripped_generators += 1
                # Section 5
                elif self.current_section == 5:
                    current_replication = int(token.value)
                    if current_replication > self.total_replications:
                        self.total_replications = current_replication
            else:
                print("Nothing")

        self.store_and_reset(current_replication)

        print(self.total_replications)
        #print(len(self.line_map))
        self.write_to_csv(filename)

    def matrix_process(self,token_collection):
        numlines = len(self.line_map)
        current_replication = 0
        self.build_matrix(self.total_replications, numlines)
        HARDEN_STATUS = 1
        INITIAL_TRIPPED = 2
        TOTAL_TRIPPED = 4
        HARDEN_AND_INITIAL_TRIPPED = HARDEN_STATUS + INITIAL_TRIPPED
        HARDENED_AND_TOTAL_TRIPPED = HARDEN_STATUS + TOTAL_TRIPPED
        TOTAL_AND_INITIAL_TRIPPED = TOTAL_TRIPPED + INITIAL_TRIPPED
        TOTAL_INITIAL_AND_HARDENED = HARDEN_STATUS + INITIAL_TRIPPED + TOTAL_TRIPPED


        print("PRIOR")
        for token in token_collection:
            #print("HERE")
            #print(token)
            if token.type == 'HEADER':
                self.isReinforced = -1
                self.current_section = self.get_section(token)
                if self.current_section == 5:
                    self.reset()
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
                        print("CURRENT REP: " + str(current_replication))
                        print("CURRENT LINE: " + str(self.current_line_number))
                        print("CURRENT VAL: " + str(self.encoded_matrix[current_replication-1][self.current_line_number-1]))
                        print("=================")
                        # if self.encoded_matrix[current_replication-1][self.current_line_number-1] == HARDEN_STATUS:
                        #     self.update_matrix(self.encoded_matrix, current_replication-1, self.current_line_number-1, HARDEN_AND_INITIAL_TRIPPED)
                        if self.encoded_matrix[current_replication-1][self.current_line_number-1] == 0:
                            self.update_matrix(self.encoded_matrix, current_replication-1, self.current_line_number-1, INITIAL_TRIPPED)
                    elif self.current_token_number == 5:             # Check if Hardened
                        self.isReinforced = int(token.value)
                        if self.isReinforced == 1:                     # If Line Hardened
                            if self.encoded_matrix[current_replication-1][self.current_line_number-1] == INITIAL_TRIPPED:
                                self.update_matrix(self.encoded_matrix, current_replication-1, self.current_line_number-1, HARDEN_AND_INITIAL_TRIPPED)
                            else:
                                self.update_matrix(self.encoded_matrix, current_replication-1, self.current_line_number-1, HARDEN_STATUS)
                        else:
                            self.total_initial_tripped += 1          # Otherwise failed
                        self.current_token_number += 1
                    else:
                        self.current_token_number += 1               # No match
                # Section 2
                elif self.current_section == 2:
                    next_line = int(token.line)
                    if next_line > self.curr_line:
                        self.curr_line = next_line
                        self.total_tripped += 1
                #     # Section 3
                # elif self.current_section == 3:
                #     if self.current_token_number > 4:
                #         self.current_token_number = 1
                #     if self.current_token_number == 4:
                #         self.total_shedding_load += float(token.value)
                #     self.current_token_number += 1
                # # Section 4
                # elif self.current_section == 4:
                #     self.total_tripped_generators += 1
                elif self.current_section == 5:
                    current_replication = int(token.value)

        self.reset()
        self.print_matrix()
        # self.encoded_matrix[1][232] = 7
        # print(self.encoded_matrix[1][232])
        # outfile = "matrix.txt"
        # dmatrix = np.matrix(smatrix)
        # outfile = "matrix.txt"
        # dmatrix.tofile(outfile, sep=" ", format="%s")

    def update_matrix(self, matrix, row, col, value):
        matrix[row][col] = value

    def print_matrix(self):
        # self.encoded_matrix[1][232] = 7
        # print(self.encoded_matrix[1][232])
        outfile = "matrix.txt"
        with open(outfile, "w") as f:
            for line in self.encoded_matrix:
                listline = str(line)
                listline = listline.replace("[", "")
                listline = listline.replace("]", "\n")
                templine = listline.replace(",", " ")
                f.write(templine)
                # f.write("end")

    def build_matrix(self, row, col):
        # Build matrix using list comprehension
        self.encoded_matrix = [[0 for x in range(col)] for y in range(row)]

    def store_and_reset(self, replication):
        plan_key = self.add_plan(self.plans, replication)
        if self.plan_length < len(plan_key):
            self.plan_length = len(plan_key)
        self.add_to_csvrows(self.hardening_plans[plan_key])
        self.reset()

    def reset(self):
        self.total_initial_tripped = 0
        self.total_tripped = 0
        self.total_tripped_generators = 0
        self.current_line_number = 0
        self.current_token_number = 1
        self.total_shedding_load = 0.0
        self.isReinforced = 0
        self.plans = [ ]

    def add_to_csvrows(self, replication_index):
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
        with open(filename, 'w') as csvfile:
            filewriter = csv.writer(csvfile)
            column_header = self.build_csv_header()
            filewriter.writerow(column_header)
            for row in self.csv_rows:
                filewriter.writerow(row)

    def build_csv_header(self):
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
        plan_key = tuple(plan)

        # if(plan_key in self.hardening_plans):
        #     temp = self.hardening_plans[plan_key]
        #     temp = replication
        #     self.hardening_plans[plan_key] = temp
        # else:
        self.hardening_plans[plan_key] = replication
        return plan_key

    def get_section(self, token):
        myregex = '=+'
        line = str(token.value)
        section_text = re.sub(myregex, '', line)
        section_text = section_text.strip()
        #print(section_text)
        return section_table[section_text]
