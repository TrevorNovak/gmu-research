import collections
import re
import csv
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
        self.hardening_plans = { }
        self.csv_rows = []
        self.plans = [ ]

    def process(self, token_collection, filename, flag):
        current_replication = 0

        for token in token_collection:
            if flag:
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
            else:
                print("Nothing")

        self.store_and_reset(current_replication)

        self.write_to_csv(filename)

    def store_and_reset(self, replication):
        plan_key = self.add_plan(self.plans, replication)
        if self.plan_length < len(plan_key):
            self.plan_length = len(plan_key)
        self.add_to_csvrows(self.hardening_plans[plan_key])
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
