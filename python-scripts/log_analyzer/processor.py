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
        self.hardening_plans = { }
        self.csv_rows = []
        self.plans = [ ]

    def process(self, token_collection):

        for token in token_collection:
            if token.type == 'HEADER':
                #print("Header: ")
                #self.current_token_number = 1
                #self.current_line_number = -1
                # self.total_initial_tripped = 0
                # self.total_tripped = 0
                self.isReinforced = -1
                self.current_section = self.get_section(token)
                #print(self.current_section)
                if self.current_section == 1:
                    #print("PLANS VAL: " + str(self.plans))
                    if len(self.plans) > 0:
                        #print("PLANS VAL AFTER: " + str(self.plans))
                        self.store_and_reset()
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
                        #print("FIRST TOKEN")
                        self.current_line_number = int(token.value)  # Save LineNum
                        self.current_token_number += 1
                    elif self.current_token_number == 5:             # Check if Hardened
                        #print("FIFTH TOKEN")
                        self.isReinforced = int(token.value)
                        #print("REINFOCED: " + str(self.isReinforced))
                        #print("CURR LINE: " + str(self.current_line_number))
                        if self.isReinforced == 1:                   # If Line Hardened
                            self.plans.append(self.current_line_number)   # Add to plan
                            #print(self.plans)
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
                    #print(token)
                    if self.current_token_number > 4:
                        self.current_token_number = 1
                    if self.current_token_number == 4:
                        self.total_shedding_load += float(token.value)
                        #print(total_shedding_load)
                    self.current_token_number += 1
                # Section 4
                elif self.current_section == 4:
                    self.total_tripped_generators += 1
                    #print(self.total_tripped_generators)

            else:
                print("Nothing")

        self.store_and_reset()

        self.write_to_csv("results.csv")
        #print(self.csv_rows)

    def store_and_reset(self):
        plan_key = self.add_plan(self.plans)
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
        #print("PLANS: " + str(self.plans))
        self.csv_rows.append([str(self.plans[0]),str(self.total_initial_tripped),str(self.total_tripped),str(self.total_shedding_load),str(self.total_tripped_generators),str(replication_index)])

    def write_to_csv(self, filename):
        with open(filename, 'w') as csvfile:
            filewriter = csv.writer(csvfile)
            column_header = self.build_csv_header()
            filewriter.writerow(column_header)
            for row in self.csv_rows:
                filewriter.writerow(row)

    def build_csv_header(self):
        return ['LineNum','Total Initial Tripped Lines','Total Tripped Lines','Total Shedding Load Amount','Total Tripped Generators','Replication Index']

    def add_plan(self, plan):
        plan_key = tuple(plan)
        #print(plan_key)

        if(plan_key in self.hardening_plans):
            temp = self.hardening_plans[plan_key]
            temp += 1
            self.hardening_plans[plan_key] = temp
        else:
            self.hardening_plans[plan_key] = 1
        #print(self.hardening_plans)
        return plan_key

    def get_section(self, token):
        myregex = '=+'
        line = str(token.value)
        section_text = re.sub(myregex, '', line)
        section_text = section_text.strip()
        #print(section_text)
        return section_table[section_text]
