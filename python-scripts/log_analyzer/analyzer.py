import re
from strings import *
from collections import defaultdict

# Find replication header
# Read replication until next header
# Determine column number
# Determine hardening plan
# Add to hashmap if not present, otherwise increment replication index and assign
# Count initial tripped lines
# Count total tripped lines
# Sum load shed
# Convert to csv and write to file
# Repeat

COLUMN_NUMBER = 5 # Optional explicit definition of hardened line column.
class Analyzer:
    output_string = ""

    def __init__(self):
        #self.hardening_plans = defaultdict(tuple)
        self.hardening_plans = { }
        self.rep_headers_list = []
        self.replication_indicies = []

    def parse(self):
        pass

    def parse_header(self, searchstring):
        myregex = r"=+" + searchstring + r"=+"
        pattern = myregex.compile(myregex)
        return pattern

    def parse_section(self):
        pass

    def check_header(self, target, string):
        myregex = r"=+\s*" + target + r"\s*=+"
        pattern = re.compile(myregex)
        return pattern.search(string)

    def search_headers(self, target, string):
        myregex = r"=+" + target + r"=+"
        pattern = re.compile(myregex)
        return pattern.finditer(string)

    def search_columns(self, target, string):
        myregex = r"" +target
        pattern = re.compile(myregex)
        return pattern.search(string)

    def determine_plan(self, replication_text, linecol, hardencol):
        ln = 0
        line_num = -1
        plan = []
        isReinforced = 0
        total_initial_tripped = 0
        cursection = 1

        for line in replication_text.splitlines():
            line_string = str(line)
            #print(str(ln) + ": " + line)
            #ln += 1
            if self.check_header(TRIP_HEADER, line_string) is not None:

            elif self.search_columns(REINFORCE, line_string) is None:
                #self.extract(line_string)
                values = line.split()
                print(values)
                isReinforced = int(values[hardencol])
                print("Harden: " + str(isReinforced))
                line_num = int(values[linecol])
                if isReinforced == 1:
                    plan.append(line_num)
                else:
                    total_initial_tripped += 1
            else:
                #parse header
                pass
        print(plan, total_initial_tripped)
        plan_key = tuple(plan)
        print(plan_key)

        if(plan_key in self.hardening_plans):
            temp = self.hardening_plans[plan_key]
            temp += 1
            self.hardening_plans[plan_key] = temp
        else:
            self.hardening_plans[plan_key] = 1
        print(self.hardening_plans)

        return plan

    def extract(self, line):
        pass

    def extract_line_num(self, column, line):
        pass

    def extract_hardening_value(self, column, line):
        pass

    def build_rep_header_indices(self, log, header):
        headers = self.search_headers(header, log)
        #print(headers)
        #print(log)
        for header in headers:
            self.rep_headers_list.append((header.start(), header.end()))

        print(self.rep_headers_list)
            #print("Start: " + str(header.start()))
            #print("End: " + str(header.end()))
            #print (replication)

    def build_replication_indices(self):
        start = 0
        end = 1
        length = len(self.rep_headers_list)

        for index, header in enumerate(self.rep_headers_list):
            if (index+1 < length):
                rep_index = (header[end]+1, self.rep_headers_list[index+1][start]-1)
                #rep_length = self.rep_headers_list[index+1][start] - header[end]
                self.replication_indicies.append(rep_index)
                print(rep_index)
            #print(header[1])
            #print(self.rep_headers_list[index+1][0])

    def analyze_replication(self, replication_text):
        plans = self.determine_plan(replication_text, 0, COLUMN_NUMBER-1)


    def read_file(self, filename):
        """
        Read line until we encounter LINE_FAILURE header. Load text until we encounter next header.
        """
        with open(filename, "r") as f:
            log = f.read()
            log_string = str(log)
            self.build_rep_header_indices(log_string, REP_HEADER)
            self.build_replication_indices()
            replication_text = log_string[self.replication_indicies[0][0] : self.replication_indicies[0][1]]
            print(replication_text)
            self.analyze_replication(replication_text)
            #replication = log_string[head_start : head_end]
            #print(str(line))
            #print(LINE_FAILURE)



            """
            replication =
            if header != None:
                #print(line)
                #print(header)
                print(header.string + " " + str(header.endpos))
                #print(line[10:87])
            else:
                print("Fail")

            self.determine_plan(4)
            """

    def aggregate(self):
        pass

    def calculate(self):
        pass

    def display(self):
        pass


analyzer = Analyzer()
analyzer.read_file("cas_results_5.txt")
#print(LINE_FAILURE)
