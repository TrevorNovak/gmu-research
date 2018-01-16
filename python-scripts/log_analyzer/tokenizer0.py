import collections
import re
import csv
from strings import *

HARDEN_COL = 5
LINENUM_COL = 1
TOTAL_COL = 6
Token = collections.namedtuple('Token', ['type', 'value', 'line', 'start', 'end'])

def tokenize(text):
    #keywords = {'IF', 'THEN', 'ENDIF', 'FOR', 'NEXT', 'GOSUB', 'RETURN'}
    token_specification = [
        ('HEADER',    r'=+\s*[a-zA-Z ]+\s*=+'),     # Header: start of a section.
        ('COLUMN',    r'[a-zA-Z ]+'),               # Column: data columns.
        ('DATA',      r'[0-9.]+'),                  # Data: data to process.
        ('NEWLINE',   r'\n'),                       # Newline character.
        ('MISC',      r'\w*'),                      # Blank lines and things.
        ('UNDEFINED', r'.'),                        # No token match, skip.
    ]
    tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
    line_num = 1
    line_start = 0
    for mo in re.finditer(tok_regex, text):
        kind = mo.lastgroup
        value = mo.group(kind)
        start = mo.start()
        end = mo.end()

        if kind == 'NEWLINE':
            line_start = mo.end()
            line_num += 1
        elif kind == 'MISC':
            pass
        elif kind == 'UNDEFINED':
            raise RuntimeError(f'{value!r} unexpected on line {line_num}')
        else:
            value = value.strip()
            if value == '' or value == None:
                pass
            else:
                column = mo.start() - line_start
                yield Token(kind, value, line_num, start, end)

def process(token_collection):
    hardening_plans = { }
    csv_rows = []
    current_section = 0
    current_line_number = -1
    current_token_number = 1
    total_initial_tripped = 0
    curr_line = 0
    total_tripped = 0
    isReinforced = -1
    plans = [ ]
    total_shedding_load = 0.0
    total_tripped_generators = 0


    for token in token_collection:
        if token.type == 'HEADER':
            print("Header: ")
            #current_token_number = 1
            #current_line_number = -1
            # total_initial_tripped = 0
            # total_tripped = 0
            isReinforced = -1
            current_section = get_section(token)
            if current_section == 1:
                print("PLANS VAL: " + str(plans))
                if len(plans) > 0:
                    store_and_reset()
        elif token.type == 'COLUMN':
            pass
        elif token.type == 'DATA':

            # # # # # # # # # # # # #
            # Process Each Section  #
            # # # # # # # # # # # # #

            # Section 1
            if current_section == 1:                        # Initial Tripped
                if current_token_number > TOTAL_COL:
                    current_token_number = 1
                if current_token_number == 1:               # LineNum token
                    print("FIRST TOKEN")
                    current_line_number = int(token.value)  # Save LineNum
                    current_token_number += 1
                elif current_token_number == 5:             # Check if Hardened
                    print("FIFTH TOKEN")
                    isReinforced = int(token.value)
                    print("REINFOCED: " + str(isReinforced))
                    print("CURR LINE: " + str(current_line_number))
                    if isReinforced == 1:                   # If Line Hardened
                        plans.append(current_line_number)   # Add to plan
                        print(plans)
                    else:
                        total_initial_tripped += 1          # Otherwise failed
                    current_token_number += 1
                else:
                    current_token_number += 1               # No match
            # Section 2
            elif current_section == 2:
                next_line = int(token.line)
                if next_line > curr_line:
                    curr_line = next_line
                    total_tripped += 1
            # Section 3
            elif current_section == 3:
                print(token)
                if current_token_number > 4:
                    current_token_number = 1
                if current_token_number == 4:
                    total_shedding_load += float(token.value)
                    #print(total_shedding_load)
                current_token_number += 1
            # Section 4
            elif current_section == 4:
                total_tripped_generators += 1
                #print(total_tripped_generators)

        else:
            print("Nothing")

    store_and_reset()

    write_to_csv("results.csv")
    print(csv_rows)
    # print(plans)
    # print(total_initial_tripped)
    # print(total_tripped)
    # print(section_table)

    def store_and_reset():
        plan_key = add_plan(plans)
        add_to_csvrows(plans, total_initial_tripped, total_tripped, total_shedding_load, total_tripped_generators, hardening_plans[plan_key])
        total_initial_tripped = 0
        total_tripped = 0
        total_tripped_generators = 0
        current_line_number = 0
        current_token_number = 1
        total_shedding_load = 0.0
        isReinforced = 0
        plans = [ ]

def add_to_csvrows(plans, total_initial_tripped, total_tripped, total_shedding_load, total_tripped_generators, replication_index):
    csv_rows.append([str(plans[0]),str(total_initial_tripped),str(total_tripped),str(total_shedding_load),str(total_tripped_generators),str(replication_index)])
    #line = str(plans[0])+","+str(total_initial_tripped)+","+str(total_tripped)+","+str(total_shedding_load)+","+str(total_tripped_generators)+","+str(replication_index)
    #print(line)
    # with open('results.csv', 'rw') as csvfile:
    #     filewriter = csv.writer(csvfile)
    #     plans_string = str(plans[0])
    #     init_tripped_string = str(total_initial_tripped)
    #     #filewriter.writerow(['this', 'that'])
    #     filewriter.writerow([str(plans[0]),str(total_initial_tripped),str(total_tripped),str(total_shedding_load),str(total_tripped_generators),str(replication_index)])

def write_to_csv(filename):
    with open(filename, 'w') as csvfile:
        filewriter = csv.writer(csvfile)
        for row in csv_rows:
        #plans_string = str(plans[0])
        #init_tripped_string = str(total_initial_tripped)
            filewriter.writerow(row)
        #filewriter.writerow([str(plans[0]),str(total_initial_tripped),str(total_tripped),str(total_shedding_load),str(total_tripped_generators),str(replication_index)])


def add_plan(plan):
    plan_key = tuple(plan)
    print(plan_key)

    if(plan_key in hardening_plans):
        temp = hardening_plans[plan_key]
        temp += 1
        hardening_plans[plan_key] = temp
    else:
        hardening_plans[plan_key] = 1
    print(hardening_plans)
    return plan_key


def get_section(token):
    #myregex = r'[a-zA-Z ]'
    myregex = '=+'
    #pattern = re.compile(myregex)
    #mo = pattern.search(token.value)
    #print(mo.string)
    line = str(token.value)
    section_text = re.sub(myregex, '', line)
    section_text = section_text.strip()
    print(section_text)
    return section_table[section_text]

def read(filename):
    with open(filename, "r") as f:
        log = f.read()
        return str(log)


statements = read("cas_results_5.txt")
process(tokenize(statements))
# for token in tokenize(statements):
#     print(token)
