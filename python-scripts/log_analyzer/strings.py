"""
Strings used in analyzer
"""

# Header Names
LINE_FAILURE = "Line Failure and Hardening in Storm Path"
TRIPPED_LINES = "Tripped Lines"
LOSS_OF_LOADS = "Loss of Loads"
TRIPPED_GENERATORS = "Tripped Generators"
REPLICATION_COUNT = "Replication count"

# Column Names
REINFORCE = "IsReinforce"

# Define Log Categories
REP_HEADER = LINE_FAILURE
TRIP_HEADER = TRIPPED_LINES
HARDEN_COLUMN = REINFORCE
REP_COUNT = REPLICATION_COUNT

# Section Map
total_sections = 5
section_table = {
      LINE_FAILURE       : 1,
      TRIPPED_LINES      : 2,
      LOSS_OF_LOADS      : 3,
      TRIPPED_GENERATORS : 4,
      REPLICATION_COUNT  : 5,
      }
