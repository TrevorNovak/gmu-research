import os
import glob
import re

def sortLines(columnValues, correspondingLines):
    sortedLines = []
    sortedColumns = []
    lineCounter = 0
    for value in columnValues:
        if len(sortedColumns) == 0: # if there are no values that are in the sorted list add the value and it is still sorted
            sortedColumns.append(value)
            sortedLines.append(correspondingLines[lineCounter].replace(",", "\t"))
        else: # values are in the sorted list. finds the correct place and inserts the line
            counter = 0
            addedLine = False
            for val in sortedColumns:
                if int(value) < int(val):
                    sortedLines.insert(counter, correspondingLines[lineCounter].replace(",", "\t"))
                    sortedColumns.insert(counter, value)
                    addedLine = True
                    break
                counter += 1
            if addedLine == False:
                sortedLines.append(correspondingLines[lineCounter].replace(",", "\t"))
                sortedColumns.append(value)
        lineCounter += 1
    return sortedLines

def aggregate():
    """aggregates LoadShedding.csv files into one called SimulateOut.txt """
    log_dir = r'simulation/'
    csvfilename='LoadShedding'
    outfile = 'SimulateOut.txt'
    regex = re.compile(r'\d+')

    temp_logs = [log for log in glob.glob(log_dir+'*.csv') if re.match(log_dir+csvfilename+'[0-9]+.csv', str(log))]

    if temp_logs:
        with open(outfile, "w+") as writer:
            lines = []
            column1Values = []
            for log in temp_logs:
                with open(log, "r+") as reader:
                    numbers = regex.findall(log)
                    line = reader.read()
                    rows = line.split("\n")
                    for row in rows:
                        if (len(row.strip()) == 0):
                            continue
                        lines.append(str(numbers[0]) + "," + str(row))
                        firstComma = row.find(",")
                        firstComma += 1
                        firstColumn = row[0:firstComma - 1]
                        column1Values.append(firstColumn)
            lines = sortLines(column1Values, lines)

            groupsOfLines = {}
            groupsOfColumn2Values = {}

            # groups the rows by first column values and records second column value
            for line in lines:
                rowValues = line.split("\t")
                column1 = rowValues[1]
                column2 = rowValues[2]
                if column1 not in groupsOfLines:
                    groupsOfLines[column1] = []
                    groupsOfColumn2Values[column1] = []
                groupsOfLines[column1].append(line)
                groupsOfColumn2Values[column1].append(column2)

            # gets the keys to the groups dictionary
            groupsOfLinesKeys = groupsOfLines.keys()

            # converts the keys from strings to numbers
            numericKeys = []
            for key in groupsOfLinesKeys:
                numericKeys.append(int(key))

            # sorts the keys numericly
            numericKeys.sort()

            # sorts the rows by the second column and converts the keys back to strings
            groupsOfLinesKeys = numericKeys
            for key in groupsOfLinesKeys:
                sortedGroup = sortLines(groupsOfColumn2Values[str(key)], groupsOfLines[str(key)])
                for line in sortedGroup:
                    writer.write(line + "\n")

        print("Files successfully read to " + outfile)
    else:
        print("No files could be located.")


if __name__ == '__main__':
    aggregate()
