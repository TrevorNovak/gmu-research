def read_config_file(conf):
    values = {}
    with open(conf, 'r') as config:
        for line in config:
            temp = line.split("=")
            print(temp[0])
            values[temp[0]] = temp[1]

    print(values)
    return values

read_config_file("input_files/config.txt")
