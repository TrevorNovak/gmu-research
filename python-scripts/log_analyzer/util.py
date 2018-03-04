def read_config_file(conf):
    values = {}
    with open(conf, 'r') as config:
        for line in config:
            temp = line.split("=")
            values[temp[0]] = temp[1].replace("\n", "")

    return values
