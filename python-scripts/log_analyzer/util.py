def read_config_file(conf):
    print(conf)
    values = {}
    try:
        with open(conf, 'r') as config:
            for line in config:
                temp = line.split("=")
                values[temp[0]] = temp[1].replace("\n", "")
    except:
        print("Unable to read config file. Please try again.")
    return values
