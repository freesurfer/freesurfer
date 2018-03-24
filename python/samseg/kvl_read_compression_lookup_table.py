from _operator import itemgetter


def kvlReadCompressionLookupTable(fileName):
    #
    # The format is: FreeSurferLabel compressedLabel name R G B A
    #
    # Things are sorted according compressedLabels
    #
    print('Reading contexts of file {0}'.format(fileName))
    table = []
    with open(fileName) as fid:
        for line in fid.readlines():
            FreeSurferLabel, compressedLabel, name, R, G, B, A = [
                data_type(value) for data_type, value in zip(
                    (int, int, str, int, int, int, int),
                    line.split())]
            #   % Add contents to output matrices
            table.append({
                'FreeSurferLabel': FreeSurferLabel,
                'compressedLabel': compressedLabel,
                'name': name,
                'color': [R, G, B, A],
            })
    # Sort output according to compressedLabel
    table = sorted(table, key=itemgetter('compressedLabel'))
    FreeSurferLabels, names, colors = [
        [entry[key] for entry in table] for key in ['FreeSurferLabel', 'name', 'color']
    ]
    return FreeSurferLabels, names, colors
