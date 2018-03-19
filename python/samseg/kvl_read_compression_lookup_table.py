# function [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( fileName )
from _operator import itemgetter


def kvlReadCompressionLookupTable(fileName):
    # %
    # % function [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( fileName )
    # %
    # % The format is: FreeSurferLabel compressedLabel name R G B A
    # %
    # % Things are sorted according compressedLabels
    # %
    #
    # disp( [ 'Reading contexts of file ' fileName ] )
    print('Reading contexts of file {0}'.format(fileName))
    # FreeSurferLabels = [];
    # compressedLabels = [];
    # names = [];
    # colors = [];
    table = []
    #
    # fid = fopen( fileName );
    with open(fileName) as fid:
        # if ( fid < 0 )
        #   error( [ 'Couldn''t read from file ' fileName ] )
        # end
        # while 1
        for line in fid.readlines():
            #   % Get the next line
            #   textLine = fgetl( fid );
            #   if ~ischar( textLine )
            #     break
            #   end
            #
            #   % Extract the contents
            #   [ result, count, errorMessage ] = sscanf( textLine, '%d %d %s %d %d %d %d' );
            #   if ( ~isempty( errorMessage ) )
            #     error( errorMessage )
            #   end
            #
            #   FreeSurferLabel = result( 1 );
            #   compressedLabel = result( 2 );
            #   name = ( char( result( 3 : end-4 ) ) )';
            #   R = result( end-3 );
            #   G = result( end-2 );
            #   B = result( end-1 );
            #   A = result( end );
            FreeSurferLabel, compressedLabel, name, R, G, B, A = [
                data_type(value) for data_type, value in zip(
                    (int, int, str, int, int, int, int),
                    line.split())]
            #
            #   % Add contents to output matrices
            #   FreeSurferLabels = [ FreeSurferLabels; FreeSurferLabel ];
            #   compressedLabels = [ compressedLabels; compressedLabel ];
            #   names = strvcat( names, name );
            #   colors = [ colors; [ R G B A ] ];
            #
            table.append({
                'FreeSurferLabel': FreeSurferLabel,
                'compressedLabel': compressedLabel,
                'name': name,
                'color': [R, G, B, A],
            })
            # end
        # fclose( fid );
    #
    #
    # % Sort output according to compressedLabel
    # [ dummy, sortedIndices ] = sort( compressedLabels );
    # FreeSurferLabels = FreeSurferLabels( sortedIndices, : );
    # names = names( sortedIndices, : );
    # colors = colors( sortedIndices, : );
    table = sorted(table, key=itemgetter('compressedLabel'))
    FreeSurferLabels, names, colors = [
        [entry[key] for entry in table] for key in ['FreeSurferLabel', 'name', 'color']
    ]
    #
    #
    #
    return FreeSurferLabels, names, colors
