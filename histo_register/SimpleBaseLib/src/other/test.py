import sbl

# display text
print( "beginning python test" )

# create a command router
cr = sbl.CommandRouter()

# run a command with args
cr.cmdtest( 1, "hello", 3.14, True )

# create a config object
conf1 = sbl.Config()
conf1.stringParam = "foo"
conf1.doubleParam = 7.77
conf1.boolParam = True
conf1.intParam = 137

# save and load
conf1.save( "test.conf" );
conf2 = sbl.Config()
conf2.load( "test.conf" );

# read config values
cr.disp( 0, "config:" )
cr.disp( 1, "intParam: %s" % conf2.intParam )
cr.disp( 1, "stringParam: %s" % conf2.stringParam )
cr.disp( 1, "doubleParam: %s" % conf2.doubleParam )
cr.disp( 1, "boolParam: %s" % conf2.boolParam )

# run a command using a config
cr.cmdtest( conf2 )

# test command cancelling
cr.statustest()
if cr.checkCommandCancel():
    cr.disp( 0, "cancelled" );

# test warning
cr.warning( "this is a warning" );

# done with test
cr.disp( 0, "done with python test" );
