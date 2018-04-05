# Pyni configuration
import configparser, os
import multiprocessing as mp

configFileDirOptions = [
    '.',
    os.path.expanduser('~'),
    '/usr/local/etc'
]

# Default configuration
config = configparser.ConfigParser()
config['pyni'] = {
    'datadir': '{}/data'.format(os.path.dirname(__file__)), #TODO change to subdir pyni when making public
    'projectdir': '/tmp',
    'threads': mp.cpu_count()
}

# Read configuration file
for configFileDir in configFileDirOptions:
    configFile = os.path.join(configFileDir,'pyni.cfg')
    if os.path.exists(configFile):
        config.read(configFile)
        break #only reads the first config file found
