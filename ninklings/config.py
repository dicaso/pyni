# Ninklings configuration
import configparser, os
configFileDirOptions = [
    '.',
    os.path.expanduser('~'),
    '/etc'
]

# Default configuration
config = configparser.ConfigParser()
config['ninklings'] = {
    'datadir': '{}/../data'.format(os.path.dirname(__file__)), #TODO change to subdir ninklings when making public
    'projectdir': '/tmp'
}

# Read configuration file
for configFileDir in configFileDirOptions:
    configFile = os.path.join(configFileDir,'ninklings.cfg')
    if os.path.exists(configFile):
        config.read(configFile)
        break #only reads the first config file found
