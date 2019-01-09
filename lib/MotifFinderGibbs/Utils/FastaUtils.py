import os

def RemoveRepeats(path):
    command = '/kb/deployment/bin/meme/bin/dust ' + path + ' > ' + path
    os.system(command)
