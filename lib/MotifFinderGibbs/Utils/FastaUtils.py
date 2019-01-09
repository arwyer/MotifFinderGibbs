import os

def RemoveRepeats(path,newpath):
    command = '/kb/deployment/bin/meme/bin/dust ' + path + ' > ' + newpath
    os.system(command)
