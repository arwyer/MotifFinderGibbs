import json
import sys
import os

#Pseudo:
#copy report into directory
#write motif set to directory
#pass info back

class GenerateReport:
  def __init__(self):
      pass
  
  def GenerateMotifReport(self,htmlDir,motifSet):
      reportPath = '/kb/module/lib/MotifFinderGibbs/Utils/Report/*'
      CopyCommand = 'cp -r ' + reportPath + ' ' + htmlDir
      os.system(CopyCommand)
      jsonFName = htmlDir + '/ReportMotif.json'
      with open(jsonFName,'w') as motifjson:
          json.dump(motifSet,motifjson)
      return
