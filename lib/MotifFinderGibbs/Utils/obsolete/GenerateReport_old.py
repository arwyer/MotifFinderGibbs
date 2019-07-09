import json
import sys
import os

#Pseudo:
#copy report into directory
#write motif set to directory
#pass info back

class GenerateReport:
  def __init__(self, htmlDir, motifSet):
      self.htmlDir = htmlDir
      self.motifSet = motifSet
  
  def GenerateReport(self):
      reportPath = '/kb/module/lib/MotifFinderGibbs/Utils/Report/*'
      CopyCommand = 'cp -r ' + reportPath + ' ' + self.htmlDir
      os.system(CopyCommand)
      jsonFName = self.htmlDir + '/ReportMotif.json'
      with open(jsonFName,'w') as motifjson:
          json.dump(self.motifSet,motifjson)
      return
