#!/usr/bin/python

class Employee:
   'Common base class for all employees'
   empCount = 0

   def __init__(self, htmlDir, motifSet):
      self.htmlDir = htmlDir
      self.motifSet = motifSet
   
   def displayCount(self):
     print ("Total Employee %d" % Employee.htmlDir)

   def displayEmployee(self):
      print ("htmlDir : ", self.htmlDir,  ", motifSet: ", self.motifSet)


