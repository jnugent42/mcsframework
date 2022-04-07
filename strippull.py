import re
import os 

filepath = '../log'

file = open('../log', 'r')
lines = file.readlines()

for index, line in enumerate(lines):
    #print("Line {}: {}".format(index, line.strip()))
    if 'png' in line.strip():
        result = re.search('textwidth]{(.*)png}', line.strip())
        if index < 35:
             os.system('cp '+result.group(1)+'png'+' paperplots')
        else: 
             os.system('cp '+result.group(1)+'png'+' paperplots2')
    if 'pdf' in line.strip():
        result = re.search('textwidth]{(.*)pdf}', line.strip())
        if index < 35:
             os.system('cp '+result.group(1)+'pdf'+' paperplots')
        else: 
             os.system('cp '+result.group(1)+'pdf'+' paperplots2')
    print(result.group(1))

    
file.close()

'''
with open(filepath) as fp:
   line = fp.readline()
   while line:
       s = line.readline()
       result = re.search('{(.*)}', s)

       #s = s.substring(s.indexOf("{") + 1);
       #s = s.substring(0, s.indexOf("}"));
       print result
'''



