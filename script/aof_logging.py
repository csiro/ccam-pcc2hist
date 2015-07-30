# Annotate functions with logging calls AOP-style given a .f90 file,
# sending result to stdout.
#
# David.Benn@csiro.au, July 2015

import re
import sys

if len(sys.argv) < 2:
    print "No F90 file specified."
    sys.exit(1)

func_start = re.compile('\s*function\s+(\w+)\s*\(')
func_end = re.compile('\s*end\s+function\s+(\w+)')

f = open(sys.argv[1], "r")
lines = f.readlines()
f.close()

for line in lines:
    start_match = func_start.search(line)
    if start_match is not None:
        print line.rstrip()
        print '        print *, ">> Entering {0}"'.format(start_match.group(1))
    else:
        end_match = func_end.search(line)
        if end_match is not None:
            print '        print *, ">> Exiting {0}"'.format(end_match.group(1))
            
        print line.rstrip()


    
