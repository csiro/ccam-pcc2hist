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
func_assign = re.compile('ncf90_\w+\s*=')
func_end = re.compile('\s*end\s+function\s+(\w+)')

f = open(sys.argv[1], "r")
lines = f.readlines()
f.close()

state = "SEARCH_FUNC_START"

for line in lines:
    line = line.rstrip()

    if state == "SEARCH_FUNC_START":
        start_match = func_start.search(line)
        if start_match is not None:
            # Skip until past declarations
            state = "SKIP_DECL"
        print line
    elif state == "SKIP_DECL":
        if func_assign.search(line) is not None or \
                line.find('present') != -1 or \
                line.find('=>') != -1: 
            print '        print *, ' \
                '">> Entering {0}"\n'.format(start_match.group(1))
            state = "SEARCH_FUNC_END"
        print line
    elif state == "SEARCH_FUNC_END":
        end_match = func_end.search(line)
        if end_match is not None:
            print '        print *, ' \
                '">> Exiting {0}"\n'.format(end_match.group(1))
            state = "SEARCH_FUNC_START"
        print line


    
