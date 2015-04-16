# Call ARF_listmode_v2 to generate multiple ARF tables
# Python3
# Jie (Laurie) Zhang
# 04/09/15
# e.g. python Gen_Tables.py lower_energy upper_energy input1.lmf [input2.lmf ...]
import sys
import re
from subprocess import call
from functools import partial

def main():
  # This command-line parsing code is provided.
  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  if len(sys.argv) < 4:
    print('usage: lower_energy upper_energy input1.lmf [input2.lmf ...]')
    sys.exit(1)

  lower_energy = sys.argv[1]
  upper_energy = sys.argv[2]
  input_files = sys.argv[3:]

  for input_file in input_files:
    output_file = 'ARF'+re.search(r'\d+',input_file).group()+'_win1_v1'
    print('Output is: ',output_file)
    call(['python','../ARF_listmode_v1.py',input_file,output_file,lower_energy,upper_energy])

if __name__ == '__main__':
  main()