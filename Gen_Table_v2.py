# Call ARF_listmode_v2 to generate multiple ARF tables
# Python3
# Jie (Laurie) Zhang
# 04/15/15
# e.g. python Gen_Tables_v2.py [lower_energy upper_energy]/[size1 size2] input1.lmf [input2.lmf ...]
import sys
import re
from subprocess import call
from functools import partial

def main():
  """
  This command-line parsing code is provided.
  Make a list of command line arguments, omitting the [0] element
  which is the script itself.
  """
  if len(sys.argv) < 2:
    print('usage: lower_energy upper_energy input1.lmf [input2.lmf ...]')
    sys.exit(1)

  input_file = sys.argv[3]
  if re.search('lmf', input_file):
    lower_energy = sys.argv[1]
    upper_energy = sys.argv[2]
    input_files = sys.argv[3:]

    for input_file in input_files:
      output_file = 'ARF'+re.search(r'\d+',input_file).group()+'_lmf_v2'
      print('Output is: ',output_file)
      call(['python','../ARF_listmode_v2.py',input_file,output_file,lower_energy,upper_energy])

  elif re.search('bim', input_file):
    size1 = sys.argv[1]
    size2 = sys.argv[2]
    input_files = sys.argv[3:]

    for input_file in input_files:
      output_file = 'ARF'+re.search(r'\d+',input_file).group()+'_bin_v2'
      print('Output is: ',output_file)
      call(['python','../ARF_listmode_v2.py',input_file,output_file,size1,size2])

if __name__ == '__main__':
  main()