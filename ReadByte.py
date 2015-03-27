# Process binary history file from SimSET byte by byte
# Jie (Laurie) Zhang
# 03/25/15

# e.g. python ReadByte.py input
import sys
import struct
import csv

def readFile(input):
	with open(input, "rb") as openfile:
		header = openfile.read(32768)
		byte = openfile.read(1)
		while byte != '':
			print struct.unpack('b', byte)
			byte = openfile.read(1)

def main():
	readFile(sys.argv[1])
		
if __name__ == "__main__":
	main()