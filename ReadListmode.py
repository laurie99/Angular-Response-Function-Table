# Generate ARF table from a single listmode file
# Jie (Laurie) Zhang
# 03/16/15

# e.g. python ReadListmode.py input listmode output
import sys
import numpy as np
import struct
import csv

def readFile(input, output):
	if sys.argv[2] == "listmode":
		# SPECT listmode data
		# 10 int16, 1 float, 1 interger*1
		# data = []
		with open(input, "rb") as openfile, open(output, 'wb') as csvfile:
			positions = openfile.read(2*10)
			writer = csv.writer(csvfile)
			writer.writerow(['X0','Y0','Z0','XPHANT','YPHANT','ZPHANT','XCRYSTAL','YCRYSTAL','ZCRYSTAL','Energy in Crystal','Photon Weight','Scatter Order'])

			while positions != "":
				locations = struct.unpack('h'*10, positions)
				energy = struct.unpack('d', openfile.read(8))
				scatter = struct.unpack('b', openfile.read(1))
				# print locations+energy+scatter
				# data.append(locations+energy+scatter)
				writer.writerow(locations+energy+scatter)
				positions = openfile.read(2*10)

	else:
		# SPECT projection data
		# 64*64*64 4 bytes float
		with open(input, "rb") as openfile:
			fileContent = openfile.read()
			projection = struct.unpack('f'*64*64*64, fileContent)
			# projection.append(struct.unpack('f'*64, openfile.read(4*64)))
		print projection, len(projection)

def main():
	photons = readFile(sys.argv[1], sys.argv[3])
		
if __name__ == "__main__":
	main()