# Read simind listmode file
# Jie (Laurie) Zhang
# 03/31/15

# e.g. python ReadSimind.py input listmode output
import sys
import struct
import csv

def read_file(input, output):
	if sys.argv[2] == "listmode":
		"""SPECT listmode data: 10 int16, 1 float, 1 interger*1"""
		with open(input, "rb") as openfile, open(output, 'w') as csvfile:
			positions = openfile.read(2*10)
			writer = csv.writer(csvfile)
			writer.writerow(['X0','Y0','Z0','XPHANT','YPHANT','ZPHANT','XCRYSTAL','YCRYSTAL','ZCRYSTAL','Energy in Crystal','Photon Weight','Scatter Order'])

			while positions != "":
				locations = struct.unpack('h'*10, positions)
				weight = struct.unpack('d', openfile.read(8))
				scatter = struct.unpack('b', openfile.read(1))
				# data.append(locations+energy+scatter)
				writer.writerow(locations+weight+scatter)
				positions = openfile.read(2*10)

	else:
		"""SPECT projection data: 64*64*64 4 bytes float"""
		with open(input, "rb") as openfile:
			fileContent = openfile.read()
			projection = struct.unpack('f'*64*64*64, fileContent)
			# projection.append(struct.unpack('f'*64, openfile.read(4*64)))
		print(projection, len(projection))
		return projection

def main():
	read_file(sys.argv[1], sys.argv[3])
		
if __name__ == "__main__":
	main()