# Process PHG history file from SimSET-PHG module
# Jie (Laurie) Zhang
# 03/31/15

# e.g. python ReadPHGhistory.py input output
import sys
import struct
import csv

def readFile(input, output):
	with open(input, "rb") as openfile, open(output, 'w') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(['photon_number','x_position','y_position','z_position','x_cosine','y_cosine','z_cosine','scatters_in_object','weight','energy'])

		header = openfile.read(32768)
		photon_number = openfile.read(1)

		while photon_number != "":
			try:
				photon_number = struct.unpack('b', photon_number)
			except Exception as e:
				print(str(e), photon_number)
				break

			positions = struct.unpack('d'*6, openfile.read(6*8))
			scatter = struct.unpack('i', openfile.read(4))
			weight = struct.unpack('d', openfile.read(8))
			energy = struct.unpack('d', openfile.read(8))
			# detector = struct.unpack('d'*4, openfile.read(4*8))
			writer.writerow(photon_number+positions+scatter+weight+energy)
			# print photon_number+positions+scatter+energy+detector
			photon_number = openfile.read(1)

def main():
	readFile(sys.argv[1], sys.argv[2])
		
if __name__ == "__main__":
	main()