import sys
import csv

# total_energy = 0
total_weights = 0
with open(sys.argv[1]) as f:
	f_csv = csv.reader(f)
	headers = next(f_csv)
	nrow = 0
	for row in f_csv:
		nrow += 1
		# total_energy += int(row[9])
		total_weights += float(row[10])
		# print total_weights, row[10]
		if float(row[10]) < 0:
			print "WARNING!!!",nrow
			print row
		# if total_weights < 0:
		# 	print nrow
		# 	print row
		# 	break
# print "Total Energy Deposited:", total_energy/10
print "Total Photon Weights:", total_weights