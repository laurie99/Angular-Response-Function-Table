# Process simind listmode file into ARF table
# Python3
# Jie (Laurie) Zhang
# 04/01/15
# e.g. python ARF_listmode_v1.py *.lmf output logfile_name
import sys
import cProfile
import csv
import struct
from math import sqrt, atan, degrees, pi, floor
import numpy as np

def read_file(filename):
	"""Read listmode data: 10 int16, 1 float, 1 interger*1"""
	data = []
	with open(filename, "rb") as openfile:
		positions = openfile.read(2*10)
		while positions != "":
			try:
				locations = struct.unpack('h'*10, positions)
				weight = struct.unpack('d', openfile.read(8))
				scatter = struct.unpack('b', openfile.read(1))
				data.append(locations+weight+scatter)
				positions = openfile.read(2*10)
			except:
				break
	return data

def angles(data):
	"""Calculate the cosine of the azimuthal angle, and tan/cot of polar angle"""
	photon_angles = []
	for photon in data:
		energy = float(photon[9])/10
		if (126. <= energy <= 154.):
			x0 = float(photon[0])/100
			y0 = float(photon[1])/100
			z0 = float(photon[2])/100
			xp = float(photon[3])/100
			yp = float(photon[4])/100
			zp = float(photon[5])/100
			xc = float(photon[6])/100
			yc = float(photon[7])/100
			zc = float(photon[8])/100
		
			weight = float(photon[10])
			scatter_order = int(photon[11])

			# print(x0,y0,z0,xp,yp,zp,xc,yc,zc,weight,scatter_order)

			x_vec = xc-x0
			y_vec = yc-y0
			z_vec = zc-z0
			quadrant = 0
			tan_phi = cot_phi = None
			travel_dist = sqrt(pow(x_vec,2) + pow(y_vec,2) + pow(z_vec,2))
			cos_theta = z_vec/travel_dist

			if (x_vec >= 0 and y_vec > 0):
				quadrant = 1
			elif (x_vec < 0 and y_vec >= 0):
				quadrant = 2
			elif (x_vec <= 0 and y_vec < 0):
				quadrant = 3
			elif (x_vec > 0 and y_vec <= 0):
				quadrant = 4

			try:
				tan_phi = y_vec/x_vec
			except ZeroDivisionError:
				tan_phi = float("inf")
			try:
				cot_phi = x_vec/y_vec
			except ZeroDivisionError:
				cot_phi = float("inf")
	
			photon_angles.append([quadrant, cos_theta, tan_phi, cot_phi, energy, weight])
	return photon_angles

def ARF_table(photon_angles):
	"""Bin the photons into a 2048*2048 matrix according to cos_theta and tan_phi/cot_phi. Then normalize the table"""
	table = np.zeros((2048, 512*4))
	cos_list = np.concatenate([np.linspace(1., 0.99, 1025), np.linspace(0.99, 0.95, 1535-1024+2)[1:], np.linspace(0.95, 0.75, 1791-1536+2)[1:], np.linspace(0.75, 0., 2047-1792+2)[1:]]) 
	tan_list13 = np.linspace(0., 1., 257)
	cot_list13 = np.linspace(1., 0., 257)
	tan_list24 = np.linspace(-1., 0., 257)
	cot_list24 = np.linspace(0., -1., 257)

	for photon in photon_angles:
		quadrant = photon[0]
		cos_theta = photon[1]
		tan_phi = photon[2]
		cot_phi = photon[3]
		weight = photon[5]
		# print(quadrant,cos_theta,tan_phi,cot_phi,weight)
		
		theta_ind = phi_ind = None

		if quadrant == 0:
			theta_ind = phi_ind = 0
		else:
			for ii in cos_list:
				if ii <= cos_theta:
					theta_ind = max(int(np.argwhere(cos_list == ii)-1), 0)
					break

			if abs(tan_phi) <= 1:
				if quadrant == 1:
					for ii in tan_list13:
						if ii >= tan_phi:
							phi_ind = max(int(np.argwhere(tan_list13 == ii)-1), 0)
							break
				elif quadrant == 2:
					for ii in tan_list24:
						if ii >= tan_phi:
							phi_ind = int(1*512 + 255 + np.argwhere(tan_list24 == ii))
							break
				elif quadrant == 3:
					for ii in tan_list13:
						if ii >= tan_phi:
							phi_ind = int(2*512 -1 + np.argwhere(tan_list13 == ii))
							break
				elif quadrant == 4:
					for ii in tan_list24:
						if ii >= tan_phi:
							phi_ind = int(3*512 + 255 + np.argwhere(tan_list24 == ii))
							break
			elif abs(cot_phi) <= 1:
				if quadrant == 1:
					for ii in cot_list13:
						if ii <= cot_phi:
							phi_ind = int(255 + np.argwhere(cot_list13 == ii))
							break
				elif quadrant == 2:
					for ii in cot_list24:
						if ii <= cot_phi:
							phi_ind = int(1*512 -1 + np.argwhere(cot_list24 == ii))
							break
				elif quadrant == 3:
					for ii in cot_list13:
						if ii <= cot_phi:
							phi_ind = int(2*512 + 255 + np.argwhere(cot_list13 == ii))
							break
				elif quadrant == 4:
					for ii in cot_list24:
						if ii <= cot_phi:
							phi_ind = int(3*512 -1 + np.argwhere(cot_list24 == ii))
							break

		# print(quadrant, int(theta_ind), int(phi_ind), tan_phi, cot_phi, weight)
		table[theta_ind, phi_ind] += weight

	# Normalize the sum of photon weights
	solid_angles = np.zeros((2048, 512*4))
	delta_phi = np.zeros(2048)

	for ii in range(len(delta_phi)):
		if ii <= 255:
			delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(tan_list13[ii+1])-atan(tan_list13[ii])))
		elif 255 < ii <= 511:
			if cot_list13[ii%256+1] != 0:
				delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(1/cot_list13[ii%256+1])-atan(1/cot_list13[ii%256])))
			else:
				delta_phi[ii] = delta_phi[ii+1024] = 90 - degrees(atan(1/cot_list13[ii%256]))
		elif 511 < ii <= 767:
			delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(tan_list24[ii%256+1])-atan(tan_list24[ii%256])))
		elif 767 < ii <= 1023:
			if cot_list24[ii%256] != 0:
				delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(1/cot_list24[ii%256+1])-atan(1/cot_list24[ii%256])))
			else:
				delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(1/cot_list24[ii%256+1]))) - 90
		elif ii >= 1024:
			break

	for ii in range(2048):
		solid_angles[ii] = abs(cos_list[ii+1]-cos_list[ii]) * delta_phi

	# table = 4*pi*table/(0.8910*1e6*solid_angles)
	table = 4*pi*table/(1e6*solid_angles)
	return abs(table)  # eliminate -0.0 entries

def main():
	# check if there are negative entries
	data = read_file(sys.argv[1])
	photon_angles = angles(data)
	np.savetxt(sys.argv[2],ARF_table(photon_angles),fmt='%.5f')
		
if __name__ == "__main__":
	cProfile.run('main()',sys.argv[3])
	# main()