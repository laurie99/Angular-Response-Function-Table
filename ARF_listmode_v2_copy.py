# Process simind listmode or projection file into ARF table in an OOP style
# Faster way to compute the index
# Python3
# Jie (Laurie) Zhang
# 04/06/15
# e.g. python ARF_listmode_v2.py input_name output_name [lower_energy upper_energy]/[size1 size2]
import sys
import cProfile
import csv
import struct
import re
from math import sqrt, atan, degrees, pi, floor
import numpy as np

class PhotonListMode(object):
    """Listmode photon: 'X0','Y0','Z0','XPHANT','YPHANT','ZPHANT','XCRYSTAL','YCRYSTAL','ZCRYSTAL','Energy in Crystal','Photon Weight','Scatter Order'"""

    def __init__(self, locations, energy, weight, scatter):
        # self.X0 = locations[0]
        # self.Y0 = locations[1]
        # self.Z0 = locations[2]

        # self.Xp = locations[3]
        # self.Yp = locations[4]
        # self.Zp = locations[5]

        # self.Xc = locations[6]
        # self.Yc = locations[7]
        # self.Zc = locations[8]

        self.energy = energy
        self.weight = weight
        self.scatter = scatter

        self.X_vec = locations[6]-locations[0]
        self.Y_vec = locations[7]-locations[1]
        self.Z_vec = locations[8]-locations[2]

    def mod(self):
        travel_dist = sqrt(self.X_vec**2 + self.Y_vec**2 + self.Z_vec**2)
        return travel_dist

    def quadrant(self):
        quadrant = 0
        if (self.X_vec > 0 and self.Y_vec >= 0):
            quadrant = 1
        elif (self.X_vec <= 0 and self.Y_vec > 0):
            quadrant = 2
        elif (self.X_vec < 0 and self.Y_vec <= 0):
            quadrant = 3
        elif (self.X_vec >= 0 and self.Y_vec < 0):
            quadrant = 4
        return quadrant
        
    def cos_theta(self):
        cos_theta = self.Z_vec/self.mod()
        return cos_theta

    def tan_phi(self):
        tan_phi = None
        try:
            tan_phi = self.Y_vec/self.X_vec
        except ZeroDivisionError:
            tan_phi = float("inf")
        return tan_phi

    def cot_phi(self):
        cot_phi = None
        try:
            cot_phi = self.X_vec/self.Y_vec
        except ZeroDivisionError:
            cot_phi = float("inf")
        return cot_phi

class PhotonBinMode(object):
    """Binned photon: 'Photon Weight', (XCRYSTAL_IND, YCRYSTAL_IND)
    """
    def __init__(self, weight, multi_index, size1, size2):
        self.weight = weight

        self.X_vec = (multi_index[0] - float(size1/2)-0.5) * 0.540
        self.Y_vec = - (multi_index[1] - float(size2/2)-0.5) * 0.540
        self.Z_vec = 200    # cm

    def mod(self):
        travel_dist = sqrt(self.X_vec**2 + self.Y_vec**2 + self.Z_vec**2)
        return travel_dist

    def quadrant(self):
        quadrant = 0
        if (self.X_vec > 0 and self.Y_vec >= 0):
            quadrant = 1
        elif (self.X_vec <= 0 and self.Y_vec > 0):
            quadrant = 2
        elif (self.X_vec < 0 and self.Y_vec <= 0):
            quadrant = 3
        elif (self.X_vec >= 0 and self.Y_vec < 0):
            quadrant = 4
        return quadrant
        
    def cos_theta(self):
        cos_theta = self.Z_vec/self.mod()
        return cos_theta

    def tan_phi(self):
        tan_phi = None
        try:
            tan_phi = self.Y_vec/self.X_vec
        except ZeroDivisionError:
            tan_phi = float("inf")*np.sign(self.Y_vec)
        return tan_phi

    def cot_phi(self):
        cot_phi = None
        try:
            cot_phi = self.X_vec/self.Y_vec
        except ZeroDivisionError:
            cot_phi = float("inf")*np.sign(self.X_vec)
        return cot_phi

def read_listmode(filename, lower_energy, upper_energy):
    """Read listmode data: 10 int16, 1 float, 1 interger*1"""
    data = []
    with open(filename, "rb") as openfile:
        positions = openfile.read(2*9)
        while positions != "":
            try:
                locations = struct.unpack('h'*9, positions)
                energy = float(struct.unpack('h', openfile.read(2))[0])/10
                weight = float(struct.unpack('d', openfile.read(8))[0])
                scatter = int(struct.unpack('b', openfile.read(1))[0])

                # 20% energy window
                if float(lower_energy) <= energy <= float(upper_energy):
                    photon = PhotonListMode(locations, energy, weight, scatter)
                    data.append(photon)

                positions = openfile.read(2*9)
            except:
                break
    return data

def read_projection(filename, size1, size2):
    """Read projection data: 4 byte floats"""
    data = []
    matrix = np.fromfile(filename, dtype = 'f4').reshape(size1, size2)
    it = np.nditer(matrix, flags=['multi_index'])
    while not it.finished:
        # print("%f <%s>" % (it[0], it.multi_index),)
        photon = PhotonBinMode(it[0], it.multi_index, size1, size2)
        data.append(photon)
        it.iternext()
    return data

def ARF_table(photon):
    """Find the correct position in the table"""
    theta_ind = phi_ind = None
    quadrant = photon.quadrant()
    cos_theta = photon.cos_theta()
    tan_phi = photon.tan_phi()
    cot_phi = photon.cot_phi()

    if quadrant == 0:
        theta_ind = phi_ind = 0
    else:
        if 1 >= cos_theta >= 0.99:
            theta_ind = floor(1023*(cos_theta-1)/(0.99-1))
        elif 0.99 >= cos_theta >= 0.95:
            theta_ind = floor(512*(cos_theta-0.99)/(0.95-0.99)) + 1023
        elif 0.95 >= cos_theta >= 0.75:
            theta_ind = floor(256*(cos_theta-0.95)/(0.75-0.95)) + 512 + 1023
        elif 0.75 >= cos_theta >= 0:
            theta_ind = floor(256*(cos_theta-0.75)/(-0.75)) + 256 + 512 + 1023

        if abs(tan_phi) <= 1:
            if quadrant == 1:
                phi_ind = floor(255*tan_phi)
            elif quadrant == 2:
                phi_ind = 768 + floor(255*(tan_phi+1))
            elif quadrant == 3:
                phi_ind = 1024 + floor(255*tan_phi)
            elif quadrant == 4:
                phi_ind = 1792 + floor(255*(tan_phi+1))
        elif abs(cot_phi) <= 1:
            if quadrant == 1:
                phi_ind = 256 + floor(-255*(cot_phi-1))
            elif quadrant == 2:
                phi_ind = 512 + floor(-255*cot_phi)
            elif quadrant == 3:
                phi_ind = 1280 + floor(-255*(cot_phi-1))
            elif quadrant == 4:
                phi_ind = 1536 + floor(-255*cot_phi)

    if theta_ind == None or phi_ind == None:
        print(quadrant, cos_theta, tan_phi, cot_phi, theta_ind, phi_ind)

    return (theta_ind, phi_ind)

def normalize_table(table):
    solid_angles = np.zeros((2048, 512*4))
    delta_phi = np.zeros(2048)

    cos_list = np.concatenate([np.linspace(1., 0.99, 1025), np.linspace(0.99, 0.95, 1535-1024+2)[1:], np.linspace(0.95, 0.75, 1791-1536+2)[1:], np.linspace(0.75, 0., 2047-1792+2)[1:]]) 
    tan_list13 = np.linspace(0., 1., 257)
    cot_list13 = np.linspace(1., 0., 257)
    tan_list24 = np.linspace(-1., 0., 257)
    cot_list24 = np.linspace(0., -1., 257)

    for ii in range(len(delta_phi)):
        if ii <= 255:
            delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(tan_list13[ii+1])-atan(tan_list13[ii])))
        elif 255 < ii <= 511:
            if cot_list13[ii%256+1] != 0:
                delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(1/cot_list13[ii%256+1])-atan(1/cot_list13[ii%256])))
            else:
                delta_phi[ii] = delta_phi[ii+1024] = 90 - degrees(atan(1/cot_list13[ii%256]))
        elif 511 < ii <= 767:
            if cot_list24[ii%256] != 0:
                delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(1/cot_list24[ii%256+1])-atan(1/cot_list24[ii%256])))
            else:
                delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(1/cot_list24[ii%256+1]))) - 90
        elif 767 < ii <= 1023:
            delta_phi[ii] = delta_phi[ii+1024] = degrees(abs(atan(tan_list24[ii%256+1])-atan(tan_list24[ii%256])))
        elif ii >= 1024:
            break

    for ii in range(2048):
        solid_angles[ii] = abs(cos_list[ii+1]-cos_list[ii]) * delta_phi

    table = 4*pi*table/(0.8910*1e6*solid_angles)
    # table = 4*pi*table/(1e6*solid_angles)
    return abs(table)  # eliminate -0.0 entries

def main():
    """Bin the photons into a 2048*2048 matrix according to cos_theta and tan_phi/cot_phi. Then normalize the table
    """
    # check the mode of the input file
    input_file = sys.argv[1]
    if re.search('lmf', input_file):
        data = read_listmode(sys.argv[1], sys.argv[3], sys.argv[4])
        table = np.zeros((2048, 512*4))

        for photon in data:
            index = ARF_table(photon)
            table[index[0], index[1]] += photon.weight
    
        table = normalize_table(table)
        np.savetxt(sys.argv[2]+'.txt',table,fmt='%.5f')

    elif re.search('bim', input_file):
        size1 = int(sys.argv[3])
        size2 = int(sys.argv[4])
        data = read_projection(sys.argv[1], size1, size2)
        table = np.zeros((2048, 512*4))

        for photon in data:
            # if photon.weight > 0:
            index = ARF_table(photon)
            table[index[0], index[1]] += photon.weight
    
        table = normalize_table(table)
        np.savetxt(sys.argv[2]+'.txt',table,fmt='%.5f')
        
if __name__ == "__main__":
    cProfile.run('main()',sys.argv[2]+'.log')
    # main()