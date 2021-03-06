# Process simind listmode file into ARF table in MapReduce non-Parallel Processing
# Python3
# Jie (Laurie) Zhang
# 04/06/15
# e.g. python ARF_listmode_v3_copy.py *.lmf output_name lower_energy upper_energy parallel_num
import sys
import cProfile
import csv
import struct
from math import sqrt, atan, degrees, pi
import numpy as np
from multiprocessing import Pool

class PhotonListMode(object):
    """
    Listmode photon: 'X0','Y0','Z0','XPHANT','YPHANT','ZPHANT','XCRYSTAL','YCRYSTAL','ZCRYSTAL','Energy in Crystal','Photon Weight','Scatter Order
    '"""

    def __init__(self, locations, energy, weight, scatter):
        self.X0 = locations[0]
        self.Y0 = locations[1]
        self.Z0 = locations[2]

        self.Xp = locations[3]
        self.Yp = locations[4]
        self.Zp = locations[5]

        self.Xc = locations[6]
        self.Yc = locations[7]
        self.Zc = locations[8]

        self.energy = energy
        self.weight = weight
        self.scatter = scatter

        self.X_vec = self.Xc-self.X0
        self.Y_vec = self.Yc-self.Y0
        self.Z_vec = self.Zc-self.Z0

    def mod(self):
        travel_dist = sqrt(self.X_vec**2 + self.Y_vec**2 + self.Z_vec**2)
        return travel_dist

    def quadrant(self):
        quadrant = 0
        if (self.X_vec >= 0 and self.Y_vec > 0):
            quadrant = 1
        elif (self.X_vec < 0 and self.Y_vec >= 0):
            quadrant = 2
        elif (self.X_vec <= 0 and self.Y_vec < 0):
            quadrant = 3
        elif (self.X_vec > 0 and self.Y_vec <= 0):
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

    def ARF_table(self, cos_list, tan_list13, cot_list13, tan_list24, cot_list24):
        """Find the correct position in the table"""
        theta_ind = phi_ind = None

        if self.quadrant() == 0:
            theta_ind = phi_ind = 0
        else:
            for ii in cos_list:
                if ii <= self.cos_theta():
                    theta_ind = max(int(np.argwhere(cos_list == ii)-1), 0)
                    break

            if abs(self.tan_phi()) <= 1:
                if self.quadrant() == 1:
                    for ii in tan_list13:
                        if ii >= self.tan_phi():
                            phi_ind = max(int(np.argwhere(tan_list13 == ii)-1), 0)
                            break
                elif self.quadrant() == 2:
                    for ii in tan_list24:
                        if ii >= self.tan_phi():
                            phi_ind = int(1*512 + 255 + np.argwhere(tan_list24 == ii))
                            break
                elif self.quadrant() == 3:
                    for ii in tan_list13:
                        if ii >= self.tan_phi():
                            phi_ind = int(2*512 -1 + np.argwhere(tan_list13 == ii))
                            break
                elif self.quadrant() == 4:
                    for ii in tan_list24:
                        if ii >= self.tan_phi():
                            phi_ind = int(3*512 + 255 + np.argwhere(tan_list24 == ii))
                            break
            elif abs(self.cot_phi()) <= 1:
                if self.quadrant() == 1:
                    for ii in cot_list13:
                        if ii <= self.cot_phi():
                            phi_ind = int(255 + np.argwhere(cot_list13 == ii))
                            break
                elif self.quadrant() == 2:
                    for ii in cot_list24:
                        if ii <= self.cot_phi():
                            phi_ind = int(1*512 -1 + np.argwhere(cot_list24 == ii))
                            break
                elif self.quadrant() == 3:
                    for ii in cot_list13:
                        if ii <= self.cot_phi():
                            phi_ind = int(2*512 + 255 + np.argwhere(cot_list13 == ii))
                            break
                elif self.quadrant() == 4:
                    for ii in cot_list24:
                        if ii <= self.cot_phi():
                            phi_ind = int(3*512 -1 + np.argwhere(cot_list24 == ii))
                            break
        if theta_ind == None or phi_ind == None:
            print(self.quadrant(), self.cos_theta(), self.tan_phi(), self.cot_phi())
        return (theta_ind, phi_ind)

def Map(photons, cos_list, tan_list13, cot_list13, tan_list24, cot_list24):
    """
    Given a list of photons, return a list of tuples containing the positions of the photon in the ARF_table, and the photon weight: (position, weight)
    """
    results = []
    for photon in photons:
        results.append([photon.ARF_table(cos_list, tan_list13, cot_list13, tan_list24, cot_list24), photon.weight])
    return results

def Partition(photon_list):
    """ Group the sublist of [position, photon weight] pairs into a position-weight-list map, so that the Reduce operation later can work on summing all the weights. The returned result is a dictionary with the structure {position: [weight1, weight2, ...] .. }
    """
    pw = {}
    for photon in photon_list:
        # print(photon)
        try:
            pw[photon[0]].append(photon[1])
        except KeyError:
            pw[photon[0]] = [photon[1]]
    return pw

def Reduce(Mapping):
    """ Given a (position: [weight1, weight2, ...]) tuple, collapse all the count tuples from the Map orientaion into a single total weight for this position, and return a final tuple (position, total_weight).
    """
    return (Mapping[0], sum(weight for weight in Mapping[1]))

def read_file(filename, lower_energy, upper_energy):
    """
    Read listmode data: 10 int16, 1 float, 1 interger*1
    """
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

def normalize_table(table, cos_list, tan_list13, cot_list13, tan_list24, cot_list24):
    """
    Normalize the ARF_table
    """
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
    data = read_file(sys.argv[1], sys.argv[3], sys.argv[4])

    """
    Bin the photons into a 2048*2048 matrix according to cos_theta and tan_phi/cot_phi. Then normalize the table
    """
    cos_list = np.concatenate([np.linspace(1., 0.99, 1025), np.linspace(0.99, 0.95, 1535-1024+2)[1:], np.linspace(0.95, 0.75, 1791-1536+2)[1:], np.linspace(0.75, 0., 2047-1792+2)[1:]]) 
    tan_list13 = np.linspace(0., 1., 257)
    cot_list13 = np.linspace(1., 0., 257)
    tan_list24 = np.linspace(-1., 0., 257)
    cot_list24 = np.linspace(0., -1., 257)

    single_weight_list = Map(data, cos_list, tan_list13, cot_list13, tan_list24, cot_list24)

    organize_photon_dict = Partition(single_weight_list)

    position_weight = list(map(Reduce, organize_photon_dict.items()))

    # Sort the position_weight into a ARF table
    table = np.zeros((2048, 512*4))
    for weight in position_weight:
        table[weight[0][0], weight[0][1]] = weight[1]
 
    table = normalize_table(table, cos_list, tan_list13, cot_list13, tan_list24, cot_list24)
    np.savetxt(sys.argv[2]+'.txt',table,fmt='%.5f')
        
if __name__ == "__main__":
    cProfile.run('main()',sys.argv[2]+'.log')
    # main()