'''
Created on Apr 15, 2015
Class containing the Data matrix
@author: Harald
'''

import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import re
from scipy.stats import chisquare
from matplotlib.widgets import Slider
from scipy.stats import binned_statistic
from scipy.stats import binom
from time import time

class Data(object):
    '''
    Class containing the Data matrix and methods to manipulate it
    '''
    header = []
    data = []  # Field for data
    subdata = []  # Field for data under Analysis
    data_p_mean = []  # Field for matrix of estimated mean allele Frequency
    x_cords_ind = 0  # Where to find x coordinates (Integer)
    y_cords_ind = 0  # Where to find y coordinates (Integer)
    color_ind = 0  # Where to find color
    id_ind = 0  # Where to find Plant ID
    sNP_inds = [2, 181]  # Which columns to find SNPs (python indexing)
    # sNP_inds = [1, 4 * 26 + 14]  # For the SNPs in Model Course
    sNP_okay = np.array([0])  # Whether a SNP is valid: 1 for good; 0 for out of range, -1 for defect marker
    p = []  # Array of allele frequencies 
    del_list = []
    del_list_SNP = []
    max_differences = 10  # Number of maximal SNP differences for individual to be deleted in double detection
    max_failed_snps = 8
    coords = []  # Placeholder for coordinates

    
    
    def __init__(self, path, snp_inds=0):
        if snp_inds: self.sNP_inds = snp_inds  # In case SNP-indices given; set them
        
        data = np.genfromtxt(path, delimiter=',', dtype=None)
        self.header = data[0, :]
        self.data = data[1:, :]
        self.clean_data()  # Clean all the not genotyped stuff and missing coordinates data
        
        self.y_cords_ind = np.where(self.header == "DistNorthofCentre")[0][0]   
        self.x_cords_ind = np.where(self.header == "DistEastofCentre")[0][0]
        
        self.sNP_okay = np.array([0 for _ in range(0, len(self.header))])  # Set default set of sNP_okay array to 0
        self.sNP_okay[self.sNP_inds[0]:(self.sNP_inds[1] + 1)] = 1  # Flag valid SNPs
        
        self.subdata = self.data
        self.p = np.array([0.5 * np.mean(self.subdata[:, i].astype('float')[self.subdata[:, i].astype('float') > -8])
                            for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)]).astype('float')  # Calculate mean allele frequencies for every working SNP
        self.data_p_mean = 2.0 * np.array([self.p for _ in self.subdata[:, 0]])  # Create Default p_mean matrix
        self.id_ind = np.where(self.header == "ID")[0][0]
        self.color_ind = np.where(self.header == "PhenoCat")[0][0]  # Extract color index
        self.coords = self.subdata[:, self.x_cords_ind:self.y_cords_ind + 1].astype('float')  # Load coordinates
        
    def set_subdata(self, xlim_down, xlim_up, ylim_down, ylim_up):
        ''' Extract data from given rectangle by using mask'''
        xcords = (self.data[:, self.x_cords_ind].astype(float) > xlim_down) & (self.data[:, self.x_cords_ind].astype(float) < xlim_up)
        ycords = (self.data[:, self.y_cords_ind].astype(float) > ylim_down) & (self.data[:, self.y_cords_ind].astype(float) < ylim_up)
        cords = np.nonzero(xcords & ycords)[0]
        self.subdata = self.data[cords, :]  # Cut out the right subdata
        self.data_p_mean = self.data_p_mean[cords, :]  # Cut out mean allele frequencies
        self.coords = self.subdata[:, self.x_cords_ind:self.y_cords_ind + 1].astype('float')  # Load  according coordinates
        
    def add_subdata(self, xlim_down, xlim_up, ylim_down, ylim_up):
        xcords = (self.data[:, self.x_cords_ind].astype(float) > xlim_down) & (self.data[:, self.x_cords_ind].astype(float) < xlim_up)
        ycords = (self.data[:, self.y_cords_ind].astype(float) > ylim_down) & (self.data[:, self.y_cords_ind].astype(float) < ylim_up)
        cords = np.nonzero(xcords & ycords)[0]
        self.subdata = np.append(self.subdata, self.data[cords, :], axis=0)
    
    def faulty_inds(self, max_failed):
        '''Shows and eliminates individuals with too many missing SNPs'''
        
        # Get a list with the number of missing snps
        missing_snps = np.array([np.sum(self.subdata[i, self.sNP_inds[0]:self.sNP_inds[1]].astype('float') == -9) for i in range(0, len(self.subdata[:, 0]))])
        print(missing_snps)
        plt.figure()
        plt.hist(missing_snps, range=[0, max(missing_snps)], bins=max(missing_snps))
        plt.show()
        
        # Delete failed inds
        ind_okay = np.where(missing_snps <= max_failed)[0]
        self.subdata = self.subdata[ind_okay, :]
        # self.data = self.subdata[ind_okay,:]
        
                    
    def SNP_analysis(self, i):
        ''' Do analysis of the ith SNP from Subdata'''
        arr_snp = self.subdata[:, i].astype("float")  # Load SNPS
        
        valid_ind = np.where(arr_snp > -8)  # Find measurements with valid SNPs.
        arr_snp = arr_snp[valid_ind] / 2.0
        arr_xcords = self.subdata[valid_ind[0], self.x_cords_ind].astype("float")
        arr_ycords = self.subdata[valid_ind[0], self.y_cords_ind].astype("float")
        
        return (arr_xcords, arr_ycords, arr_snp)
        
    def plot_SNP_slider(self):
        '''Does a Plot of SNPs in subdata with slider'''  
        fig = plt.figure()
        ax = plt.subplot(311)
        ax1 = plt.subplot(312)
        ax2 = plt.subplot(313)
        
        fig.subplots_adjust(left=0.25, bottom=0.25)
 
        # Do first Image:
        [arr_xcords, arr_ycords, arr_snp] = self.SNP_analysis(self.sNP_inds[0])  # Get SNP-Data for individuals with valid SNP-Entries
        distance_mean, _, _ = binned_statistic(arr_xcords, arr_snp, bins=10, statistic='mean')  # Extract data for Allele Frequency Plot
        
        ax.set_xlabel('x-Coord')
        ax.set_ylabel('y_Coord')
        ax.scatter(arr_xcords, arr_ycords, c=arr_snp, alpha=0.6)

        ax1.plot(np.arange(10), distance_mean, 'ro-')
        ax1.set_ylim([0, 1])
        
        ax2.hist(arr_snp, bins=3, range=[-0.01, 1.01], alpha=0.6)

        # Define slider
        axcolor = 'lightgoldenrodyellow'
        bx = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        
        slider = Slider(bx, 'SNP: ', self.sNP_inds[0], self.sNP_inds[1], valinit=0, valfmt='%i')
        
     
        def update(val):
            [arr_xcords, arr_ycords, arr_snp] = self.SNP_analysis([int(val)])  # Exctract Data
            distance_mean, _, _ = binned_statistic(arr_xcords, arr_snp, bins=10, statistic='mean')  # Extract data for Allele Frequency Plot
            counts, _ = np.histogram(arr_snp, bins=3, range=[-0.001, 1.001])  # Calculate counts for 0,1,2 states
            counts_freq = counts / float(len(arr_snp))
            p = np.mean(arr_snp)  # Calculate Allele Frequency
            hw_freq = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])  # Calculate Expected Frequencies from HW
            
            _, pval = chisquare(counts, hw_freq * len(arr_snp), ddof=1)
            
            ax.cla()
            ax.scatter(arr_xcords, arr_ycords, c=arr_snp, alpha=0.6)
            ax.set_xlabel('x-Coord')
            ax.set_ylabel('y_Coord')
            
            ax1.cla()
            ax1.plot(np.arange(10), distance_mean, 'ro-')
            ax1.set_ylim([0, 1])
            ax1.text(0.1, 0.8, "Mean Allele Frequency: %.2f" % p)
            
            ax2.cla()
            weights = [1.0 / len(arr_snp) for i in range(0, len(arr_snp))]  # @UnusedVariable
            ax2.bar([0, 1, 2], counts_freq, 0.4, alpha=0.6, color='r', label="Measured Genotypes")
            ax2.bar([0.4, 1.4, 2.4], hw_freq, 0.4, alpha=0.6, color='b', label="Expected in HW")
            ax2.set_ylim([0, 1])
            ax2.legend()
            ax2.set_xticks([0.4, 1.4, 2.4], ('00', '01', '11'))
            ax2.text(0.1, 0.8, "Xi^2 P-Value: %.4f " % pval)
            plt.title("SNP: " + self.header[val])
            
            fig.canvas.draw()
            
        slider.on_changed(update)
     
        plt.show()

    def SNP_cleaning(self, min_maf=0.05, pval=0.01):
        '''Do quality analysis of SNPs and flag non suitable ones'''
        
        p = np.array([np.mean(self.SNP_analysis(i)[2]) for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)])  # Get mean frequencies
        # p_vals = np.array([self.chi2(i) for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)])  # Get p-Values
        p_vals = np.array([self.chi2_geo_struc(i) for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)])  # Get p-Values WITH EXPLICIT GEOGRAPHIC STRUCTURE
        
        print("SNPs with p-Val<%.5f:" % pval)
        print(np.where(p_vals < pval)[0] + self.sNP_inds[0])
        # Visualize results
        plt.figure()
        ax = plt.subplot(211)
        ax1 = plt.subplot(223)
        ax2 = plt.subplot(224)
        
        ax.hist(p, range=[0, 1], bins=10, alpha=0.8, color='g')  # Allele Frequency spectrum
        ax.set_title("Allele Frequency spectrum")
        ax1.hist(p_vals, range=[0, 1], bins=20, alpha=0.8, color='g')  # P-value Distribution
        ax1.set_title("P-Values")
        ax2.plot(p_vals, 'ro')
        plt.show()
        
        # Flag the failed SNPs
        fail1 = np.where(p_vals < pval)[0] + self.sNP_inds[0]  # Failed for HW
        fail2 = np.where((p < min_maf) | (p > (1 - min_maf)))[0] + self.sNP_inds[0]  # Failed for too low MAF
        self.sNP_okay[fail2] = -2  # Flag Entries with too low MAF
        self.sNP_okay[fail1] = -1  # Flag SNPs which cannot be trusted because HW
        print("Failed SNPs for HW: ")
        print(self.header[self.sNP_okay == -1])
        print("Failed SNPs for MAF: ")
        print(self.header[self.sNP_okay == -2])
        
        self.del_list_SNP = np.where((self.sNP_okay == -1) | (self.sNP_okay == -2))[0]
        print("Total SNPs flagged: %.1f" % len(self.del_list_SNP))
        
        
    def chi2(self, SNP_ind):
        ''' Does chi2 analysis of SNP_indth SNP and return p-Value'''
        arr_snp = self.SNP_analysis(SNP_ind)[2]
        counts, _ = np.histogram(arr_snp, bins=3, range=[-0.001, 1.001])  # Calculate counts for 0,1,2 states
        p = np.mean(arr_snp)  # Calculate Allele Frequency
        hw_freq = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])  # Calculate Expected Frequencies from HW    
        _, pval = chisquare(counts, hw_freq * len(arr_snp), ddof=1)
        return pval
    
    def chi2_geo_struc(self, SNP_ind):
        ''' Do chi2 analysis of SNP-indth SNP and return p-Value
        this time take geographic structure into account (via p_mean)'''
        arr_snp = self.subdata[:, SNP_ind].astype("float")  # Load SNPS
        valid_ind = np.where(arr_snp > -8)[0]  # Find measurements with valid SNPs.
        print(len(valid_ind))
        arr_snp = arr_snp[valid_ind] / 2.0  # Go to actual allele frequency and valid SNPs
        p_mean = self.data_p_mean[valid_ind, SNP_ind - self.sNP_inds[0]] / 2.0  # Extract the right p_mean values
        
    
        counts, _ = np.histogram(arr_snp, bins=3, range=[-0.001, 1.001])  # Calculate counts for 0,1,2 states
        
        hw_freq = np.array([[(1 - p) ** 2, 2 * p * (1 - p), p ** 2] for p in p_mean])  # Calculate expected numbers assuming HW   
        hw_freq = np.sum(hw_freq, 0)
        _, pval = chisquare(counts, hw_freq, ddof=1)
        return pval
        
        
        
    def double_elemination(self):
        '''Tries to eliminate all potential double entries'''
        # Flag and delete double entries:
        
        for i in range(1, len(self.subdata[:, 0])):
            other_SNPs = self.subdata[0:i, self.sNP_inds[0]:self.sNP_inds[1]].astype('float')  # Matrix of sNPs to compare
            self_SNPs = self.subdata[i, self.sNP_inds[0]:self.sNP_inds[1]].astype('float')  # Own SNPs
            
            diffs = other_SNPs - self_SNPs  # Vector whether SNPs are different
            bool_diffs = (diffs != 0)  # Create Vector with boolean entries whether stuff is different 
            diffs = np.sum(bool_diffs, axis=1)  # Gives array of pairwise differences to individual i
            sig_diffs = np.where(diffs < (2 * self.max_failed_snps + self.max_differences))  # Where there are big enough differences warranting further investigation

            for j in sig_diffs[0]:
                self.double_check(i, j, self.max_differences)  # For every candidate do the actual double_check
            # for j in range(0, i):
                # self.double_check(i, j, max_differences=10)
        
        print("Deleting %.1f entries" % len(self.del_list))
        self.del_entries(self.del_list)  # Delete double entries
        self.del_list = []  # Empty deletion list.
                
    def double_check(self, i, j, max_differences=10):
        ''' Return 1 if ith and jth entry are identical, 0 else'''
        differences = [self.subdata[i, k] != self.subdata[j, k] for k in range(self.sNP_inds[0], self.sNP_inds[1])]  # Vector wether entries different
        opposing_homozygotes = [((float(self.subdata[i, k]) - float(self.subdata[j, k])) in [-2, 2]) for k in range(self.sNP_inds[0], self.sNP_inds[1])]  # Vector with number of opposing homozygotes
        
        trusted_entry1 = [self.subdata[i, k] != "-9" for k in range(self.sNP_inds[0], self.sNP_inds[1])]  # Binary vector whether entry1 trusted
        trusted_entry2 = [self.subdata[j, k] != "-9" for k in range(self.sNP_inds[0], self.sNP_inds[1])]  # Binary vector whether entry2 trusted

        true_differences = np.array(trusted_entry1) * np.array(trusted_entry2) * np.array(differences)

        if sum(true_differences) < max_differences:  # 10 seems to be a good value
            print("\nSubject %.1f and subject %.1f share %.1f different genotypes" % (i, j, sum(true_differences)))
            print("Opposing Homozygotes: %.1f" % sum(opposing_homozygotes))
            delta = [self.subdata[i, self.x_cords_ind].astype(float) - self.subdata[j, self.x_cords_ind].astype(float), self.subdata[i, self.y_cords_ind].astype(float) - self.subdata[j, self.y_cords_ind].astype(float)]
            print("Geographic distance: %.3f" % np.linalg.norm(delta))  # Give out Euclidean Norm   
            print("Individual 1:" + str(self.subdata[i, 0]))
            print("Individual 2: " + str(self.subdata[j, 0])) 
            self.del_list.append(j)  # Append to deletion-list
            
    def del_entries(self, del_list):
        '''Deletes given entries from subdata'''
        self.subdata = np.delete(self.subdata, del_list, axis=0)
        
    def del_SNPs(self):
        '''Deletes the given SNPS from the analysis and corrects necessary fields'''
        print("Deleting %.1f SNPs." % len(self.del_list_SNP))
        self.subdata = np.delete(self.subdata, self.del_list_SNP, axis=1)
        self.header = np.delete(self.header, self.del_list_SNP, axis=0)
        self.sNP_okay = np.delete(self.sNP_okay, self.del_list_SNP, axis=0)
        self.p = np.delete(self.p, self.del_list_SNP, axis=0)
        self.sNP_inds = [self.sNP_inds[0], self.sNP_inds[1] - len(self.del_list_SNP)]
        self.x_cords_ind -= len(self.del_list_SNP)
        self.y_cords_ind -= len(self.del_list_SNP)
        self.del_list_SNP = []
    
    def create_halfsibs(self, sNPs1, sNPs2, sNPs3):
        '''Pairs ind1 with ind2 and ind3 and returns SNPs'''
        new_SNPs = [-9 for _ in range(0, len(sNPs1))]  # Initialize new SNP, default is error
        new_SNPs1 = [-9 for _ in range(0, len(sNPs1))]
        
        for k in range(0, len(sNPs1)):
            if sNPs1[k] + sNPs2[k] + sNPs3[k] > -1:
                new_SNPs[k] = (np.random.binomial(1, sNPs1[k] / 2.0) + np.random.binomial(1, sNPs2[k] / 2.0))
                new_SNPs1[k] = (np.random.binomial(1, sNPs1[k] / 2.0) + np.random.binomial(1, sNPs3[k] / 2.0))
        return np.array[new_SNPs, new_SNPs1]
                           
    def update_p(self):
        '''Calculates new Allele frequencies for subdata'''
        a = self.p
        self.p = np.array([0.5 * np.mean(self.subdata[:, i].astype('float')[self.subdata[:, i].astype('float') > -8]) 
                           for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)]).astype('float')  # Calculate mean allele frequencies for every working SNP
        
        self.data_p_mean = np.array([2.0 * self.p for _ in self.subdata[:, 0]])
        plt.figure()
        plt.plot(a, 'ro')
        plt.plot(self.p, 'bo')
        plt.show()
        
    def compare_pass_p(self):
        ''' Compares Allele Frequencies from over the pass with the ones on this side of the pass'''
        self.set_subdata(-20000, -8000, -10000, 10000)
        p_pass = np.array([0.5 * np.mean(self.subdata[:, i].astype('float')[self.subdata[:, i].astype('float') > -8]) 
                           for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)]).astype('float')  # Calculate mean allele frequencies for every working SNP
        self.set_subdata(-5000, 3000, -10000, 2000)
        p_main = np.array([0.5 * np.mean(self.subdata[:, i].astype('float')[self.subdata[:, i].astype('float') > -8]) 
                           for i in range(self.sNP_inds[0], self.sNP_inds[1] + 1)]).astype('float')  
        
        # Plot:
        plt.figure()
        plt.plot(p_pass, 'ro', label="Pass Frequencies")
        plt.plot(p_main, 'bo', label="Main core Frequencies")
        plt.legend()
        plt.grid()
        plt.show()
    
    def ld_check(self, r2_score):
        '''Checks for LD via correlation'''
        good_SNPs = np.where(self.sNP_okay == 1)[0]
        data = self.subdata[:, good_SNPs].astype('float')  # Load data
        names = self.header[good_SNPs]
        
        masked_data = np.ma.array(data, mask=data == -9)  # Mask the missing data values
        r = np.ma.corrcoef(masked_data, rowvar=0)
        print("Mean Correlation: %.4f: " % np.mean(r))
        print("Standard Deviation: %.4f:" % np.std(r))
        
        for i in range(0, len(data[0, :])):
            r[i, i] = 0.15
        
        plt.figure()
        plt.pcolor(r)
        plt.xlim(0, len(names))
        plt.ylim(0, len(names))
        plt.colorbar()
        plt.xticks(np.arange(len(names)) + 0.5, names, rotation='vertical')
        plt.tick_params(labelsize=6)
        plt.show()
        
        paircorr = np.square([r[i, i + 1] for i in range(0, len(data[0, :]) - 1)])
        paircorr1 = np.square([r[i, i + 2] for i in range(0, len(data[0, :]) - 2)] + [0])  # 2 Individuals ahead
        paircorr2 = np.square([r[i, i + 3] for i in range(0, len(data[0, :]) - 3)] + [0, 0])  # 3 Individuals ahead
        paircorr = np.fmax(paircorr, paircorr1)  # Calculate Maximum
        paircorr = np.fmax(paircorr, paircorr2)
        plt.plot(paircorr, 'ro')
        plt.show()
        
        print("%.1f SNPs flagged." % np.sum(paircorr > r2_score))
        sNP_okay = np.where(self.sNP_okay == 1)[0]
        flag_ind = sNP_okay[paircorr > r2_score]
        self.sNP_okay[flag_ind] = -3  # Flag Subset of good SNPs with high correlation
        print(self.sNP_okay)
    
    def color_correlation(self, r2_score):
        '''Gives the correlation of all loci with the color.
        Flags loci with too high correlation'''
        # Define a proper color vector...
        color_ind = np.where(self.header == "Red")[0][0]
        color = self.subdata[:, color_ind]
        print(color)
        
        # Check where only color values
        is_float = np.array([1 if re.match("^\d+(\.\d+)*$", i) else 0 for i in color]).astype(bool)  # Only entries without minus and numbers or dots
        

    
        color_cor = np.array([0 for _ in range(self.sNP_inds[0], self.sNP_inds[1] + 1)]).astype('float')
        for i in range(self.sNP_inds[0], self.sNP_inds[1]):
            inds = np.where((self.subdata[:, i] > -1) & (is_float == 1))[0]  # Where the SNPs are okay
            
            p1 = self.subdata[inds, i].astype('float')
            c1 = color[inds].astype('float')
            color_cor[i - self.sNP_inds[0]] = np.corrcoef(p1, c1)[0, 1]
            
        
        plt.figure()
        plt.plot(color_cor, 'ro')
        plt.xlabel("Locus")
        plt.ylabel("Pearson R")
        plt.show()  
        
        color_cor = np.square(color_cor)
        
        plt.hist(color_cor, range=[0, 0.1])
        plt.ylabel("R^2")
        plt.show()
        
        print("%.1f SNPs flagged." % np.sum(color_cor > r2_score))
        self.sNP_okay[np.where(color_cor > r2_score)[0] + self.sNP_inds[0]] = -4  # Flag Subset of good SNPs with high correlation to color
        
    def extract_good_SNPs(self):
        '''Extracts good SNPs: Returns SNP-Matrix, their allele frequency, their coordinates and name as numpy float arrays'''
        good_SNPs = np.where(self.sNP_okay == 1)[0]
        data = self.subdata[:, good_SNPs].astype('float')  # Load data
        p = self.p[good_SNPs - self.sNP_inds[0]].astype('float')  # Load allele frequencies
        coords = self.subdata[:, self.x_cords_ind:self.y_cords_ind + 1].astype('float')
        names = self.header[good_SNPs]
        color = self.subdata[:, self.sNP_inds[1] + 1].astype(float)
        return(data, p, coords, names, color)    
        
    def kernel_estimation(self, sigma=500):
        '''Given sigma, estimates the mean allele frequencies per individual'''
        self.coords = self.subdata[:, self.x_cords_ind:self.y_cords_ind + 1].astype('float')  # Load coordinates    # WORKAROUND - ADDRESS
        data, p, coords = self.subdata[:, self.sNP_inds[0]:(self.sNP_inds[1] + 1)].astype('float'), self.p, self.coords  # Load necessary Data
        
        # Impute missing genotypes
        for (ind, sNP), value in np.ndenumerate(data):  # Iterate over all data points
            if value == -9:  # In case of missing data
                data[ind, sNP] = np.random.binomial(2, p[sNP])  # Draw random genotype for missing data
        
        # Calculate Pairwise distance (Overkill - but only by a factor of 2 - I cannot be bothered.
        # dist_mat1 = np.array([[np.linalg.norm(coords[i, :] - coords[j, :]) for i in range(n)] for j in range(n)]) # The old way
        dist_mat = self.calc_distance_matrix(coords)
        # np.fill_diagonal(dist_mat, 10000)  # Set self distance to very far away (comp. trick instead of deleting it)
        
        
        # Calculate likelihood for various sigmas
#         y=[]
#         sigmas=[i for i in range(150,801,75)]
#         for sigma in sigmas:
#             print("Doing Sigma: %.0f " % sigma)
#             p_mean = np.array([[self.calc_p_mean(dist_mat[i, :], data[:, j], sigma) for j in range(len(p))] for i in range(n)])
#             y.append(self.calc_ll_pmean(data, p_mean/2.0))
#             
#         plt.figure()
#         plt.plot(sigmas,y,'r-')
#         plt.xlabel("Sigma")
#         plt.ylabel("Log Likelihood")
#         plt.show() 
        
        self.picture_kernel_smoothing(dist_mat, data, coords, p)  # Call picture function for smoothing the kernel

        p_mean = self.calc_p_mean(dist_mat, data, sigma)
        self.data_p_mean = p_mean
        
    def kernel_estimation_rare(self, sigma=500, rare_factor=20):
        '''Given sigma, estimates the mean allele frequencies per individual - from rarified data'''
        self.coords = self.subdata[:, self.x_cords_ind:self.y_cords_ind + 1].astype('float')  # Load coordinates
        data, p, coords = self.subdata[:, self.sNP_inds[0]:(self.sNP_inds[1] + 1)].astype('float'), self.p, self.coords  # Load necessary Data
        
        # Impute missing genotypes
        for (ind, sNP), value in np.ndenumerate(data):  # Iterate over all data points
            if value == -9:  # In case of missing data
                data[ind, sNP] = np.random.binomial(2, p[sNP])  # Draw random genotype for missing data
        
        dist_mat = self.calc_distance_matrix(coords)
        keep_inds = self.rarefy(dist_mat, rare_factor)
        
        plt.figure()
        plt.scatter(coords[keep_inds, 0], coords[keep_inds, 1])
        plt.title("Rarefied inds")
        plt.show()
        
        
        
        p_mean = self.calc_p_mean(dist_mat[:, keep_inds], data[keep_inds, :], sigma)
        self.data_p_mean = p_mean
        
    def rarefy(self, dist_mat, rare_dist, factor = 10):
        '''Rarefies according to distance_matrix; such that on average 1 ind in radius rare_dist
        Return list of rarefied individuals'''
        print("Rarefying...")
        near_dist_mat = dist_mat < rare_dist*factor
        near_nbrs = np.sum(near_dist_mat, axis=1)
        chance = np.random.random(len(near_nbrs))
        keep_inds = np.where(chance < (factor**2 / near_nbrs))[0]  # Print where the nearest neighbours are to be found
        print("Original inds: %i " % len(near_nbrs))
        print("Rarefied inds: %i " % len(keep_inds))
        return keep_inds

        
    def picture_kernel_smoothing(self, dist_mat, data, coords, p):
        '''Calculates Kernel smoothing with 3 kernel distances and give back a picture for every SNP'''
        # Calculate the mean matrix
        
        print("Calculate Mean allele frequs for every individual")
        
        
        sigma = 200
        p_mean200 = self.calc_p_mean(dist_mat, data, sigma)
        print(self.calc_ll_pmean(data, p_mean200 / 2.0))
        
        sigma = 500
        p_mean = self.calc_p_mean(dist_mat, data, sigma)
        print(self.calc_ll_pmean(data, p_mean / 2.0))

        sigma = 1000
        p_mean1000 = self.calc_p_mean(dist_mat, data, sigma)
        print(self.calc_ll_pmean(data, p_mean1000 / 2.0))
        
        print("Size of empirical Data: %i %i" % (len(data[0, :]), len(data[:, 0])))
        print("Size of estimated all-freq Data: %i %i" % (len(p_mean[0, :]), len(p_mean[:, 0])))
        
        fig = plt.figure()
        # Do first Image:
        l, = plt.plot(coords[:, 0], data[:, 0], 'ko', label='Real Allele Frequency')
        l200, = plt.plot(coords[:, 0], p_mean200[:, 0], 'yo', label='200m')
        l500, = plt.plot(coords[:, 0], p_mean[:, 0], 'ro', label='500m')
        l1000, = plt.plot(coords[:, 0], p_mean1000[:, 0], 'bo', label='1000m')
        plt.legend()
        
        plt.xlabel('x-Coord')
        plt.ylabel('Allele Freq')
        plt.ylim([-0.3, 2.3])

        # Define slider
        axcolor = 'lightgoldenrodyellow'
        bx = plt.axes([0.25, 0, 0.65, 0.03], axisbg=axcolor)
        plt.axes
        slider = Slider(bx, 'SNP: ', 0, len(p), valinit=0, valfmt='%i')
     
        def update(val):
            l.set_ydata(data[:, val])
            l200.set_ydata(p_mean200[:, val])
            l500.set_ydata(p_mean[:, val])
            l1000.set_ydata(p_mean1000[:, val])
            fig.canvas.draw_idle()
            
        slider.on_changed(update)
        plt.show()
        
#     def calc_p_mean_or(self, dists, p, sigma):
#         '''Given a distance matrix and the matrix of allele frequencies,
#         calculate the mean allele frequency'''
#         p, dists = np.array(p), np.array(dists)  # Just in case that not numpy array
#         
#         weights = 1 / (2.0 * np.pi * sigma ** 2) * np.exp(-dists ** 2 / (2.0 * sigma ** 2))  # Calculate the Gaussian weights
#         p_mean = np.dot(weights * p) / np.sum(weights)  # Calculate weighted mean
#         return(p_mean)
    
    def calc_p_mean(self, dist_mat, p_mat, sigma):
        '''Given a distance matrix and the matrix of allele frequencies,
        calculate the mean allele frequency'''
        print("Smoothing out...")
        start = time()
        p_mat, dist_mat = np.array(p_mat), np.array(dist_mat)  # Just in case that not numpy array
        
        weights = 1 / (2.0 * np.pi * sigma ** 2) * np.exp(-dist_mat ** 2 / (2.0 * sigma ** 2))  # Calculate the Gaussian weights
        p_mean = np.dot(weights, p_mat) / np.sum(weights, axis=1)[:, None]  # Calculate weighted mean
        print("Time taken %.2f" % (time() - start))
        return(p_mean)
    
    
    def calc_ll_pmean(self, x, p_mean):
        '''Calculates the likelihood of having this exact p-mean. x is empirical data - p_mean estimated mean'''
        print("Calculating likelihood")
        l_matrix = binom.pmf(x, 2, p_mean)
        ll = np.sum(np.log(l_matrix))
        return ll
    
    def calc_distance_matrix(self, coords):
        '''Calculate the distance matrix between all coords. Requires numpy array of 2d coordinates'''
        print("Calculating distance matrix")
        start = time()
        dist_mat = np.linalg.norm(coords[:, None] - coords, axis=2)
        print("Time taken %.2f" % (time() - start))
        return dist_mat
        
    def clean_data(self):
        '''Method that cleans data. Gets rid of "NA GPS Values
        and NGY genotypes'''
        data, header = self.data, self.header
        y_cords_ind = np.where(self.header == "DistNorthofCentre")[0][0]   
        x_cords_ind = np.where(self.header == "DistEastofCentre")[0][0]
        print("Raw data Nr. individuals: %i " % len(data[:, 0]))
        
        gps_good_ind = np.where((data[:, x_cords_ind] != 'NA') * (data[:, y_cords_ind] != 'NA'))[0]
        data = data[gps_good_ind, :]  # Delete Entries without GPS
        

        # Now to NGY errors
        
        no_gen = np.array([np.sum(data[:, i] == 'NGY') for i in np.arange(self.sNP_inds[0], self.sNP_inds[1] + 1)])
        
        print("Number of loci with <10 NGY SNPs: % i" % np.sum(no_gen < 10))
        bad_snps = np.where(no_gen > 10)[0]  # Detect bad loci; i.e. loci with more than 10 non-genotyped individuals
        
        data = np.delete(data, bad_snps + self.sNP_inds[0], axis=1)  # Delete bad SNPs
        header = np.delete(header, bad_snps + self.sNP_inds[0])  # Delete header SNPs
        self.sNP_inds[1] -= len(bad_snps)
        
        data[data[:, :] == 'NGY'] = -9  # Replace the few not genotyped individuals with failed SNPs
        self.data, self.header = data, header  # Save what was done
        
    def save_data(self, path):
        '''Method that pickles the data.'''
        temp = self.subdata
        self.data = self.subdata  # Save subdata
        self.subdata = []  # Temporary delete subdata
        pickle.dump(self, open(path, "wb"), protocol=2)  # Pickle the data
        self.subdata = temp  # Restore subdata
        
    def save_p(self):
        '''Method that saves allele-frequencies. Written to extract for Nick'''
        np.savetxt("allele_freqs.csv", self.p, delimiter=",")
        
        
        
        
        
        

    
    
    

        



        

        
