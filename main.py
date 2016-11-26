'''
Various methods to analyse the SNP-Data
@author: Harald
'''

from Data import Data  # @UnusedImport
from FlaggedSNPs import FlaggedSNPs
from extract_data import Extractor
import numpy as np  # @UnusedImport
import matplotlib.pyplot as plt
import cPickle as pickle  # @UnusedImport

# path, snP_inds = '../Data/SNP_data.csv', [2, 181]    # Where to load from
# path, snP_inds ='../Data/SNPsModelCourse.csv', [1, 4 * 26 + 14]
# path, snP_inds = '../Data/HZ2014.csv', [3, 9 * 26 + 13]

# Where to save the pickle file to
# path1 = '../Data/SNP_data.p'        # This is the SNP-data based on the cline-dataset    
# path1 = '../Data/SNP_data2013.p'    # This is the SNP-data based on the 2013 data-set
path1 = '../Data/SNP_data2014.p'    # SNP-data from 2014

min_maf = 0.1  # Defines the threshold for IBD-analysis
p_val = 0.0001  # Defines the critical p-Value for flagging loci according to HW 0.0005
r2_score = 0.03  # Max r2 for LD-Flagging
max_failed_snps = 8  # Max Nr of failed SNPs per individual
sigma = 3000  # Kernel width used for smoothing

# subregion = [-300, 500, -150, 150]  # Core Hybrid Zone
# subregion = [-1084, -300, 840, 1075]  # Yellow Flank 1
# subregion = [285, 1200, 200, 480]  # Magenta Flank 1
# subregion=[-20000,-8000,-10000,10000] #Over the pass
# subregion = [-50000, 100000, -100000, 11000]  # Everything
subregion = [-3000, 2000, -10000, 1100]  # Big Core
#subregion = [-2000,1500, -10000, 1100]   # Slight less than Big Core
# subregion = [-4200, 3000, -10000, 2000]  # Very Big Core
# subregion=[-2000,-500,-10000,1100]    # Core left half
# subregion=[500,2000,-10000,1100]    # Core right half
# subregion=[-500,500,-10000,600]    # Vertical core
    
def main():
    print("Welcome back buddy! \n Loading data...")
    #data = Data(path, snP_inds)  # Load original data!! (UNCOMMENT)
    data = pickle.load(open(path1, "rb"))  # Load data from pickle file!
    #data.data_p_mean = np.array([2.0 * data.p for _ in data.data[:, 0]])  # Create Default p_mean matrix NORMALLY DONE IN CONSTRUCTOR - HOTFIX
    data.set_subdata(-100000, 100000, -100000, 100000)  # Set everything as subdata
    flagged_SNPs = FlaggedSNPs(data)  # Load default flagged SNPs
    flagged_SNPs.maf = min_maf  # A bit of a workaround to set the MAF
    
    print("Nr of Columns: %i" % len(data.data[0, :]))
    
    while True:
        inp = int(input("\n What you wanna do? \n(1) Print geographic data \n(2) Extract sub-data \n(3) Data Cleaning  \n(4) Data Analysis "
                        "\n(5) Data Extraction \n(6) Exit\n"))
        
        if inp == 1:
            flagged_SNPs.plot_geography()
        
        elif inp == 2:
            data.set_subdata(subregion[0], subregion[1], subregion[2], subregion[3])  # Extract subdata
            flagged_SNPs = FlaggedSNPs(data)  # Update flagged SNPs
            data.save_p()
            # data.add_subdata(subregion1[0], subregion1[1], subregion1[2], subregion1[3])
            # data.add_subdata(subregion2[0], subregion2[1], subregion2[2], subregion2[3])
            print("Nr of subsamples %.2f:" % len(data.subdata[:, 0]))
        
            plt.figure()
            plt.scatter(data.subdata[:, data.x_cords_ind].astype(float), data.subdata[:, data.y_cords_ind].astype(float), alpha=0.5)
            plt.axis("equal")
            plt.grid()
            plt.show()
        
        elif inp == 3:
            while True:
                inp1 = int(input("\nWhat do you want to do?\n(1) Flag failed SNPs \n(2) Delete doubles \n(3) Update allele Frequencies "
                                 "\n(4) Detect faulty individuals \n(5) LD-Check \n(6) Save Data "
                                 "\n(7) Color Correlation \n(8) Kernel Estimation \n(9) Go back\n"))
                if inp1 == 1:
                    data.SNP_cleaning(min_maf, p_val)  
                if inp1 == 2:
                    data.double_elemination()
                if inp1 == 3:
                    data.update_p()       
                if inp1 == 4:
                    data.faulty_inds(max_failed_snps)
                if inp1 == 5:
                    data.ld_check(r2_score)
                if inp1 == 6:
                    data.save_data(path1)
#                     data.data = data.subdata  # Save subdata
#                     data.subdata = []
#                     pickle.dump(data, open(path1, "wb"))  # Pickle the data
                if inp1 == 7:
                    data.color_correlation(r2_score)  
                if inp1 == 8:
                    data.kernel_estimation(sigma)
                    data.kernel_estimation_rare(sigma=sigma)  
                if inp1 == 9:
                    flagged_SNPs = FlaggedSNPs(data)
                    break
                                
        
        elif inp == 4:
            while True:
                inp2 = int(input("\n(1) SNP-analysis \n(2) IBD-analysis \n(3) Y-M comparison \n(4) PC-Analysis \n(5) Permute colors "
                "\n(6) Local Covariances \n(7) Geographic Comparison \n(8) Extract Coancestry-Txts \n(9) Back\n"))
                if inp2 == 1:
                    data.plot_SNP_slider()
                if inp2 == 2:
                    geo_scale = int(input("\nOn what geographic scale? (m)\n"))
                    flagged_SNPs.analyze_correlations(geo_scale)
                if inp2 == 3:
                    geo_scale = int(input("\nOn what geographic scale? (m)\n"))
                    flagged_SNPs.ym_comparison(geo_scale)
                if inp2 == 4:
                    flagged_SNPs.principal_components()
                if inp2 == 5:
                    # data.compare_pass_p()   # Secret option
                    flagged_SNPs.permute_colors()
                if inp2 == 6:
                    # flagged_SNPs.homozygosity_analysis()
                    # flagged_SNPs.nbh_analysis(nbh_dist=75, min_nbh=5) 
                    flagged_SNPs.loc_all_freq_covariances()
                    flagged_SNPs.x_dist_cov()
                if inp2 == 7:  # Bonus Line
                    geo_scale = int(input("\nOn what geographic scale? (m)\n"))
                    flagged_SNPs.geo_comparison(geo_scale)
                if inp2 == 8: 
                    flagged_SNPs.save_coancestry_data()
                if inp2 == 9:
                    break
        
        elif inp == 5:
            extractor = Extractor(data)
            while True:
                inp2 = int(input("\n(1) Extract Coancestry-Txts\n(2) Extract simulated related pairs"
                                 "\n(3) Extract highly related pairs\n(4) Extract Data for Alex \n(5) Back \n"))
                if inp2 == 1:
                    extractor.save_coancestry_data()
                elif inp2 == 2:
                    extractor.extract_pairs(200, 200, 200, 200)  # nr PO, fullsibs, half_sibs, unrelated
                elif inp2 == 3:
                    extractor.extract_high_rel_pairs()
                elif inp2 == 4:
                    extractor.extract_snp_data()
                elif inp2 == 5:
                    break
                                       
        elif inp == 6:
            print("See u!!")
            break
    
    
if __name__ == '__main__':
    main()
    
    
