# File: DecodeRNASequences.py
# Decoding RNA sequences. Translated from matlab
# Petter Ranefall 2014

# Import neccessary modules:
import sys
import numpy as np
import time
import os
import re
import copy

class DecodeRNASequences(object):
    """
    Read, register and tile large images from e.g. slide scanner, and create a .csv file as input for CellProfiler
    """
    
    def __init__(self):
        pass

    def run(self, Output_path, \
            CPcsv_name, \
            taglist_file, \
            output_prefix, \
            num_hybs = 4, \
            quality_Th_D = 0.5, \
            DO_intensity_Th_D = 0.005, \
            verbose = False):
        """
        Input parameters:
        num_hybs: Number of hybs (DEFAULT = 4)
        quality_Th_D: Quality threshold (DEFAULT = 0.5) # minimum relevant is 1/num_hybs
        DO_intensity_Th_D: (DEFAULT = 0.005) # minimum relevant is 0.05 as this is the threshold value applied at detection
        Input_path_Csv_Img: (DEFAULT = '.') #folder with CP Csv file and images used as input
        Output_path
        CPcsv_name = 'DefaultOUT_blobs.csv'
        output_prefix = 'OT_130913_slideD'
        verbose: True if progress strings should be written. (DEFAULT = True)

        Output:
        csv file

        """

        start = time.clock()
        prev = start
        Output_path = Output_path.replace('\\', '/')
        if not Output_path.endswith('/'):
            Output_path = Output_path + '/'
        self.decode(num_hybs, quality_Th_D, DO_intensity_Th_D, CPcsv_name, Output_path, output_prefix, taglist_file)
        
        # Done!
        if verbose:
            stop = time.clock()
            print("Total time %s s" % str(stop-start))
        pass

    def decode(self, num_hybs, quality_Th_D, DO_intensity_Th_D, CPcsv_name, Output_path, output_prefix, taglist_file):
        #-----------------------------------------------
        #read the taglist file
        taglist = np.genfromtxt(taglist_file, delimiter=',', dtype=None)
        
        #-----------------------------------------------
        #read the csv file
        seq = np.genfromtxt(CPcsv_name, delimiter=',')
        lenSeq = len(seq[:,0])
        seq = seq[1:lenSeq,:] # Remove header
        num_blobs = len(seq[:,1])/num_hybs # 2nd column show object number, but the counter is set to zero for each new tile
        tile_ID = seq[:,5] # 3rd column shows tile ID
        tiles = np.unique(tile_ID)
        parent_cell = seq[:,14]

        x_pos = np.zeros(num_blobs)
        y_pos = np.zeros(num_blobs)
        seq_res = np.zeros((num_blobs, num_hybs))
        seq_strength = np.zeros((num_blobs, num_hybs))
        seq_quality = np.zeros((num_blobs, num_hybs))
        seq_cell_ID = np.zeros((num_blobs, 1))

        #------------------------------------------
        # set thresholds 
        #------------------------------------------

        DO_intensity_T = DO_intensity_Th_D # minimum relevant is 0.05 as this is the threshold value applied at detection

        quality_T = quality_Th_D # minimum relevant is 1/num_hybs

        #-------------------------------------------

        start_blob_ID=0
        start_cell_ID=0

        for t in range(0,len(tiles)):
            tt = tiles[t]
            tile_seq_data = seq[seq[:,5] == tt,:]
            tile_num_blobs = int(max(tile_seq_data[:,1]))
            tile_num_cells = int(max(tile_seq_data[:,14]))

            for i in range(0, tile_num_blobs):                
                pp = tile_seq_data[tile_seq_data[:,1] == i+1, 8:12] # collects data from all hybsteps of given blob tile
                if len(pp) == num_hybs:
                    hyb0row = tile_seq_data[tile_seq_data[:,1] == i+1, :]
                    hyb0row = hyb0row[hyb0row[:,2] == "hyb0", :]
                    x = hyb0row[3] + hyb0row[12]
                    y = hyb0row[4] + hyb0row[13]
                    x_pos[i+start_blob_ID] = x
                    y_pos[i+start_blob_ID] = y
                    #remove NaN
                    f = np.isnan(pp)
                    pp[f] = 0
                    #remove negative_values
                    f = np.where(pp<0)
                    pp[f] = 0
                    a = []
                    b = []
                    for j in range(0, num_hybs):
                        maxInd = np.where(pp[j,:] == max(pp[j,:]))[0]
                        a.append(pp[j,maxInd[0]])
                        b.append(maxInd[0])
                    # [a b] = max(pp(1:num_hybs,:),[],2);   %strongest letter at each hyb step

                    for j in range(0, num_hybs):
                        seq_res[i+start_blob_ID,j] = b[j]
                        seq_strength[i+start_blob_ID,j] = tile_seq_data[tile_seq_data[:,1] == i+1, 6][0] # DO intensity
                        sumpp = sum(pp[j, 0:num_hybs])
                        if sumpp < 1e-6:
                            seq_quality[i+start_blob_ID,j] = 0
                        else:
                            seq_quality[i+start_blob_ID,j] = a[j]/sumpp
                        
                    if tile_seq_data[i,14]>0:
                        seq_cell_ID[i+start_blob_ID] = start_cell_ID+tile_seq_data[i,14]
            
            start_blob_ID = start_blob_ID+tile_num_blobs
            start_cell_ID = start_cell_ID+tile_num_cells

            if num_hybs==1:
                seq_num=seq_res[:,0] + 1
            elif num_hybs==2:
                seq_num=seq_res[:,0]*10+seq_res[:,1] + 11
            elif num_hybs==3:
                seq_num=seq_res[:,0]*100+seq_res[:,1]*10+seq_res[:,2] + 111
            elif num_hybs==4:
                seq_num=seq_res[:,0]*1000+seq_res[:,1]*100+seq_res[:,2]*10+seq_res[:,3] + 1111
            elif num_hybs==5:
                seq_num=seq_res[:,0]*10000+seq_res[:,1]*1000+seq_res[:,2]*100+seq_res[:,3]*10+seq_res[:,4] + 11111
            elif num_hybs==6:
                seq_num=seq_res[:,0]*100000+seq_res[:,1]*10000+seq_res[:,2]*1000+seq_res[:,3]*100+seq_res[:,4]*10+seq_res[:,5] + 111111
            else:
                print('uncompatible number of hybridizations')


        f = np.isnan(seq_quality)
        seq_quality[f] = 0
        seq_quality_min = np.zeros(num_blobs)
        seq_strength_max = np.zeros(num_blobs)
        for j in range(0, num_blobs):
            seq_quality_min[j] = min(seq_quality[j, :])
            seq_strength_max[j] = max(seq_strength[j, :])

        #---------------------------------
        # all found transcripts before QT

        ffbt = []
        for j in range(0, num_blobs):
            if seq_quality_min[j]>0 and seq_strength_max[j]>0:
                ffbt.append(j)
        allbt = seq_num[ffbt] 
        x_allbt = x_pos[ffbt]
        y_allbt = y_pos[ffbt]
        par_cell_allbt = parent_cell[ffbt]
        seq_quality_min_allbt = seq_quality_min[ffbt]
        all_lettersbt = []
        for uubt in range(0,len(allbt)):
            all_lettersbt.append(self.num2letters(allbt[uubt], num_hybs))

        # all found transcripts after QT

        ff = []
        for j in range(0, num_blobs):
            if seq_quality_min[j]>quality_T and seq_strength_max[j]>DO_intensity_T:
                ff.append(j) # filter out sequences with low quality 
         
        all = seq_num[ff] 
        x_all = x_pos[ff]
        y_all = y_pos[ff]
        par_cell_all = parent_cell[ff]
        seq_quality_min_all = seq_quality_min[ff]
        all_letters = []
        for uu in range(0,len(all)):
            all_letters.append(self.num2letters(all[uu], num_hybs))

        # ------------------------------------------

        # unique found transcripts 
        hh = np.unique(seq_num[ff])
        seq_letters = []
        for u in range(0,len(hh)):
            seq_letters.append(self.num2letters(hh[u], num_hybs))
        binEdges = np.append(hh, hh[-1] + 0) # add rightmost edge
        [ab, b] = np.histogram(seq_num[ff], bins=binEdges)

        #----------------------------------------------
        # extract transcript name 
        #----------------------------------------------
        nExpected = len(taglist)
        expected_list = []
        exp_tags = []
        for n in range(0,nExpected):
            expected_list.append(taglist[n][0])
            exp_tags.append(taglist[n][1])

        # --------------------------------------------
        # if unexpected write NNNN
        # ---------------------------------------------
        name_tag = []
        wildCardsInExpectedList = []
        wildCardTags = []
        for j in range(0, len(expected_list)):
            findNs = re.finditer('N', expected_list[j])
            c = 0
            for i in findNs:
                c = c + 1
            if c > 0 and c < num_hybs:
                wildCardsInExpectedList.append(expected_list[j])
                wildCardTags.append(exp_tags[j])

        for kk in range(0, len(b)):
            letters = self.num2letters(b[kk], num_hybs)
            if letters in expected_list:
                posit = expected_list.index(letters)
                name_tag.append(exp_tags[posit])
            else:
                ##check against wildcard N in taglist
                found = False
                foundStr = ""
                for j in range(0, len(wildCardsInExpectedList)):
                    findNs = re.finditer('N', wildCardsInExpectedList[j])
                    ll = list(letters)
                    for i in findNs:
                        ll[i.span()[0]] = 'N'
                    modLetters = "".join(ll)
                    if modLetters == wildCardsInExpectedList[j]:
                        found = True
                        foundStr = wildCardTags[j]
                        break
                if found:
                    name_tag.append(foundStr)
                else:
                    name_tag.append('NNNN')

        name_tag_all = []
        for kkk in range(0, len(all)):
            letters = self.num2letters(all[kkk], num_hybs)
            if letters in expected_list:
                posit = expected_list.index(letters)
                name_tag_all.append(exp_tags[posit])
            else:
                ##check against wildcard N in taglist
                found = False
                foundStr = ""
                for j in range(0, len(wildCardsInExpectedList)):
                    findNs = re.finditer('N', wildCardsInExpectedList[j])
                    ll = list(letters)
                    for i in findNs:
                        ll[i.span()[0]] = 'N'
                    modLetters = "".join(ll)
                    if modLetters == wildCardsInExpectedList[j]:
                        found = True
                        foundStr = wildCardTags[j]
                        break
                if found:
                    name_tag_all.append(foundStr)
                else:
                    name_tag_all.append('NNNN')
                
        name_tag_allbt = []
        for kkkk in range(0, len(allbt)):
            letters = self.num2letters(allbt[kkkk], num_hybs)
            if letters in expected_list:
                posit = expected_list.index(letters)
                name_tag_allbt.append(exp_tags[posit])
            else:
                ##check against wildcard N in taglist
                found = False
                foundStr = ""
                for j in range(0, len(wildCardsInExpectedList)):
                    findNs = re.finditer('N', wildCardsInExpectedList[j])
                    ll = list(letters)
                    for i in findNs:
                        ll[i.span()[0]] = 'N'
                    modLetters = "".join(ll)
                    if modLetters == wildCardsInExpectedList[j]:
                        found = True
                        foundStr = wildCardTags[j]
                        break
                if found:
                    name_tag_allbt.append(foundStr)
                else:
                    name_tag_allbt.append('NNNN')

        # ________________________________________________________________________________________
        # Write text file
        # -------------
           
        # open a file for writing
        fid = open(Output_path + output_prefix + '_codes_n_counts.csv', 'w')
        # print a title, followed by a blank line
        fid.write('Code Count\n\n');
        # print values in column order
        for row in range(0, len(hh)):
            fid.write('%s %d %s\n' % (seq_letters[row],ab[row], name_tag[row]))
        fid.close()
                
        # ________________________________________________________________________________________
        # Write text file with x and y positions of all found transcripts 
        # -------------
   
        # open a file for writing
        fid = open(Output_path + output_prefix + '_xy_position_all_found_afterQT.csv', 'w')
        # print a title, followed by a blank line
        fid.write('letters,name,X_pos,Y_pos,parent_cell,seq_quality\n\n')
        # print values in column order
        for row in range(0, len(all)):
            fid.write('%s,%s,%d,%d,%d,%.7f\n' %(all_letters[row], name_tag_all[row], np.floor(x_all[row]), np.floor(y_all[row]), par_cell_all[row], seq_quality_min_all[row]))
        fid.close()
        
        # ________________________________________________________________________________________
        # Write text file with x and y positions of all found transcripts before QT
        # -------------
   
        # open a file for writing
        fid = open(Output_path + output_prefix + '_xy_position_all_found_beforeQT.csv', 'w')
        # print a title, followed by a blank line
        fid.write('letters,name,X_pos,Y_pos,parent_cell,seq_quality\n\n')
        # print values in column order
        for row in range(0, len(allbt)):
            fid.write('%s,%s,%d,%d,%d,%.7f\n' %(all_lettersbt[row], name_tag_allbt[row], np.floor(x_allbt[row]), np.floor(y_allbt[row]), par_cell_allbt[row], seq_quality_min_allbt[row]))
        fid.close()

    def num2letters(self, num, num_hybs):
        #print("num: %s" %(str(num)))
        if num_hybs > 0:
            letters = ['A', 'C', 'G', 'T']
            numDiv10 = int(np.floor(num/10))
            rest = int(num - numDiv10 * 10)
            #print("  rest: %s" %(str(rest)))
            lastLetter = letters[rest-1]
            res = self.num2letters(numDiv10, num_hybs-1) + lastLetter
        else:
            res = ''
        return res
