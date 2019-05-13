f1=open("/Users/ahmedelewa/Documents/Postdoc/2-Blastema/2-padlock/Ch038_QT_0.4_0.001_details.csv","r")
from collections import defaultdict


## variables   ##
radius =100                                             # maximum radius analyzed from an target
cell_call_limit =2                                      # number of self targets to call a cell type
abort_limit =1                                          # number of other targets to rule-out cell type
##    lists    ##
housek =['40S','60S','atp5','eif2']			                # housekeeping genes
group2 =['DMBT1','MARCS','AMPN','Copia','Gypsy','LINE.2']               # markers for group 2 - expressing TEs (Nuria)
group6 =['NDUA7','PELO','unkG6','ANGL2']                                # markers for group 6 - expressing ECM genes (Ahmed)
group7 =['SRGN','CYTF','LYSC','unkG7','CCR4', 'RORC', 'CD11b','MPEG1']  # markers for group 7 - probably immunecells  (Connie)
litgen =['KAZD','CIRBP','Sox9']                                         # genes targeted based on literature
unwant =['NNNN']					                # for filtering out (like NNNN)
target_group = group2                                   # this is the group of interest
## dictionaries ##
dict_map	= defaultdict(str) 			# dictionary where target (serial number) is key and x,y position is value
dict_id		= defaultdict(str) 			# dictionary where target (serial number) is key and ID (gene name) is value
dict_coor	= defaultdict(str)			# dictionary where coordinate (x,y) is key and ID (gene name) is value. 
dict_coor_sanity= defaultdict(list)			# make dict_coor value a list in case more than one gene in same location
dict_coor_mirror= defaultdict(str)

#########################
## Step 1a: import data ##

f1.readline()						# skip header
c = 1							# for serial number
for a in f1:
	a		= a.strip()
	serial		= str(c).zfill(8)
	c+=1
	gene		= a.split(",")[1]
	posx		= str(int(float(a.split(",")[2])))	# convert float to integer and then to string
	posy		= str(int(float(a.split(",")[3])))	# convert float to integer and then to string

## Step 1b: filter data ##
	
	if not gene in unwant:					# remove entries in the unwant list (like NNNN)
		dict_map[serial]=posx+","+posy
		dict_id[serial]=gene
		dict_coor_sanity[posx+","+posy].append(gene)		
		dict_coor[posx+","+posy]=gene
		dict_coor_mirror[posx+","+posy]=gene

###########################
## Step 2a: sanity check ##					# make sure no position has more than 1 gene

conflicts = 0
for coor, gene in dict_coor_sanity.iteritems():
	if len(gene) > 1: 					
		conflicts+=1
print "Number of overlaping genes =" ,conflicts
		

## Step 2b: start walking ##

if conflicts == 0:
	#dict_coor_mirror = dict_coor				# duplicate dictionary to iterate of one and extract from the second
	for coor, gene in dict_coor_mirror.iteritems():
                neighbor=[]                                     # general list for all neighbors
                count_target_group=0                            # this is for Step 3.
                count_other_group=0                             # this is for Step 3.
                count_housek=0                                  # this is for Step 3.
                status = 'in progress'                          # this switch turns to 'complete' during evaluation (Step 3b)
                report = ''                                     # print out of the cell type sucess report
		posx	= coor.split(",")[0]
		posy	= coor.split(",")[1]
		if gene in target_group:                        # only process cell-types of interest (should make the process faster)
                        l=1                                     # start with level 1 (one degree away from target)
                        while l<radius:                         # i.e. while level is less than maximum defined at begining
                                exec("neighbor" +str(l)+ "=[]") # create list for level
                                exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(0)+ ")+" + '","+str(int(posy)+' +str(l)+ ")])") 
                                exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(0)+ ")+" + '","+str(int(posy)-' +str(l)+ ")])") 
                                exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(l)+ ")+" + '","+str(int(posy)+' +str(0)+ ")])")
                                exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)-" +str(l)+ ")+" + '","+str(int(posy)+' +str(0)+ ")])")
                                s=1                             # continue with step 1 (step 0 was the previous command)
                                while s< l  :                   # i.e. while steps are within level range minus 1
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(s)+ ")+" + '","+str(int(posy)+' +str(l)+ ")])") #+s +l (s loop 0 to l)
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(s)+ ")+" + '","+str(int(posy)-' +str(l)+ ")])") #+s -l (s loop 0 to l)
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)-" +str(s)+ ")+" + '","+str(int(posy)+' +str(l)+ ")])") #-s +l (s loop 0 to l)
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)-" +str(s)+ ")+" + '","+str(int(posy)-' +str(l)+ ")])") #-s -l (s loop 0 to l)
                                        s+=1

                                s=1                             # reset to step 1 to search the perimeters
                                while s< l+1  :                 # i.e. while steps are within level range
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(l)+ ")+" + '","+str(int(posy)+' +str(s)+ ")])") #+s +l (s loop 0 to l)
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)+" +str(l)+ ")+" + '","+str(int(posy)-' +str(s)+ ")])") #+s -l (s loop 0 to l)
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)-" +str(l)+ ")+" + '","+str(int(posy)+' +str(s)+ ")])") #-s +l (s loop 0 to l)
                                        exec("neighbor" +str(l)+ ".append(dict_coor[str(int(posx)-" +str(l)+ ")+" + '","+str(int(posy)-' +str(s)+ ")])") #-s -l (s loop 0 to l)
                                        s+=1
                                exec("neighbor+=neighbor" +str(l)) # append to general list

## Step 3a: evaluate after each level ##
                                
                                exec("evaluation=neighbor" +str(l))     # pass the list of the current level into another list called 'evaluation'
                                for n in evaluation:
                                        if n != '':
                                                if n in housek:
                                                        count_housek+=1
                                                if n in target_group:   # add score to target group
                                                        count_target_group+=1
                                                if not n in target_group and not n in housek and not n in litgen:      # if it's not in target group
                                                        count_other_group+=1                                            # or housek or litgene then it's another group
                                
## Step 3b: the breaks for aborting after failed cell-calling ##
                                
                                if count_other_group < abort_limit:             # continue to evaluate until the number of other targets exceeds limit
                                        if count_target_group >= cell_call_limit and status == 'in progress': # if number of self targets exceeds limit call cell and compelte evaluation
                                                report = gene, posx, posy, set(neighbor),"level ",l, count_target_group, count_other_group, count_housek
                                                status = 'complete'
                                l+=1
                        if not report =='':
                                print report,count_housek
		
