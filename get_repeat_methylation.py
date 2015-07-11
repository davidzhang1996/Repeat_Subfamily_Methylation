"""
README
Overview

Inputs: This program takes in three input files: A sam file that is outputted by BsMap, a repeat-masker file that 
contains all the repeat occurences in the human genome and the alignments to their subfamily consensus sequences, and an 
embl file that contains the consensus sequences for the repeat subfamilies. 

Outputs: This program outputs the name of the subfamily, the locations of each cytosine in the consensus sequence, 
and the methylation level for each cytosine. 

Procedure: The program first loads the consensus sequences and the sam file into memory. It then parses the
repeat-masker file such that it analyzes each repeat occurence one by one. It collects all the sam reads whose location 
in the reference intersects with the location of the repeat occurence, and then computes the ratio of thymines to 
thymines plus cytosines in order to estimate a methylation level for each cytosine in the consensus sequence. Bases 
in the sam reads that are not c/t but corresbond to a cytosine in the consensus sequence are treated as either SNPs/reading 
errors and are not factored into the calculation. 


"""
import math 
import time
start_time = time.time()

rmsk_type=""
check_white_space=True 

#Input: None
#Output: None 

#Procedure: 
# The main method drives the whole program. The get_subfam_consensus_seqs method first creates the methylation dictionary 
# that contains the locations of the cytosines in the repeat consensus sequences. The create_sam method parses the Sam 
# file for multimapped reads. The parse_rmsk method takes each repeat occurence from the rmsk file one by one
# obtains all the sam files that overlap with that repeat, and compares the sam reads to the repeat consensus sequence 
# in order to estimate a methylation level for the cytosines. It then prints the position of each cytosine 
# in the consensus sequence and the percent methylation level. 


def main():
	methylation_dict={}
	
	methylation_dict=get_subfam_consensus_seqs(methylation_dict)
	sam_dict=create_sam()
	methylation_dict=parse_rmsk(sam_dict, methylation_dict)
	convert_methylation_dict(methylation_dict)

#==================================================
#1
#Input: Emty methylation dictionary 
#Output: Methylation dictionary that has the name of the repeat subfamily as the key and a list as the value. 
#The entries of the list are of the format, [location_of_cytosine_along_consensus|number_converted_t's/total_number_c/t]

#Procedure: 
#This method parses the repeat consensus sequence file in order to obtain the consensus sequences. The repeat consensus 
#sequence file is in embl format, which means that each line is annotated with a 2 letter ID. The method identifies
#which lines contain the actual sequence by the ID and appends the locations of the cytosines in that sequence to the 
#methylation dictionary. 

def get_subfam_consensus_seqs(methyl_dict):
	methylation_dict=methyl_dict

	with open ("repeat_consensus_seqs_987.txt", "r") as repeats_file:
		subfam_name, cons_position = "", 1  
		sequence=False, 1
		for line in repeats_file: 
			line=line.rstrip("\n")
			line_list=line.split()

			if (line_list[0]=="ID"):
				sequence=False
				subfam_name=line_list[1]
				cons_position=1
				methylation_dict[subfam_name]=[[],[]]

			if (sequence==True):
				for bp_block in xrange(len(line_list)-1):
					for bp in line_list[bp_block]: 
						if (bp=="C" or bp=="c"):
							methylation_dict[subfam_name][0].append("%s|0/0" % (cons_position))
						if (bp=="G" or bp=="g"):
							methylation_dict[subfam_name][1].append("%s|0/0" % (cons_position))
						cons_position+=1

			if (line_list[0]=="//"):
				reverse=[]
				for cons_location in xrange(len(methylation_dict[subfam_name][1])-1, -1, -1): 
					cons_bp=int(methylation_dict[subfam_name][1][cons_location][0:methylation_dict[subfam_name][1][cons_location].index("|")])
					reverse.append("%s|0/0" %(cons_position-cons_bp+1))
				methylation_dict[subfam_name][1]= reverse
				reverse= []

			if (line_list[0]=="SQ"): 
				sequence=True

	return methylation_dict 

#==================================================
# 2
# Inputs: None 

# Outputs: A dictionary that contains the locations of all the cytosines and thymines relative to the start of the 
# chromosome. The keys of the dictionary are of the format: chromsome|start_position~stop_position.repeater. 
# 
# The repeater is added in case there are multiple sam reads that have the sam start site. In that case,  
# their keys would be the same and only one of those reads would be recorded. Prelimary analysis suggests this  
# occurs around once every 800 multimapped reads. 

# Procedure: 
# This method block contains two methods, the top level create_sam method and the auxiliary 
# reverse_complements method. The create_sam method scans through the sam file and filters for reads that have a flag
# that is greater than or equal to 256 and less than or equal to 511. Flags that are greater than 256 are multimapped. 
# Flags that are greater than 272 are sequences that are mapped to the reverse complement. They are processed with the 
# reverse_complement method that takes the reverse complement of the sam sequence so it can directly be compared to the
# repeat masker reference sequence. 

def create_sam():
	sam_dict={}
	with open ("adrenal_gland_100000000_filtered_forward.sam", "r") as sam_file:
		repeater=1
		for line in sam_file: 
			line=line.rstrip("\n")
			line_list=line.split()
			key, location, sam_seq, chrm, chrm_start, chrm_stop = "", "", "", "", "", ""
			value=[]

			if(line_list[0][0]!="@" and (int(line_list[1])==256 or int(line_list[1])>=272)):
				sam_seq, chrm, chrm_start= line_list[9], line_list[2][3:len(line_list[2])], line_list[3]

				chrm_stop= "%s" %(int(chrm_start)+len(sam_seq)-1)

				if (len(chrm)==1): 
					chrm="0"+chrm

				key=chrm+"|"+chrm_start+"~"+chrm_stop
				if (sam_dict.has_key(key+"."+"%s" %(repeater))):
					repeater+=1
				key=key+"."+"%s" %(repeater)
				sam_dict[key]=[]

				if (int(line_list[1]==256)):
					sam_dict[key].append("F")

				elif (int(line_list[1])==272): 
					sam_dict[key].append("C")

				for index in xrange(len(sam_seq)):
					location=""

					location= "%s|%s" %(index+int(chrm_start), sam_seq[index].upper())

					if (location!=""):
						sam_dict[key+"."+"%s" %(repeater)].append(location)

				if (len(sam_dict[key+"."+"%s" %(repeater)])==0):
					sam_dict[key+"."+"%s" %(repeater)].append('-1')

	return sam_dict

#==================================================
# 3
# Inputs:		1: The sam_dict that contains the chrosome and chromosome starting site for each read
#		    	2: The methylation dictionary that contains the name of all the repeat subfamilies 
#			   and the locations of all the cytosines in the consensus sequence of the subfamily
#
# Output: 	    The methylation dictionary that contains the updated methylation levels for all the repeat subfamilies
#				That are found in the repeat masker file 
#
# Procedure:	The parse_rmsk method is a fairly top-level method that loops through the repeat-masker file to analyze 
#				each repeat occurence one by one. It creates a variable rmsk_dict, which the get_repeats_rmsk()  method 
#				modifies to contain information regarding the location of repeat, the locations of the cytosines relative 
#				to the start of the chromosome and the locations the cytosines relative to the start of the 
#				consensus sequence. Once the the complete repeat subfamily is retrieved fromt the file, it is passed on 
#				to other methods and the rmsk dictionary is cleared. 				


#				The get_repeats_sam() method gets all the sam reads that intersect with the repeat 
#				occurence and appends them to the end of the rmsk_dict value. 

#				The calculate_methylation_levels() method compares the locations of the cytosines in the repeat consensus
#				sequence with the intersected sam reads, and updates the methylation levels for the cytosines based on 
#				the ratio of converted cytosines to non-converted cytosines.  
#				
#				The update_methylation_dictionary() method then uses the methylation data in the rmsk_dict to 
#				update the methylation_dict dictionary. By the end, the methylation_dict dictionary should reflect 
#				the methylation pattern for that particular repeat occurence added to all the previous methylation patterns
#				for that repeat subfamily. 

rmsk_dict={}

def parse_rmsk(sam_dict, methyl_dict):
	rmsk_type, line_count, complement, complete_subfamily= "Header", 1, False, False
	rmsk_dict={}

	with open ("hg19_154.fa.align", "r") as rmsk_file: 
		for line in rmsk_file: 
			rsmk_dict, rmsk_type, line_count, complement, complete_subfamily = get_repeats_rmsk(rmsk_dict, line, rmsk_type, line_count, complement, complete_subfamily) #3a

			if (complete_subfamily):
				rmsk_dict=get_repeats_sam(rmsk_dict, sam_dict) #3b
				rmsk_dict=calculate_methylation_levels(rmsk_dict) #3c
				methylation_dict=update_methylation_dictionary(rmsk_dict, methylation_dict) #3d

		return methylation_dict

#==================================================
# 4
# Inputs:		1: The methylation_dict that contains the complete methylation patters for all the repeat subfamilies 
#				   found in the repeat-masker file 

# Outputs:		Prints to the screen the name of the subfamily, the positions of all the cytosines, and the methylation 
#				levels for those cytosines

# Procedure: 	This is the main output method. It just loops through the methylation dictionary and ouputs the info 
#				in the appropriate format. 	

def convert_methylation_dict(methyl_dict):
	methylation_dict=methyl_dict
	for key in methylation_dict: 
		methyl=methylation_dict[key]
		
		for methyl_bp in xrange(len(methyl)):
			methyl_loc_separator=methyl[methyl_bp].index("|")
			methyl_tot_separator=methyl[methyl_bp].index("/")

			position=methyl[methyl_bp][0:methyl_loc_separator]
			methyl_number=float(methyl[methyl_bp][methyl_loc_separator+1: methyl_tot_separator])
			total=float(methyl[methyl_bp][methyl_tot_separator+1:])

			if (total!=0.0):
				print key
				print position + "\t" + "%s" %(methyl_number/total)

	print("--- %s seconds ---" % (time.time() - start_time))

#=========================
# 2a
# Inputs: Sam sequence that maps to reverse complement 
# Outputs: Sam sequence that is reverse complemented 

# Procedure: Start from the end of the original sequence, go to start of sequence, and return the complement 
# of that sequence 

def reverse_complement(sequence): 
	sam_seq=sequence 
	reversed_sam_seq=""
	for bp in xrange(len(sam_seq)-1,-1,-1):
		if (sam_seq[bp]=="C" or sam_seq[bp]=="c"): 
			reversed_sam_seq+="G"
		elif (sam_seq[bp]=="T" or sam_seq[bp]=="t"):
			reversed_sam_seq+="A"
		elif (sam_seq[bp]=="A" or sam_seq[bp]=="a"):
			reversed_sam_seq+="T"
		elif (sam_seq[bp]=="G" or sam_seq[bp]=="g"):
			reversed_sam_seq+="C"
	return reversed_sam_seq

#=========================
# 3a
# Inputs:		1: The rmsk_dictionary that contains the incomplete reference and consensus sequence information. 
#				2: A line from the repeat-masker file. 

# Outputs:		The rmsk dictionary with updated information from the line. This typically includes information regarding
#				the locations of the cytosines. 

#Procedure: 	This method basically takes in the line from the repeat masker file and classifies it as either a 
#				Header, Sequence, or Footer line based on the first couple of characters in the line. If it is a 
#				Header line, then it creates an entry in the rmsk_dict based off of the subfamily name with information
#				regarding the repeats starting and stopping site with regards to both the chromosome and the consensus 
#				sequence. If it is a Sequence line, then it updates the rmsk dictionary on the locations of the cytosines
#				with respect to the chromosome and consensus sequence. If it is a Footer line, then it tells the 
#				parse_rmsk() method that called it that a complete repeat occurence has occured. 
				
#				This method calls upon the get_type() method, which basically is the method that actually classifies 
#				the line based on the first couple of characters in the file. 


def get_repeats_rmsk(rmsk_dictionary, line_string, rmsk_type_string, line_count_int, complement_bool, complete_subfamily_bool, line_number_int):
	rmsk_dict, line, rmsk_type, line_count, complement, complete_subfamily, line_number = rmsk_dictionary, line_string, rmsk_type_string, line_count_int, complement_bool, complete_subfamily_bool, line_number_int

	line=line.rstrip("\n")
	line_list=line.split()

	if (rmsk_type=="Sequence"): 
		cons_start, cons_stop, chrm_seq, cons_seq_reverse, cons_seq_forward = "","","","",""
		subfamily=rmsk_dict.keys()[0]

		if (line_count%4==2):
			if (line_list[0]=="Matrix"):
				complete_subfamily=True
				rmsk_type="Footer"
				
			else: 
				chrm_seq=rmsk_dict[subfamily][1]+line_list[2]
				rmsk_dict[subfamily][1]=chrm_seq

		elif(line_count%4==0 and rmsk_type!="Footer"):
			if (len(line_list[0])==1):
				line_list.remove(line_list[0])

			if(line_count==4):
				cons_start=line_list[1]

			if (complement): 
				cons_seq_reverse= rmsk_dict[subfamily][3]+line_list[2]
				rmsk_dict[subfamily][3]=cons_seq_reverse
			else: 
				cons_seq_forward=  rmsk_dict[subfamily][2]+line_list[2]
				rmsk_dict[subfamily][2]=cons_seq_forward
			
			cons_stop=line_list[3]	

		if (len(cons_start)!=0):
			rmsk_dict[subfamily][0][3] = cons_start

		if (len(cons_stop)!=0):
			rmsk_dict[subfamily][0][4] = cons_stop
		line_count+=1

		return rmsk_dict, rmsk_type, line_count, complement, complete_subfamily

	if (rmsk_type=="Footer" and len(line)!=0):
		if (line_list[0]=="Transitions"):
			subfamily=rmsk_dict.keys()[0]
			ref_seq, cons_seq, cons_start, cons_remaining, ref_start, complement = rmsk_dict[subfamily][1], rmsk_dict[subfamily][2], int(rmsk_dict[subfamily][0][3]), int (rmsk_dict[subfamily][0][5]), int(rmsk_dict[subfamily][0][1]), rmsk_dict[subfamily][0][6]
			cons_positions_list, ref_positions_list=[], []

			rmsk_dict[subfamily][1]=[]
			rmsk_dict[subfamily][2]=[]  

			cons_position=cons_start
			ref_position=ref_start

			for index in xrange(len(ref_seq)):
				if (ref_seq[index]==cons_seq[index]):
					if (cons_seq[index]=="C" or cons_seq[index]=="c"):
						cons_positions_list.append("%s|C" %(cons_positions))
						ref_positions_list.append("%s|C" %(ref_positions))

					elif (cons_seq[index]=="G" or cons_seq[index]=="g"):
						cons_positions_list.append("%s|G" % (cons_positions))
						ref_positions_list.append("%s|G" %(ref_positions))

				if (complement): 
					cons_position-=1 
				else: 
					cons_position+=1

			rmsk_dict[subfamily][1], rmsk_dict[subfamily][2], = ref_positions_list, cons_positions_list
			return rmsk_dict, rmsk_type, line_count, complement, complete_subfamily

		elif (line_list[0]=="Gap_init"):
			rmsk_type, line_count, complement, complete_subfamily="Header", 1, False, False 	
			rmsk_dict.clear()
			return rmsk_dict, rmsk_type, line_count, complement, complete_subfamily
		else: 
			return rmsk_dict, rmsk_type, line_count, complement, complete_subfamily		


	if(rmsk_type=="Header" and len(line)!=0):
		if (line_list[8]!="C"):
			subfam_name = line_list[8]
			cons_remaining= line_list[11][1:len(line_list[11])-1]
			
		else: 
			complement=True 
			subfam_name = line_list[9]
			cons_remaining= line_list[10][1:len(line_list[10])-1]

		subfamily, chrm, chrm_start, chrm_stop, cons_start, cons_stop, chrm_seq, cons_seq = subfam_name[:subfam_name.index("#")], line_list[4][3:len(line_list[4])], line_list[5], line_list[6], "", "", "", ""

		if (len(chrm)==1):
			chrm="0"+chrm

		info=[chrm, chrm_start, chrm_stop, cons_start, cons_stop, cons_remaining, complement]

		rmsk_dict[subfamily] = [info, chrm_seq, cons_seq]
		rmsk_type="Sequence"
		return rmsk_dict, rmsk_type, line_count, complement, complete_subfamily

	elif(rmsk_type=="Header" and len(line)==0):
		return rmsk_dict, rmsk_type, line_count, complement, complete_subfamily

#=========================
# 3b
# Inputs:		1: The rmsk_dictionary that contains the complete information of the locations of the cytosines 
#				2: The sam_dictionary that contains information regarding the chromosome and chromosome start site
#				   of the read, and the locations of the cytosines and thymines with regard to the start of the chromsome
# 
# Output:		The rmsk_dictionary with the methylation levels for each cytosine appended to the value of the dictionary 
#
# Procedure:	This method is a medium-level method that dictates how the intersecting sam reads are actually retrieved
#				from the sam_dictionary. This method creates a list that contains the keys of all the sam reads in sam_dict,  	
#	 			the keys of the reads contains information on where the sam read is located on the reference genome. 
#				Because lists are inherently ordered, you can sort the entries of the list alphanumerically, and then do 
#				a binary search through the list by comparing the locaiton of your repeat occurence to the locations in 
#				the sam read. Then you can isolate the sam reads that intersect with your repeat masker file, and then 
#				call the cytosine locations from the sam_dict using the entries in sam_list. 
#
#				The sort_sam_list() method sorts the sam list alphanumerically.
#
#				The valid_sam_repeats() method searches through the sorted sam list for the interesecting sam reads 
#				and appends the cytosine and thymine locations to the end of the rmsk_dict value.  

def get_repeats_sam(rmsk_dictionary, sam_dictionary):
	rmsk_dict=rmsk_dictionary
	sam_dict=sam_dictionary 

	sam_list=[]

	for item in sam_dict: 
		sam_list.append(item)

	sam_list=sort_sam_list(sam_list) #3b i
	rmsk_dict=valid_sam_repeats(rmsk_dict, sam_dict, sam_list) #3b ii
	return rmsk_dict

#========================
# 3c
# Inputs:		1: The complete rmsk_dictionary containing the rmsk information as well as the interseting sam sequences. 

#
# Output: 		This method returns the rmsk_dict with the methylation levels for each cytosine in the rmsk consensus 
#				analyzed. 
#
# Procedure: 	This is a low level method that basically loops through the sam reads in the rmsk dict and checks to see
#				if the reads contain thymines or cytosines in the chromosome positions where there are cytosines in the 
#				reference genome. 


def calculate_methylation_levels(rmsk_dictionary):
	rmsk_dict=rmsk_dictionary
	key=rmsk_dict.keys()[0]
	total, methyl, location=0, 0, 0
	methylation_sequence_forward, methylation_sequence_reverse=[], []
	ref_pos_list=rmsk_dict[key][1]
	cons_remaining, ref_complement=rmsk_dict[0][6], rmsk_dict[0][5]
	
	for ref_pos in xrange(len(ref_pos_list)):
		for read in xrange(3, len(rmsk_dict[key])): 
			sam_pos_list=rmsk_dict[key][read]
			sam_strang=sam_pos_list[0]
			for sam_pos in xrange(2, len(sam_pos_list)): 
				separatorIndexSam, separatorIndexRef = sam_pos_list[sam_pos].index("|"), res_pos_list[re_pos].index("|")
				sam_chrm_pos, ref_chrm_pos, sam_chrm_bp, ref_chrm_bp = sam_pos_list[sam_pos][0:separatorIndexSam], ref_pos_list[ref_pos][0:separatorIndexRef], sam_pos_list[sam_pos][separatorIndexSam+1:], ref_pos_list[ref_pos][separatorIndexRef]
				
				if (complement and sam_strand="C"):
					if ((sam_chrm_bp=="A" or sam_chrm_bp=="G") and (ref_chrm_bp=="A") and ref_chrm_pos==sam_chrm_pos):
						total+=1 
						if (sam_chrm_bp=="G"): 
							methyl+=1
					0

				elif (complement and sam_strand="F"):
				elif (not complement and sam_strand="F"):
				elif (not complement and sam_strand="C"):
				if (ref_pos_list[ref_pos]==int(sam_pos_list[sam_pos][0:separatorIndex]) and (sam_pos_list[sam_pos][separatorIndex+1:]=="C" or sam_pos_list[sam_pos][separatorIndex+1:]=="T")): 
					total+=1 
					if (sam_pos_list[sam_pos][separatorIndex+1:]=="T"):
						methyl+=1

		if (total!=0):		
			methylation_sequence.append("%s|%s/%s" %(rmsk_dict[key][2][ref_pos], methyl, total))	

		total=0
		methyl=0
			
	rmsk_dict[key]=methylation_sequence
	return rmsk_dict

#========================
# 3d 
# Inputs:		1: The rmsk dict that contiains the methylation levels for all the cytosines along the consensus sequence
#				2: The methyl dict that contains the locations of all the cytosines along the subfamily consensus sequence
#				   and the methylation levels of those cytosines 
#
# Outputs: 		3: The methyl dict that contains updated methylation levels for the subfamily contained in the rmsk dict 
#
# Procedure: 	This method is a fairly low level method that just loops through the rmsk dict and for each cytosine 
#				location in the rmsk dict, it updates the information for that cytosine in the methylation dict 

def update_methylation_dictionary(rmsk, methyl_dict):
	rmsk_dict, key, methylation_dict= rmsk, rmsk.keys()[0], methyl_dict
	methyl, total= 0, 0

	for rmsk_bp in xrange(len(rmsk_dict[key])): 
		rmsk=rmsk_dict[key]
		rmsk_loc_separator, rmsk_tot_separator=rmsk[rmsk_bp].index("|"), rmsk[rmsk_bp].index("/")
		rmsk_loc, rmsk_methyl, rmsk_total= rmsk[rmsk_bp][0:rmsk_loc_separator], rmsk[rmsk_bp][rmsk_loc_separator+1:rmsk_tot_separator], rmsk[rmsk_bp][rmsk_tot_separator+1:]

		for methyl_bp in xrange(len(methyl_dict[key])):
			methyl=methylation_dict[key]
			methyl_loc_separator, methyl_tot_separator=methyl[methyl_bp].index("|"), methyl[methyl_bp].index("/")
			methyl_loc, methyl_methyl, methyl_total = methyl[methyl_bp][0:methyl_loc_separator], methyl[methyl_bp][methyl_loc_separator+1:methyl_tot_separator], methyl[methyl_bp][methyl_tot_separator+1:]
			
			if (rmsk_loc==methyl_loc):
				methyl= int(rmsk_methyl)+int (methyl_methyl)
				total=int(rmsk_total)+ int(methyl_total)
				methylation_dict[key][methyl_bp]=methyl_loc+ "|"+ "%s/%s" %(methyl,total)	

	return methylation_dict

#============
# 3a i
# Inputs:		1: An input list you want to reverse the elements of
#
# Outputs: 		The input list with its elements reversed
#
# Procedure:	This is a low level method that takes the last element of the input list and appends it as the first 
#				element of the output list, and so on and so forth

def reverse (input_list):
	output_list=[]
	for index in xrange(len(input_list)-1, -1, -1): 
		output_list.append(input_list[index])
	return output_list



#============
#3a ii 
# Inputs:		1: An input list that contains the positions of the cytosines in a sequence

# Output: 		A list that contains the list of the cytosines in the reverse complement of a sequence 

# Procedure: 	This is a low level method that takes a position of a cytosine in a sequence, and subtracts the position
#				from the length of the sequence plus 1. 

def reverse_complement(input_list, end_positions_int):
	output_list=[]
	for index in xrange (len(input_list)-1, -1, -1): 
		output_list.append(input_list[index]-end_positions_int+1)
	output_list=reverse(output_list)
	return output_list

#============
# 3b i
# Inputs:		The sam list that contains the keys for the sam dictionary.  
#
# Outputs: 		The sam list ordered in alphanumerical order 
#
# Procedure: 	This is a low level method that uses a merge sort to order the elements of the sam list 

def sort_sam_list(sam):
	chrm, chrmPos, lhSepInd, rhSepInd = "", "", "", ""
	aList=sam
	if (len(aList)>1):  
		midList=len(aList)/2
		rightHalf=sam[midList:]
		leftHalf=sam[:midList]

		rhLength=len(rightHalf) 
		lhLength=len(leftHalf)

		rightHalf=sort_sam_list(rightHalf)
		leftHalf=sort_sam_list(leftHalf)

		i=0
		j=0
		k=0

		while (i<lhLength and j<rhLength): 
			lhSepInd=leftHalf[i].index("|")
			rhSepInd=rightHalf[j].index("|")

			if (leftHalf[i][:lhSepInd]<rightHalf[j][:rhSepInd]):
				aList[k]=leftHalf[i]
				i+=1
			elif (leftHalf[i][:lhSepInd]>rightHalf[j][:rhSepInd]): 
				aList[k]=rightHalf[j]
				j+=1

			elif (leftHalf[i][:lhSepInd]==rightHalf[j][:rhSepInd]):
				if (leftHalf[i][lhSepInd+1:]<rightHalf[j][rhSepInd+1:]):
					aList[k]=leftHalf[i]
					i+=1
				else: 
					aList[k]=rightHalf[j]
					j+=1
			k+=1 

		while (i<lhLength):

			aList[k]=leftHalf[i]
			k+=1
			i+=1

		while (j<rhLength): 
			aList[k]=rightHalf[j]
			k+=1 
			j+=1 

	return aList

#============
# 3b ii 
# Inputs:		1: The rmsk dictionary that contains the positions of the cytosines of the repeat along the chromosome 
#				2: The sam dictionary that contains all the sam reads 
#				3: The sam list that is ordered with respect to the position of the reads along the chromosome
# 
# Outputs:		The rmsk_dict with the locations of the cytosines in all the overlapping sam reads appended to the rmsk
#				dictionary. 
#
# Procedure: 	This is a mid level method that loops through all the overlapping sam reads in the sam list, and uses the
#				names of the sam reads to call upon the information stored in the sam dictionary. It determines where to 
#				start and stop looping through the sam list by calling upon the search_sam method to find the sam reads 
#				start and stop overlapping with the repeat. 

def valid_sam_repeats(rmsk_dictionary, sam_dictionary, sam): 
	rmsk_dict, sam_dict, sam_list = rmsk_dictionary, sam_dictionary, sam

	rmsk_key=rmsk_dict.keys()[0]
	chrm, chrm_start, chrm_stop = rmsk_dict[rmsk_key][0][0], rmsk_dict[rmsk_key][0][1], rmsk_dict[rmsk_key][0][2]
	
	start_position=search_sam(sam_list, chrm, chrm_start,0, "Start") # 3b ii a 
	stop_position=search_sam(sam_list, chrm,chrm_stop,0, "Stop")

	for index in xrange(start_position, stop_position+1): 
		chrm=sam_list[index][sam_list[index].index("|")+1:] 
		sequence=sam_dict[sam_list[index]]
		rmsk_dict[rmsk_key].append(sequence)
	return rmsk_dict

#======
# 3b ii a
# Inputs:		1: The ordered sam list containing the locations of the sam reads in alphanumerical order
#				2: The chromsome and chromosome position of the rmsk read 
#				3: The boundary position of where the sam reads start intersecting or stop intersecting the rmsk 
#				   repeat, based on the chromosome_position_type 
#				4: The chromosome position type, which details whether the function should look for where the sam reads 
#				   start or stop intersecting the rmsk repeat 
# 
# Outputs:		The location along the sam list of where the rmsk repeat starts or stops interesecting the sam reads
#
# Procedure: 	This method is a low level method that utilizes that fact that the sam_list is ordered with regards 
#				to chromosomal position. It uses a binary search to find where the rmsk repeat overlaps with the sam reads
	


def search_sam(searchList, chromosome, chromosome_position, boundary_position, chromosome_position_type):
	sam_list, chrm, chrm_pos, bound_pos, chrm_pos_type= searchList, chromosome, int(chromosome_position), boundary_position, chromosome_position_type
	midList=len(sam_list)/2

	separatorIndex=sam_list[midList].index("|")
	chromosomeSeparatorIndex=sam_list[midList].index("~")
	repeaterSeparatorIndex=sam_list[midList].index(".")
	sam_chrm=sam_list[midList][0:separatorIndex]
	sam_chrm_start_pos=int(sam_list[midList][separatorIndex+1:chromosomeSeparatorIndex])
	sam_chrm_stop_pos=int(sam_list[midList][chromosomeSeparatorIndex+1:repeaterSeparatorIndex])

	rightHalf=sam_list[midList:]
	leftHalf=sam_list[:midList]

	if (len(sam_list)>1): 
		if (chrm<sam_chrm):
			bound_pos+=search_sam(leftHalf, chrm, chrm_pos, bound_pos, chrm_pos_type)

		elif (chrm>sam_chrm):
			bound_pos+=midList+search_sam(rightHalf, chrm, chrm_pos, bound_pos, chrm_pos_type)

		elif (chrm==sam_chrm): 
			if (chrm_pos_type=="Start"):
				if (chrm_pos<sam_chrm_stop_pos):
					bound_pos+=search_sam(leftHalf, chrm, chrm_pos, bound_pos, chrm_pos_type)

				elif (chrm_pos>sam_chrm_stop_pos):
					bound_pos+=midList+search_sam(rightHalf, chrm, chrm_pos, bound_pos, chrm_pos_type)
			elif (chrm_pos_type=="Stop"):
				if (chrm_pos>sam_chrm_start_pos):
					bound_pos+=midList+search_sam(leftHalf, chrm, chrm_pos, bound_pos, chrm_pos_type)

				elif (chrm_pos<sam_chrm_start_pos):
					bound_pos+=search_sam(rightHalf, chrm, chrm_pos, bound_pos, chrm_pos_type)

	else:
		if (chrm<sam_chrm):
			bound_pos+=0

		elif (chrm>sam_chrm):
			bound_pos+=1

		elif (chrm==sam_chrm): 
			if (chrm_pos_type=="Start"): 
				if (chrm_pos<sam_chrm_stop_pos):
					bound_pos+=0

				elif (chrm_pos>sam_chrm_stop_pos):
					bound_pos+=1
			elif (chrm_pos_type=="Stop"):
				if (chrm_pos>sam_chrm_start_pos):
					bound_pos+=1

				elif (chrm_pos<sam_chrm_stop_pos):
					bound_pos+=0

	return bound_pos

main()

