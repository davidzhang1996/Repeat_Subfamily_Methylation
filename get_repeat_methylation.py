"""
V1: 6/16/15 

"""
import math 

methylation_dict={}
rmsk_type=""
check_white_space=True 

# The main method drives the whole program. It first creates dictionaries that contains locations 
# of the cytosines in the repeat consensus sequences (get_subfam_consensus_seqs), and parses the Sam 
# file for multimapped reads (create_sam). The program then takes each repeat occurence from the rmsk file 
# one by one, obtains all the sam files that overlap with that repeat, and updates the methylation dictionary
# based on comparing the sam reads to the reference (parse_rmsk). It it then prints the position of each cytosine 
# in the consensus sequence and the percent methylation. 

def main():
	methylation_dict={}

	methylation_dict=get_subfam_consensus_seqs(methylation_dict)
	sam_dict=create_sam()
	methylation_dict=parse_rmsk(sam_dict, methylation_dict)
	convert_methylation_dict(methylation_dict)

def get_subfam_consensus_seqs(methyl_dict):
	methylation_dict=methyl_dict

	with open ("repeat_consensus_seqs_1000.txt", "r") as repeats_file:
		subfam_name, cons_position = "", 1  
		sequence=False
		for line in repeats_file: 
			line=line.rstrip("\n")
			line_list=line.split()

			if (line_list[0]=="ID"):
				sequence=False
				subfam_name=line_list[1]
				cons_position=1
				methylation_dict[subfam_name]=[]

			if (sequence==True):
				for bp_block in xrange(len(line_list)-1):
					for bp in line_list[bp_block]: 
						if (bp=="C" or bp=="c"):
							methylation_dict[subfam_name].append("%s|0/0" % (cons_position))
						cons_position+=1

			if (line_list[0]=="SQ"): 
				sequence=True
	return methylation_dict


def create_sam():
	sam_dict={}
	with open ("adipose_head_300.sam", "r") as sam_file:
		repeater=1
		for line in sam_file: 
			line=line.rstrip("\n")
			line_list=line.split()
			key, location, sam_seq, chrm, chrm_start, chrm_stop = "", "", "", "", "", ""
			value=[]

			if(line_list[0][0]!="@" and int(line_list[1])<=511 and int(line_list[1])>=256):
				sam_seq, chrm, chrm_start= line_list[9], line_list[2][3:len(line_list[2])], line_list[3]

				if (int(line_list[1])<=511 and int(line_list[1])>=272): 
					sam_seq=reverse_complement(line_list[9])

				chrm_stop= "%s" %(int(chrm_start)+len(sam_seq))

				if (len(chrm)==1): 
					chrm="0"+chrm

				key=chrm+"|"+chrm_start+"~"+chrm_stop
				if (sam_dict.has_key(key)):
					repeater+=1
				sam_dict[key+"."+"%s" %(repeater)]=[]

				for index in xrange(len(sam_seq)):
					location=""

					if (sam_seq[index]=="C" or sam_seq[index]=="c" or sam_seq[index]=="T" or sam_seq[index]=="t"):
						location= "%s|%s" %(index+int(chrm_start), sam_seq[index].upper())

					if (location!=""):
						sam_dict[key+"."+"%s" %(repeater)].append(location)

				if (len(sam_dict[key+"."+"%s" %(repeater)])==0):
					sam_dict[key+"."+"%s" %(repeater)].append('-1')
	#print sam_dict
	return sam_dict

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


rmsk_dict={}

def parse_rmsk(sam_dict, methyl_dict):
	sam_dict, methylation_dict=sam_dict, methyl_dict
	#print methylation_dict
	global check_white_space, complete_subfamily, line_count
	check_white_space, complete_subfamily, line_count=True, False, 1
	rmsk_dict={}

	with open ("hg19_328.fa.align", "r") as rmsk_file: 
		for line in rmsk_file: 
			rsmk_dict=get_repeats_rmsk(rmsk_dict, line)

			if (complete_subfamily):
				rmsk_dict=get_repeats_sam(rmsk_dict, sam_dict) 
				rmsk_dict=calculate_methylation_levels(rmsk_dict)
				#print methylation_dict
				methylation_dict=update_methylation_dictionary(rmsk_dict, methylation_dict)

				#for item in rmsk_dict: 
				#	print item
				#	print rmsk_dict[item]
				rmsk_dict.clear()
				line_count=1
				complete_subfamily=False
	#print methylation_dict
	return methylation_dict

def get_repeats_rmsk(dictionary, string):
	global check_white_space, complete_subfamily, line_count
	rmsk_dict, line = dictionary, string

	line=line.rstrip("\n")
	line_list=line.split()

	rmsk_type, subfamily, chrm, chrm_start, chrm_start, cons_start, cons_stop, chrm_seq, cons_seq = "", "", "", "", "", "", "", "", ""
	info_list, cons_positions_list, ref_positions_list=[],[],[]

	if (check_white_space==True):
		rmsk_type=get_type(line_list, line)

	if(check_white_space==False):
		subfamily=rmsk_dict.keys()[0]

		if (line_count%4==2):
			if (line_list[0]=="Matrix"):
				check_white_space=True

			else:
				chrm_seq=rmsk_dict[subfamily][1]+line_list[2]
				#print chrm_seq
				rmsk_dict[subfamily][1]=chrm_seq

		elif(line_count%4==0):
			if (len(line_list[0])==1):
				line_list.remove(line_list[0])

			if(line_count==4):
				cons_start=line_list[1]

			cons_stop, cons_seq=line_list[3], rmsk_dict[subfamily][2]+line_list[2]
			rmsk_dict[subfamily][2]= cons_seq	

		if (len(cons_start)!=0):
			rmsk_dict[subfamily][0][3] = cons_start

		if (len(cons_stop)!=0):
			rmsk_dict[subfamily][0][4] = cons_stop
		line_count+=1
		return rmsk_dict

	if(rmsk_type=="Footer" and line_list[0]=="Gap_init"):
		subfamily=rmsk_dict.keys()[0]
		cons_seq, ref_seq, cons_start, ref_start = rmsk_dict[subfamily][1], rmsk_dict[subfamily][2], int(rmsk_dict[subfamily][0][3]), int(rmsk_dict[subfamily][0][1])
		
		rmsk_dict[subfamily][1]=[]
		rmsk_dict[subfamily][2]=[]  

		cons_position=cons_start
		ref_position=ref_start

		for index in xrange(len(cons_seq)):
			if (cons_seq[index]==ref_seq[index] and (cons_seq[index]=="C" or cons_seq[index]=="c")):
				cons_positions_list.append(cons_position)
				ref_positions_list.append(ref_position)

			if (cons_seq[index]!="-"):
				cons_position+=1

			if (ref_seq[index]!="-"):
				ref_position+=1


		complete_subfamily=True
		rmsk_dict[subfamily][2]=ref_positions_list
		rmsk_dict[subfamily][1]=cons_positions_list

		return rmsk_dict

	if(rmsk_type=="Header"):
		if (line_list[8]!="C"):
			subfam_name = line_list[8]
		else: 
			subfam_name = line_list[9]

		subfamily, chrm, chrm_start, chrm_stop, cons_start, cons_stop = subfam_name[:subfam_name.index("#")], line_list[4][3:len(line_list[4])], line_list[5], line_list[6], "", ""

		if (len(chrm)==1):
			chrm="0"+chrm

		info=[chrm, chrm_start, chrm_stop, cons_start, cons_stop]

		

		rmsk_dict[subfamily] = [info, chrm_seq, cons_seq]
		check_white_space=False
		return rmsk_dict

def get_type(line_list, line):
	global rmsk_type

	if (len(line_list)!=0 and (line_list[0]=="Matrix" or line_list[0]=="Kimura" or line_list[0]=="Transitions" or line_list[0]=="Gap_init")):
		return "Footer"

	elif (len(line_list)!=0 and line[0]!=" " and line_list[0]!="Matrix" and line[0]!="C"): 
		return "Header"



def update_methylation_dictionary(rmsk, methyl_dict):
	rmsk_dict, methylation_dict= rmsk, methyl_dict
	methyl, total= 0, 0

	for key in rmsk_dict: 
		for rmsk_bp in xrange(len(rmsk_dict[key])): 
			rmsk=rmsk_dict[key]
			rmsk_loc_separator=rmsk[rmsk_bp].index("|")
			rmsk_tot_separator=rmsk[rmsk_bp].index("/")
			for methyl_bp in xrange(len(methyl_dict[key])):
				methyl=methylation_dict[key]
				methyl_loc_separator=methyl[methyl_bp].index("|")
				methyl_tot_separator=methyl[methyl_bp].index("/")

				if (rmsk[rmsk_bp][0:rmsk_loc_separator]==methyl_bp[methyl_bp][0:methyl_loc_separator]):
					methyl= int(rmsk[rmsk_bp][rmsk_loc_separator+1: rmsk_tot_separator])+int(methyl[methyl_bp][methyl_loc_separator+1: methyl_tot_separator])
					total=int(rmsk[rmsk_bp][rmsk_tot_separator+1:])+ int(methyl[methyl_bp][methyl_tot_separator+1])	
					methylation_dict[key][methyl_bp]=rmsk[rmsk_bp][0:rmsk_loc_separator]+ "%s/%s" %(methyl/total)

	return methylation_dict


def calculate_methylation_levels(rmsk_dictionary):
	rmsk_dict=rmsk_dictionary
	key=rmsk_dict.keys()[0]
	total, methyl=0, 0
	methylation_sequence=[]
	ref_pos_list=rmsk_dict[key][2]

	for ref_pos in xrange(len(ref_pos_list)):
		for read in xrange(3, len(rmsk_dict)): 
			sam_pos_list=rmsk_dict[key][read]
			for sam_pos in xrange(len(sam_pos_list)): 
				separatorIndex=sam_pos_list[sam_pos.index("|")]
				if (ref_pos_list[ref_pos]==sam_pos[sam_pos[0:separatorIndex]]): 
					total+=1 
					if (sam_pos[separatorIndex+1:]=="T"):
						methyl+=1

		if (total!=0):		
			methylation_sequence.append("%s|%s/%s" %(rmsk_dict[key][1][ref_pos], methyl, total))	

		total=0
		methyl=0
			
	rmsk_dict[key]=methylation_sequence
	#print rmsk_dict
	return rmsk_dict

def get_repeats_sam(rmsk_dictionary, sam_dictionary):
	rmsk_dict=rmsk_dictionary
	#print rmsk_dict
	sam_dict=sam_dictionary 

	sam_list=[]

	for item in sam_dict: 
		sam_list.append(item)

	sam_list=sort_sam_list(sam_list)
	rmsk_dict=valid_sam_repeats(rmsk_dict, sam_dict, sam_list)
	return rmsk_dict

def valid_sam_repeats(rmsk_dictionary, sam_dictionary, sam): 
	rmsk_dict, sam_dict, sam_list = rmsk_dictionary, sam_dictionary, sam

	rmsk_key=rmsk_dict.keys()[0]
	chrm, chrm_start, chrm_stop = rmsk_dict[rmsk_key][0], rmsk_dict[rmsk_key][1], rmsk_dict[rmsk_key][2]
	
	start_position=search_sam(sam_list, chrm, chrm_start,0)
	stop_position=search_sam(sam_list, chrm,chrm_stop,0)

	for index in xrange(start_position, stop_position): 
		chrm=sam_list[index][sam_list[index].index("|")+1:] 
		sequence=sam_dict[sam_list[index]]
		rmsk_dict[rmsk_key].append(sequence)
	rmsk_dict=reformat_appended_sam(rmsk_dict) 

	return rmsk_dict

def reformat_appended_sam(rmsk_dictionary):
	rmsk_dict=rmsk_dictionary 
	#print rmsk_dict
	key=rmsk_dict.keys()[0]
	for read in xrange(7, len(rmsk_dict[key])): 
		positionSamRead=0 
		reformattedSequence=""
		separatorIndex=rmsk_dict[key][read].index("|")
		chrm_start, chrm_stop, chrm_seq=int(rmsk_dict[key][1]), int(rmsk_dict[key][2]), rmsk_dict[key][6]
		sam_start= int(rmsk_dict[key][read][:separatorIndex])
		sam_stop= int(sam_start+len(rmsk_dict[key][read][separatorIndex+1:]))
		sam_sequence = rmsk_dict[key][read][separatorIndex+1:]
		sam_length=len(sam_sequence)

		for bp in xrange(chrm_stop-chrm_start):
			if (bp<(sam_start-chrm_start) or chrm_seq[bp]=="-"):
				reformattedSequence+="-"
			elif (chrm_seq[bp]!="-" and positionSamRead<sam_length):
				reformattedSequence+=sam_sequence[positionSamRead]
				positionSamRead+=1
			elif (bp> (sam_start-chrm_start+sam_length)): 
				reformattedSequence+="-"

		rmsk_dict[key][read]=reformattedSequence

	return rmsk_dict

#print reformat_appended_sam(rmsk_dict)

def search_sam(searchList, chromosome, chromosome_position, boundary_position):
	sam_list, chrm, chrm_pos, bound_pos = searchList, chromosome, chromosome_position, boundary_position
	midList=len(sam_list)/2

	separatorIndex=sam_list[midList].index("|")
	chromosomeSeparatorIndex=sam_list[midList].index("~")
	sam_chrm=sam_list[midList][0:separatorIndex]
	sam_chrm_pos=int(sam_list[midList][separatorIndex+1:chromosomeSeparatorIndex])

	rightHalf=sam_list[midList:]
	leftHalf=sam_list[:midList]

	if (len(sam_list)>1): 
		if (chrm<sam_chrm):
			bound_pos+=search_sam(leftHalf, chrm, chrm_pos, bound_pos)

		if (chrm>sam_chrm):
			bound_pos+=midList+search_sam(rightHalf, chrm, chrm_pos, bound_pos)

		if (chrm==sam_chrm): 
			if (int(chrm_pos)<sam_chrm_pos):
				bound_pos+=search_sam(leftHalf, chrm, chrm_pos, bound_pos)

			if (int(chrm_pos)>sam_chrm_pos):
				bound_pos+=midList+search_sam(rightHalf, chrm, chrm_pos, bound_pos)

	else:
		if (chrm<sam_chrm):
			bound_pos+=0

		if (chrm>sam_chrm):
			bound_pos+=1

		if (chrm==sam_chrm): 
			if (int(chrm_pos)<sam_chrm_pos):
				bound_pos+=0

			if (int(chrm_pos)>sam_chrm_pos):
				bound_pos+=1

	return bound_pos

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
				#print "equal"
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

def convert_methylation_dict(methyl_dict):
	methylation_dict=methyl_dict
	for key in methylation_dict: 
		methyl=methylation_dict[key]
		for methyl_bp in xrange(len(methyl)):
			methyl_loc_separator=methyl[methyl_bp].index("|")
			methyl_tot_separator=methyl[methyl_bp].index("/")

			print methyl[methyl_bp]
			position=methyl[methyl_bp][0:methyl_loc_separator]
			methyl_number=float(methyl[methyl_bp][methyl_loc_separator+1: methyl_tot_separator])
			total=float(methyl[methyl_bp][methyl_tot_separator+1:])

		print key 
		#print position + "\t" + "%s" %(methyl_number/total)



# Input conditions: rmsk_dict that has, how each subfamily, the consensus sequence 
# appended to the item list.
# 
# Output conditions: rmsk_dict with methylation patterns for all C's  

sam_dict=create_sam()
rmsk_dict={}
#print sam_dict
#print get_repeats_sam(rmsk_dict, sam_dict)	
main()
#analyze_rmsk()
