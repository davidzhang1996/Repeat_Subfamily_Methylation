"""
V1: 6/16/15 

"""
import math 

global methylation_dict
methylation_dict={}
rmsk_type=""
check_white_space=True 

def main():
	get_subfam_consensus_seqs()
	sam_dict=create_sam()
	parse_rmsk(sam_dict)
	#convert_methylation_dict()

def get_subfam_consensus_seqs():
	global methylation_dict 
	methylation_dict["test"]="Yes"
	with open ("repeat_consensus_seqs_1000.txt", "r") as repeats_file:
		subfam_name,seq = "","" 
		addToDict=False
		sequence=False
		for line in repeats_file: 
			line=line.rstrip("\n")
			line_list=line.split()

			if (line_list[0]=="ID"):
				sequence=False
				subfam_name=line_list[1]
				if (addToDict): 
					methylation_dict[subfam_name]=[seq]
					methylation_sequence=""
					for iterations in xrange(len(seq)):
						if (seq[iterations]!="C" or seq[iteration]!="c"):
							methylation_sequence+="-"
						else: 
							methylation_sequence+="0/0"
					methylation_dict[subfam_name].append(methylation_sequence)
					subfam_name, seq, methylation_sequence = "", "", ""

				addToDict=True
				
			if (sequence==True):
				for index in xrange(len(line_list)-1):
					seq+=line_list[index]

			if (line_list[0]=="SQ"): 
				sequence=True


def create_sam():
	sam_dict={}
	with open ("adipose_head.sam", "r") as sam_file:
		repeater=1
		for line in sam_file: 
			line=line.rstrip("\n")
			line_list=line.split()
			id=""

			if(line_list[0][0]!="@" and int(line_list[1])<=511 and int(line_list[1])>=256):
				chrm=line_list[2][3:len(line_list[2])]
				if (len(chrm)==1): 
					chrm="0"+chrm
				chrm_pos=line_list[3]

				id=line_list[2][3:len(line_list[2])]+"|" +line_list[3]
				if (sam_dict.has_key(id)):
					repeater+=1
				sam_dict[id+"."+"%s" %(repeater)]=[str(int(line_list[3])+len(line_list[9])), line_list[9]]

	return sam_dict


def parse_rmsk(dictionary):
	sam_dict=dictionary
	global check_white_space, complete_subfamily, line_count
	check_white_space, complete_subfamily, line_count=True, False, 1
	rmsk_dict={}

	with open ("hg19_1000000.fa.align", "r") as rmsk_file: 
		for line in rmsk_file: 
			rsmk_dict=get_repeats_rmsk(rmsk_dict, line)

			if (complete_subfamily):
				rmsk_dict=get_repeats_sam(rmsk_dict, sam_dict) 
				rmsk_dict=calculate_methylation_levels(rmsk_dict)

				subfam_name, cons_start, cons_stop, methyl_seq = rmsk_dict.keys()[0], rmsk_dict[3], rmsk_dict[4], rmsk_dict[7]
				update_methylation_dictionary(subfam_name, cons_start, cons_stop, methyl_seq)

				rmsk_dict.clear()
				line_count=1
				complete_subfamily=False

def update_methylation_dictionary(name, start, stop, seq):
	subfam_name, cons_start, cons_stop, methyl_seq = name, start, stop, seq
	repeat_seq=methylation_dict[subfam_name][0]

	for bp in xrange(len(methyl_seq)):
		if (methyl_seq[bp]==int and methyl_seq[bp]==repeat_seq[bp+cons_start]):
			methylation_dict[subfam_name][1][bp+cons_start]="%s" % (int(repeat_seq[bp+cons_start])+int(methyl_seq[bp])



def calculate_methylation_levels(rmsk_dictionary):
	rmsk_dict=rmsk_dictionary
	key=rmsk_dict.keys()[0]
	total, methyl=0, 0
	methylation_sequence=""
	cons_seq, ref_seq = rmsk_dict[key][5], rmsk_dict[key][6]

	for bp in xrange(len(cons_seq)): 
		if (cons_seq[bp]=="C" or cons_seq=="c" and cons_seq[bp]==ref_seq[bp]):
			for read in xrange(7,len(rmsk_dict[key])):
				if (rmsk_dict[key][read][bp]=="C" or rmsk_dict[key][read][bp]=="c" or rmsk_dict[key][read][bp]=="T" or rmsk_dict[key][read][bp]=="t"):
					total+=1 
					if (rmsk_dict[key][read][bp]=="T" or rmsk_dict[key][read][bp]=="t"):
						methyl+=1 

			methylation_sequence+="%s/%s" %(methyl, total)			
			total=0
			methyl=0
			
		else: 
			methylation_sequence+="-"
	rmsk_dict[key]=rmsk_dict[key][0:7]
	rmsk_dict[key].append(methylation_sequence)
	return rmsk_dict


def get_repeats_sam(rmsk_dictionary, sam_dictionary):
	rmsk_dict=rmsk_dictionary
	sam_dict=sam_dictionary 

	sam_list=[]

	for item in sam_dict: 
		sam_list.append(item)

	sam_list=sort_sam_list(sam_list)
	rmsk_dict=valid_sam_repeats(rmsk_dict, sam_dict, sam_list)
	return rmsk_dict

def valid_sam_repeats(rmsk_dictionary, sam_dictionary, sam): 
	rmsk_dict, sam_dict, sam_list = rmsk_dictionary, sam_dictionary, sam

	print rmsk_dict

	rmsk_key=rmsk_dict.keys()[0]
	chrm, chrm_start, chrm_stop = rmsk_dict[rmsk_key][0], rmsk_dict[rmsk_key][1], rmsk_dict[rmsk_key][2]
	
	start_position=search_sam(sam_list, chrm, chrm_start,0)
	stop_position=search_sam(sam_list, chrm,chrm_stop,0)

	for index in xrange(start_position, stop_position+1): 
		chrm=sam_list[index][sam_list[index].index("|")+1:] 
		sequence=chrm+"|"+sam_dict[sam_list[index][1]]
		rmsk_dict.append(sequence)
	rmsk_dict=reformat_appended_sam(rmsk_dict) 

	return rmsk_dict

def reformat_appended_sam(rmsk_dictionary):
	rmsk_dict=rmsk_dictionary 
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

print reformat_appended_sam(rmsk_dict)

def search_sam(searchList, chromosome, chromosome_position, boundary_position):
	sam_list, chrm, chrm_pos, bound_pos = searchList, chromosome, chromosome_position, boundary_position
	midList=len(sam_list)/2
	separatorIndex=sam_list[midList].index("|")
	sam_chrm=sam_list[midList][0:separatorIndex]
	sam_chrm_pos=int(sam_list[midList][separatorIndex+1:len(sam_list[midList])-2])

	rightHalf=sam_list[midList:]
	leftHalf=sam_list[:midList]

	print "Total list: ", sam_list, len(sam_list)
	print "Right half: ", rightHalf
	print "Left half: ", leftHalf

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

		#print "Right half:", rightHalf
		#print "Left half:", leftHalf

		i=0
		j=0
		k=0

		while (i<lhLength and j<rhLength): 
			lhSepInd=leftHalf[i].index("|")
			rhSepInd=rightHalf[j].index("|")

			#print leftHalf[i][:lhSepInd]
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

	#print "Merged List:", aList
	return aList

def get_repeats_rmsk(dictionary, string):
	global check_white_space, complete_subfamily, line_count
	rmsk_dict, line = dictionary, string

	line=line.rstrip("\n")
	line_list=line.split()

	rmsk_type, subfam_name, chrm, chrm_start, chrm_start, cons_start, cons_stop, chrm_seq, cons_seq = "", "", "", "", "", "", "", "", ""

	if (check_white_space==True):
		rmsk_type=get_type(line_list, line)

	if(check_white_space==False):
		subfamily=rmsk_dict.keys()[0]

		if (line_count%4==2):
			if (line_list[0]=="Matrix"):
				check_white_space=True

			else:
				chrm_seq=line_list[2]
				rmsk_dict[subfamily][6]+=chrm_seq

		elif(line_count%4==0):
			if (len(line_list[0])==1):
				line_list.remove(line_list[0])

			if(line_count==4):
				cons_start=line_list[1]
			cons_stop=line_list[3]	
			cons_seq=line_list[2]
			rmsk_dict[subfamily][5]+=cons_seq

		if (len(cons_start)!=0):
			rmsk_dict[subfamily][3] = cons_start

		if (len(cons_stop)!=0):
			rmsk_dict[subfamily][4] = cons_stop
		line_count+=1
		return rmsk_dict

	if(rmsk_type=="Footer" and line_list[0]=="Gap_init"):
		complete_subfamily=True
		return rmsk_dict

	if(rmsk_type=="Header"):

		if (line_list[8]!="C"):
			subfam_name = line_list[8]

		else: 
			subfam_name = line_list[9]

		subfam_name, chrm, chrm_start, chrm_stop = subfam_name[:subfam_name.index("#")], line_list[4][3:len(line_list[4])], line_list[5], line_list[6]
		if (len(chrm)==1):
			chrm="0"+chrm
		rmsk_dict[subfam_name] = [chrm, chrm_start, chrm_stop, cons_start, cons_stop, cons_seq, chrm_seq]
		check_white_space=False
		return rmsk_dict

def get_type(line_list, line):
	global rmsk_type

	if (len(line_list)!=0 and (line_list[0]=="Matrix" or line_list[0]=="Kimura" or line_list[0]=="Transitions" or line_list[0]=="Gap_init")):
		return "Footer"

	elif (len(line_list)!=0 and line[0]!=" " and line_list[0]!="Matrix" and line[0]!="C"): 
		return "Header"



# Input conditions: rmsk_dict that has, how each subfamily, the consensus sequence 
# appended to the item list.
# 
# Output conditions: rmsk_dict with methylation patterns for all C's  

sam_dict=create_sam()
rmsk_dict={}
#print sam_dict
print get_repeats_sam(rmsk_dict, sam_dict)	
#main()
#analyze_rmsk()
