"""
V1: 6/16/15 

"""
import math 

sam_dict, methylation_dict={}, {}
rmsk_type=""
check_white_space=True 

def main():
	create_sam()
	create_rmsk()
	#analyze_rmsk()

def create_sam():
	with open ("adipose_head_400.sam", "r") as sam_file:
		repeater=1
		for line in sam_file: 
			line=line.rstrip("\n")
			line_list=line.split()
			id=""

			if(line_list[0][0]!="@" and int(line_list[1])<=511 and int(line_list[1])>=256):
				id=line_list[2][3:len(line_list[2])]+line_list[3]
				if (sam_dict.has_key(id)):
					repeater+=1
				sam_dict[id+"."+"%s" %(repeater)]=[str(int(line_list[3])+len(line_list[9])), line_list[9]]

def create_rmsk(): 
	rmsk_dict
	subfam_name, chrm, chrm_start, chrm_start, cons_start, cons_stop, chrm_seq, cons_seq = "", "", "", "", "", "", "", ""
	methylation_sequence=""

	with open ("hg19_test4.fa.align", "r") as rmsk_file: 
		global rmsk_type
		global check_white_space

		line_count=1
		repeat=1
		for line in rmsk_file: 
			line=line.rstrip("\n")
			line_list=line.split()

			if (check_white_space==True):
				rmsk_type=get_type(line_list, line)

			if(check_white_space==False):
				if (line_count%4==2):
					if (line_list[0]=="Matrix"):
						check_white_space=True
					chrm_seq+=line_list[2]

				elif(line_count%4==0):
					if (len(line_list[0])==1):
						line_list.remove(line_list[0])

					if (abs(int(line_list[3])-int(cons_stop))<=1):
						check_white_space=True
					cons_seq+=line_list[2]
				line_count+=1

			if(rmsk_type=="Footer" and line_list[0]=="Matrix"):
				if (rmsk_dict.has_key(subfam_name)):
					repeat+=1
				rmsk_dict[subfam_name] = [chrm, chrm_start, chrm_stop, cons_start, cons_stop, cons_seq, chrm_seq]
				subfam_name, chrm, chrm_start, chrm_start, cons_start, cons_stop, chrm_seq, cons_seq= "", "", "", "", "", "", "", ""
				rmsk_type=""

			if(rmsk_type=="Header"):
				line_count=1

				if (line_list[8]!="C"):
					subfam_name, cons_start, cons_stop = line_list[8], line_list[9], line_list[10]
				else: 
					subfam_name, cons_start, cons_stop = line_list[9], line_list[11], line_list[12]

				subfam_name, chrm, chrm_start, chrm_stop = subfam_name[:subfam_name.index("#")], line_list[4], line_list[5], line_list[6]
				check_white_space=False
				rmsk_type=""
			#line_number+=1


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

"""
rmsk_dict={
			"example1": ['0', '0', '0', '0', '0', "ACTGCA-CGC", "AGTCCAAGTC", "GTT--AATGT", "---GTAATGT", "ATTCC------"],
			"example2": ['0', '0', '0', '0', '0', "ATTC-CACC-GA", "CTACATGTC---", "--ATTGCCTGAC", "ATTACGATTACG", "TATCA-GGCGTA"]
}
"""
"""

def analyze_rmsk():
	for key in rmsk_dict:
		methyl_seq_list, methylation_level_list = [], []
		total_methylation, total_reads, methylation_composite = 0.0, 0.0, 0.0 

		for index in xrange(5, len(rmsk_dict[key])):
			methyl_seq_list.append(rmsk_dict[key][index])

		#print methyl_seq_list
		for index1 in xrange(len(methyl_seq_list[0])): 
			#print index1
			if (methyl_seq_list[0][index1]!="-"): 	
				if (methyl_seq_list[0][index1]=="C" or methyl_seq_list[0][index1]=="c"):
					for index2 in xrange(1,len(methyl_seq_list)):
						if (methyl_seq_list[index2][index1]=="C" or methyl_seq_list[index2][index1]=="c" or methyl_seq_list[index2][index1]=="T" or methyl_seq_list[index2][index1]=="t"):
							total_reads+=1
						if (methyl_seq_list[index2][index1]=="T" or methyl_seq_list[index2][index1]=="t"):
							total_methylation+=1 

					if (total_reads!=0.0):
						methylation_composite=total_methylation/total_reads
						methylation_level_list.append("%.4s"%(methylation_composite))
						total_reads, total_methylation = 0.0, 0.0
				else: 
					methylation_level_list.append("-")


		methylation_dict[key]=methylation_level_list
		methylation_sequence=""

"""	
main()
#analyze_rmsk()
