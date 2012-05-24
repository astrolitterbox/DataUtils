#!/usr/bin/env python
# -*- coding: utf-8 -*-

def print_usage():
	print """NED batch parser is to parse output files from NED batch script. It works only 
for "compact" output type!!! See http://nedwww.ipac.caltech.edu/help/batch.html"""
	print """Min-Su Shin (Princeton University) developed the code in 2008. If you have any questions 
and suggestions, feel free to contact: msshin (at) umich.edu (a new email address)."""
	print " "
	print "NED_batch_compact_parser.py [NED result file] [output columns] [-s] [-o output filename]"
	print """output columns:	Type = 0
		Object Name = 1
		Position = 2
		Velocity of z = 3
		Reference number = 4
		Field 1 = 5
		Position uncertainty = 6
		Position angle = 7
		Uncertainty = 8
		Photometry = 9
		Field 2 = 10
		Field 3 = 11
		Field 4 = 12"""
	print "-s : print only the closest object. A default setting is printing out all objects."
	print "-o : output filename. If it's not provided, standard output is default."
	print " "
	print "Example : NED_batch_compact_parser.py 3arcsec.ned 0,1,2,3 -s -o test.out"


import sys
import getopt
import re
import string


try:
	result_file = sys.argv[1]
	output_cols = string.split(sys.argv[2], ',')
except:
	print_usage()
	sys.exit()
output_check = []
for ind in range(0,13):
	output_check.append(0)
for ind in output_cols:
	i = int(ind)
	output_check[i] = 1

single_output = 0
try:
	opts, args = getopt.getopt(sys.argv[3:], "so:v", ["single match", "output="])
	for o, a in opts:
		if o == "-s":
			single_output = 1
except:
	single_output = 0

standard_output = 1
try:
	opts, args = getopt.getopt(sys.argv[3:], "so:v", ["single match", "output="])
	for o, a in opts:
		if o == "-o":
			standard_output = 0
			output_fn = a
except:
	standard_output = 1

print "Input file : ",result_file
fd = open(result_file,'r')
result_begin=0
line_num=0
obj_num=0
# Counting objects
while 1:
	oneline = fd.readline()
	if not oneline : break
	line_num += 1
	return_val = re.findall('NEARPOSN', oneline)
	if len(return_val) > 0 :
		obj_num += 1
	return_val = re.findall('SEARCH RESULTS', oneline)
	if len(return_val) > 0 :
		result_begin = line_num
		break
print obj_num, "objects in the file..."

if standard_output < 1 :
	out_fd = open(output_fn, 'w')
	print "Output file : ",output_fn

print "Parsing the result section..."
obj_count = 0
while 1:
	oneline = fd.readline() # skip the first empty line
	return_val = re.findall('Names for all the objects found in this request', oneline)
	if len(return_val) > 0 : break
	oneline = fd.readline() # Position information
	obj_count += 1
	temp = string.split(oneline)
	temp = string.join(temp[1:])
	temp = string.split(temp, ';')
	ra = temp[0]
	ra = string.strip(ra)
	ra = ra[0:len(ra)-1]
	dec = string.strip(temp[1])
	search_radius = string.strip(temp[2])
	oneline = fd.readline() # Number of match
	#print oneline, "*************"
	return_val = re.findall('Object does not exist in the data base', oneline)
	if len(return_val) < 1:
		#print string.split(oneline)[0]
		num_of_objects = int(string.split(oneline)[0])
	else :
		num_of_objects = 0
	if num_of_objects > 0 :
		list_result=[]
		output_result=[]
		output_str = string.join([str(obj_count),ra,dec,search_radius],',')
		oneline = fd.readline() # first header line
		oneline = fd.readline() # second header line
		for i in range(0,num_of_objects):
			result1 = fd.readline() # first result line
			obj_type = result1[0:6]
			list_result.append(obj_type)
			obj_name = result1[7:23]
			list_result.append(obj_name)
			obj_pos = result1[24:46]
			list_result.append(obj_pos)
			obj_z = result1[47:55]
			list_result.append(obj_z)
			obj_ref = result1[56:59]
			list_result.append(obj_ref)
			obj_fd1 = string.strip(result1[61:])
			result2 = fd.readline() # second result line
			if result2[7] :
				obj_fd1 = obj_fd1 + result2[7:23]
			list_result.append(obj_fd1)
			obj_pos_unc = result2[24:42]
			list_result.append(obj_pos_unc)
			obj_posang = result2[43:46]
			list_result.append(obj_posang)
			obj_unc = result2[47:55]
			list_result.append(obj_unc)
			obj_phot = result2[56:59]
			list_result.append(obj_phot)
			obj_fd2 = result2[60:67]
			list_result.append(obj_fd2)
			obj_fd3 = result2[68:74]
			list_result.append(obj_fd3)
			obj_fd4 = string.strip(result2[75:])
			list_result.append(obj_fd4)
			for j in range(0, 13) :
				if output_check[j] > 0 :
					temp = string.strip(list_result[j])
					if len(temp) < 1 :
						temp = 'None'
					output_result.append(temp)
			if single_output > 0 :
				if i == 0:
					output_str = output_str + "," + string.join(output_result,',')
			else :
				output_str = output_str + "," + string.join(output_result,',')
		if standard_output < 1:
			out_fd.write(output_str+"\n")
		else:
			print output_str
	oneline = fd.readline() # skip the empty line
	oneline = fd.readline() # next search line
	if obj_count == obj_num :
		break
if standard_output < 1 :
	out_fd.close()
