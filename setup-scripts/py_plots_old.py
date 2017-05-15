#TYPE: python_command
# script for plotting multiple xvg files with x vs y format (or any text file with this format)
# includes plotting function, then function to read command line args and call plotting function
 
import sys, getopt, os
from pdb import *
import numpy as np
from matplotlib.pyplot import *
#import quick_ploto

############### Syntax (including optional flags): 

#py_plot.py -f <file_names> -e <error_data_names (opt)> -xlim <xlims_text_file (opt)> -ylim <ylims_text_file (opt)> -edit (opt) -norm (opt) -i <info file (opt)> -l <legend file (opt)> -o <fig_file_name_out)

############## notes: 
# requires separate -f for each file
# if data file names include "mn" and error bar file names include "err" or "std" and are in same dir as data files, can simply use -e err or -e std and it will automatically plot errorbars
# default number of error bars is 10 (evenly spaced) , if use -errs flag will plot all available error bars
# x/y axis labels and title specified in a text file, otherwise will look for labels in xvg files, format: xlabel <label_here>, title <title_here>, etc, all on separate lines
# x/y axis limits can be specified either in info file or on command line
# legend entries specified in text file (each entry on a new line)
# edit flag triggers breakpoint before plotting
# norm flag normalizes data by first entry
# clus flag not used
# if include output file, will save fig without showing, else, will show fig

############## tips:
# if plotting a bunch of similar files, it can be helpful to write a script to call py_plot.py with the various input files 
# use latex notation to get mathy labels or legend entries, ex '$r^{2}=3\sigma$'


def py_plot1(f_names_in, f_names_err_in, fig_name1, xlab1, ylab1, title1, leg1, xlims, ylims, opt_all_error, opt_edit, opt_points,opt_norm):

	# for labels: if specified, use labels, if not, try to find in xvg, if not, no labels
	
	n_files = 0
	for i in range(0,len(f_names_in)):
		if len(f_names_in[i]) > 0:
			n_files =n_files+1		
	opt_err =1
	if len(f_names_err_in) == 0:
		opt_err = 0
	
	colors1 = ['b','g','r','c','m']
	
	xs = []; ys = []
	for n in range(0,n_files):
		
		f_name_in = f_names_in[n]
		
		dr = open(f_name_in,'r')
		
		x = np.zeros(1000000)
		y = np.zeros(1000000)
		c1 = 0
		
		xvg_title1 = 'na'; xvg_xlab1 = 'na'; xvg_ylab1 = 'na'
		for line in dr:
			dl = line.split()
			if len(dl) > 1 and dl[1] == 'title':
				i1 = line.find('"');i2 = line.rfind('"')
				xvg_title1 = line[i1+1:i2]
				
			if len(dl) > 1 and dl[1] == 'xaxis':
				i1 = line.find('"');i2 = line.rfind('"')
				xvg_xlab1 = line[i1+1:i2]
				
			if len(dl) > 1 and dl[1] == 'yaxis':
				i1 = line.find('"');i2 = line.rfind('"')
				xvg_ylab1 = line[i1+1:i2]	
			#if len(dl) > 1 and dl[0] == ids[0] and dl[1] == ids[1]:
				#break
			if len(dl) == 2:
				try:
					x[c1] = float(dl[0])
					y[c1] = float(dl[1])
					c1 = c1+1
					break
				except:
					foo = 1
			
		for line in dr:
			dl = line.split()
			x[c1] = float(dl[0])
			y[c1] = float(dl[1])
			c1 = c1+1

		x = x[0:c1]
		y = y[0:c1]
		xs.append(x)
		ys.append(y)
		dr.close()
		if opt_norm == 1:
			y_plot = y/np.sum(y)
		else:
			y_plot = y
		if opt_points ==1:
			plot(x,y_plot,'.')
		else:
			plot(x,y_plot)
		
		hold
	
	#if leg1 != 'na':
		#legend(leg1,loc=0)	
	
	for n in range(0,n_files):
		n_err_bar = 10
		x0 = xs[n]; y0 = ys[n]
		#set_trace()
		if opt_all_error == 1:
			i0 = range(0,len(x0))
		else:
			i0 = map(int, (np.round(np.linspace(1,len(x0)-1,n_err_bar))) )
		x = x0[i0]; y = y0[i0]

		# load in errors, if exitst
		if opt_err == 1:	
			f_name_err_in = f_names_err_in[n]
			dr = open(f_name_err_in,'r')

			y_err = np.zeros(1000000)
			c1 = 0
			for line in dr:
				dl = line.split()
				
				#if len(dl) > 1 and dl[0] == ids[2] and dl[1] == ids[3]:
					#break
				if len(dl) == 2:
					try:
						y_err[c1] = float(dl[1])
						c1 = c1+1
						break
					except:
						foo = 1
			for line in dr:
				dl = line.split()
				#x[c1] = float(dl[0])
				y_err[c1] = float(dl[1])
				c1 = c1+1

			y_err = y_err[0:c1]
			y_err = y_err[i0]
			dr.close()
		#set_trace()	
		if opt_err==1:   	# time series
			errorbar(x,y,yerr=y_err,fmt=None,ecolor=colors1[n]),	
		
		hold
		
		if n == n_files -1:
			if opt_edit == 1:
				set_trace()
				## errorbar(9,2.42,yerr=0.04,fmt='bx',ecolor='b')
				##  errorbar(9,2.09,yerr=0.045,fmt='go',ecolor='g') 
				## errorbar(9,1.69,yerr=0.053,fmt='ro',ecolor='r')
				print "edit..."
			
			if xlab1 != 'na':
				xlabel(xlab1)
			elif xvg_xlab1 != 'na':
				xlabel(xvg_xlab1)
			
			if ylab1 != 'na':
				ylabel(ylab1)
			elif xvg_ylab1 != 'na':
				ylabel(xvg_ylab1)	
			
			if title1 != 'na':
				title(title1)
			elif xvg_title1 != 'na':
				title(xvg_title1)
			
			if xlims != 'na':
				xlim(xlims[0],xlims[1])		
			if ylims != 'na':
				ylim(ylims[0],ylims[1])
			if leg1 != 'na':
				legend(leg1,loc=0)	
						
			if fig_name1 != 'na':
				savefig(fig_name1)
				print fig_name1+' ... saved'
			else:
				show()
				print "not saving..."
			#set_trace()

def main(argv):
	f_name_in = ''
	f_name_out = ''
	try:
		
		opts, args = getopt.getopt(argv,"hf:e:o:i:l:p",["xlim=","ylim=","errs","edit","clus=","norm"])
		
		# help, files_in, errs_in, ouptput, 
		# info: xlabel, ylabel, title, xlims, ylims
		# legend, save(01)
	except getopt.GetoptError:
		print 'py_plot.py -f# <inputfiles> -e# <errorfiles> -o <outputfile> -i <info_file> -l <legend_file>'
		sys.exit(2)
		
	# Defaults
	path1 = os.getcwd()
	
	f_name_info = 'na'
	xlab1 = 'na'
	ylab1 = 'na'
	title1 = 'na'

	xlims = 'na'; ylims = 'na'
	
	f_name_leg = 'na'
	leg1 = 'na'
	
	fig_name = 'na'
	opt_all_error = 0
	opt_edit = 0 
	opt_points = 0
	opt_norm = 0
	
	f_names_in = []; e_names_in = []
	for opt, arg in opts:
		if opt == '-h':
			print 'py_plot.py -f <inputfiles> -e <errorfiles/same> -o <outputfile.png>  -i <info_file> -l <legend_file>'
			sys.exit()
			
		elif opt in ("-f"):
			f_name_ino = arg
			f_names_in.append(f_name_ino)
		
		elif opt in ("-e"):
			e_name_ino = arg
			e_names_in.append(e_name_ino)
			
		elif opt in ("-o"):
			fig_name = arg
		elif opt in ("-i"):
			f_name_info = arg
		elif opt in ("-l"):
			f_name_leg = arg
		elif opt in ("-p"):
			opt_points = 1
		
	if f_name_info != 'na':
		dr = open(f_name_info,'r')
		for line in dr:
			dl = line.split()
			if len(dl) > 0 and dl[0] == 'xlabel':
				i1 = line.find('"');i2 = line.rfind('"')
				xlab1 = line[i1+1:i2]
			if len(dl) > 0 and dl[0] == 'ylabel':
				i1 = line.find('"');i2 = line.rfind('"')
				ylab1 = line[i1+1:i2]
			if len(dl) > 0 and dl[0] == 'title':
				i1 = line.find('"');i2 = line.rfind('"')
				title1 = line[i1+1:i2]
			if len(dl) > 0 and dl[0] == 'xlims':
				xlims = [float(dl[1]), float(dl[2])]
				print "xlims: "+str(xlims[0])+' '+str(xlims[1])
			if len(dl) > 0 and dl[0] == 'ylims':
				ylims = [float(dl[1]), float(dl[2])]
				print "ylims: "+str(ylims[0])+' '+str(ylims[1])	
		dr.close()
	
	for opt, arg in opts:
		if opt == '--xlim':
			str1 = arg
			i1 = str1.find(',')
			
			xlims = [float(str1[0:i1]), float(str1[i1+1:len(str1)])]
		if opt == '--ylim':
			str1 = arg
			i1 = str1.find(',')
			ylims = [float(str1[0:i1]), float(str1[i1+1:len(str1)])]
		if opt == '--errs':
			opt_all_error = 1
		if opt == '--edit':
			opt_edit = 1
		if opt == '--norm':
			opt_norm = 1
		if opt == '--clus':	# NOT USED FOR NOW
			p_name0 = arg  #'/bigtmp/jb2qm/therm_stat/'
			dw = open('/mnt/hds/tmp_plot/get_figs.py','w')
			dw.write('import os\n')
			cmd1 = 'cp '
			for j in range(0,len(f_names_in)):
				cmd1 = cmd1+' '+p_name0+f_names_in[j]
			cmd1 = cmd1+' ./'
			dw.write('os.system(cmd1)\n')
			dw.close()
	
	if f_name_leg != 'na':
		leg1 = []
		dr = open(f_name_leg,'r')
		for line in dr:
			leg1.append(line[0:-1])	
		dr.close()
	
	if len(e_names_in) > 0 and e_names_in[0] == 'std':
		e_names_in = []
		for n in range(0,len(f_names_in)):
			f_name0 = f_names_in[n]
			e_name0 = f_name0.replace('mn','std')
			e_names_in.append(e_name0)
			print "reading error from "+e_name0+"..."
	if len(e_names_in) > 0 and e_names_in[0] == 'err':
		e_names_in = []
		for n in range(0,len(f_names_in)):
			f_name0 = f_names_in[n]
			e_name0 = f_name0.replace('mn','err')
			e_names_in.append(e_name0)
			print "reading error from "+e_name0+"..."
	
	# call function
	py_plot1(f_names_in, e_names_in, fig_name, xlab1, ylab1, title1, leg1, xlims, ylims, opt_all_error, opt_edit, opt_points,opt_norm)	
	print "...done"

if __name__ == "__main__":
   main(sys.argv[1:])
