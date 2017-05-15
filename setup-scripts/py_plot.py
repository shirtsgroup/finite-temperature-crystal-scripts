#TYPE: python_command
# script for plotting multiple xvg files with x vs y format (or any text file with this format)
# includes plotting function, then function to read command line args and call plotting function
 
import sys, getopt, os
from pdb import *
import numpy as np
import string

#sys.path.append("/mnt/hds/funcs/")
sys.path.append("/home/jb2qm/hds/funcs/")
import read_file as rf
#import write_file as wf


opt_plot_load = 1
if opt_plot_load == 1:
	from matplotlib.pyplot import *
	import plot_funcs

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


def py_plot1(f_names_in, f_names_err_in, f_names_xerr_in, n_colx, n_cols_in, n_col_ex, n_col_ey, fig_name1, xlab1, ylab1, label_size, leg_size, num_size, size_line, title1, leg1, xlims, ylims, dlims, opt_n_err, opt_edit, opt_points, opt_stats, t_stats0, opt_save_data, y_facs, y_norms, x_facs, x_norms, opt_slog, opt_plot, opt_csv, x_csv, opt_xvg_title, opt_inv, xvert, yhorz, colors1, str_cmd, opt_tacf, opt_hist, f_avg, opt_stop, opt_custom, leg_loc_in, opt_rev):

	opt_add_zero = 0 
	#ax_square = 0
	ax_square = 1 
	
	if opt_plot_load == 1:	
		fig = figure()
		ax = fig.add_subplot(111)

	# FILES
	n_files = 0
	for i in range(0,len(f_names_in)):
		if len(f_names_in[i]) > 0:
			n_files =n_files+1		
	N_cols = len(n_cols_in)
	if N_cols == 1 and n_files > 1:
		n_cols_in = n_cols_in*n_files
	N_cols = len(n_cols_in)

	if opt_rev == 1:
		f_names_in.reverse()
		f_names_err_in.reverse()
		f_names_xerr_in.reverse()

	# DISPLAY	
	opt_dash = 0 
	#opt_dash = 1  # turn on for dashes
	if n_files > 3:
		opt_dash = 0 

	if colors1 == 'na':
		
		colors1 = ['b','g','r','c','m','y','k','0.9','0.8','0.7','.6','.5','.4','.3','.2','.1']	
		symbols1 = ['.','s','o','^','v','d','*','+','.','.','.','.','.','.','.','.']
		#symbols1 = [u'\u2B14',u'\u2B15']	# black ll, ur

		mfcs = ['b','g','r','c','m','y','k','0.9','0.8','0.7','.6','.5','.4','.3','.2','.1']	
		#mfcs = ['b','w','g','w','r','w']
		size_lines = len(colors1)*[size_line]	#[2.5,2.5,2.5,1,1,1,1]
		style_lines = len(colors1)*['-']#['-','-','-','--','--','--','-']
	
		if opt_dash == 1:
			print "DASHES..."
			# isotherm scheme
			colors1 = 5*['b']
			colors1 = 5*['k']
			style_lines = ['-','--',':']##,'-','--','--','--']
			
			mfcs = ['b','w','g','w','goldenrod','w','r','w','c','c']
			symbols1 = ['o','o','o','o','o','o']
	elif colors1[0]  == 'iso':
		print "ISOTHERM COLORS..."
		# isotherm scheme
		#colors1 = ['b','b','g','g','goldenrod','goldenrod','r','r','c','c']
		#mfcs = ['b','w','g','w','goldenrod','w','r','w','c','c']
		#style_lines = ['-','--','-','--','-','--','--','--']
		## NEW
		#colors1 = ['b','g','goldenrod','b','g','goldenrod','r','r','c','c']
		#mfcs = ['b','g','goldenrod','w','w','w','c','c']
		#style_lines = ['-','-','-','--','--','--','-','--','--','--']

		#symbols1 = ['o','o','o','o','o','o','o','o']
		#symbols1 = ['o','s','^','o','s','^']
		#size_lines = len(colors1)*[size_line]	#[2.5,2.5,2.5,1,1,1,1]
			
		## NEW B/W
		## 3 epss
		#colors1 = 6*['k','k']
		#mfcs = ['k','k','k','w','w','w']
		#style_lines = ['-','-','-','--','--','--','-','--','--','--']
		#style_lines = ['-','-','-',':',':',':','-','--','--','--']
		#symbols1 = ['o','s','^','o','s','^']
		
		## PAPER1
		#coo = []
		#coo.append('steelblue')
		#coo.append('darkseagreen')
		#coo.append('darkkhaki')
		#
		#mfcs = coo
		#mfcs = mfcs+ 3*['w']
		#style_lines = ['-','-','-','--','--','--','-','--','--','--']
		#symbols1 = ['o','s','^','o','s','^']
		#symbols1 = ['o','s','^','o','s','^']

		##POLY PROPS
		coo = []
		mfcs = []
		coo.append('b')
		coo.append('g')
		coo.append('r')
		#coo.append('c')
		#coo = 4*['k']
		#mfcs = 4*['w']
		mfcs = coo
		style_lines = 4*['-']
		symbols1 = ['o','s','D']#o','s','^']
		

		#mfcs = 4*['w']
		#style_lines = 4*['-']

		colors1 = []
		for i in range(0,2):
			colors1 = colors1+coo

		size_lines = len(colors1)*[size_line]	#[2.5,2.5,2.5,1,1,1,1]

	# Opt to add letter (for publication)
	let_loc = 0
	#let_loc = [0.02, 0.92]

	let_add = '(a)'	
	#let_add = '(b)'	
	let_font = 20 

	#elif colors1[0] == 'lines':
	#	colors1 = ['b','g','r','c','m','y','k--','0.25','0.5','0.75']	# supports up to 6 curves

	# Error
	opt_err = 1	
	opt_xerr = 1	

	#opt_n_err = 1
	if len(f_names_err_in) == 0:
		opt_err = 0
	if len(f_names_xerr_in) == 0:
		opt_xerr = 0

	if opt_err == 1:
		if f_names_err_in[0] == 'std2err':
			opt_std2err = 1
			del f_names_err_in[0]
		else:
			opt_std2err = 0
		if f_names_err_in[0] == 'std':
			del f_names_err_in[0]
		n_err_bar = 10

	# CSV for labels: if specified, use labels, if not, try to find in xvg, if not, no labels
	if opt_csv == 1:
		dw_csv = open('tmp_csv.csv','w')
		c1 = 0
	# SCALING	
	if len(y_facs) == 1:
		#y_facs =  n_files*y_facs
		y_facs =  N_cols*y_facs
	if len(y_norms) == 1:
		#y_norms = n_files*y_norms
		y_norms = N_cols*y_norms
	y_facs = map(float,y_facs)
	y_norms = map(float,y_norms)
	
	if len(x_facs) == 1:
		x_facs = x_facs*N_cols#*n_files
	if len(x_norms) == 1:
		x_norms = x_norms*N_cols#*n_files
	x_facs = map(float,x_facs)
	x_norms = map(float,x_norms)
	
	xs = []; ys = []
	i0s = []

	# LOOP OVER DATA TO PLOT
	for n in range(0,N_cols):	# number of columns to plot
		if N_cols > n_files:
			f_name_in = f_names_in[0]
			y_facs = y_facs*N_cols
		#elif N_cols == n_files:
		else:
			f_name_in = f_names_in[n]
		n_cols = n_cols_in[n]

		#if N_cols < n_files:
		#	f_name_in = f_names_in[n]
		#n_cols = n_cols_in[0]*n_files

		i0s.append([])		
	
		for j in range(0,len(n_cols)):
	
			j_col = n_cols[j]
			xy_col_in = [n_colx, j_col]
	
			x,y, xvg_title1, xvg_xlab1, xvg_ylab1 = rf.gen_read_col_info(f_name_in, xy_col_in)	# read in data
	
		########### ERROR DATA, if exists
		if opt_err == 1:
			f_name_err_in = f_names_err_in[n]
			print "reading error from "+f_name_err_in

			if f_name_err_in == f_names_in[n]: # err in same file

				j_col = n_col_ey[0]

				#if f_names_err_in != f_names_xerr_in:
				#	j_col = 'last'
				#else:
				#	j_col = 'last'
				
				x_err0,y_err0 = rf.gen_read_col(f_name_err_in, j_col)	# read in file

			else:
				j_col = n_col_ey[n]
				n_cols = n_cols_in[n]
				
				for j in range(0,len(n_cols)):
				
					j_col = n_cols[j]
					#xy_col_in = [n_colx, j_col]
					x_err0,y_err0 = rf.gen_read_col(f_name_err_in, j_col)	# read in file
		if opt_xerr == 1:	
			f_name_xerr_in = f_names_xerr_in[n]
			print "reading x-error from "+f_name_xerr_in

			if f_name_xerr_in == f_names_in[n]:
				#if f_names_err_in != f_names_xerr_in:
				#	j_col = 'last'
				#else:
				#	j_col = 'slast'
				#	#j_col = 2

				#if n_col_ex > 0:
				#	j_col = n_col_ex

				x_xerr0,y_xerr0 = rf.gen_read_col(f_name_xerr_in, n_col_ex)	# read in file
			else:
				n_cols = n_cols_in[n]
				
				for j in range(0,len(n_cols)):
				
					j_col = n_cols[j]
					#xy_col_in = [n_colx, j_col]

					x_xerr0,y_xerr0 = rf.gen_read_col(f_name_xerr_in, j_col)	# read in file

		# SCALING
		if opt_inv == 1:
			y = np.divide(1.,y)
		
		y = y*y_facs[n]
		if y_norms[n] == 0:
			y = y/y[0]	# normalize by y(x=0)
		else:	
			y = y/y_norms[n]
		x = x*x_facs[n]
		if x_norms[n] == 0:
			x = x/x[0]	# normalize by y(x=0)
		else:	
			x = x/x_norms[n]


		x_err = 'na'
		y_err = 'na'
		x_xerr = 'na'
		y_xerr = 'na'
		if opt_err == 1:	
			y_err0 = y_err0*y_facs[n]	
			if opt_std2err == 1:
				# try to find N_err
				N_err = 0
				dr0 = open(f_name_err_in,'r')
				for line in dr0:
					if line[0] == 'N':
						N_err = int(line[2:-1])
				dr0.close()
	
				if N_err == 0:
					N_err = float(raw_input('Enter number of data points for error...'))
				else:
					print "Found N_err: "+str(N_err)
				y_err0 = np.divide(y_err0, np.sqrt(N_err-1))

		# Determine how much data to consider
		if dlims == 'na':
			x_plot = np.copy(x)
			y_plot = np.copy(y)
			if opt_err == 1:
				x_err = np.copy(x_err0)
				y_err = np.copy(y_err0)
			if opt_xerr == 1:
				x_xerr = np.copy(x_xerr0)
				y_xerr = np.copy(y_xerr0)
		else:
			if dlims[0] == -1:
				foo = 0 #reduces number of input data points, USE WITH CAUTION!
				#i_kp = np.arange(0,len(x),3)	
				i_kp = ((x >= dlims[0]) & (x <= dlims[1]))
			else:
				i_kp = ((x >= dlims[0]) & (x <= dlims[1]))

			x_plot = x[i_kp]
			y_plot = y[i_kp]
			if opt_err == 1:
				x_err  = x_err0[i_kp]
				y_err  = y_err0[i_kp]
			if opt_xerr == 1:
				x_xerr  = x_xerr0[i_kp]
				y_xerr  = y_xerr0[i_kp]
		#xs.append(x)
		#ys.append(y)

		# PLOT
		if opt_plot_load == 1 and opt_plot == 1:	
			coln = colors1[n]
			symn = symbols1[n]
			lsn = style_lines[n]
			lwn = size_lines[n]
			mfcn = mfcs[n]
			legn = 'na'
			if leg1 != 'na':
				if leg1 == 'fnames':
					legn = f_name_in
				else:
		
					legn = leg1[n]
			if opt_add_zero:
				ia0 = 0
				x_plot = np.insert(x_plot, ia0, 0.0)
				y_plot = np.insert(y_plot, ia0, 0.0)
				#x_err = np.insert(x_err, ia0, 0.0)
				#y_err = np.insert(y_err, ia0, 0.0)

			################CUSTOMIZE PLOT
			# LOOK FOR OSURF OR OSYS
			if opt_custom == 'qhat':
				print "###############qhat ISO FIT"

				opt_points[1] = 'custom1'	

			if opt_custom == 'cp':
				print "############### cp ISO FIT"
				#cp	
				i_osys = f_name_in.find('osys')
				i_10 = f_name_in.find('10')

				#if i_osys >= 0 and i_10 >= 10:
				#	opt_points[1] = 'langmuir_lin'	
				#else:
				#	opt_points[1] = 'lin0'	
				if i_osys >= 0:
					opt_points[1] = 'langmuir'	
				else:
					opt_points[1] = 'lin0'	


			ax = plot_funcs.plot_pyplot(f_name_in, ax, x_plot, y_plot, x_err, y_err, x_xerr, y_xerr, coln, symn, mfcn, lsn, lwn, opt_points, opt_n_err, opt_hist, opt_slog, legn, opt_stop)
			#hold


	#for n in range(0,n_files):
	#	try:
	#		x0 = xs[n]; y0 = ys[n]
	#	except:
	#		set_trace()
	#		foo = 0		

	#	if opt_all_error == 1:
	#		i0 = range(0,len(x0))
	#	else:
	#		if len(i0s) == 0 or len(i0s) == len(f_names_in):
	#			i0 = map(int, (np.round(np.linspace(1,len(x0)-1,n_err_bar))) )
	#			if len(x0) <= n_err_bar:
	#				i0 = range(0,len(x0))
	#		else:
	#			i0 = []
	#			i0n = i0s[n]
	#			for ni in range(0,len(i0n)):
	#				i00 = np.where(x0 == i0n[ni])



	#		ax.errorbar(x,y,yerr=y_err[i0],fmt=None,ecolor=colors1[n]),	
		if opt_plot_load == 1:
			#hold

			if opt_stats != 0:
				save_stats = 0
				if t_stats0 == 'na' or t_stats0 == 'save':
					t_lim0 = x[0]
					t_lim1 = x[-1]
					save_stats = 1
				else:
					# could either be st/end or st/end/inc (with sv option)
					t_lims = t_stats0.split(',')

					if len(t_lims) >= 3 and t_lims[-1] == 'save' or t_lims[-1] == 'sv':
						save_stats = 1
					
					if t_lims[0][0] == 'l': # take last x amount of t-series
						t_lim0 = x[-1] - float(t_lims[0][1:])
						t_lim1 = x[-1]
					else:
						t_lim0 = t_lims[0]
						t_lim1 = t_lims[1]
						
						if t_lim0 == ':':
							t_lim0 = x[0]
						else:
							t_lim0 = float(t_lim0)
						
						if t_lim1 == ':':
							t_lim1 = x[-1]
						else:
							t_lim1 = float(t_lim1)
				
				ix =  (x >= t_lim0) & (x <= t_lim1) 
				if len(ix) == 0:
					print "no data..."
				else:
					x_stats = x[ix]
					y_stats = y[ix]
					fit1 = np.polyfit(x_stats,y_stats, 1)

					p1a = fit1[0]	# slope
					p1b = fit1[1]	# yint			
					#print (p1a,p1b)
					#set_trace()			
				
					##HC
					#y_stats = y_stats[y_stats<0.01]	
					#x_stats = x_stats[y_stats<0.01]	
						
					y_mn = np.mean(y_stats)
					y_std = np.std(y_stats)
					y_err = y_std/np.sqrt(len(y_stats))
					y_std_per = y_std/y_mn
					y_area = np.trapz(y_stats,x_stats)
					
					if n == 0:
						y_mns = []
					y_mns.append(y_mn)

				if opt_stats == 'tacf':# == 1:
					t_stats = map(float,t_stats0.split(',')) #t1,t2,t_st_mn
					fig_acf = figure()
					import general_py
					#y_acf_in = y_stats-y_mn
					y_acf_in = y_stats/y_stats[0]
					y_mn_end = np.mean(y_stats[x_stats>t_stats[-1]])
					y_acf_in = y_stats-y_mn_end

					# INPUTS
					#c_ads
					n_pts_acf = 100#190
					dt_d = 2*1200 
					#tlims_acf = [0, n_pts_acf*dt_d]
					#tlims_acf = [0, 0.75*x[-1]]
					tlims_acf = [0, 360000]
					tlims_acf = [0, 3E5]
					tlims_acf = [x[0], x[-1]]


					#dt_d = 3 
					#tlims_acf = [0, 3000]
					## prot U
					#dt_d = 0.6
					#tlims_acf = [0., 60000.]
										
					tlims_d = [x[0], x[-1]]
					## if x is not time
					#tlims_d = [x[0]*dt_d, x[-1]*dt_d]

					del_to = 1
					opt_acf_norm = 0 
					 
					
					
					f_name_out = f_name_in+'acf'
					n_acf_ref = 1; opt_acf_norm = 0

					t_acf,y_acf = general_py.acf_direct1(y_acf_in, dt_d, tlims_acf[0], tlims_acf[1], tlims_d[0], tlims_d[1], del_to, opt_acf_norm)

					out1 = general_py.fit_model(t_acf[t_acf<60000], y_acf[t_acf<60000], [100, -0.01], 'exp0', [])

					
					p_fit = out1[0]
					print p_fit
					y_fit = p_fit[0]*np.exp(p_fit[1]*t_acf)
					tau = -1./p_fit[1]
					print tau

					print "plotting acf..."
					plot(t_acf, y_acf, t_acf, y_fit); show()
					set_trace()				
					wf.gen_write(f_name_out, t_acf, y_acf, [])
					print "ACF done..."

				print f_name_in
				print "x_range: "+str(t_lim0)+'  '+str(t_lim1)
				print "########## MEAN/ERR: "
				print "y_mean:     "+str(y_mn)
				print "y_stdev_raw:    "+str(y_std)
				print "y_err_raw:    "+str(y_err)
				print "y_stdev_per:    "+str(y_std_per)
				print "y_area:     "+str(y_area)+'\n'
				print "########## LIN FIT: "
				print "slope:     "+str(p1a)+'\n'
				print "y-int:     "+str(p1b)+'\n'
				print "slope/6 (cm2/s):     "+str(p1a/6*0.01)+'\n'
				print "slope/4 (cm2/s):     "+str(p1a/4*0.01)+'\n'
				if save_stats == 1:
					df_name_in = f_name_in.split('.')
					f_name_stats = 'stats_'+f_name_in.replace(df_name_in[-1],'txt')
					dw_stats = open(f_name_stats,'w')
					dw_stats.write('## Stats from \n')
					dw_stats.write('## t1: '+str(t_lim0)+'   t2:'+str(t_lim1)+'\n')
					dw_stats.write('y_mn:        '+str(y_mn)+'\n')
					dw_stats.write('y_std:       '+str(y_std)+'\n')
					dw_stats.write('y_stdev_per: '+str(y_std_per)+'\n')
					dw_stats.write('y_area:      '+str(y_area)+'\n')
					dw_stats.close()
				if opt_csv == 1:
					#l1 = "%8s %8.8f\n" %(x_csv[c1], y_mean)]
					opt_mn_diff = 0
					
					if n > 0 and opt_mn_diff == 1:
						y_mn_diff = y_mns[n]-y_mns[n-1]
						l1 = "%8s %1s %8.8f %1s %8.8f %1s %8.8f\n" %(x_csv[c1], ',', y_mn, ',', y_std, ',', y_mn_diff)
						y_mn_diffs.append(y_mn_diff)
						
					else:
						y_mn_diffs = []
						l1 = "%8s %1s %8.8f %1s %8.8f\n" %(x_csv[c1], ',', y_mn, ',', y_std)
					dw_csv.write(l1)
					
					c1+=1
					if n == n_files-1 and opt_mn_diff == 1:
						l1 = "%8s %1s %8.8f %1s %8.8f %1s %8.8f %1s %8.8f\n" %('Mn(std)', ',', np.mean(y_mns), ',', np.std(y_mns), ',', np.mean(y_mn_diffs),',',np.std(y_mn_diffs))
						dw_csv.write(l1)

			if f_avg != '0' and n == N_cols-1:
				ys_all = np.array(ys)
				y_avg = np.mean(ys_all,axis=0)						

				if f_avg == '1':
					plot(x,y_avg); show()
				else:
					wf.gen_write(f_avg,x,y_avg,[])	
			if opt_save_data != 0:
				if opt_save_data == 'same':
					f_name_sv = f_name_in+'sv'
				elif opt_save_data == 'err':
					f_name_sv = f_name_in+'err'
				else:
					f_name_sv = opt_save_data#f_name_in+'sv'
				f_name_sv = wf.gen_write(f_name_sv, x, y_plot, [str_cmd])

		# Add labels/edit	
		#if n == n_files -1 and opt_plot_load == 1:
		if n == N_cols -1 and opt_plot_load == 1:
			if opt_edit == 1:
				set_trace()
				## errorbar(9,2.42,yerr=0.04,fmt='bx',ecolor='b')
				##  errorbar(9,2.09,yerr=0.045,fmt='go',ecolor='g') 
				## errorbar(9,1.69,yerr=0.053,fmt='ro',ecolor='r')
				print "edit..."

			if xlab1 != 'na':
	
				xlabel(xlab1,fontsize=label_size)
				gcf().subplots_adjust(bottom=0.15)
					
			elif xvg_xlab1 != 'na':
				xlabel(xvg_xlab1,fontsize=label_size)
			
	
			if ylab1 != 'na':
				
				ylabel(ylab1,fontsize=label_size)
			elif xvg_ylab1 != 'na':
				ylabel(xvg_ylab1,fontsize=label_size)	
			
			if title1 != 'na':
				title(title1,fontsize=label_size)
			elif xvg_title1 != 'na' and opt_xvg_title == 1:
				title(xvg_title1,fontsize=label_size)
			
			if xlims != 'na':
				xlim(xlims[0],xlims[1])	
			#else:
				#x_lmin = min(x_plot)
				#x_lmax = max(x_plot)
				#x_lrange = x_lmax-x_lmin
				#xlim(x_lmin-0.1*x_lrange, x_lmax+0.1*x_lrange)
				
				
			if ylims != 'na':
				ylim(ylims[0],ylims[1])	
			#else:
				#y_lmin = min(y_plot)
				#y_lmax = max(y_plot)
				#y_lrange = y_lmax-y_lmin
				#ylim(y_lmin-0.1*y_lrange, y_lmax+0.1*y_lrange)

			# Get x/y lims
			xlim_in = xlim()
			ylim_in = ylim()

			if ax_square == 1:

				v2 = xlim_in
				v1 = ylim_in
				fac1 = (v2[-1]-v2[0])/(v1[-1]-v1[0])
				ax.set_aspect(fac1)


			# Plot horizontal or vertical line
			if len(xvert) > 0:
				for i in range(0,len(xvert)):
					
					ax.plot([xvert[i], xvert[i]], [ylim_in[0], ylim_in[1]],'k')
			if len(yhorz) > 0:
				for i in range(0,len(yhorz)):
					
					ax.plot([xlim_in[0], xlim_in[1]], [yhorz[i], yhorz[i]], 'k')		

			if num_size != 12:
				for tick in ax.xaxis.get_major_ticks():
                        		tick.label1.set_fontsize(num_size)			
				for tick in ax.yaxis.get_major_ticks():
                        		tick.label1.set_fontsize(num_size)			

			
			#LEGLOC
			#ll1 =   'ur'
			#ll1 = 'ul'
			#ll1 = 'll'
			#ll1 =   'lr'
			if leg_loc_in == 'na':
				ll1 =   'lr'
			else:
				ll1 = leg_loc_in	

			if ll1 == 'ur':
				leg_loc = 1
			elif ll1 == 'ul':
				leg_loc = 2
			elif ll1 == 'll':
				leg_loc = 3
			elif ll1 == 'lr':
				leg_loc = 4

			# KEY: 	1 2 3 4 ur ul ll lr
			# KEY: 	8 9 lc uc
			#leg_loc = 1 
			#leg_loc = 9
			if leg1 == 'fnames':
				legend(prop={'size':10},loc=leg_loc)

			elif leg1 != 'na':
				
				ncol_leg = 1 
				#ncol_leg = 3 

				if opt_points[0] == '0':
					foo = 1	

					if foo==1:

						lega = legend(loc=leg_loc,columnspacing=0.01,borderpad=0.1,ncol=ncol_leg,prop={'size':leg_size}, handletextpad=0.1 )
						lega.draw_frame(False)
					else:
						legend(leg1,loc=leg_loc,prop={'size':leg_size})	
				else:
				#if opt_points == 3:
					#legend(loc=leg_loc,ncol=ncol_leg,prop={'size':leg_size})
					#legend(loc=leg_loc,ncol=ncol_leg,prop={'size':leg_size}, columnspacing=0.3,borderpad=0.3, handletextpad=0.3 )
					lega = legend(loc=leg_loc,ncol=ncol_leg,prop={'size':leg_size}, numpoints=1, columnspacing=1,borderpad=0.3, handletextpad=0.1 )
					lega.draw_frame(False)

			if let_loc != 0 and let_loc != 'na':
				let_x = let_loc[0]*np.diff(xlim_in)[0]
				let_y = let_loc[1]*np.diff(ylim_in)[0]

				ax.text(let_x, let_y, let_add, fontsize=let_font)
						
			##HC-TEXT
			#ax.text(0.35,30, '0.0D$_s$')	
			#ax.text(0.5,60, '0.2D$_s$')	
			#ax.text(0.58,90, '0.4D$_s$')	
			#ax.text(0.65,110, '0.6D$_s$')	
			#ax.text(0.71,140, '0.8D$_s$')	
			#ax.text(0.75,170, '1.0D$_s$')	


			if fig_name1 != 'na':

				dpi_in = 160
				#savefig(fig_name1,pad_inches=0.25,dpi=dpi_in)
				savefig(fig_name1,bbox_inches='tight',pad_inches=0.1,dpi=dpi_in)
				print "saving..."
				print fig_name1
			else:
				
				if opt_plot == 1:
					#HC
					#ax.set_xlim([0.1,5])
					#ax.set_ylim([0,35])
					#ax.axes.get_xaxis().set_visible(False)
					#ax.axes.get_yaxis().set_visible(False)
					#ax.axis('off')

					show()
					print "not saving..."


	if opt_csv == 1:
		dw_csv.close()
		
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hf:x:y:o:i:p:j:I:",["fex=","fey=","nex=","ney=","xlim=","ylim=","dlim=","errs","edit","clus=","vert=","horz=","all","mean=","slope=","hist=","logx=","noplot","csv=","colors=","save_data=","Fy=","Ny=","Fx=","Nx=","tacf=","avg=","lleg","stop","custom=",'ll=',"rev"])
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

	dlims = 'na'
	xlims = 'na'
	ylims = 'na'
	
	leg1 = 'na'
	
	fig_name = 'na'
	opt_n_err = 10	# number of errorbars
	opt_edit = 0 
	opt_points = ['0']
	opt_stats =0; t_stats0 = 'na'
	opt_hist = [0]
	opt_save_data = 0
	opt_slog =0
	opt_plot = 1
	opt_csv = 0; x_csv = 'na'
	opt_inv = 0
	opt_tacf = 0
	opt_rev = 0   # reverse the order of file names
	leg_loc = 'na'

	f_avg = '0'
	y_fac = ['1']
	y_norm = ['1']
	x_fac = ['1']
	x_norm = ['1']

	label_size = 14
	num_size = 12
	leg_size = 14#12
	size_line = 1.5
	size_line = 2.5

	xvert = []
	yhorz = []
	colors = 'na'
	opt_stop = 0
	opt_stop = 'na'
	opt_custom = 'na'
	
	f_name_in = ''
	f_name_out = ''
	f_names_in = []
	ey_names_in = []
	ex_names_in = []
	i_x = 0
	i_y = ['na']
	i_ex = 0
	i_ey = 0
	
	for opt, arg in opts:
		
		if opt == '-h':
			print 'py_plot.py -f <inputfiles> -e <errorfiles/same> -o <outputfile.png>  -i <info_file> -l <legend_file>'
			sys.exit()
				
		if opt == '--all':	
			f_names_in = []
			for j in range(0,len(argv)):	# find xvg or txt file with wildcard
				i_ext = argv[j].find('xvg')
				if i_ext != -1:
					f_names_in.append(argv[j])
				else:
					i_ext = argv[j].find('txt')
					if i_ext != -1:
						f_names_in.append(argv[j])
			
			if len(f_names_in) == 0:
				print "No data files..."
				set_trace()
			
		if opt in ("-f"):
			f_name_ino = arg
			f_names_in.append(f_name_ino)
							
			
		if opt in ("-y"):	# Column indeces to plot
			if arg == 'qs':
				n_col_ino = [2]
			elif arg == 'cp':
				n_col_ino = [3]
			else:	
				n_col_ino = map(int, arg.split(',')) # new list
			for i0 in n_col_ino:
				i_y.append([i0])
		if opt in ("-x"):
			if arg == 'cp':
				i_x = 3
			elif arg == 'C':
				i_x = 0
			else:
				i_x = int(arg)	
		if opt in ("-o"):
			fig_name = arg
		if opt in ("-i"):
			f_name_info = arg

		if opt in ("-I"):
			opt_inv = 1
		if opt in ("-p"):
			opt_points = arg.split(',')
	
		if opt in ("--Fy"):		# scaling factor
			if os.path.isfile(arg):		# if input file
				dr = open(arg,'r')
				line = dr.readline()
				y_fac = line[0:len(line)-1]
				dr.close()
			else:
				y_fac = map(float, arg.split(','))	# string
		if opt in ("--Ny"):		# scaling factor
			if os.path.isfile(arg):		# if input file
				dr = open(arg,'r')
				line = dr.readline()
				y_norm = line[0:len(line)-1]
				dr.close()
			else:
				y_norm = map(float, arg.split(','))	# string
				
		
		if opt in ("--Fx"):		# scaling factor
			if os.path.isfile(arg):		# if input file
				dr = open(arg,'r')
				line = dr.readline()
				x_fac = line[0:len(line)-1]
				dr.close()
			else:
				x_fac = map(float, arg.split(','))	# string
		if opt in ("--Nx"):		# scaling factor
			if os.path.isfile(arg):		# if input file
				dr = open(arg,'r')
				line = dr.readline()
				x_norm = line[0:len(line)-1]
				dr.close()
			else:
				x_norm = map(float, arg.split(','))	# string
						
						
		if opt == "--fey":	# Error files
			ey_names_in = arg.split(',')
		if opt == "--fex":	# x Error files
			ex_names_in = arg.split(',')
		if opt == '--dlim':
			dlims = map(float, arg.split(','))
		if opt == '--xlim':
			xlims = map(float, arg.split(','))
		if opt == '--ylim':
			ylims = map(float, arg.split(','))
		
		if opt == '--nex':	# Error indecs
			if arg == 'cp':
				i_ex = 7
			elif arg == 'C':
				i_ex = 4
			else:
				i_ex = int(arg)	

		if opt == '--ney':	# Error indecs
			if arg == 'cp':
				i_ey = [7]
			elif arg == 'qs':
				i_ey = [6]
			else:
				i_ey = map(int, arg.split(','))

		if opt == '--errs':	# plot all error bars
			opt_n_err = 'all'
		if opt == '--edit':
			opt_edit = 1
		if opt == '--lleg':
			leg1 = 'fnames'
		if opt == '--rev':
			opt_rev = 1
		if opt == '--tacf':
			opt_stats = 'tacf'
			t_stats0 = arg
			#opt_tacf = 1
		if opt == '--vert':
			xvert = map(float, arg.split(','))	# list
		if opt == '--horz':
			yhorz = map(float, arg.split(','))	# list
		if opt == '--colors':
			colors = arg.split(',')	# list
		
		if opt == '--hist':
			opt_hist = arg.split(',')
			opt_hist = [1]+opt_hist
		if opt == '--mean':
			opt_stats = 'mean'
			t_stats0 = arg
		if opt == '--slope':
			opt_stats = 'slope'
			t_stats0 = arg
		if opt == '--avg':
			f_avg = arg			
		if opt == '--ll':
			leg_loc = arg
		if opt == '--logx':
			opt_slog = arg
			try:
				opt_slog = int(opt_slog)
			except:
				foo = 0
		if opt == '--stop':
			opt_stop = 1
		if opt == '--custom':
			opt_custom = arg
		if opt == '--noplot':
			opt_plot = 0
		if opt == '--save_data':
			opt_save_data = arg
			
		if opt == '--csv':
			opt_csv = 1	
			# input x point
			if os.path.isfile(arg):
				x_csv = []
				dr = open(arg,'r')
				for line in dr:
					dl = line.split()
					x_csv.append(dl[0])
				dr.close()
			else:
				
				x_csv = arg.split(',')
				
	for opt, arg in opts:	# catch for journal option error
		if opt in ("-j"):
			if fig_name == 'na':
				print "Specify figure name..."
				set_trace()		
						
	# PLOT INFO
	opt_xvg_title = 1
	
	if f_name_info != 'na':
		go_ylab = 0; go_title = 0
		dr = open(f_name_info,'r')
		for line in dr:
			dl = line.split()
			if len(dl) > 0 and dl[0] == 'xlabel':
				i1 = line.find('"');i2 = line.rfind('"')
				xlab1 = line[i1+1:i2]
			if len(dl) > 0 and dl[0] == 'ylabel':
				i1 = line.find('"');i2 = line.rfind('"')
				ylab1 = line[i1+1:i2]
				go_ylab = 1
			if len(dl) > 0 and dl[0] == 'title':
				i1 = line.find('"');i2 = line.rfind('"')
				title1 = line[i1+1:i2]
				go_title = 1
			if len(dl) > 0 and dl[0] == 'xlims':
				xlims = [float(dl[1]), float(dl[2])]
				print "xlims: "+str(xlims[0])+' '+str(xlims[1])
			if len(dl) > 0 and dl[0] == 'ylims':
				ylims = [float(dl[1]), float(dl[2])]
				print "ylims: "+str(ylims[0])+' '+str(ylims[1])	
			if len(dl) > 0 and dl[0] == 'fontsize':
				label_size = float(dl[1])
			if len(dl) > 0 and dl[0] == 'legsize':
				leg_size = float(dl[1])
			if len(dl) > 0 and dl[0] == 'numsize':
				num_size = float(dl[1])
			if len(dl) > 0 and dl[0] == 'linesize':
				size_line = float(dl[1])
			if len(dl) > 0 and dl[0] == 'leg':
				if leg1 == 'na':
					leg1 = []
				leg1.append(line[len('leg'):-1])
				
		dr.close()
		if go_ylab == 1 and go_title == 0:
			opt_xvg_title = 0
	
	if len(i_y) > 1:
		del i_y[0]
	else:	
		i_y = len(f_names_in)*[[1]]

	if len(ey_names_in) > 0:
		if ey_names_in[0] == 'same':
			ey_names_in = []
			for n in range(0,len(f_names_in)):
				f_name0 = f_names_in[n]
				ey_names_in.append(f_name0)
		elif ey_names_in[0] == 'std' or ey_names_in[0] == 'std2err':
			for n in range(0,len(f_names_in)):
				f_name0 = f_names_in[n]
				e_name0 = f_name0.replace('mn','std')
				ey_names_in.append(e_name0)
		elif ey_names_in[0] == 'err':
			ey_names_in = []
			for n in range(0,len(f_names_in)):
				f_name0 = f_names_in[n]
				e_name0 = f_name0.replace('mn','err')
				ey_names_in.append(e_name0)

	if len(ex_names_in) > 0:
		if ex_names_in[0] == 'same':
			ex_names_in = []
			for n in range(0,len(f_names_in)):
				f_name0 = f_names_in[n]
				ex_names_in.append(f_name0)

	
	str_cmd = ''
	for i in range(0,len(argv)):
		str_cmd = str_cmd+argv[i]+' '

	# call function	
	py_plot1(f_names_in, ey_names_in, ex_names_in, i_x, i_y, i_ex, i_ey, fig_name, xlab1, ylab1, label_size, leg_size, num_size, size_line, title1, leg1, xlims, ylims, dlims, opt_n_err, opt_edit, opt_points, opt_stats, t_stats0, opt_save_data, y_fac, y_norm, x_fac, x_norm, opt_slog, opt_plot, opt_csv, x_csv, opt_xvg_title, opt_inv, xvert, yhorz, colors, str_cmd, opt_tacf, opt_hist, f_avg, opt_stop, opt_custom, leg_loc, opt_rev)
	
	# save off command
	cmd_out = 'n'
	
	if fig_name != 'na':
		cmd_out = 'y'
	#else:
		#cmd_out = raw_input("Save desc? (y/n)")	 
	
	if cmd_out == 'y':
		
		fig_name0 = fig_name.split('/')
		fig_name0[-1] = 'desc_'+fig_name0[-1][0:-4]
		f_name_cmd = '/'.join(fig_name0)
		
		dw = open(f_name_cmd,'w')
		dw.write(str_cmd+'\n')
		dw.close()
	
	# copy files over to journal
	for opt, arg in opts:
		if opt in ("-j"):

			#jour_opts = arg.split(',')
			#if jour_opts[0] == 'def':
				#p_name_jour = '/home/jb2qm/jour/'
			#else:
				#p_name_jour = jour_opts[0]
			#n_dir_jour = jour_opts[1]
			p_name_jour = '/home/jb2qm/jour/'
			n_dir_jour = arg
			
			n_last = -1	# find number so far
			dr = open(p_name_jour+'title_list.txt','r')
			if n_dir_jour == 'n':	# just next available number
				for line in dr:
					dl = line.split()
					if len(dl) > 0:
						n_last = int(dl[0]) 
				n_dir_jour = n_last+1
			else:
				n_dir_jour = int(n_dir_jour)
			dr.close()
			
			p_name_plot = p_name_jour+'plot'+str(n_dir_jour)
			if os.path.lexists(p_name_plot):
				print("Re-writing data in "+p_name_plot+' ???')
				set_trace()
				foo = 0
			else:
				os.system('mkdir '+p_name_plot)

			dw = open(p_name_jour+'tmp_list.txt','w')
			dr = open(p_name_jour+'title_list.txt','r')
			if n_dir_jour > n_last:	# just next available number
				for line in dr:
					dw.write(line)
				dw.write(str(n_dir_jour)+'   '+title1+'\n')	
			else:
				for line in dr:
					dl = line.split()
					if dl[0] == int(n_dir_jour):
						dw.write(str(n_dir_jour)+'   '+title1+'\n')
			dw.close()
			dr.close()
			os.system('mv '+p_name_jour+'tmp_list.txt '+p_name_jour+'title_list.txt')
			
			for j in range(0,len(f_names_in)):
				os.system('cp '+f_names_in[j]+' '+p_name_plot)
			for j in range(0,len(e_names_in)):
				os.system('cp '+f_names_in[j]+' '+p_name_plot)		
			if f_name_info != 'na':
				os.system('cp '+f_name_info+' '+p_name_plot)		
			if f_name_leg != 'na':
				os.system('cp '+f_name_leg+' '+p_name_plot)	
			if fig_name != 'na':
				os.system('cp '+fig_name+' '+p_name_plot)
				os.system('cp '+fig_name+' '+p_name_jour+'fig'+str(n_dir_jour)+'.png')
	#print "...exiting"

if __name__ == "__main__":
   main(sys.argv[1:])
