#MODULE general_plot.py
from pdb import *
import os
import sys
sys.path.append("/mnt/hds/funcs/")
#import general_py
from matplotlib.pyplot import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import random
import scipy.signal

#from matplotlib import scale as mscale
#from matplotlib import transforms as mtransforms
#from matplotlib import rc
matplotlib.rcParams.update({'font.size': 14})


#plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def add_text(ax, x_in, y_in, let_in, font_in):
	xl1 = xlim() 
	yl1 = ylim() 
	xlet = xl1[1]*x_in
	ylet = yl1[1]*y_in

	print "Label loc (x,y): "
	print xlet
	print ylet
	ax.text(xlet,ylet,let_in, fontsize=20)

	return ax

def plot_pyplot(f_name, ax, x_plot, y_plot, x_err, y_err, x_xerr, y_xerr, coln, symn, mfcn, lsn, lwn, opt_points, opt_n_err, opt_hist, opt_slog, legn, opt_stop):
	#import general_py

	#opt_points: 0 (points/fit_opt), 1/2 (fit_bnds)
	opt_fit = 0
	try:
		int(opt_points[0])
	except:
		opt_fit = 1

	if opt_hist[0] == 0:
		#if opt_points[0] == 'fit':
		if opt_fit == 1:
			# FITTING DATA
			#import general_py

			lsn2 = lsn
			lsn = 'None'
			n_fit_pts = 1000
			
			# Find fitting bounds	
			try:
				fit_bnds = [float(opt_points[1]), float(opt_points[2])]

			except:
				try:
					fit_bnds = [float(opt_points[2]), float(opt_points[3])]
				except:
					print "Prob with fit_bnds"
					set_trace()

			if len(x_plot) > 1:
				print "######### FITTING: "+f_name

				# CHECK FOR EXTRA PARAMS (for iso fitting)
				ex_params = find_ex_params(opt_points, f_name)

				if opt_points[1] == 'pore':
					opt_hf_fit = 0 
					
					fit_lim = 0
					fit_lim = 0.99
					#fit_lim = 1.0 
					if fit_lim == 0:
						opt_hf_fit = 1					
		

					if opt_hf_fit == 1:
		        	                i_ft_st = len(x_plot)-1
						ex_params = ['opt_hf_fit',1]	
		        	        else:
		        	                i_ft_st = np.where((y_plot>=fit_lim*max(y_plot)))[0][0]
						ex_params = ['opt_hf_fit',0]	

					fit_bnds = [0, x_plot[i_ft_st]]

				x_fit, y_fit, p_fit = general_py.fits(opt_points, fit_bnds, x_plot, y_plot, n_fit_pts, ex_params)

				##HC write out slope
				#f_outo = f_name+'lin'
				#dw = open(f_outo,'w')
				#dw.write(str(p_fit[0])+'\n')
				#dw.close()

				# OPTION TO CHANGE LEGEND
				#if n < 2:
				#	qm = y_plot[-1]
				#else:
				#	qm = p_fit[0]
				#qm = p_fit[0]
				#legn = legn+', $q_m$= '+"%3.0f"%(qm)


		mke = 'None'	# default, mark every point

		if opt_points[0] == '0':	# just lines, no points
			symn = ' '
		elif opt_points[0] == '1':	# just points
			lsn = 'None'

		elif opt_points[0] == '3':	# stagger points
			mke = 'None'
			if opt_slog == 0:
				fr0 = 400
				i00 = random.randint(0,fr0)
			elif opt_slog == 'x':# or opt_slog == 'y':
				t_dens = [[0,10], [10,100],[100,1000]]
				fr0s = [10, 100, 400]
				for nt in range(0,len(t_dens)):
					fr0 = fr0s[nt]
					i_kp = np.where((x_plot >= t_dens[nt][0]) & (x_plot <= t_dens[nt][1]+1))
					x_kp = x_plot[i_kp]

					i00 = random.randint(0,fr0)
					mark_tup = [i00,fr0]
					mark_range = np.arange(i00, i_kp[0][-1], fr0)
					mark_range = mark_range[mark_range < len(i_kp[0])]
					mark_x = x_kp[mark_range]
					for ni in range(0,len(mark_range)):
						i0s[n].append(mark_x[ni])

					if nt == 0:
						ax.semilogx(x_plot[i_kp], y_plot[i_kp], color=coln,marker=symn,markevery=mark_tup,label=leg1[n],markersize=8,linestyle=ls1,linewidth=size_lines[n])
					else:
						ax.semilogx(x_plot[i_kp], y_plot[i_kp], color=coln,marker=symn,markevery=mark_tup,markersize=8,linestyle=ls1,linewidth=size_lines[n])


		if opt_slog == 0:
			#MARKERSIZE 
			msz = 8 
			msz = 12 
			#msz = 4 
			mew = 2
			#mew = 1.5
			mew = 1 
			mke = [8,30]   # mark every 
			mke = 'None'   # mark every 
			#mke = 2   # mark every 
			mec = coln
			mec = 'w' 

			if opt_fit == 1 and len(x_plot) > 1:
				ax.plot(x_fit,y_fit, color=coln, linestyle=lsn2, linewidth=lwn)

			if mke == 'None':
				if type(symn) != unicode:
					######MAIN PLOTTING

					ax.plot(x_plot, y_plot, color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=mec, label=legn, markersize=msz, markeredgewidth=mew)#, markevery=mke)
					if opt_stop == 1:
						x_stop = 2*[x_plot[-1]]
						y_stop = [min(y_plot), 1.25*max(y_plot)]
						y_stop = [min(y_plot), 8]
						ax.plot(x_stop,y_stop,color=coln,linestyle='--',linewidth=2)

				else:
					ax.plot(x_plot, y_plot, linestyle=lsn, linewidth=lwn, label=legn)
					for x0, y0 in zip(x_plot, y_plot):	
						text(x0, y0, symn, size=30,fontname='STIXGeneral', va='center', ha='center', clip_on=True)						
				
#, fontname='STIXGeneral')
			elif type(mke) == list:  # different marker frequencies
				x_bk = 0.95
				i1 = np.where((y_plot<x_bk*max(y_plot)))
				i2 = np.where((y_plot>=x_bk*max(y_plot)))

				ax.plot(x_plot[i1], y_plot[i1], color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=coln, label=legn, markevery=mke[0],markersize=msz,markeredgewidth=mew)
				ax.plot(x_plot[i2], y_plot[i2], color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=coln, markevery=mke[1],markersize=msz,markeredgewidth=mew)
			else:
				ax.plot(x_plot, y_plot, color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=coln, label=legn, markevery=mke,markersize=msz,markeredgewidth=mew)




		elif opt_slog == 'x':	
			ax.semilogx(x_plot, y_plot, color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=coln, label=legn)#, markevery=mke)
			#if opt_points[0] == 'fit':
			if opt_fit == 1:
				ax.semilogx(x_fit,y_fit, color=coln, linestyle=lsn2, linewidth=lwn)
		elif opt_slog == 'y':	
			ax.semilogy(x_plot, y_plot, color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=coln, label=legn)#, markevery=mke)
			#if opt_points[0] == 'fit':
			if opt_fit == 1:
				ax.semilogy(x_fit,y_fit, color=coln, linestyle=lsn2, linewidth=lwn)
		elif opt_slog == 'xy':	
			ax.loglog(x_plot, y_plot, color=coln, linestyle=lsn, linewidth=lwn, marker=symn, markerfacecolor=mfcn, markeredgecolor=coln, label=legn)#, markevery=mke)
			#if opt_points[0] == 'fit':
			if opt_fit == 1:
				ax.loglog(x_fit,y_fit, color=coln, linestyle=lsn2, linewidth=lwn)
		
		if x_err != 'na' or x_xerr != 'na':
			if opt_n_err != 'all':
				if opt_points != '3':

					i0 = general_py.int_linspace_x(x_plot, opt_n_err)
					#i0 = map(int, (np.round(np.linspace(1, len(x_plot)-1, opt_n_err))) )
					if len(x_plot) <= opt_n_err:
						i0 = range(0, len(x_plot))
				else:	# staggered
					foo = 0
		#               else:
        	#                       i0 = []
        	#                       i0n = i0s[n]
        	#                       for ni in range(0,len(i0n)):
        	#                               i00 = np.where(x0 == i0n[ni])
        	#                               i0.append(i00[0][0])
        	#       x = x0[i0]; y = y0[i0]

			else:	# all errorbars
				i0 = range(0, len(x_plot))


			# PLOT ERROR

		elw = 2
		if x_err != 'na' and x_xerr != 'na':

			ax.errorbar(x_plot[i0], y_plot[i0], yerr=y_err[i0], xerr=y_xerr[i0], fmt=None, ecolor=coln, elinewidth=elw)
		elif x_err != 'na':
			ax.errorbar(x_plot[i0], y_plot[i0], yerr=y_err[i0], fmt=None, ecolor=coln, elinewidth=elw)
		elif x_xerr != 'na':
			ax.errorbar(x_plot[i0], y_plot[i0], xerr=y_xerr[i0], fmt=None, ecolor=coln, elinewidth=elw)

	elif opt_hist[0] == 1:
		import general_py
		n_bins = int(opt_hist[1])
		opt_norm = 0
		if len(opt_hist) > 2 and opt_hist[2] == 'norm':
			opt_norm = 1
		b1, hist1 = general_py.quick_hist(y_plot, n_bins, opt_norm)
		if opt_points[0] == '0':
			#ax.plot(b1, hist1)#, coln,linewidth=lsn)
			ax.plot(b1, hist1, label=legn)#, coln,linewidth=lsn)
			#ax.plot(b1, hist1, coln,linewidth=lsn)
		elif opt_points[0] == '1':
			ax.plot(b1, hist1,'.'+coln)
		elif opt_points[0] == '2':
			ax.plot(b1, hist1,'.',x,y_plot,coln)
		

	# DONE
	return ax

		
def sub_plot_hist(fig_name, d_in, de_in, y_lims, fig_rows, fig_cols, title1, x_lab1, x_lab2, y_lab1, leg_in, colors_in, hatch_in, opt_leg, size_num, size_lab, size_leg, ax_ticks):
	#TYPE: plotting

	# optional ticker edits
	n_num1 = 1
	n_num2 = 1
	if ax_ticks != 'na':
		opt_ax_edit = 1

	majorLocator = MultipleLocator(ax_ticks[0])
	minorLocator = MultipleLocator(ax_ticks[1])
	majorFormatter = FormatStrFormatter('%'+str(n_num1)+'.'+str(n_num2)+'f')	# EDIT 

	dim1 = np.shape(d_in)
	n_x = dim1[0]
	n_y = dim1[1]
	n_figs = fig_rows*fig_cols

	fig = figure()
	for n in range(0,n_figs):
		str1 = str(fig_rows)+str(fig_cols)+str(n+1)

		d1 = d_in[n,:]
		p1 = subplot(str1)
		ind0 = range(0,len(d1))
		ind1 = np.array(ind0)+0.1
		wid1 = 1.0

		if de_in == 'na':	# NO ERROR BARS
			for i in range(0, n_y):
				if leg_in != 'na':
					legi = leg_in[i]
				else:
					legi = ''
				colori = colors_in[i]
				hatchi = hatch_in[i]
				if hatchi == '':
					p1.bar(ind1[i],d1[i], label=legi, color=colori)
				else:
					p1.bar(ind1[1],d1[i], label=legi, color=colori, hatch=hatchi)

		else: 			# ERROR BARS
			ecol1 = 'red'
			ecap1 = 4

			de1 = de_in[n,:]
			for i in range(0, n_y):
				if leg_in != 'na':
					legi = leg_in[i]
				else:
					legi = ''
				colori = colors_in[i]
				hatchi = hatch_in[i]

				if hatchi == '':
					p1.bar(ind1[i], d1[i], yerr=de1[i], label=legi, color=colori, ecolor=ecol1, capsize=ecap1)
				else:
					p1.bar(ind1[i], d1[i], yerr=de1[i], label=legi, color=colori, hatch=hatchi, ecolor=ecol1, capsize=ecap1)
		
		xticklabels = p1.get_xticklabels()
		#yticklabels = p1.get_yticklabels()
		#if n == 0 or n == 3:
		#	foo = 0
		#else:
		#	setp(yticklabels,visible=False)


		for tick in p1.yaxis.get_major_ticks():
			tick.label1.set_fontsize(size_num)
			#set_trace()
	
		## OPTION1	
		#setp(xticklabels,visible=False)
		#xlabel(x_lab1[n],fontsize=size_lab)#, **font)

		# OPTION2	
		tk_range = np.array([2.,6.])
		#tk_range = np.array([1.5,4.5])
		tk_labs = x_lab1
		xticks( tk_range, tk_labs, fontsize=26)

		if opt_ax_edit == 1:
			p1.yaxis.set_major_locator(majorLocator)
			p1.yaxis.set_major_formatter(majorFormatter)
			p1.yaxis.set_minor_locator(minorLocator)

		ylim(y_lims[0],y_lims[1])
	
		# LEGEND PLACEMENT	
		if n == 0:
			if y_lab1 != 'na':
				#ylabel('r'+"'"+y_lab1+"'")
				#ylabel('r'+'"'+y_lab1+'"')
				n_sp_ylab = 35 # waterD
				n_sp_ylab = 0 # water tau

				ylab_in = y_lab1+n_sp_ylab*' '
				ylabel(ylab_in, position=(0.5,0.5),fontsize=size_lab)
				ylabel(ylab_in, fontsize=size_lab)
		#if n == 1 or n==0:
		if n == 0:
			if opt_leg != 0: # Lang. Dyn
				#p1.legend(loc=9,prop={'size':size_leg},bbox_to_anchor=(0.4,1.20),ncol = 7,labelspacing=0.9) #, **font0)
				p1.legend(loc=9,prop={'size':size_leg},bbox_to_anchor=(0.5,1.1),ncol = 7,labelspacing=0.9) #, **font0)


	subplots_adjust(hspace=0.45)
	#subplots_adjust(hspace=0.45, bottom=0.12,right=0.85)
	if fig_name == 'na':
	
		show()
		set_trace()	
		foo = 1
	else:
		#f_out1 = fig_name+'.png'
		f_out1 = fig_name+'.eps'
		f_out2 = fig_name+'.png'
		print "Saving... "+f_out1
		dpi_in = 160
		#savefig(f_out1,bbox_inches='tight',pad_inches=0.15,dpi=dpi_in)
		savefig(f_out1,bbox_inches='tight',pad_inches=0.02)
		savefig(f_out2,bbox_inches='tight',pad_inches=0.02)

		#savefig(f_out2,bbox_inches='tight',pad_inches=0.15)
		close()	
	print "done..."
	

def comp_dists(fig_name_ext, opt_f_ext, x_ref, y_ref, f_names_in, t_lims,  n_bins, col_in, sym_in, leg_in, xlab, ylab, x_lim, y_lim,dpi_in):

	if x_ref != 'na':
		fig = figure()
		ax = fig.add_subplot(111)
		ax.plot(x_ref, y_ref, sym_in[-1], color=str(col_in[-1]), label=leg_in[-1]); hold
		foo = 0	

	bs = []
	hists = []
	n_sims = 15
	p_sim0 = 10
	#n_sims = 1
	#p_sim0 = 1
	ext2 = '.xvg'
	
	labs1 = ['NH     ','mol-NH  ','atom-NH']

	for n in range(0, len(f_names_in)):
		f_name_in0 = f_names_in[n]
		if opt_f_ext == 0:
			f_name_ns = [f_name_in0]
		else:
			f_name_ns = []
			for m in range(0,n_sims):
				f_name_m =  f_name_in0+str(p_sim0+m)+ext2 
				f_name_ns.append( f_name_m )

		x_in = []
		y_in = []
		for m in range(0,len(f_name_ns)):
			f_name_in = f_name_ns[m]
			print f_name_in		
			x0, y0 = rf.gen_read2(f_name_in)

			i_kp = np.where((x0 > t_lims[0]) & (x0 < t_lims[1]))
	
			y1 = y0[i_kp]

			### FFT (optional)
			#t1 = 300.
			#t2 = t_lims[1]#699.92
			#dt = 0.004
			#Fs = 250	
			#fq0, pw0 = general_py.quick_fft(y1, Fs, t1, t2, dt)

			#fq1 = fq0[1:]
			#pw1 = pw0[1:]

			#n_filt = 25
			#pw2 = scipy.signal.medfilt(pw1,51)
			#mpw2 = max(pw2)
			#imax = np.where((pw2 >= 0.95*mpw2))
			#mfqs = fq1[imax]
			#mtaus = 1./mfqs
			#mxtau = mtaus[0]
			#mxtau = np.median(mtaus)
			#print mfqs	
			#print mtaus	
	
			##plot(fq1, pw1)
			#lab1 = labs1[n]+'  '+"%1.3f"%(mxtau)+'ps'
			#plot(fq1, pw2,label=lab1)
			#
			#hold
	
			##set_trace()

			y_mn = np.mean(y1)
			y_std = np.std(y1)
			print "Mean: "+str(y_mn)
			print "Std:  "+str(y_std)

			for i in range(0,len(y1)):
				y_in.append(y1[i])
			#print len(y_in)

		print "**********Mean: "+str(np.mean(y_in))
		print "**********Err: "+str(np.std(y_in)/np.sqrt(len(y_in)))
		b1, hist = general_py.quick_hist(y_in, n_bins)
		bs.append(b1)
		hists.append(hist)
		
		ax.plot(b1, hist, sym_in[n], color=str(col_in[n]), label=leg_in[n], markersize=8)
	#xlabel('Freq (1/ps)')
	#ylabel('Power')
	#legend()

	#show()
	#set_trace()

	num_size=16
	lab_size=20
	leg_size=16

	#legend(prop={'size':leg_size})
	xlabel(xlab,fontsize=lab_size)
	ylabel(ylab,fontsize=lab_size)
	ylim(y_lim[0], y_lim[1])
	xlim(x_lim[0], x_lim[1])
	for tick in ax.xaxis.get_major_ticks():
        	tick.label1.set_fontsize(num_size)			
	for tick in ax.yaxis.get_major_ticks():
        	tick.label1.set_fontsize(num_size)			

	if fig_name_ext != 'na':	
		savefig(fig_name_ext+'.eps',bbox_inches='tight',pad_inches=0.02)
		#savefig(fig_name_ext+'.png',bbox_inches='tight',pad_inches=0.15,dpi=dpi_in)
		#savefig(fig_name_ext+'.pdf',bbox_inches='tight',pad_inches=0.15)
	else:
		show()

	
def find_ex_params(opt_points, f_name):
	ex_params = []

	if opt_points[1] == 'custom1':  # for qhat with linear 
		print "CUSTOM1"
	
		ids = ['surf:','poly_lang:','poly_lin0:']

		for j in range(0,len(ids)):
			idj = ids[j]	
	
			dr = open(f_name, 'r')
			for line in dr:
				dl = line.split()
				if dl[0] == idj:
					ex_params.append(dl)
			dr.close()

		ex_params.append( f_name )

		#dr = open(f_name, 'r')
		#for line in dr:
		#	dl = line.split()
		#	if dl[0] == 'surf:' or dl[0] == 'poly_lang:' or dl[0] == 'poly_lin0:':
		#		ex_params.append(dl)
		#dr.close()
		##ex_params.append( 0.9 )
	
		##if f_name.find('osurf') >= 0:
		##	ex_params.append('osurf')	
		##else:
		##	ex_params.append('osys')	
		#ex_params.append( f_name )

	return ex_params
