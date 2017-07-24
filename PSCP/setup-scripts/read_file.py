# read_file.py
# read_file module
import string
import numpy as np 
import scipy
from pdb import *



def gen_read1(f_name_in):
		
	dr = open(f_name_in,'r')
	c1=0
	x = np.zeros(1000000)
		
	for line in dr:
		dl = line.split()	
		try:
			x[c1] = float(dl[0])
			c1 = c1+1
			break
		except:
			foo = 1	
	for line in dr:
		dl = line.split()
		x[c1] = float(dl[0])
		c1 = c1+1

	x = x[0:c1]
	return (x)

def gen_read1l(f_name_in):
		
	dr = open(f_name_in,'r')
	c1=0
	x = np.zeros(10000000)
		
	for line in dr:
		dl = line.split()	
		try:
			x[c1] = float(dl[0])
			c1 = c1+1
			break
		except:
			foo = 1	
	for line in dr:
		dl = line.split()
		x[c1] = float(dl[0])
		c1 = c1+1

	x = x[0:c1]
	return (x)

def gen_read2(f_name_in):
		
	dr = open(f_name_in,'r')
	c1=0
	n_pts = 3000000
	n_pts = 1000000
	x = np.zeros(n_pts)
	y = np.zeros(n_pts)

	for line in dr:
		dl = line.split()	
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
	return (x,y)

def gen_read_any_num(f_name_in):
		
	dr = open(f_name_in,'r')
	c1=0
	n_pts = 3000000
	#n_pts = 25000

	for line in dr:
		dl = line.split()
		try:

			float(dl[0])
			float(dl[1])
			x = np.zeros(n_pts)
			y = np.zeros((n_pts,len(dl)-1))

			x[c1] = float(dl[0])
			len_dl = len(dl)
			for j in range(1,len_dl):
				y[c1,j-1] = float(dl[j])			
			c1 = c1+1
			break
		except:
			foo = 1	
	for line in dr:
		dl = line.split()
		try:
			x[c1] = float(dl[0])
		except:
			if type(dl[0]) == str:
				break

		for j in range(1,len_dl):
			y[c1,j-1] = float(dl[j])			
		c1 = c1+1

	try:
		x = x[0:c1]
		y = y[0:c1,:]
	except:
		print f_name_in
		set_trace()
	
	return (x,y)


def gen_read3(f_name_in):
		
	dr = open(f_name_in,'r')
	c1=0
	x = np.zeros(1000000)
	y = np.zeros(1000000)
	z = np.zeros(1000000)

	for line in dr:
		dl = line.split()	
		try:
			x[c1] = float(dl[0])
			y[c1] = float(dl[1])
			z[c1] = float(dl[2])
			c1 = c1+1
			break
		except:
			foo = 1	
	for line in dr:
		dl = line.split()
		x[c1] = float(dl[0])
		y[c1] = float(dl[1])
		z[c1] = float(dl[2])
		c1 = c1+1

	x = x[0:c1]
	y = y[0:c1]
	z = z[0:c1]
	return (x,y,z)

def gen_read4(f_name_in,opt_out):
		
	dr = open(f_name_in,'r')
	c1=0
	t = np.zeros(1000000)
	x = np.zeros(1000000)
	y = np.zeros(1000000)
	z = np.zeros(1000000)

	for line in dr:
		dl = line.split()	
		try:
			t[c1] = float(dl[0])
			x[c1] = float(dl[1])
			y[c1] = float(dl[2])
			z[c1] = float(dl[3])
			c1 = c1+1
			break
		except:
			foo = 1	
	for line in dr:
		dl = line.split()
		t[c1] = float(dl[0])
		x[c1] = float(dl[1])
		y[c1] = float(dl[2])
		z[c1] = float(dl[3])
		c1 = c1+1

	t = t[0:c1]
	x = x[0:c1]
	y = y[0:c1]
	z = z[0:c1]
	if opt_out == 1:   # to make output y a single array
		x = np.array([x,y,z])
		#return (t,x,y,z)
		return (t,x)
	else:
		return (t,x,y,z)
	

def gen_read_y2(f_name_in):
		
	dr = open(f_name_in,'r')
	c1=0
	x = np.zeros(1000000)
	y = np.zeros((1000000,2))

	for line in dr:
		dl = line.split()	
		try:
			x[c1] = float(dl[0])
			y[c1,0] = float(dl[1])
			y[c1,1] = float(dl[2])
			c1 = c1+1
			break
		except:
			foo = 1	
	for line in dr:
		dl = line.split()
		x[c1] = float(dl[0])
		y[c1,0] = float(dl[1])
		y[c1,1] = float(dl[2])
		c1 = c1+1

	x = x[0:c1]
	y = y[0:c1,:]
	return (x,y)

def gen_read_cat(f_names_in):

	n_fs = len(f_names_in)
		
	dr = open(f_name_in,'r')
	c1=0
	x = np.zeros(1000000)
	y = np.zeros(1000000)

	for line in dr:
		dl = line.split()	
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
	return (x,y)


def gen_readn(f_name_in, i_in):
		
	dr = open(f_name_in,'r')
	c1=0
	y = np.zeros((1000000, len(i_in) ) )

	for line in dr:
		dl = line.split()	
		try:
			for i in range(0,len(i_in)):
				y[c1,i] = float( dl[i_in[i]])	
			c1 = c1+1
			break
		except:
			foo = 1
	for line in dr:
		dl = line.split()	
		for i in range(0,len(i_in)):
			y[c1,i] = float( dl[i_in[i]])	
		c1 = c1+1

	y = y[0:c1,:]

	return (y)

def gen_read_any_rowst(f_name_in, opt_float, row_st):
	# reads numbers or text screens
	print f_name_in	

	dr = open(f_name_in,'r')
	c1=0
	x = []
		
	for line in dr:
		if c1 >= row_st: 

			dl = line.split()	
			if c1 == row_st:
				for i in range(0,len(dl)):
					x.append([])
			
			for i in range(0,len(dl)):
				if opt_float == 1:
					try:
						x[i].append(float(dl[i]))
					except:
						print line
						set_trace()
				else:
					x[i].append(dl[i])

		c1 = c1+1

	return (x)

def gen_read_any(f_name_in, opt_float):
	# reads numbers or text screens
	print "reading... "+f_name_in	

	dr = open(f_name_in,'r')
	c1=0
	x = []
		
	for line in dr:
		dl = line.split()	
		if c1 == 0:
			for i in range(0,len(dl)):
				x.append([])
		
		for i in range(0,len(dl)):
			if opt_float == 1:
				try:
					x[i].append(float(dl[i]))
				except:
					print line
					set_trace()
			else:
				x[i].append(dl[i])

		c1 = c1+1

	return (x)
	
def gen_read_col(f_name_in,i_col):
	print f_name_in
	print i_col	
	dr = open(f_name_in,'r')
	c1=0
	x = np.zeros(1000000)
	y = np.zeros(1000000)
		
	for line in dr:
		dl = line.split()	
		try:	
			x[c1] = float(dl[0])
			if i_col == 'last':
				i_col1 = -1
			elif i_col == 'slast':
				i_col1 = -2
			else:
				i_col1 = i_col
			y[c1] = float(dl[i_col1])
			c1 = c1+1
			break
		except:
			foo = 1	

	for line in dr:
		dl = line.split()
		try:
			x[c1] = float(dl[0])
			if i_col == 'last':
				i_col1 = -1
			elif i_col == 'slast':
				i_col1 = -2
			else:
				i_col1 = i_col
			y[c1] = float(dl[i_col1])
			c1 = c1+1
		except:
			if type(dl[0]) == str:
				break
			else:

				x[c1] = float(dl[0])
				y[c1] = float(dl[i_col+1])
				c1 = c1+1

	x = x[0:c1]
	y = y[0:c1]
	return (x,y)


def gen_read_field_id(f_name_in, ids, inds):
	#ids:  len(N) list of len(2) lists, each with [0] col_id of id_str, [1] id_str to check for
	#inds: len(N) list with each col_id to get the value for
	opt_try_float = 1 

	dr = open(f_name_in,'r')
	out = []

	for line in dr:
		dl = line.split()
		for i in range(0,len(ids)):
			id1 = ids[i]
			if len(dl) > id1[0]:
				
				for j in range(1,len(id1)):
					if dl[id1[0]] == id1[j]:
						try:
							if opt_try_float == 1:
								out1 =  float(dl[inds[i]])
							else:
								out1 =  dl[inds[i]]
						except:
							out1 =  dl[inds[i]]

						out.append( out1)

		if len(out) == len(ids):
			break
	return out

	#for line in dr:
	#	dl = line.split()
	#	if len(dl) > 0:
	#		for i in range(0,len(ids)):
	#			id1 = ids[i]
	#			for j in range(1,len(id1)):
	#				if dl[id1[0]] == id1[j]:
	#					try:
	#						out1 =  float(dl[inds[i]])
	#					except:
	#						out1 =  dl[inds[i]]

	#					out.append( out1)

	#	if len(out) == len(ids):
	#		break
	#return out
			



def gen_read_col_info(f_name_in,i_cols):
		
	dr = open(f_name_in,'r')
	c1=0
	n_pts = 3000000
	n_pts = 1000000
	x = np.zeros(n_pts)
	y = np.zeros(n_pts)
	
	xvg_title1 = 'na'; xvg_xlab1 = 'na'; xvg_ylab1 = 'na'
	
	for line in dr:		# find first data point
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
					
		if len(dl) >= 2:			
			
			try:
				x0 = float(dl[0])
				x[c1] = float(dl[i_cols[0]])
				y[c1] = float(dl[i_cols[1]])
				c1 = c1+1
				break		
			except:
				try: 
					x0 = float(dl[0])
					x[c1] = float(dl[i_cols[0]])
                                	y[c1] = float(dl[i_cols[1]+1])	# for comma
                                	c1 = c1+1
					break
				except: 
					foo = 1		
	for line in dr:
		dl = line.split()
		if len(dl) > 1:
			try:
				x[c1] = float(dl[i_cols[0]])
				y[c1] = float(dl[i_cols[1]])
				c1 = c1+1
			except:
				try:
					x[c1] = float(dl[i_cols[0]])
					y[c1] = float(dl[i_cols[1]+1])
					c1 = c1+1
				except:
					break

						
	x = x[0:c1]
	y = y[0:c1]
	
	return (x,y,xvg_title1, xvg_xlab1, xvg_ylab1)


def gro(f_name, n_row, n_col, pos1):
# n_row = num tsteps, n_col = num atoms

	
	cx=np.zeros(shape=(n_row,n_col))  #initialize time by atom arrays 
	cy=np.zeros(shape=(n_row,n_col))
	cz=np.zeros(shape=(n_row,n_col))
	vx=np.zeros(shape=(n_row,n_col))
	vy=np.zeros(shape=(n_row,n_col))
	vz=np.zeros(shape=(n_row,n_col))
	
	C = [cx, cy, cz, vx, vy, vz]
	
	n_col_pre = 3  # number of columns w/o data
	i_row = 0     # init index for tstep
	i_col = -1    # init index for atom 
	
	#pdb.set_trace()
	d = open(f_name, 'r')
	print "file opened"
	c = 0


	for line in d:
		if len(line)>pos1:  # if line with data
			
			c = c+1
			
			if c%(n_col*100) == 0:
				print "Reading, ", 100*i_row/n_row,"%"
				#pdb.set_trace()
			
			d1 = line.split()
			i_col = int(d1[2])-1  # read atom index from line
			
			cx[i_row,i_col] = float(d1[3])
			cy[i_row,i_col] = float(d1[4])
			cz[i_row,i_col] = float(d1[5])
			vx[i_row,i_col] = float(d1[6])
			vy[i_row,i_col] = float(d1[7])
			vz[i_row,i_col] = float(d1[8])
			#pdb.set_trace()
			if i_col == n_col-1:   # if last atom, increment row 	
				#pdb.set_trace()
				i_row = i_row+1
		# end if 
	# end for
	d.close()
	
	return (cx,cy,cz,vx,vy,vz)

def gro_pos(f_name, n_row, n_col, pos1):
# n_row = num tsteps, n_col = num atoms
	
	cx=np.zeros(shape=(n_row,n_col))  #initialize time by atom arrays 
	cy=np.zeros(shape=(n_row,n_col))
	cz=np.zeros(shape=(n_row,n_col))
	
	n_col_pre = 3  # number of columns w/o data
	i_row = 0     # init index for tstep
	i_col = -1    # init index for atom 
	
	#pdb.set_trace()
	d = open(f_name, 'r')
	print "file opened"
	c = 0

	for line in d:
		d1 = line.split()
		if len(d1)==pos1:  # if line with data
			
			c = c+1
			
			if c%(n_col*100) == 0:
				print "Reading, ", 100*i_row/n_row,"%"
				#pdb.set_trace()

			i_col = int(d1[2])-1  # read atom index from line
			
			cx[i_row,i_col] = float(d1[3])
			cy[i_row,i_col] = float(d1[4])
			cz[i_row,i_col] = float(d1[5])

			#pdb.set_trace()
			if i_col == n_col-1:   # if last atom, increment row 	
				#pdb.set_trace()
				i_row = i_row+1
		# end if 
	# end for
	d.close()
	
	return (cx,cy,cz)
			
def gro_to_comtxt(f_name_in, f_name_out, n_row, n_col, pos1, pos2, pos3):
	# arg in: filenames in/out, number of time steps, number of atoms, min_length, initial time (ps), end time (ps)
	
	# Time params
	dt = 0.02
	t_wr = round(pos2,5); t_wro = round(pos2,5); t_len = round(pos3,5) - t_wro + dt
	
	dr = open(f_name_in, 'r')   # file to read in
	dw = open(f_name_out, 'w')  # file to write to 
	
	n_col_pre = 3  # number of columns w/o data
	i_row = 0     # init index for tstep
	i_col = -1    # init index for atom 
	len_dline = 9   # expected line length 
	str_n_col = str(n_col)
	
	c = 0

	c3 = 0 
	spc = ' '
	spc1 = '  '
	
	# coordinates
	c_o = np.zeros( (1,6), dtype = float)
	c_h1 = np.zeros( (1,6), dtype = float)
	c_h2 = np.zeros( (1,6), dtype = float)
	c_mol = np.zeros( (1,6), dtype = float)
	i_c = range(0,6)
	MW_H = 1.00794
	MW_O = 15.9994
	MW_sum = 2*MW_H + MW_O
	box_str_id = 3 # box size
	n_atom_str = '895 \n'
	
	for line in dr:
		len_line = len(line)
		dline = line.split()
		if len_line<pos1:  # if line with out data
			c2 = 0
			#pdb.set_trace()
			if len(dline) == box_str_id:  # box dimensions
				line1 = '   '+line
				if c > 1:
					dw.write(line1)
				c3 = c3+1
				if c3%500 == 0:
					print "Calc..."
				if t_wr == pos3:	
					break
				t_wr = round((t_wr + dt),5) # increment time 
				#if t_wr == 816.26:
					#pdb.set_trace()
			elif len_line < 10:   # atom number
				line1 = n_atom_str
				dw.write(line1)
			else:    		# time header string
				line1 = line
				dw.write(line1)
			
		else:  # line with data
			c = c+1
			#pdb.set_trace()
			
			# if coordinates lumped together, separate
			if len(dline) < len_dline:  

				tmp = [0]
				for nc in range(3,len(dline)):
					i_neg = dline[nc].rfind('-')
					if i_neg > 0:
						tmp.append([nc, dline[nc][0:i_neg], dline[nc][i_neg:len(dline[nc])] ] )
				del tmp[0]
				for nc in range(0,len(tmp)):
					tmp[nc][0] = tmp[nc][0]+nc
					i1 = tmp[nc][0]
					dline.insert(i1,tmp[nc][1])
					dline.insert(i1+1,tmp[nc][2])
					del dline[i1+2]
			
			# COM calculation
			if dline[1] == "OW":
				for n in i_c:
					c_o[0,n] = float(dline[n+n_col_pre])
			if dline[1] == "HW1":
				for n in i_c:
					c_h1[0,n] = float(dline[n+n_col_pre])
			if dline[1] == "HW2":
				c2 = c2+1

				for n in i_c:
					c_h2[0,n] = float(dline[n+n_col_pre])
			
				sc1 = np.multiply(MW_H, np.add(c_h1[0,0:3], c_h2[0,0:3]))
				sc = np.multiply((1/MW_sum), np.add(sc1, np.multiply(MW_O, c_o[0,0:3])) )
			
				sv1 = np.multiply(MW_H, np.add(c_h1[0,3:6], c_h2[0,3:6]))
				sv = np.multiply((1/MW_sum), np.add(sv1, np.multiply(MW_O, c_o[0,3:6])) )

				str1 = dline[0]+'  '
				str2 = 'OW'
				#pdb.set_trace()
				dw.write("%10s %4s %4.0f %9.5f %9.5f %9.5f %10.6f %10.6f %10.6f\n" %(str1, str2, c2, sc[0],sc[1],sc[2],sv[0],sv[1],sv[2]))
				# formating for mol num,mol name, atom name, atom num, 3 coords, 3 vels
	dr.close()
	dw.close()							

	foo = 1
	return (foo)
				
				

				
	
def xvg1(f_name, id_st1, id_st2):
# file name, number of rows data, number of y_vars, identifier for start of data
	
	n_row = 1000000
	x = np.zeros(n_row)
	y = np.zeros(n_row)
	
	d = open(f_name, 'r')
	print f_name+" opened"

	for line in d:
		d1 = line.split()
		if len(d1) > 1 and d1[0] == id_st1 and d1[1] == id_st2:
			break
	
	c = 0
	for line in d:

		d1 = line.split()
		
		try:
			x[c] = float(d1[0])
			y[c] = float(d1[1])
		except: 
			print "may have exceeded data size..."
			pdb.set_trace()
		c = c+1
	x_out = x[0:c]
	y_out = y[0:c]
	
	d.close()		

	return (x_out,y_out)	


def avg_xvg(p_names_in, id_st1, id_st2):
        N = len(p_names_in)
        xo,yo = xvg1(p_names_in[0], id_st1, id_st2)

        ys = np.zeros((len(yo),N))

        for n in range(0,N):

                xo,yo = xvg1(p_names_in[n], id_st1, id_st2)
                ys[:,n] = yo
        y_mn = np.mean(ys,axis=1)
        y_std = np.std(ys,axis=1)
        y_err = np.divide(y_std,np.sqrt(N-1))
        return (xo,y_mn,y_std,y_err)




def xvg01(f_name, n_row, n_y, id_st1, id_st2, id_end1):
# file name, number of rows data, number of y_vars, identifier for start of data
	
	x = np.zeros(shape=(n_row,1))
	y = np.zeros(shape=(n_row,n_y))
	
	d = open(f_name, 'r')
	print "file opened"
	c = 0

	for line in d:
		c = c+1
		d1 = line.split()

		if len(d1) > 1 and d1[0] == id_st1 and d1[1] == id_st2:

			break
		
	c = 0
	
	if n_y == 1:
		for line in d:

			if line == id_end1:
				break
			if int(n_row) == c:
				break
			d1 = line.split()
			
			try:
				x[c] = float(d1[0])
				y[c] = float(d1[1])
			except: 
				pdb.set_trace()
			c = c+1
			#if c == 500:
				#pdb.set_trace()
		#pdb.set_trace()			
			
	else:
		for line in d:
			#pdb.set_trace()
			if line == id_end1:
				break
			d1 = line.split()
			try:
				x[0,c] = float(d1[0])
			except:
				pdb.set_trace()	
			for n in range(1,n_y+1):
				#pdb.set_trace()
				y[n-1,c] = float(d1[n])
			c = c+1
			
	d.close()		

	return (x,y)

def xvg2(f_name, n_row, n_y, id1_len, id1_sts, id_end1):
# file name, number of rows data, number of y_vars
# identifier for start of data: number of fields and either first 1 or 2 fields
	
	x = np.zeros(shape=(n_row,1))
	y = np.zeros(shape=(n_row,n_y))
	
	d = open(f_name, 'r')
	print "file opened"
	c = 0

	for line in d:
		c = c+1
		d1 = line.split()
		
		if len(d1) == id1_len:
                        if id1_len == 1 and d1[0] == id1_sts[0]:
                                break
                        elif id1_len > 1 and d1[0] == id1_sts[0] and d1[1] == id1_sts[1]:
                                break
		
	c = 0
	
	if n_y == 1:
		for line in d:

			if line == id_end1:
				break
			if int(n_row) == c:
				break
			d1 = line.split()
			
			try:
				x[c] = float(d1[0])
				y[c] = float(d1[1])
			except: 
				pdb.set_trace()
			c = c+1
			#if c == 500:
				#pdb.set_trace()
		#pdb.set_trace()			
			
	else:
		for line in d:
			#pdb.set_trace()
			if line == id_end1:
				break
			d1 = line.split()
			try:
				x[0,c] = float(d1[0])
			except:
				pdb.set_trace()	
			for n in range(1,n_y+1):
				#pdb.set_trace()
				y[n-1,c] = float(d1[n])
			c = c+1
			
	d.close()		

	return (x,y)
			
def xvg_text(f_name, row_id, col_id, char_id):
# finds a particular item

	d = open(f_name, 'r')
	#print "file opened"
	c = 0
	d_out = [0]
	
	for line in d:
		
		for i in range(0,len(row_id)):
			
			if c == row_id[i]:
				d1 = line.split()

				if char_id[i] == 'na':
					
					d_val1 = d1[col_id[i]]
				else:
					d_val1o = d1[col_id[i]]
					ind_char = d_val1o.rfind(char_id[i])
					d_val1 = d_val1o[0:ind_char]
				#pdb.set_trace()
				d_out.append(float(d_val1))
				
		if len(d_out) == len(row_id)+1:
			break		
				
				
		c = c+1	
	
	del d_out[0]
	
	return (d_out)
			
		
def count_occ(f_name, str1):
	ct = 0
	dr = open(f_name,'r')
	for line in dr:	
		ct0 = string.count(line, str1)
		ct += ct0


	dr.close()
	return ct
		
