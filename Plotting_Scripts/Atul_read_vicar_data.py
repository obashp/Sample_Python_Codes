import numpy as np
import struct
import re
import pandas as pd

def read_binary(file):
    
    with open(file, 'rb') as f:
        text = f.read()

    pos=44
    NUMVAR = struct.unpack('i',text[pos:pos+4])[0]
    pos+=4;
    print('NUMVAR = ',NUMVAR)
    
    i=0
    varname=['']*NUMVAR
    while (i < NUMVAR):
        temp = struct.unpack('i',text[pos:pos+4])[0]
        pos +=4
        
        if (temp == 0):
            i += 1
        else:
            varname[i] += chr(temp)
    
    pos=pos+40
    
    NI = int.from_bytes(text[pos:pos+4],byteorder='little',signed=True)
    pos+=4; print('NI = ',NI)
    NJ = int.from_bytes(text[pos:pos+4],byteorder='little',signed=True)
    pos+=4; print('NJ = ',NJ)
    NK = int.from_bytes(text[pos:pos+4],byteorder='little',signed=True)
    pos+=4; print('NK = ',NK)
    
    if NK == 1:
        ndim = 2
    else:
        ndim = 3
    
    DATAMARKER = struct.unpack('f',text[pos:pos+4])[0]
    print('DATAMARKER = ',DATAMARKER); pos+=4; 
    ZONEMARKER = struct.unpack('f',text[pos:pos+4])[0]
    print('ZONEMARKER = ',ZONEMARKER); pos+=4; 
    NUMRP = struct.unpack('i',text[pos:pos+4])[0]
    pos+=4; print('NUMRP = ',NUMRP)
    
    for i in range(NUMVAR):
        DATATYPE = struct.unpack('i',text[pos:pos+4])[0]
        pos+=4;
    
    #print('pos = ', pos)
    if DATATYPE ==1:
        size = 4
        fmt = 'f'
    else:
        size = 8
        fmt = 'd'
        
    xc = np.zeros((NI,1),dtype=float)
    yc = np.zeros((NJ,1),dtype=float)
    qq = np.zeros((NJ,NI,NK,NUMVAR),dtype=float)
    
    for i in range(NI):
        str1 = text[pos:pos+size]
        xc[i] = struct.unpack(fmt,str1)[0]
        pos+=size;
    #    print('xc(i) = ',xc[i][0], ', binary = ',str1);
        
    for j in range(NJ-1):
        for i in range(NI):
            pos+=size;
            
    for j in range(NJ):
        str1 = text[pos:pos+size]
        yc[j] = struct.unpack(fmt,str1)[0]
        pos+=size; 
    #    print('yc(i) = ',xc[j][0], ', binary = ',str1)
        for i in range(NI-1):    
            pos+=size;
            
    #print('pos = ', pos)
    
    for kk in range(NUMVAR-ndim):
        for j in range(NJ):
            for i in range(NI):
                str1 = text[pos:pos+size]
                qq[j,i,0,kk] = struct.unpack(fmt,str1)[0]
                pos+=size

    return(xc,yc,qq)


def read_marker2d(file):
	nline = 2
	pos = 0
	N = 4
	file1  = open(file,'r').read()
	nbody = file1.count('ZONE')
	fp = open(file)
	ibody = 0
	for i, line in enumerate(fp):
		if i == 1:
			varname=(re.findall(r'".*"', line)[0].split(','))
			nb_var = len(varname)
		if (line.find('ZONE') != -1):
			ibody += 1
			N1 = int(re.findall(r'N=.*E=', line)[0].strip('N= * E'))
			N1 = N1//3
			N = max(N,N1)
	fp.close()
	print('No of body', nbody)
	data = np.zeros((N,nb_var,nbody), dtype=float)
	N = np.zeros((2,1),dtype = int)
	ibody = 0
	fp = open(file)
	for i, line in enumerate(fp):
		if (line.find('ZONE') != -1):
			N1 = int(re.findall(r'N=.*E=', line)[0].strip('N= * E'))
			N1 = N1//3
			N[ibody] = N1
			print('number of markers in body ',ibody+1,'= ',N1)
			nline = i
			found = 1
		if pos == N1:
			ibody += 1
			pos = 0
			found = 0
		if ibody == nbody:
			break
		if i > nline and pos < N1 and found == 1:
			temp = np.fromstring(line, sep=' ',count=nb_var)
			data[pos,:,ibody] = temp
			pos += 1
	fp.close()
	return(nbody, N, data)
	
def read_Nudata(file):
	line=open(file).readlines()
	varname=(re.findall(r'".*"', line[0])[0].split(','))
	nb_var = len(varname)
	N = int(re.findall(r'I =.*', line[1])[0].strip('I ='))
	data = np.zeros((N,nb_var), dtype=float)
	i = 0
	pos = 2
	ist = 0
	iend = 0
	while (i < N):
		temp = np.fromstring( line[pos], dtype=np.float, sep=' ' )
		iend = ist + temp.size
		data[i,ist:iend] = temp
		ist = iend
		pos += 1
		if iend == nb_var:
			i += 1
			ist = 0
	return(data)
    