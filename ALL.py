#By CUMT stu
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)  #保证numpy矩阵显示完全
#定义全局变量:
GM=3.986005e+14   #GM在WGS-84坐标系中的地球引力常数
we=7.29211567e-5  #地球自转角速度（GPS）
#我的学号：07152660，按要求需计算14时37分40秒的卫星位置，即14 37 30，瞬时时：52650
LV=299792458 #光速
'''----------------------------------------------------------------------------------------------------------------------------'''
 #读N文件
'''----------------------------------------------------------------------------------------------------------------------------'''
def readN(ephemerisfile,outputfile):
	fide = open(ephemerisfile,'r+')
	#line=fide.readlines	#读整个文件
	lines = fide.readlines()#.splitlines(True)
	fide.seek(0)

	for s in lines:
		fide.write(s.replace('D','e'))# replace是替换，write是写入
	
	fide.close()
	fide = open(ephemerisfile,'r')
	line = fide.readlines()#.splitlines(True)
	ALPHA = np.zeros((4))
	BETA = np.zeros((4))
	for i,content in enumerate(line):#i为行数，content为内容
		if "ION ALPHA" in content:
			ALPHA[0]=line[i][3:14]
			ALPHA[1]=line[i][15:26]
			ALPHA[2]=line[i][27:38]
			ALPHA[3]=line[i][39:50]
			i=i+1
			BETA[0]=line[i][3:14]
			BETA[1]=line[i][15:26]
			BETA[2]=line[i][27:38]
			BETA[3]=line[i][39:50]
			break

	for i,content in enumerate(line):#i为行数，content为内容
		if "ENe OF HEAeER" in content:
			i=i+1
			break

	#print(line[20][:20])#检查是否把D转成e

	noeph = int((len(line)-i)/8)#卫星数
	#print(noeph)
	fide.close()

	year=np.zeros((noeph))
	month=np.zeros((noeph))
	day=np.zeros((noeph))
	hour=np.zeros((noeph))
	minute=np.zeros((noeph))
	second=np.zeros((noeph))
	svprn=np.zeros((noeph))
	weekno=np.zeros((noeph))
	t0c=np.zeros((noeph))
	tgd=np.zeros((noeph))
	aodc=np.zeros((noeph))
	toe=np.zeros((noeph))
	af2=np.zeros((noeph))
	af1=np.zeros((noeph))
	af0=np.zeros((noeph))
	aode=np.zeros((noeph))
	deltan=np.zeros((noeph))
	M0=np.zeros((noeph))
	ecc=np.zeros((noeph))
	roota=np.zeros((noeph))
	toe=np.zeros((noeph))
	cic=np.zeros((noeph))
	crc=np.zeros((noeph))
	cis=np.zeros((noeph))
	crs=np.zeros((noeph))
	cuc=np.zeros((noeph))
	cus=np.zeros((noeph))
	Omega0=np.zeros((noeph))
	omega=np.zeros((noeph))
	i0=np.zeros((noeph))
	Omegadot=np.zeros((noeph))
	idot=np.zeros((noeph))
	accuracy=np.zeros((noeph))
	health=np.zeros((noeph))
	fit=np.zeros((noeph))
	IODE=np.zeros((noeph))
	codes=np.zeros((noeph))
	L2flag=np.zeros((noeph))
	svaccur=np.zeros((noeph))
	svhealth=np.zeros((noeph))
	iodc=np.zeros((noeph))
	tom=np.zeros((noeph))
	spare=np.zeros((noeph))
	eph=np.zeros((21,noeph))
	
	j=0
	for x in range(noeph):
		svprn[j]=(line[i][0:2])
		year[j] = line[i][3:5]
		month[j] = line[i][6:8]
		day[j] = line[i][9:11]
		hour[j] = line[i][12:14]
		minute[j] = line[i][15:17]
		second[j] = line[i][19:22]
		af0[j] = (line[i][22:41])
		af1[j] = line[i][41:60]
		af2[j] = line[i][60:79]
		i=i+1
		IODE[j] = line[i][3:22]
		crs[j] = line[i][22:41]
		deltan[j] = line[i][41:60]
		M0[j] = (line[i][60:79])
		i=i+1
		cuc[j] = (line[i][3:22])
		ecc[j] = (line[i][22:41])
		cus[j] = (line[i][41:60])
		roota[j] = (line[i][60:79])
		i=i+1
		toe[j] = (line[i][3:22])
		cic[j] = (line[i][22:41])
		Omega0[j] = (line[i][41:60])
		cis[j] = (line[i][60:79])
		i=i+1
		i0[j] =  (line[i][3:22])
		crc[j] = (line[i][22:41])
		omega[j] = (line[i][41:60])
		Omegadot[j] = (line[i][60:79])
		i=i+1
		idot[j] = (line[i][3:22])
		codes[j] = (line[i][22:41])
		weekno[j] = (line[i][41:60])
		L2flag[j] = (line[i][60:79])
		i=i+1
		svaccur[j] = (line[i][3:22])
		svhealth[j] = (line[i][22:41])
		tgd[j] = (line[i][41:60])
		iodc[j] = line[i][60:79]
		i=i+1
		tom[j] = ((line[i][3:22]))
		spare[j] = line[i][22:79]
		i=i+1
		j=j+1
	
	eph[0][:noeph]  = svprn  #卫星prn号
	eph[1][:noeph]  = af2 #钟漂
	eph[2][:noeph]  = M0 # 参考时刻的平近点角（开普勒参数）
	eph[3][:noeph]  = roota #参考时刻的轨道半径平方根（开普勒参数）
	eph[4][:noeph]  = deltan #卫星平均角速度的改正值
	eph[5][:noeph]  = ecc #轨道偏心率（开普勒参数）
	eph[6][:noeph]  = omega #近地角距（开普勒参数）
	eph[7][:noeph]  = cuc #升交角距的正弦余值
	eph[8][:noeph]  = cus #升角距的改正项振幅
	eph[9][:noeph]  = crc #轨道向径的正余弦调
	eph[10][:noeph] = crs #轨道向径的改正项振幅
	eph[11][:noeph] = i0 #轨道倾角
	eph[12][:noeph] = idot #卫星轨道倾角变化率
	eph[13][:noeph] = cic #轨道倾角的正余弦调
	eph[14][:noeph] = cis #轨道倾角的改正项振幅
	eph[15][:noeph] = Omega0 #升交点赤径（开普勒参数）
	eph[16][:noeph] = Omegadot #升交点赤径变化率
	eph[17][:noeph] = toe #星历参考时间
	eph[18][:noeph] = af0 #钟差
	eph[19][:noeph] = af1 #钟速
	eph[20][:noeph] = tom #信息传输时间
		
	with open(outputfile,'w') as fidu:
		fidu.write(str(eph))
		#print(svprn[2])
	return eph,noeph,ALPHA,BETA
'''----------------------------------------------------------------------------------------------------------------------------'''
#读O文件
'''----------------------------------------------------------------------------------------------------------------------------'''
def readO(ephemerisfile,outputfile):
	fide= open(ephemerisfile,'r') 
	lines = fide.readlines()
	fide.close()
	for i,content in enumerate(lines):#i为行数，content为内容
		if "APPROX POSITION XYZ " in content:#大致坐标
			nx=float(lines[i][:15])
			ny=float(lines[i][15:28])
			nz=float(lines[i][30:42])	
			break
	i=i+1
	TXH=lines[i][:15] #天线高
	
	for i,content in enumerate(lines):#i为行数，content为内容
		if "INTERVAL" in content:
			interval=int(float(lines[i][4:10])) #时间间隔
			break
	i=i+1
	year=int(lines[i][2:6]) #用于时间转换
	month=int(lines[i][10:12])
	day=int(lines[i][16:18])
	first_clock=lines[i][22:25]#为了计算循环参数K
	first_minute=lines[i][28:31]
	first_second=lines[i][33:35]
	i=i+1
	last_clock=lines[i][22:25]
	last_minute=lines[i][28:31]
	last_second=lines[i][33:35]
	
	K=int(((int(last_clock)-int(first_clock))*3600+(int(last_minute)-int(first_clock))*60+(int(last_second)-int(first_second))*1)/interval)
	#观测数(历元数)
	#print(K)
	for i,content in enumerate(lines):#i为行数，content为内容
		if "END OF HEADER" in content:
			break

	j=i+1#用于确定各个观测历元中观测到的卫星最大值，从而确定Ocontent的维数
	sates=[] #每个历元观测的卫星数

	for x in range(K):
		sa=0
		sa=int(lines[j][33:35])
		sates.append(lines[j][33:35])
		j=j+sa+1
		

	M=int(max(sates)) #观测卫星数最多的历元

	Ocontent=np.zeros((8,K,M))
	
	for x in range(K):   #大循环
		i=i+1
		sate_num=lines[i][33:35] #卫星数，非数组
		
		for n in range(int(sate_num)):
			i=i+1
			Ocontent[0][x][n]=int(lines[i][1:3])#卫星的prn号

			if lines[i][5:17].isspace():
				Ocontent[1][x][n]=None
			else:
				Ocontent[1][x][n]=lines[i][5:17]
			
			if lines[i][20:33].isspace():
				Ocontent[2][x][n]=None
			else:
				Ocontent[2][x][n]=lines[i][20:33]

			if lines[i][36:49].isspace():
				Ocontent[3][x][n]=None
			else:
				Ocontent[3][x][n]=lines[i][36:51] #可能出现空列表的情况

			if lines[i][53:65].isspace():
				Ocontent[4][x][n]=None
			else:
				Ocontent[4][x][n]=lines[i][53:65]

			if lines[i][69:81].isspace():
				Ocontent[5][x][n]=None
			else:
				Ocontent[5][x][n]=lines[i][69:81]

			if lines[i][90:97].isspace():
				Ocontent[6][x][n]=None
			else:
				Ocontent[6][x][n]=lines[i][90:97]

			if (lines[i][107:114].isspace())or(not lines[i][107:114]): #可能存在全是空格或为空两种情况
				Ocontent[7][x][n]=None
			else:
				Ocontent[7][x][n]=lines[i][107:114]
			#print(Ocontent[6][x][n])
			
	with open(outputfile,'w') as fidu:
		fidu.write(str(Ocontent))
	
	return Ocontent,sates,year,month,day,nx,ny,nz
'''----------------------------------------------------------------------------------------------------------------------------'''
#年月日时转为儒略日，再将儒略日转为GPS时
'''----------------------------------------------------------------------------------------------------------------------------'''
def TIMEchange(year,month,day,ST):
	#年月日时转为儒略日
	clock=ST/3600
	if(month<=2):
		year=year-1
		month=month+12
	else:
		jd = math.floor(365.25*(year+4716))+math.floor(30.6001*(month+1))+day+clock/24-1537.5
	#儒略日转为GPS时
	a = math.floor(jd+0.5)
	b = a+1537
	c = math.floor((b-122.1)/365.25)
	e = math.floor(365.25*c)
	f = math.floor((b-e)/30.6001)
	d = b-e-math.floor(30.6001*f)+(jd+0.5)%1
	day_of_week = math.floor(jd+0.5)%7
	week = math.floor((jd-2444244.5)/7)
	sec_of_week = ((d%1)+day_of_week+1)*86400
	if(day_of_week==6):#如果当天就是周日
		sec_of_week=sec_of_week-7*24*3600
	#print(day_of_week)
	#print(sec_of_week)#周内秒
	return sec_of_week
'''----------------------------------------------------------------------------------------------------------------------------'''
#确定需要计算用的卫星及数据准备
'''----------------------------------------------------------------------------------------------------------------------------'''
def getsat(eph,noeph,Ocontent,sec_of_week,sates,ST):
	i=int(ST/30)     #瞬时时间在Ocontent矩阵中的位置
	M=np.zeros((noeph))	
	j=[]		#eph矩阵离瞬时矩阵最近的时间的列号
	for x in range(noeph): #185是noeph的值
		M[x]=abs(sec_of_week-eph[17][x])#取绝对值
	MIN=min(M)
	#print(MIN)
	for x in range(noeph):
		if(sec_of_week-eph[17][x]==MIN):
			j.append(x) #看了下矩阵应该有12个卫星,但是瞬时只观测到9个卫星，应当把12个卫星位置都计算出来,j存了N文件卫星prn号的列号
	prn=np.zeros((len(j)))
	sateeph=np.zeros((21,int(sates[i])))
	for x in range(len(j)):
		prn[x]=eph[0][j[0]+x]
	for x in range(int(sates[i])):
		for n in range(len(j)):
			if(Ocontent[0][i][x]==prn[n]):
				for w in range(21):
					sateeph[w][x]=eph[w][j[0]+n] #让该历元N文件数据按O文件卫星顺序排列
					#print(sateeph[w][x])
	inde=[]

	for index,value in enumerate(sateeph[2][:]):#找出该数据块中在O文件但不在N文件的卫星行号
		if(value==0):
			inde.append(index)
	#print(inde)
	ture_O_satenum=int(sates[i])-len(inde)
	sateeph_new=np.zeros((21,ture_O_satenum))
	#print(prn[0])
	
	for w in range(21):
		sateeph_new[w][:]=np.delete(sateeph[w][0:int(sates[i])],inde)
	#以上完成了N文件数据sateeph转换为O文件卫星排列对应格式,得到sateeph_new
	#以下完成了O文件删除无对应卫星数据的工作，并存到Ocontent_new矩阵
	Ocontent_new=np.zeros((8,ture_O_satenum))	#prn
	#print(sateeph_new)
	for w in range(8):
		Ocontent_new[w][:]=np.delete(Ocontent[w][i][0:int(sates[i])],inde)
	#print(Ocontent_new[0][:])
	return sateeph_new,Ocontent_new,ture_O_satenum
'''----------------------------------------------------------------------------------------------------------------------------'''
#卫星钟差改正
'''----------------------------------------------------------------------------------------------------------------------------'''
def SateTimecorrection(ture_O_satenum,sateeph_new,Ocontent_new,sec_of_week):
	#tx_RAW=np.zeros((len(o1)))
	tx_RAW = sec_of_week - Ocontent_new[4][:]/LV #考虑信号传播过程中的卫星运动
	toe=sateeph_new[17][:]
	tt=tx_RAW-toe#卫星钟的参考时刻
	#GPS时间超限或下溢的修复:
	half_week = 302400
	for x in range(ture_O_satenum):
		if(tt[x] > half_week):
			tt[x] = tt-2*half_week
		if(tt[x] < -half_week):
			tt[x] = tt+2*half_week
	#钟差改正：
	tcorr = (sateeph_new[1][:]*tt + sateeph_new[19][:])*tt + sateeph_new[18][:]
	tx_GPS = tx_RAW-tcorr
	#print(tcorr*LV)
	return tx_GPS,tcorr
'''----------------------------------------------------------------------------------------------------------------------------'''
#卫星位置计算
'''----------------------------------------------------------------------------------------------------------------------------'''
def satelliteposition(sateeph_new,Ocontent_new,sates,ture_O_satenum,tx_GPS):

	prn=np.zeros((ture_O_satenum))
	toe=np.zeros((ture_O_satenum))
	n0=np.zeros((ture_O_satenum)) #未改正的平均角速度
	n =np.zeros((ture_O_satenum)) #改正后的平均角速度
	PJD0=np.zeros((ture_O_satenum))#参考历元平近点角
	PJD=np.zeros((ture_O_satenum)) #瞬时历元平近点角
	E=np.zeros((ture_O_satenum))  #偏近点角
	ecc=np.zeros((ture_O_satenum))#偏心率
	f=np.zeros((ture_O_satenum))  #真近点角
	u0=np.zeros((ture_O_satenum))  #升交距角(未考虑摄动改正)
	u=np.zeros((ture_O_satenum))   #升交距角（已改正）
	r0=np.zeros((ture_O_satenum)) #卫星矢径（未改正）
	r=np.zeros((ture_O_satenum))  #卫星矢径（已改正）
	i0=np.zeros((ture_O_satenum)) #卫星轨道倾角（未改正）
	ii=np.zeros((ture_O_satenum)) #卫星轨道倾角（已改正）
	idot=np.zeros((ture_O_satenum))#卫星轨道倾角变化率
	omega=np.zeros((ture_O_satenum)) #近地点角距
	cuc=np.zeros((ture_O_satenum))#摄动改正参数
	cus=np.zeros((ture_O_satenum))
	crc=np.zeros((ture_O_satenum))
	crs=np.zeros((ture_O_satenum))
	cic=np.zeros((ture_O_satenum))
	cis=np.zeros((ture_O_satenum))
	GZu=np.zeros((ture_O_satenum))#摄动改正项
	GZr=np.zeros((ture_O_satenum))
	GZi=np.zeros((ture_O_satenum))
	for x in range(ture_O_satenum):
		n0[x]=math.sqrt(GM)/(math.pow(sateeph_new[3][x],3)) #eph[3][:]为参考时刻的轨道半径平方根

	n=n0+sateeph_new[4][:]
	PJD0=sateeph_new[2][:]
	ecc=sateeph_new[5][:]
	toe=sateeph_new[17][:]
	prn=sateeph_new[0][:]

	tk=tx_GPS-toe
	#print(tk)
	#GPS时间超限或下溢的修复:
	half_week = 302400
	for x in range(ture_O_satenum):
		if(tk[x] > half_week):
			tk[x] = tk-2*half_week
		if(tk[x] < -half_week):
			tk[x] = tk+2*half_week
	#print(tk)
	PJD=PJD0+np.multiply(n,tk) #矩阵点乘
	PJD=(PJD+2*math.pi)%(2*math.pi)
	#print(PJD)
	#迭代求E，偏近点角,这里要注意不用把E化为角度，因为在python使用sin（）内部默认为弧度
	for x in range(7):#发现7次迭代基本收敛
		E=PJD+np.multiply(ecc,np.sin(E))
		#print(E)
	E=(E+2*math.pi)%(2*math.pi)
	for x in range(ture_O_satenum):
		f[x]=math.atan2((math.sqrt(1-ecc[x]*ecc[x])*math.sin(E[x])),(math.cos(E[x])-ecc[x]))
	
	omega=sateeph_new[6][:]
	u0=f+omega #升交距角
	#u0=(u0+2*math.pi)%(2*math.pi)
	#计算摄动改正项
	
	cuc=sateeph_new[7][:]
	cus=sateeph_new[8][:]
	crc=sateeph_new[9][:]
	crs=sateeph_new[10][:]
	cic=sateeph_new[13][:]
	cis=sateeph_new[14][:]
	GZu=np.multiply(cuc,np.cos(2*u0))+np.multiply(cus,np.sin(2*u0))
	GZr=np.multiply(crc,np.cos(2*u0))+np.multiply(crs,np.sin(2*u0))
	GZi=np.multiply(cic,np.cos(2*u0))+np.multiply(cis,np.sin(2*u0))
	#计算受摄卫星矢径r0（改正计算中要用到）
	for x in range(ture_O_satenum):
		r0[x]=math.pow(sateeph_new[3][x],2)*(1-ecc[x]*math.cos(E[x]))
	i0=sateeph_new[11][:]
	idot=sateeph_new[12][:]
	#进行改正计算
	u=u0+GZu
	r=r0+GZr
	ii=i0+GZi+np.multiply(idot,tk)
	
	#卫星在轨道平面坐标系的坐标计算：
	xk=np.zeros((ture_O_satenum))
	yk=np.zeros((ture_O_satenum))
	xk=np.multiply(r,np.cos(u))
	yk=np.multiply(r,np.sin(u))

	#计算观测瞬间升交点的经度L：
	Omega0=np.zeros((ture_O_satenum))#升交点赤径（开普勒参数）
	#注意：广播星历给出的不是Omegatoe，而是该值与本周起始时刻的格林尼治恒星时GASTweek之差，Omega0
	Omegadot=np.zeros((ture_O_satenum))#升交点赤径变化率
	#Omegas=np.zeros((len(j)))#观测瞬时的升径交点赤径（待计算）
	#for x in range(len(j)):#导出Omega0和Omegadot数据
	Omega0=sateeph_new[15][:]
	Omegadot=sateeph_new[16][:]
	L=Omega0+(Omegadot-we)*tk-we*toe#注意，书上的公式是错的
	L=(L+2*math.pi)%(2*math.pi) #保证在小周期中
	#print(L)
	#计算在地固坐标系下的坐标
	Xdg=np.zeros((ture_O_satenum))
	Ydg=np.zeros((ture_O_satenum))
	Zdg=np.zeros((ture_O_satenum))
	Xdg=np.multiply(xk,np.cos(L))-np.multiply(yk,(np.multiply(np.cos(ii),np.sin(L))))
	Ydg=np.multiply(xk,np.sin(L))+np.multiply(yk,(np.multiply(np.cos(ii),np.cos(L))))
	Zdg=np.multiply(yk,np.sin(ii))

	#print(L)
	#print(Xdg,Ydg,Zdg)
	return Xdg,Ydg,Zdg
'''----------------------------------------------------------------------------------------------------------------------------'''
#消除地球自转误差
'''----------------------------------------------------------------------------------------------------------------------------'''
def earthcorr(ture_O_satenum,Xdg,Ydg,Zdg,nx,ny,nz):
	D0=np.zeros((ture_O_satenum))
	for x in range(ture_O_satenum):
		D0[x]=math.sqrt((Xdg[x]-nx)*(Xdg[x]-nx)+(Ydg[x]-ny)*(Ydg[x]-ny)+(Zdg[x]-nz)*(Zdg[x]-nz))
	#消除地球自转影响：
	traveltime=D0/LV
	omegatau=traveltime*we
	X0=Xdg*np.cos(omegatau)+Ydg*np.sin(omegatau)
	Y0=-Xdg*np.sin(omegatau)+Ydg*np.cos(omegatau)
	Z0=Zdg
	xjsj=nx
	yjsj=ny
	zjsj=nz
	return xjsj,yjsj,zjsj,X0,Y0,Z0
'''----------------------------------------------------------------------------------------------------------------------------'''
#用给出的笛卡尔坐标X，Y，Z计算大地坐标求对应参考椭球下的大地坐标（经纬度和大地高）
'''----------------------------------------------------------------------------------------------------------------------------'''
def togeod(nx,ny,nz):
	a=6378137 #椭圆长半径
	finv=298.257223563 #扁率
	h = 0
	tolsq = 1.e-10
	maxit = 10
	rtd = 180/math.pi
	esq=(2-1/finv)/finv
	oneesq=1-esq
	P=math.sqrt(nx*nx+ny*ny)
	dlambda=math.atan2(ny,nx)*rtd
	if(dlambda<0):
		dlambda=dlambda+360
	r=math.sqrt(P*P+nz*nz)
	sinphi=nz/r
	dphi=math.asin(sinphi)
	h=r-a*(1-sinphi*sinphi/finv)
	for i in range(maxit):
		sinphi=math.sin(dphi)
		cosphi=math.cos(dphi)
		N_phi=a/math.sqrt(1-esq*sinphi*sinphi)
		dP=P-(N_phi+h)*cosphi
		dZ=nz-(N_phi*oneesq+h)*sinphi
		h=h+(sinphi*dZ+cosphi*dP)
		dphi = dphi+(cosphi*dZ-sinphi*dP)/(N_phi + h)
		if(dP*dP+dZ*dZ<tolsq):
			break
	dphi=dphi*rtd
	return dphi,dlambda,h
'''----------------------------------------------------------------------------------------------------------------------------'''
#计算卫星高度角
'''----------------------------------------------------------------------------------------------------------------------------'''
def getel(nx,ny,nz,X0,Y0,Z0,ture_O_satenum):
	dx=X0-nx
	dy=Y0-ny
	dz=Z0-nz
	dtr=math.pi/180
	phi,dlambda,h=togeod(nx,ny,nz)
	cl = math.cos(dlambda*dtr)
	sl = math.sin(dlambda*dtr)
	cb = math.cos(phi*dtr)
	sb = math.sin(phi*dtr)
	E=-sl*dx+(-sb)*cl*dy+cb*cl*dz
	N=cl*dx+(-sb)*sl*dy+cb*sl*dz
	U=cb*dy+sb*dz
	hor_dis=np.zeros((ture_O_satenum))
	Az=np.zeros((ture_O_satenum))
	EI=np.zeros((ture_O_satenum))
	for i in range(ture_O_satenum):
		hor_dis[i]=math.sqrt(E[i]*E[i]+N[i]*N[i])
		Az[i]=math.atan2(E[i],N[i])/dtr
		EI[i]=math.atan2(U[i],hor_dis[i])/dtr
		if(Az[i]<0):
			Az[i]=Az[i]+360
	#D=math.sqrt(dx*dx+dy*dy+dz*dz)
	return EI
'''----------------------------------------------------------------------------------------------------------------------------'''
#电离层改正（双频改正）
'''----------------------------------------------------------------------------------------------------------------------------'''
def Klobuchar(Ocontent_new):
	#el=getel(nx,ny,nz,X0,Y0,Z0,ture_O_satenum)
	#sida=445/(el+20)-4#地心夹角
	Vion=-(Ocontent_new[4][:]-Ocontent_new[5][:])*1.54573
	return Vion
'''----------------------------------------------------------------------------------------------------------------------------'''
#伪距单点定位
'''----------------------------------------------------------------------------------------------------------------------------'''
def SPLocation(tcorr,sateeph_new,Ocontent_new,ture_O_satenum,Xdg,Ydg,Zdg,nx,ny,nz):

	vT=0#接收机钟差
	xjsj,yjsj,zjsj,X0,Y0,Z0=earthcorr(ture_O_satenum,Xdg,Ydg,Zdg,nx,ny,nz)
	el=getel(nx,ny,nz,X0,Y0,Z0,ture_O_satenum)
	#删除高度角小于10度的卫星:
	el_dele=[]
	for i in range(ture_O_satenum):
		if(abs(el[i])<10):
			el_dele.append(i)
	Vion=Klobuchar(Ocontent_new)
	Vion=np.delete(Vion[0:ture_O_satenum],el_dele)
	#print(Vion)
	X0=np.delete(X0[0:ture_O_satenum],el_dele)
	Y0=np.delete(Y0[0:ture_O_satenum],el_dele)
	Z0=np.delete(Z0[0:ture_O_satenum],el_dele)
	tcorr=np.delete(tcorr[0:ture_O_satenum],el_dele)
	el=np.delete(el[0:ture_O_satenum],el_dele)

	ture_O_satenum=ture_O_satenum-len(el_dele)
	if(ture_O_satenum<4):
		print('卫星数不足，无法计算接收机位置')
	else:
		Ocontent_P1_fin=np.zeros((ture_O_satenum))
		Ocontent_P1_fin=np.delete(Ocontent_new[4][:],el_dele)#多维矩阵非要让左边矩阵和右边删除后的矩阵大小一致，因此必须定义，而一维矩阵不存在这个问题

		low0=np.zeros((ture_O_satenum))
		cs=np.zeros((ture_O_satenum))

		P=np.eye(ture_O_satenum) #暂且让其为1
		#print(el)
			#利用卫星高度角进行定权：
		el=math.pi*el/180#角度化弧度
		#print(el)
		for i in range(ture_O_satenum):
			P[i][i]=math.sin(el[i])*math.sin(el[i])

		#print(Ocontent_P1_fin)
		for e in range(6):#迭代6次
			for x in range(ture_O_satenum):
				low0[x]=math.sqrt((X0[x]-xjsj)*(X0[x]-xjsj)+(Y0[x]-yjsj)*(Y0[x]-yjsj)+(Z0[x]-zjsj)*(Z0[x]-zjsj))
			#print(low0)
			#准备位置计算矩阵数据：
			l=((X0-xjsj)/low0).T
			l.shape=(ture_O_satenum,1) #一维数组转置的时候有个坑，光transpose没有用，需要指定shape参数
			m=((Y0-yjsj)/low0).T
			m.shape=(ture_O_satenum,1)
			n=((Z0-zjsj)/low0).T
			n.shape=(ture_O_satenum,1)
			for x in range(ture_O_satenum):
				cs[x]=1
			cs.shape=(ture_O_satenum,1)
			B=np.hstack((-l,-m))#在行上合并
			B=np.hstack((B,-n))
			B=np.hstack((B,cs))
			
			
			L=(Ocontent_P1_fin-low0).T+tcorr*LV-Vion
			L.shape=(ture_O_satenum,1)

			x=np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(B),P),B)),B.T),P),L)
			xjsj=xjsj+x[0]
			yjsj=yjsj+x[1]
			zjsj=zjsj+x[2]
			vT=vT+x[3]
		#精度评定：
		jindu=math.sqrt((nx-xjsj)*(nx-xjsj)+(ny-yjsj)*(ny-yjsj)+(nz-zjsj)*(nz-zjsj))
		print(jindu)
		#print(vT) #LV*pos机钟差
	return xjsj,yjsj,zjsj,vT,jindu
'''----------------------------------------------------------------------------------------------------------------------------'''
#误差可视化（类似GPSeasy3出图，分别对比x，y，z方向的误差）
'''----------------------------------------------------------------------------------------------------------------------------'''
def picture(error,x_error,y_error,z_error):
	time = np.linspace(1,20,20)
	#print(time)
	plt.figure(1)
	plt.xlabel("time")  
	plt.ylabel("error/m")  
	plt.title("Location error") 
	plt.plot(time,x_error,"b-",label='x_direction') 
	plt.plot(time,y_error,"r-",label='y_direction') 
	plt.plot(time,z_error,"g-",label='z_direction') 
	plt.plot(time,x_error, 'bo')
	plt.plot(time,y_error, 'ro')
	plt.plot(time,z_error, 'go')
	plt.grid(True) #添加格子
	#plt.legend(loc = 0) #图例位置自动
	plt.figure(2)
	plt.xlabel("time")  
	plt.ylabel("distance error/m")  
	plt.title("Location distance error") 
	plt.plot(time,error,"m-",lw = 1.5)
	plt.plot(time,error, 'mo')
	plt.show() 
'''----------------------------------------------------------------------------------------------------------------------------'''
#主函数
'''----------------------------------------------------------------------------------------------------------------------------'''
if __name__ == '__main__':
	eph,noeph,ALPHA,BETA=readN('D:/pygps/gold2950.17n','D:/pygps/eph.txt')
	Ocontent,sates,year,month,day,nx,ny,nz=readO('D:/pygps/gold2950.171o','D:/pygps/Ocontent.txt')
	error=np.zeros((20))
	x_error=np.zeros((20))
	y_error=np.zeros((20))
	z_error=np.zeros((20))
	for time_num in range(20):
		ST=52650+30*time_num
		sec_of_week=TIMEchange(year,month,day,ST)
		sateeph_new,Ocontent_new,ture_O_satenum=getsat(eph,noeph,Ocontent,sec_of_week,sates,ST)
		tx_GPS,tcorr=SateTimecorrection(ture_O_satenum,sateeph_new,Ocontent_new,sec_of_week)
		Xdg,Ydg,Zdg=satelliteposition(sateeph_new,Ocontent_new,sates,ture_O_satenum,tx_GPS)
		xjsj,yjsj,zjsj,vT,jindu=SPLocation(tcorr,sateeph_new,Ocontent_new,ture_O_satenum,Xdg,Ydg,Zdg,nx,ny,nz)
		error[time_num]=jindu
		x_error[time_num]=xjsj-nx
		y_error[time_num]=yjsj-ny
		z_error[time_num]=zjsj-nz
	average_error=np.mean(error)
	print('20个历元的定位平均误差：',average_error)
	picture(error,x_error,y_error,z_error)