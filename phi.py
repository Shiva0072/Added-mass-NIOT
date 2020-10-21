import numpy as np
import sys

def matrix(npanels,face_coordn,face_norml,mode,face_ar):
	a_Fc=np.zeros((npanels,npanels))
	b_Fc=np.zeros(npanels)
	for i in range(0,npanels):
		r_fp=face_coordn[i]
		b_j=0
		for j in range(0,npanels):
			if j==i:
				a_Fc[i,j]=2*3.1414
				continue;
			r_sp=face_coordn[j]
			r=r_fp-r_sp
			n=face_norml[j]
			a_Fc[i,j]=(1*r[mode]*n[mode]*face_ar)/((sum(r**2)**1.5)*(sum(n**2)**0.5))
			b_j=b_j+((n[mode]*face_ar)/(sum(r**2)**0.5))
		b_Fc[i]=1*b_j
	return(a_Fc,b_Fc)

if __name__=="__main__":
#extracting the no. in between
	a=int(input('Enter Value of length : '))
	b=int(input('Enter value of breadth : '))
	c=int(input("Enter value of height : "))
	print("Density of water used : 1.025tonnes/m^3");
	print("....NOTE....")
	print("X is along front and back face")
	print("Y is along left and right face")
	print("Z is along bottom and top face")
	print("values of modes, heave=2, sway=1.surge=0");
	mode_k = int(input('Enter the mode : '))

	if mode_k >= 3:
		print("Wrong input of mode !, exiting")
		exit()

	x=np.linspace(0,1*a, (2*a+1))                    
	y=np.linspace(0,1*b, (2*b+1))
	z=np.linspace(0,1*c, (2*c+1))

	xx=np.array([x[i] for i in range(0,len(x)) if(i%2!=0)])
	yy=np.array([y[i] for i in range(0,len(y)) if(i%2!=0)])
	zz=np.array([z[i] for i in range(0,len(z)) if(i%2!=0)])
	###panels centers coordinates storing
	Fc1=[]     #bottom
	Fc2=[]     #top
	Fc3=[]     #front
	Fc4=[]	   #back
	Fc5=[]	   #right
	Fc6=[]     #left

	for i in range(len(xx)):
		for j in range(len(yy)):
				Fc1.append(tuple((xx[i],yy[j],0*c)))	#bottom
				Fc2.append(tuple((xx[i],yy[j],1*c)))    #top
	for i in range(len(yy)):
		for j in range(len(zz)):
				Fc3.append(tuple((1*a,yy[i],zz[j])))	#front
				Fc4.append(tuple((0*a,yy[i],zz[j])))    #back
	for i in range(len(xx)):
		for j in range(len(zz)):
				Fc5.append(tuple((xx[i],1*b,zz[j])))	#right
				Fc6.append(tuple((xx[i],0*b,zz[j])))    #left

	Fc1=np.array(Fc1)       #bottom
	Fc2=np.array(Fc2)       #top
	Fc3=np.array(Fc3)       #front
	Fc4=np.array(Fc4)	 	#back
	Fc5=np.array(Fc5)	 	#right
	Fc6=np.array(Fc6)		#left
	##panels centers end

	#panels normals making
	n_Fc1=np.array([(0,0,-1) for i in range(0,len(Fc1))])   #bottom
	n_Fc2=np.array([(0,0,1) for i in range(0,len(Fc2))])	#top
	n_Fc3=np.array([(1,0,0) for i in range(0,len(Fc3))])	#front
	n_Fc4=np.array([(-1,0,0) for i in range(0,len(Fc4))])	#back
	n_Fc5=np.array([(0,1,0) for i in range(0,len(Fc5))])	#right
	n_Fc6=np.array([(0,-1,0) for i in range(0,len(Fc6))])	#left
	#panels normals end

	mode=mode_k           # 0 for surge,1 for sway, 2 for heave; 
	m=len(Fc1)       # no. of pannels m is same for all the faces
	dA=1;

	# for Face 1(bottom)
	a_Fc1,b_Fc1=matrix(len(Fc1), Fc1, n_Fc1, mode, dA)
	#for face 2(top)
	a_Fc2,b_Fc2=matrix(len(Fc2), Fc2, n_Fc2, mode, dA)
	#for face 3(front)
	a_Fc3,b_Fc3=matrix(len(Fc3), Fc3, n_Fc3, mode, dA)
	#for face 4(back)
	a_Fc4,b_Fc4=matrix(len(Fc4), Fc4, n_Fc4, mode, dA)
	#for face 5(right)
	a_Fc5,b_Fc5=matrix(len(Fc5), Fc5, n_Fc5, mode, dA)
	# for face 6(left)
	a_Fc6,b_Fc6=matrix(len(Fc6), Fc6, n_Fc6, mode, dA)
	
	# finding phi
	phi_Fc1=np.dot(np.linalg.inv(a_Fc1),b_Fc1)
	phi_Fc2=np.dot(np.linalg.inv(a_Fc2),b_Fc2)
	phi_Fc3=np.dot(np.linalg.inv(a_Fc3),b_Fc3)
	phi_Fc4=np.dot(np.linalg.inv(a_Fc4),b_Fc4)
	phi_Fc5=np.dot(np.linalg.inv(a_Fc5),b_Fc5)
	phi_Fc6=np.dot(np.linalg.inv(a_Fc6),b_Fc6)

	
	print("Phi of bottom Face(face_1) : ",phi_Fc1)
	print('\r')
	print("Phi of Top Face(face_2) : ",phi_Fc2)
	print('\r')
	print("Phi of Front Face(face_3) : ",phi_Fc3)
	print('\r')
	print("Phi of Back Face(face_4) : ",phi_Fc4)
	print('\r')
	print("Phi of Right Face(face_5): ",phi_Fc5)
	print('\r')
	print("Phi of Left Face(face_6) : ",phi_Fc6)
	print('\r')	



	print('Added Mass........')
	m_1_1=0    #bottom face added mass
	m_2_2=0    #top face added mass
	m_3_3=0    #front
	m_4_4=0    #back
	m_5_5=0    #right
	m_6_6=0		#left		

	rho=1.025

	if (mode==2):                           #heave (top-bottom)
		for i in range(len(phi_Fc1)):
			m_1_1=m_1_1+phi_Fc1[i]*(n_Fc1[i][mode])
		for i in range(len(phi_Fc2)):
			m_2_2=m_2_2+phi_Fc2[i]*(n_Fc2[i][mode])
		m_2_2=dA*rho*m_2_2
		m_1_1=dA*rho*m_1_1
		print('Heave Added mass : ',(m_1_1+m_2_2))
	elif (mode==1):                           #sway (left-right)
		for i in range(len(phi_Fc5)):
			m_5_5=m_5_5+phi_Fc5[i]*(n_Fc5[i][mode])
		for i in range(len(phi_Fc6)):
			m_6_6=m_6_6+phi_Fc6[i]*(n_Fc6[i][mode])
		m_6_6=dA*rho*m_6_6
		m_5_5=dA*rho*m_5_5
		print('Sway Added mass : ',(m_6_6+m_5_5))
	elif (mode==0):                           #surge (front-back)
		for i in range(len(phi_Fc3)):
			m_3_3=m_3_3+phi_Fc3[i]*(n_Fc3[i][mode])
		for i in range(len(phi_Fc4)):
			m_4_4=m_4_4+phi_Fc4[i]*(n_Fc4[i][mode])
		m_4_4=dA*rho*m_4_4
		m_3_3=dA*rho*m_3_3
		print('Surge Added mass : ',(m_3_3+m_4_4))
