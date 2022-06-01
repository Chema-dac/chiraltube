import sys
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import random as rd
import statistics as stat

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

ERROR=0.025

class punto3d:
    def __init__(self, x, y, z, ele):
        self.x=x #coordenadas en X, Y y Z
        self.y=y
        self.z=z
        self.ele=ele #elemento
    def __str__(self):
        return "{}, {}, {}, {}".format(self.x, self.y, self.z, self.ele)
    def __repr__(self):
        return "({}, {}, {}, {})".format(self.x, self.y, self.z, self.ele)
    def __eq__(self, other): #igualdad
        return self.x==other.x and self.y==other.y and self.z==other.z
    def __add__(self, other): #suma
        return punto3d(self.x+other.x, self.y+other.y, self.z+other.z, other.ele)
    def __mul__(self, other): #multiplicacion
        if isinstance(self, punto3d):
            if isinstance(other, punto3d):
                return punto3d(self.x*other.x, self.y*other.y, self.z*other.z, self.ele)
            else:
                return punto3d(self.x*other, self.y*other, self.z*other, self.ele)
        elif isinstance(other, punto3d):
            return punto3d(self*other.x, self*other.y, self*other.z, other.ele)
        else:
            return punto3d(self*other.x, self*other.y, self*other.z, other.ele)
    def __sub__(self,other): #resta
        return punto3d(self.x-other.x, self.y-other.y, self.z-other.z, other.ele)
    def mag(self):#magnitud del vector
        return math.sqrt(self.x**2+self.y**2+self.z**2)


def distancia(p1, p2):
    return math.sqrt((p2.x-p1.x)**2 + (p2.y-p1.y)**2 + (p2.z-p1.z)**2)

def testgraf(Arr):
    x=[ele.x for ele in Arr]
    y=[ele.y for ele in Arr]
    z=[ele.z for ele in Arr]
    C=[]
    for ele in Arr:
        if ele.ele=='S':
            C.append('b')
        elif ele.ele=='C':
            C.append('g')
        else:
            C.append('r')
    fig=plt.figure()
    ax=plt.axes(projection='3d')
    ax.scatter3D(x,y,z, c=C)
    plt.show()
    

def arr_inicial(n, m, x, y):
    '''Devuelve el arreglo inicial'''
    y_tot=n+x+1
    x_inf=-y-1
    x_sup=m+1
    Arr=[]
    py=0
    pz=0
    for i in range(y_tot+1):
        px=a2*x_inf+a1*i
        for _ in range(x_sup-x_inf+1):
            for atom in UnitCell:
                Arr.append(px+atom)            
            px=px+a2
    return Arr

def robtenerxy(n, m, error=ERROR, MAX=100):
    '''Devuelve el ángulo theta, el arreglo de pares X, Y.'''
    Res=[]
    if n==0:
        theta=0
        for x in range(1,MAX):
            y=x*a1.mag()*math.cos(phi)/a2.mag()
            if abs(round(y)-y)<error:
                Res.append((x, y, abs(round(y)-y), (abs(round(y)-y)/y)*100))
        for y in range(1,MAX):
            x=y*a2.mag()/(a1.mag()*math.cos(phi))
            if abs(round(x)-x)<error:
                Res.append((x, y, abs(round(x)-x), (abs(round(x)-x)/x)*100))
    elif m==0:
        theta=phi
        for x in range(1,MAX):
            y=x*a1.mag()/(a2.mag()*math.cos(phi))
            if abs(round(y)-y)<error:
                Res.append((x, y, abs(round(y)-y), (abs(round(y)-y)/y)*100))
        for y in range(1,MAX):
            x=y*a2.mag()*math.cos(phi)/a1.mag()
            if abs(round(x)-x)<error:
                Res.append((x, y, abs(round(x)-x), (abs(round(x)-x)/x)*100))
    else:
        theta= math.acos((a2.x*(a2*m+a1*n).x + a2.y*(a2*m+a1*n).y+a2.z*(a2*m+a1*n).z)/(a2.mag()*(a2*m+a1*n).mag()))
        for x in range(1, MAX):
            #y=((n*a1.mag()*a1.mag())/(m*a2.mag()*a2.mag()))*(math.sin(math.pi-2*phi+theta)/math.sin(theta))*x
            y=((n*a1.mag()*a1.mag())/(m*a2.mag()*a2.mag()))*(math.sin(phi-theta)*math.cos(phi-theta))*x/(math.sin(theta)*math.cos(theta))
            if abs(round(y)-y)<error:
                Res.append((x, y, abs(round(y)-y), (abs(round(y)-y)/y)*100))
        for y in range(1, MAX):
            #x=((m*a2.mag()*a2.mag())/(n*a1.mag()*a1.mag()))*(math.sin(theta)/math.sin(math.pi-2*phi+theta))*y
            x=((m*a2.mag()*a2.mag())/(n*a1.mag()*a1.mag()))*(math.sin(theta)*math.cos(theta))*y/(math.sin(phi-theta)*math.cos(phi-theta))
            if abs(round(x)-x)<error:
                Res.append((x, y, abs(round(x)-x), (abs(round(x)-x)/x)*100))
    Res.sort(key=sortfuncT)
    return theta, Res

def obtenerxy(n, m, error=ERROR, MAX=100):
    '''Imprime el angulo correspondiente a n y m, y encuentra los pares X y Y que cumplen la relacion con un error menor a "error".
        Busca numeros enteros Y y Y que sean menores a MAX.'''
    theta, Res= robtenerxy(n,m, error, MAX)
    print("Ángulo: ", math.degrees(theta),  "\n\nResultados a un error <", error)
    print("X,    Y,    Error absoluto, Error porcentual")
    for ele in Res:
        x,y,err,por=ele
        print("{}, {}, {}, {}%".format(round(x,4),round(y,4), err, por))
        
def obtenerpq(n, m, x,y, MAX=100):
    '''Devuelve los indices del vecotr R: P y Q'''
    N=x*m + y*n
    for p in range(MAX):
        q=(1-y*p)/x
        if int(q)==q and (0< m*p - n*q <= N):
            break
    return p, q

def rotar(Arr, theta):
    ''' Rota todos los puntos del arreglo Arr, phi-theta radianes al contrario de las manecillas del reloj. Devuelve el arreglo rotado'''
    for ele in Arr:
        tmpx=ele.x
        tmpy=ele.y
        ele.x=tmpx*math.cos(phi-theta)-tmpy*math.sin(phi-theta)
        ele.y=tmpx*math.sin(phi-theta)+tmpy*math.cos(phi-theta)
    return Arr

def eliminar(Arr, n, m, x, y):
    '''Elimina los puntos en Arr que no son parte del "listón". Devuelve el arreglo, la distancia en Y, la distancia en X.'''
    nuevo_Arr=[]
    disy= distancia(punto3d(0,0,0,0), (a1*x-a2*y))
    disx= distancia(punto3d(0,0,0,0), (a1*n+a2*m))
    for ele in Arr:
        if ele.y>=-0.000001 and ele.y<(disy-0.000001) and ele.x>=0 and ele.x<disx:
            nuevo_Arr.append(ele)
    return nuevo_Arr, disy, disx

def coordenadas(n, m, x, y, prt=1, arch_out=False, grf=False):
    if not arch_out:
        arch_out= "./Ribbon({},{}).xyz".format(n,m)
    theta, Res = robtenerxy(n, m)
    if prt==1:
        print("\nÁngulo obtenido correctamente: Theta=", math.degrees(theta), "\n")
    Arr=arr_inicial(n, m, x, y)
    if prt==1:
        print("Arreglo inicial creado correctamente\n")
    Arr= rotar(Arr, theta)
    if prt==1:
        print("Rotación completada correctamente\n")
    Arr, disy, disx = eliminar(Arr, n, m, x, y)
    if prt==1:
        print("Listón creado correctamente, {} átomos en total.".format(len(Arr)))
    with open(arch_out, mode = 'w') as outfile:
        print(len(Arr), file=outfile)
        print("Listón sin doblar con índices n={} y m={}. Parámetro en y: {}. Parámetro en x: {}".format(n, m, disy, disx), file=outfile)
        for ele in Arr:
            print("{}    {}    {}    {}".format(ele.ele, ele.x, ele.y, ele.z), file=outfile)
    if grf:
        testgraf(Arr)
        
def nanotubo(Arr, n, m, z_med=False):
    '''Crea el nanotubo, devuelve el arreglo, el radio.'''
    x_max= distancia(punto3d(0,0,0,0), (a1*n+a2*m))
    radio= x_max/(2*math.pi)
    if z_med==False:
        z_med= (a3.z)/2        
    for ele in Arr:
        gamma= ele.x/radio
        tmpy=ele.y
        radio_rel= ele.z-z_med
        ele.x= (radio+radio_rel)*math.cos(gamma)
        ele.y= (radio+radio_rel)*math.sin(gamma)
        ele.z=tmpy
    return Arr, radio

def coordenadasNT(n, m, x, y, repeat=0,prt=1, arch_out=False, grf=False):
    extra=""
    theta, Res = robtenerxy(n, m)
    if prt==1:
        print("\nÁngulo obtenido correctamente: Theta=", math.degrees(theta), "\n")
    Arr=arr_inicial(n, m, x, y)
    if prt==1:
        print("Arreglo inicial creado correctamente\n")
    Arr= rotar(Arr, theta)
    if prt==1:
        print("Rotación completada correctamente\n")
    Arr, disy, disx = eliminar(Arr, n, m, x, y)
    if prt==1:
        print("Listón creado correctamente, {} átomos en total.\n".format(len(Arr)))
    Arr, radio=nanotubo(Arr, n, m)
    if prt==1:
        print("Nanotubo creado correctamente, {} átomos en total.\n".format(len(Arr)))
    if repeat!=0:
        Arr=repetir(Arr, disy, repeat)
        extra="x{}".format(repeat+1)
        print("Nanotubo repetido correctamente {} veces, {} átomos en total.\n".format(repeat, len(Arr)))
    if not arch_out:
        arch_out= "Nanotubo({},{}){}.xyz".format(n,m, extra)    
    with open(arch_out, mode = 'w') as outfile:
        print(len(Arr), file=outfile)
        print("Nanotube with chiral indices (n,m)=({},{}). Parameters: z={}, radius={}, X & Y={}. The translational indices used were (x,y)=({},{})".format(n, m, disy*(repeat+1), radio, 2*radio, x,y), file=outfile)
        for ele in Arr:
            print("{}    {}    {}    {}".format(ele.ele, ele.x, ele.y, ele.z), file=outfile)
    if grf:
        testgraf(Arr)
        
def repetir(Arr, disy, veces):
    Copia=[]
    for ele in Arr:
        Copia.append(ele)
    py=0
    for _ in range(veces):
        for ele in Copia:
            Arr.append(punto3d(0,0,py+disy,'na')+ele)
        py=py+disy
    return Arr

def buscar(tope, factor=1,  r=0):
    nRes=[]
    for n in range(1,tope):
        for m in range(1,tope):
            t, Res = robtenerxy(n, m)
            for ele in Res:
                x, y, err, errp=ele
                if x<=factor*n and y<=factor*m and n!=m:
                    nRes.append((n, m, x, y, err, errp))
    print("n, m, X,    Y,    Error absoluto, Error porcentual")
    for ele in nRes:
        n, m, x, y, err, por=ele
        print("{}, {}, {}, {}, {}, {}".format(n, m, x,y, err, por))
    if r==1:
        return nRes

def buscarnum(limite, tope=10, r=0):
    nRes=[]
    for n in range(1,tope):
        for m in range(1,tope):
            t, Res = robtenerxy(n, m)
            x,y,err,perr=Res[0]
            x=round(x)
            y=round(y)
            Arr, disy,disx,rad =hacer_todo(n,m,x,y)
            if len(Arr)<=limite:
                nRes.append((len(Arr),n,m,x,y,err,perr))
    print("no., n, m, X,    Y,    Error absoluto, Error porcentual")
    for ele in nRes:
        num, n, m, x, y, err, por=ele
        print("{}, {}, {}, {}, {}, {}, {}".format(num, n, m, x,y, err, por))
    if r==1:
        return nRes


def printcoords(Arr, file_name, comment):
    with open(file_name, mode = 'w') as outfile:
        print(len(Arr), file=outfile)
        print(comment, file=outfile)
        for ele in Arr:
            print("{}    {}    {}    {}".format(ele.ele, ele.x, ele.y, ele.z), file=outfile)


def leer_arch(arch):
    A=[]
    UnitCell=[]
    with open(arch, mode='r') as infile:
        for line in infile:
            line=line.strip()
            if line.startswith("nat"):
                nat= line.replace("nat","").replace("=","")
                nat=int(nat.strip().strip(','))
            elif line.startswith("CELL_PARAMETERS"):
                for _ in range(3):
                    valores= infile.readline().strip().split()
                    A.append(punto3d(float(valores[0]), float(valores[1]), float(valores[2]), 'NA'))
            elif line.startswith("ATOMIC_POSITIONS"):
                units = line.replace("ATOMIC_POSITIONS", "").strip().strip('{}()')
                if units=='angstrom' or units=='Angstrom':
                    for _ in range(nat):
                        valores = infile.readline().strip().split()
                        UnitCell.append(punto3d(float(valores[1]), float(valores[2]), float(valores[3]), valores[0]))
                elif units=='crystal' or units =='Crystal':
                    print("\n***ERROR*** No se ha implementado la lectura de archivos de entrada con unidades en 'crystal'. Porfavor convierta las unidades a angstroms.\n")
                    sys.exit()
                else:
                    print("\n***ERROR*** No se conoce o no se ha implementado el uso de unidades '{}'.".format(units))
                    sys.exit()
    return nat, A, UnitCell

def shift_cell(A, UnitCell):
    Alturas=[]
    for ele in UnitCell:
        if ele.z not in Alturas:
            Alturas.append(ele.z)
    rango= max(Alturas)-min(Alturas)
    #prom=stat.mean(Alturas)
    punto_medio=min(Alturas)+rango/2
    altura_tot=A[2].z
    medio= altura_tot/2
    #distancia=medio-prom
    distancia=medio-punto_medio
    for ele in UnitCell:
        ele.z=ele.z+distancia
    for ele in UnitCell:
        if ele.z>altura_tot or ele.z<0:
            print("\n***ERROR*** Problemas con el eje z de la celda unitaria (Vector a3)")
            sys.exit()
    return UnitCell

def sortfunc(ele):
    x, y, err, perr=ele
    return math.sqrt(x*x + y*y)

def sortfuncT(ele):
    x,y,err,perr=ele
    return distancia(punto3d(0,0,0,0), a1*x+ a2*y)

def buscarnm(radio, error=0.1):
    Res = []
    for n in range(100):
        for m in range(100):
            rad_prueba= distancia(punto3d(0,0,0,0), (a1*n+a2*m))/(2*math.pi)
            if abs(rad_prueba-radio)<error:
                Res.append((n, m, abs(rad_prueba-radio), abs(rad_prueba-radio)/radio))
    Res.sort(key=sortfuncT)    
    return Res

def rbuscarnm(radio, error=0.1):
    Res = []
    for n in range(100):
        for m in range(100):
            rad_prueba= distancia(punto3d(0,0,0,0), (a1*n+a2*m))/(2*math.pi)
            if abs(rad_prueba-radio)<error:
                Res.append((n, m, rad_prueba,abs(rad_prueba-radio), abs(rad_prueba-radio)/radio))
    return Res
    
def hacer_todo(n,m,x,y, nt=True):
    '''Devuelve el arreglo, la altura y disx, y el radio'''
    theta, Res = robtenerxy(n,m)
    Arr, disy, disx = eliminar(rotar(arr_inicial(n,m,x,y), theta), n,m,x,y)
    radio=0
    if nt:
        Arr, radio = nanotubo(Arr, n,m)
    return Arr, disy, disx, radio       
    
def escalar(Arr, disy_vieja, disy_nueva):
    for ele in Arr:
        rel= ele.z/disy_vieja
        ele.z = rel*disy_nueva
    return Arr

def prueba(max_n, max_m, rad_min=2, disy_max=15):
    Resultados=[]
    for n in range(max_n):
        for m in range(max_m):
            print(n, " , ", m)
            if n==0 and m==0:
                continue
            if distancia(punto3d(0,0,0,0), (a1*n+a2*m))/(2*math.pi)<rad_min:
                continue
            t, Res = robtenerxy(n,m)
            x,y,err,perr= Res[0]
            x=round(x)
            y=round(y)
            Arr, Arr_int, Arr_ext, radio, disy, disy_int, disy_ext = ajustar(n,m,x,y,0)
            if disy<=disy_max and disy_int<=disy_max and disy_ext<=disy_max:
                Resultados.append((n,m,round(disy,4), round(disy_int,4), round(disy_ext,4)))
                print("append")
    return Resultados
    

def ajustar(n,m,x,y, repeat):
    theta, Res= robtenerxy(n,m)
    Arr, disy, disx = eliminar(rotar(arr_inicial(n,m,x,y), theta), n,m,x,y)
    
    ntArr, ntdisy, ntdisx = eliminar(rotar(arr_inicial(n,m,x,y), theta), n,m,x,y)
    ntArr, radio = nanotubo(ntArr, n,m)
    ntArr= repetir(ntArr, ntdisy, repeat)
##    ntArr2, ntdisy2, ntdisx2 = eliminar(rotar(arr_inicial(n,m,x,y), theta), n,m,x,y)
##    ntArr2, radio2 = nanotubo(ntArr2, n, m, z_med=h)
##    ntArr2= repetir(ntArr2, ntdisy2, repeat)

    ##Anillo interior:
    radio_int= radio + (3.434413167-5.0)
    error=0.01
    Res= buscarnm(radio_int, error=error)
    while len(Res)==0:
        error=error+0.01
        Res= buscarnm(radio_int, error=error)
    if len(Res)>=2:
        n_int, m_int, err_int, perr_int = Res[0]
        n_int2, m_int2, err_int2, perr_int2 = Res[1]
        theta_int, tmpArr =robtenerxy(n_int, m_int)
        x1, y1, tmp, ptmp=tmpArr[0]
        x1=round(x1)
        y1=round(y1)
        theta_int2, tmpArr2 =robtenerxy(n_int2, m_int2)
        x2, y2, tmp, ptmp=tmpArr2[0]
        x2=round(x2)
        y2=round(y2)
        if (x1*x1 + y1*y1) > (x2*x2 + y2*y2):
            n_int=n_int2
            m_int=m_int2
            x1=x2
            y1=y2
            theta_int=theta_int2
    else:
        n_int, m_int, err_int, perr_int = Res[0]
        theta_int, tmpArr= robtenerxy(n_int, m_int)
        x1, y1, tmp, ptmp =tmpArr[0]
        x1=round(x1)
        y1=round(y1)
    Arr_int, disy_int, disx_int= eliminar(rotar(arr_inicial(n_int,m_int,x1,y1), theta_int), n_int,m_int,x1,y1)
    difradio_int=radio_int-distancia(punto3d(0,0,0,0), (a1*n_int+a2*m_int))/(2*math.pi)
    Arr_int, rad_int = nanotubo(Arr_int, n_int, m_int, z_med=3.434413167-difradio_int)
    print("Interior:" , n_int, ' , ', m_int, ' , ', x1, ' , ', y1)
    ##Anillo exterior:

    radio_ext= radio + (6.565586833-5.0)
    error=0.01
    Res= buscarnm(radio_ext, error=error)
    while len(Res)==0:
        error=error+0.01
        Res= buscarnm(radio_ext, error=error)
    if len(Res)>=2:
        n_int, m_int, err_int, perr_int = Res[0]
        n_int2, m_int2, err_int2, perr_int2 = Res[1]
        theta_int, tmpArr =robtenerxy(n_int, m_int)
        x1, y1, tmp, ptmp=tmpArr[0]
        x1=round(x1)
        y1=round(y1)
        theta_int2, tmpArr2 =robtenerxy(n_int2, m_int2)
        x2, y2, tmp, ptmp=tmpArr2[0]
        x2=round(x2)
        y2=round(y2)
        if (x1*x1 + y1*y1) > (x2*x2 + y2*y2):
            n_int=n_int2
            m_int=m_int2
            x1=x2
            y1=y2
            theta_int=theta_int2
    else:
        n_int, m_int, err_int, perr_int = Res[0]
        theta_int, tmpArr= robtenerxy(n_int, m_int)
        x1, y1, tmp, ptmp =tmpArr[0]
        x1=round(x1)
        y1=round(y1)
    
    Arr_ext, disy_ext, disx_ext= eliminar(rotar(arr_inicial(n_int,m_int,x1,y1), theta_int), n_int,m_int,x1,y1)
    difradio_ext=radio_ext-distancia(punto3d(0,0,0,0), (a1*n_int+a2*m_int))/(2*math.pi)
    Arr_ext, rad_ext = nanotubo(Arr_ext, n_int, m_int, z_med=6.565586833-difradio_ext)
    print("Exterior:" , n_int, ' , ', m_int, ' , ', x1, ' , ', y1)
    
    return ntArr, Arr_int, Arr_ext, radio, ntdisy, disy_int, disy_ext

def juntar(Medio, Interior, Exterior, radio, disy):
    delta=0.1
    radio_int = radio + (3.434413167-5.0)
    radio_ext = radio + (6.565586833-5.0)
    NTfinal=[]
    NTcortado=[]
    for ele in Medio:
        magnitud = math.sqrt(ele.x*ele.x + ele.y*ele.y)
        if magnitud>=(radio-delta) and magnitud<=(radio+delta):
            NTfinal.append(ele)
    for ele in Interior:
        magnitud = math.sqrt(ele.x*ele.x + ele.y*ele.y)
        if magnitud>=(radio_int-delta) and magnitud<=(radio_int+delta):
            NTfinal.append(ele)
    for ele in Exterior:
        magnitud = math.sqrt(ele.x*ele.x + ele.y*ele.y)
        if magnitud>=(radio_ext-delta) and magnitud<=(radio_ext+delta):
            NTfinal.append(ele)
    for ele in NTfinal:
        if ele.z<=disy:
            NTcortado.append(ele)
    return NTfinal, NTcortado
        
def dividir(n,m,x,y):
    theta, Res= robtenerxy(n,m)
    Arr, disy, disx = eliminar(rotar(arr_inicial(n,m,x,y), theta), n,m,x,y)
    Alturas=[]
    for ele in Arr:
        if ele.z not in Alturas:
            Alturas.append(ele.z)
    Alturas.sort()
    capas=len(Alturas)
    Capas=[]
    for i in range(capas):
        Tmp=[]
        Capas.append(Tmp)  
    for ele in Arr:
        for i in range(capas):
            if ele.z==Alturas[i]:
                Capas[i].append(ele)
    return Capas, disy, disx

def ajustar_dif(Capas, disx, suma, sub):
    capas=len(Capas)
    X=[]
    NT=[]
    for i in range(capas):
        Tmp=[]
        X.append(Tmp)
    for i in range(capas):
        for ele in Capas[i]:
            X[i].append(ele.x)
    for i in range(len(X)):
        X[i].sort()
        
    radio=disx/(2*math.pi)
    z_med=(a3.z)/2
    
    #anillo interior:
    for i in range(sub):
        x_last=X[0].pop()
        Capas[0].pop()
        disx_int= x_last-X[0][0]
        
    rad_int=disx_int/(2*math.pi)
    for ele in Capas[0]:
        gamma= ele.x/rad_int
        tmpz=ele.y
        radio_rel=ele.z-z_med
        tmpx= (radio+radio_rel)*math.cos(gamma)
        tmpy= (radio+radio_rel)*math.sin(gamma)
        NT.append(punto3d(tmpx, tmpy, tmpz, ele.ele))
        
    #anillo medio:
    for ele in Capas[1]:
        gamma= ele.x/radio
        tmpz=ele.y
        radio_rel=ele.z-z_med
        tmpx= (radio+radio_rel)*math.cos(gamma)
        tmpy= (radio+radio_rel)*math.sin(gamma)
        NT.append(punto3d(tmpx, tmpy, tmpz, ele.ele))
        
    #anillo exterior:
    for i in range(suma):
        dist_entre_atomos= X[2][1]-X[2][0]
        X[2].append(X[2][-1]+dist_entre_atomos)
        Capas[2].append(punto3d(dist_entre_atomos*2,0,0,'NA')+Capas[2][-2])
        disx_ext= disx+ dist_entre_atomos

    rad_ext=disx_ext/(2*math.pi)
    for ele in Capas[2]:
        gamma= ele.x/rad_ext
        tmpz=ele.y
        radio_rel=ele.z-z_med
        tmpx= (radio+radio_rel)*math.cos(gamma)
        tmpy= (radio+radio_rel)*math.sin(gamma)
        NT.append(punto3d(tmpx, tmpy, tmpz, ele.ele))
    
    
    return X, Capas, NT, radio
    

def main():
    ribbon=False
    r=0
    graph=False
    Error=ERROR
    
    if "-l" in opts:
        ribbon=True
        
    if "-g" in opts:
        graph=True
        
    for word in opts:
        if word.startswith("-r"):
            r= int(word.strip().replace("-r", ""))
        if word.startswith("-E"):
            Error= float(word.strip().replace("-E", ""))
            
    n= int(input("\nEscriba el valor de n: "))
    m= int(input("\nEscriba el valor de m: "))
    
    if "-o" in opts:
        obtenerxy(n,m,error=Error)
        sys.exit()
        
    try:
        arch_out=args[1]
    except IndexError:
        print("\n***No se especificó archivo de salida. Se usará un nombre genérico***\n\n")
        arch_out=False
    
    theta, Res = robtenerxy(n, m, error=Error)
    try:
        x, y, err, perr=Res[0]
    except IndexError:
        print("\n***ERROR*** No se encontraron soluciones enteras para los indices del vector traslacional en un rango de 0 a 100\nTrate con otro par de índices quirales.")
        sys.exit
    x= round(x)
    y= round(y)
    
    if ribbon:
        coordenadas(n,m,x,y,arch_out=arch_out, grf=graph)
    else:
        coordenadasNT(n,m,x,y,repeat=r, arch_out=arch_out, grf=graph)
        



if __name__== "__main__":
    if ("-h" in opts) or ("--h" in opts) or ("-help" in opts):
        print("\nEste es un programa que imprime las coordenadas de diferentes estructuras en archivos .xyz para su utilización en visualizadores u otros programas.")
        print("Creado por José María de Albornoz Caratozzolo.")
        print("\nUtilización: \n\t\t python3 ~/nanotubos.py <comandos opcionales> <archivo de entrada> <archivo de salida>")
        print("\nSi no se especifíca un archivo de salida, se utilizará un nombre generado automáticamente.")
        print("El archivo de entrada debe contener únicamente la celda unitaria a partir de la cual se quiere crear el nanotubo o listón correspondiente. Este archivo debe estar en formato .in de Quantum Espresso, y únicamente utiliza los parámetros 'nat', 'CELL_PARAMETERS' y 'ATOMIC_POSITIONS'.")
        print("NOTA: Por el momento el programa sólo acepta ATOMIC_POSITIONS en Angstroms.")
        print("\nCOMANDOS OPCIONALES: (por orden de ejecución en el programa)")
        print("\tComando \t\tResultado\n")
        print("\t-h / -help / --h \tImprime este mensaje de ayuda.\n")
        print("\t-p \t\tImprime los parámetros de red y las posiciones atómicas\n\t\t\tde la celda unitaria. Utilizalo para checar que el\n\t\t\tprograma leyó correctamente el archivo de entrada.\n")
        print("\t-l \t\tDa las coordenadas del listón en vez de las coordenadas\n\t\t\tdel nanotubo.\n")
        print("\t-g \t\tAdemás de dar las coordenadas, imprime una gráfica del\n\t\t\tresultado. Requiere que esté instalada la librería\n\t\t\tde python 'matplotlib'.\n")
        print("\t-r<num> \tRepite el nanotubo <num> veces en z. <num> debe\n\t\t\tser un numero entero.\n\t\t\tEjemplo: '-r3' va a repetir el nanotubo 3 veces en z.\n")
        print("\t-E<num> \tCambia el error aceptable al buscar soluciones\n\t\t\tenteras para el vector traslacional a <num>\n\t\t\tEjemplo: '-E0.05' cambia el error aceptable a 0.05 Angstroms.\n")
        print("\t-o \t\tObtiene e imprime en pantalla los valores X y Y\n\t\t\tcorrespondientes a los indices n y m proporcionados.\n")
        sys.exit()
    try:
        arch_in=args[0]
    except IndexError:
        print("\n***ERROR*** No se especificó archivo de entrada.\nEscribe 'python3 nanotubos.py -h' para obtener ayuda.")
        sys.exit()
    try: nat, A, UnitCell =leer_arch(arch_in)
    except FileNotFoundError:
        print("\n***ERROR*** No se encontró el archivo '{}'\n".format(arch_in))
        sys.exit()
    UnitCell=shift_cell(A, UnitCell)
    print("\nCelda unitaria leída del archivo de entrada {}\n".format(arch_in))
    if "-p" in opts:    
        print("\nnat=", nat)
        print("\nCell Parameters:\n a1={}\n a2={}\n a3={}".format(A[0], A[1], A[2]))
        print("\nUnit Cell:")
        for ele in UnitCell:
            print("", ele)
        
        
    a1= A[0]
    a2= A[1]
    a3= A[2]
    a2= a2*(-1)
    for i in range(len(UnitCell)):
        UnitCell[i]=a2+UnitCell[i]
    phi= math.acos((a1.x*a2.x+a1.y*a2.y+a1.z*a2.z)/(a1.mag()*a2.mag()))


    if "-p" in opts:
        sys.exit()
    main()

    
        
