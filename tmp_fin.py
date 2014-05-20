# -*- charset: UTF-8 -*-
import sys
from math import sin, pi
def gauss(dim, h):
    """
    """

    u = [[0 for i in range(dim)] for i in range(dim)]
    rho = [[0 for i in range(dim)] for i in range(dim)]
    
    for i in range(dim):
        u[0][i] = cont_u
        u[dim-1][i] = cont_d
        u[i][0] = cont_l
        u[i][dim-1] = cont_r
    
    #rho[int(0.23*dim)][int(0.23*dim)] = char_f
    #rho[int(0.67*dim)][int(0.67*dim)] = char_s
    rho[int(0.63*dim)][int(0.37*dim)] = char_f
    rho[int(0.37*dim)][int(0.63*dim)] = char_f
    rho[int(0.37*dim)][int(0.37*dim)] = char_s
    rho[int(0.63*dim)][int(0.63*dim)] = char_s    


    cont_gauss = 0
    cont_converg = 0
    error = 0


    while cont_converg < 31:
        for i in range(1, dim-1):
            for j in range(1, dim-1):
                u_old = u[i][j]
                u[i][j] = 1./4.*(u[i-1][j] + u[i+1][j] + u [i][j-1] + u[i][j+1] \
                                    + h**2*rho[i][j])
                mean = (u[i][j] + u_old)/2
                if mean != 0:
                    error += abs((u_old - u[i][j])/mean)


        cont_gauss += 1
        converg = error <= 1.E-3
        if converg:
            cont_converg += 1
        sys.stdout.write("\r")
        sys.stdout.write("[%s] %s %s" % (converg, error, cont_gauss))
        sys.stdout.flush()
        error = 0
        

    if operation == "1":

        file_gauss = open('gauss.dat', 'w')

        for i in range(dim):
            for j in range(dim):
                file_gauss.write(str(i))
                file_gauss.write("\t")
                file_gauss.write(str(j))
                file_gauss.write("\t")
                file_gauss.write(str(u[i][j]))
                file_gauss.write("\n")

        file_gauss.close()


    return cont_gauss



def sor(dim, h, omega):
    """
    """

    u = [[0 for i in range(dim)] for i in range(dim)]
    rho = [[0 for i in range(dim)] for i in range(dim)]
    
    for i in range(dim):
        u[0][i] = cont_u
        u[dim-1][i] = cont_d
        u[i][0] = cont_l
        u[i][dim-1] = cont_r

    #rho[int(0.23*dim)][int(0.23*dim)] = char_f
    #rho[int(0.67*dim)][int(0.67*dim)] = char_s
    rho[int(0.63*dim)][int(0.37*dim)] = char_f
    rho[int(0.37*dim)][int(0.63*dim)] = char_f
    rho[int(0.37*dim)][int(0.37*dim)] = char_s
    rho[int(0.63*dim)][int(0.63*dim)] = char_s  


    cont_sor = 0
    cont_converg = 0
    error = error_old =  0


    while cont_converg < 31:
        for i in range(1, dim-1):
            for j in range(1, dim-1):
                u_old = u[i][j]
                u[i][j] = omega/4.*(u[i-1][j] + u[i+1][j] + u[i][j-1]\
                                   + u[i][j+1] + h**2*rho[i][j])\
                                   + (1 - omega)*u[i][j]
                mean = (u[i][j] + u_old)/2
                if mean != 0:
                    error += abs((u_old - u[i][j])/mean)


        cont_sor += 1
        converg = error <= 1.E-3
        if converg:# or error_old == error:
            cont_converg += 1

        
        sys.stdout.write("\r")
        sys.stdout.write("[%s] %s %s" % (converg, error, cont_sor))
        sys.stdout.flush()
        #error_old = error
        error = 0

    if operation == "1":

        file_sor = open('sor.dat', 'w')

        for i in range(dim):
            for j in range(dim):
                file_sor.write(str(i))
                file_sor.write("\t")
                file_sor.write(str(j))
                file_sor.write("\t")
                file_sor.write(str(u[i][j]))
                file_sor.write("\n")

        file_sor.close()


    return cont_sor



def calore(dim_s, dim_t):
    """DA FINIRE!!!"""
    temp = [[(10.) for i in range(dim_t)] for i in range(dim_s)]

    A = 10.
    B = 12.
    #for i in range(dim_t):
     #   temp[0][i] = 11.
      #  temp[dim_s-1][i] = A + B*sin(20*pi*float(i)/dim_t)
    for i in range(int(0.3*dim_s), int(0.7*dim_s)):
        temp[i][0] = t_0


    for time in range(1, dim_t):
        for space in range(1, dim_s-1):

            temp[space][time] = temp[space][time-1] + tau/(2.*t_s)*(temp[space+1][time-1] + temp[space-1][time-1]\
                                                    - 2*temp[space][time-1])

    file_calore = open('calore.dat', 'w')
    for i in range(dim_s):
        for j in range(int(dim_t-(9*dim_t/10.))):
            file_calore.write(str(i))
            file_calore.write("\t")
            file_calore.write(str(j))
            file_calore.write("\t")
            file_calore.write(str(temp[i][j]))
            file_calore.write("\n")
    file_calore.close()
    
    file_stag = open('stag.dat', 'w')
    for j in range(200,int(10*(dim_t/10.-1)), 91):
        for i in range(dim_s):


            file_stag.write(str(i))
            file_stag.write("\t")
            file_stag.write(str(temp[i][j]))
            file_stag.write("\n")
    file_stag.close()

def crost_terr():
    global dim_t, dim_s, h, t_i, t_s, tau

    dim_t = 365*10
    dim_s = 20
    h = 1
    K = 1.0
    t_i = t_0 =10
    t_s = h**2/(2*K)
    tau = 0.2*t_s
    calore(dim_s, dim_t)

def confr():
    global cont_u, cont_d, cont_l, cont_r
    global char_f, char_s

    cont_u = cont_d = 0.0
    cont_l = cont_r = -0.0

    char_f = 5.
    char_s = -5.

    print """Opzioni:
             1) Confronto metodo Gauss - SOR al variare della grandezza della matrice
             2) Calcolo omega ottimale metodo SOR"""
    ans = raw_input(">>>")


    while True:

        if ans == "1":

            file_conf = open('conf_tmp.dat', 'w')

            for dim in range(30, 200, 10):

                print "Dimensione matrice:\t", dim
                cont_gauss = gauss(dim, 1./dim)
                print "Numero pass. Gauss:\t", cont_gauss
                cont_sor = sor(dim, 1./dim, 1.9)
                print "Numero pss. S.O.R.:\t", cont_sor                
                print


                file_conf.write(str(dim))
                file_conf.write("\t")
                file_conf.write(str(cont_gauss))
                file_conf.write("\t")
                file_conf.write(str(cont_sor))
                file_conf.write("\n")

            file_conf.close()

            break


        elif ans == "2":

            file_omega = open('omega.dat', 'w')
            #pass_gauss = sor(50, 1./50, 1)

            for k in range (1000, 2000, 50):

                omega = k/1000.
                cont_sor = sor(100, 1./100, omega)
                print "Omega:\t", omega, "\tNumero pass.:\t", cont_sor

                file_omega.write(str(omega))
                file_omega.write("\t")
                file_omega.write(str(cont_sor))
                file_omega.write("\n")

            file_omega.close()

            break


        else:
            print "Opzione non valida"


            
print """Balblabla
1) Studio distribuzione di carica (Poisson)
2) Studio eq. calore
3) Crosta terr.
4) Ottimiz."""

while True:

    operation = raw_input(">>>")

    if operation == "1":

        print "Inserire dimensione matrice:"
        dim = int(raw_input(">>>"))

        # print "Inserire h etc etc etc"

        print "Inserire condizioni a contorno:"
        cont_u = float(raw_input("\t- Lato superiore:\t"))
        cont_d = float(raw_input("\t- Lato inferiore:\t"))
        cont_l = float(raw_input("\t- Lato sinistro:\t"))
        cont_r = float(raw_input("\t- Lato destro:\t"))
        
        print "Inserire cariche e posizione:"
        char_f = float(raw_input("\t- Prima carica:\t"))
        char_s = float(raw_input("\t- Seconda carica:\t"))
        
        print "Scegliere metodo:"
        print "\t1) Gauss-Seidel"
        print "\t2) SOR"

        while True:

            ans = raw_input(">>>")

            if ans == "1":
                gauss(dim, h)
                break


            elif ans == "2":
                print "Inserire omega:"
                omega = float(raw_input(">>>"))
                sor(dim, 1./dim, omega)
                break


            else:
                print "Opzione non valida"
            
        break


    elif operation == "2":

        print "Inserire dimensione griglia:"
        dim_s = int(raw_input("\t- spaziale:\t"))
        dim_t = int(raw_input("\t- temporale:\t"))
        h = float(raw_input("Inserire passo temporale:\t"))
        k = float(raw_input("Inserire coeff. diffusione:\t"))
        t_s = h**2/(2*k)
        print "Massimo passo temporale consigliato %s" %t_s
        tau = float(raw_input("Inserire passo temporale:\t"))

        t_0 = float(raw_input("Inserire temperatura iniziale:\t"))
        if tau <= t_s:
            print "Sistema stabile."
            calore(dim_s, dim_t)
        else:
            print "Sistema instabile."

        break


    elif operation == "3":
        crost_terr()
        break


    elif operation == "4":
        confr()
        break


    else:
        print "Opzione non valida"
