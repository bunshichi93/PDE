# -*- charset: UTF-8 -*-
from pylab import sum

def gauss(dim, h):
    """
    """

    u = [[1 for i in range(dim)] for i in range(dim)]
    rho = [[0 for i in range(dim)] for i in range(dim)]
    
    for i in range(dim):
        u[0][i] = cont_u
        u[dim-1][i] = cont_d
        u[i][0] = cont_l
        u[i][dim-1] = cont_r
    
    rho[int(0.87*dim)][int(0.87*dim)] = char_f
    rho[int(0.19*dim)][int(0.19*dim)] = char_s


    cont_gauss = 0
    cont_converg = 0
    error = 0


    while cont_converg < 11:
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

    
    rho[int(0.87*dim)][int(0.87*dim)] = char_f
    rho[int(0.19*dim)][int(0.19*dim)] = char_s


    cont_sor = 0
    cont_converg = 0
    error = 0


    while cont_converg < 11:
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
        if converg:
            cont_converg += 1
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



def calore():
    """DA FINIRE!!!"""

    for i in range(1, dim-1):
        for j in range(1, dim-1):
            temp[i+1][j] = temp[i][j] + tau/(2*ts)*(temp[i][j+1] + temp[i][j-1]\
                                                    - 2*temp[i][j])



def confr():
    global cont_u, cont_d, cont_l, cont_r
    global char_f, char_s

    cont_u = cont_d = 0.1
    cont_l = cont_r = -0.1

    char_f = 500.
    char_s = -500.

    print """Opzioni:
             1) Confronto metodo Gauss - SOR al variare della grandezza della matrice
             2) Calcolo omega ottimale metodo SOR"""
    ans = raw_input(">>>")


    while True:

        if ans == "1":

            file_conf = open('conf_tmp.dat', 'w')

            for dim in range(30, 200):

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

            for k in range (1000, 2005, 5):

                omega = k/1000.
                cont_sor = sor(100, 1./50, omega)
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
        #calore()
        break


    elif operation == "3":
        #ins
        break


    elif operation == "4":
        confr()
        break


    else:
        print "Opzione non valida"
