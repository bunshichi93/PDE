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
    
    rho[int(0.8*dim)][int(0.8*dim)] = char_f
    rho[int(0.2*dim)][int(0.2*dim)] = char_s

    sum_old = sum_new = 0.
    cont_gauss = 0
    converg = False


    while not converg:
        for i in range(1, dim-1):
            for j in range(1, dim-1):
                u[i][j] = 1./4.*(u[i-1][j] + u[i+1][j] + u [i][j-1] + u[i][j+1] \
                                    + h**2*rho[i][j])
                sum_new +=u[i][j]

        cont_gauss += 1
        flut = sum_new - sum_old
        perc = abs(flut/sum_new)*100
        converg = perc <= 1.E-03
        #converg = abs(flut) == 0
        #print flut, sum_new, sum_old
        sum_old = sum_new#(u)
        sum_new = 0

        #print flut, "\t", converg, "\t#", cont_gauss, "\r"


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

    u = [[1 for i in range(dim)] for i in range(dim)]
    rho = [[0 for i in range(dim)] for i in range(dim)]
    
    for i in range(dim):
        u[0][i] = cont_u
        u[dim-1][i] = cont_d
        u[i][0] = cont_l
        u[i][dim-1] = cont_r

    
    rho[int(0.8*dim)][int(0.8*dim)] = char_f
    rho[int(0.2*dim)][int(0.2*dim)] = char_s

    sum_old = sum_new = 10.
    cont_sor = 0
    converg = False


    while not converg:
        for i in range(1, dim-1):
            for j in range(1, dim-1):
                u[i][j] = omega/4.*(u[i-1][j] + u[i+1][j] + u[i][j-1]\
                                   + u[i][j+1] + h**2*rho[i][j])\
                                   + (1 - omega)*u[i][j]
                sum_new +=u[i][j]

        cont_sor += 1
        flut = sum_new - sum_old
        perc = abs(flut/sum_new)*100
        converg = perc <= 1.E-03
        #converg = abs(flut) == 0
        #print flut, sum_new, sum_old
        sum_old = sum_new#(u)
        sum_new = 0

        #print flut, "\t", converg, "\t#", cont_sor, "\r"
        

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

    char_f = 5.# 5
    char_s = -5.# -5

    print """Opzioni:
             1) Confronto metodo Gauss - SOR al variare della grandezza della matrice
             2) Calcolo omega ottimale metodo SOR"""
    ans = raw_input(">>>")


    while True:

        if ans == "1":

            file_conf = open('conf.dat', 'w')

            for dim in range(30, 100):
                cont_sor = sor(dim, 0.1, 1.9)
                cont_gauss = gauss(dim, 0.1)
                
                file_conf.write(str(dim))
                file_conf.write("\t")
                file_conf.write(str(cont_gauss))
                file_conf.write("\t")
                file_conf.write(str(cont_sor))
                file_conf.write("\n")
                
                print "Dimensione matrice:\t", dim
                print "Numero pass. Gauss:\t", cont_gauss
                print "Numero pss. S.O.R.:\t", cont_sor
                print

            file_conf.close()

            break


        elif ans == "2":

            file_omega = open('omega.dat', 'w')

            for omega in range (1000, 1999, 1):
                cont_sor = sor(50, 0.1, omega/1000.)

                file_omega.write(str(omega))
                file_omega.write("\t")
                file_omega.write(str(cont_sor))
                file_omega.write("\n")
                
                print "Omega:\t", omega/1000., "\tNumero pass.:\t", cont_sor

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
        dim = raw_input(">>>")

        # print "Inserire h etc etc etc"

        print "Inserire condizioni a contorno:"
        cont_u = raw_input("\t- Lato superiore:\t")
        cont_d = raw_input("\t- Lato inferiore:\t")
        cont_l = raw_input("\t- Lato sinistro:\t")
        cont_r = raw_input("\t- Lato destro:\t")
        
        print "Inserire cariche:"
        char_f = raw_input("\t- Prima carica:\t")
        char_s = raw_input("\t- Seconda carica:\t")
        
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
                omega = raw_input(">>>")
                sor(dim, h, omega)
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
